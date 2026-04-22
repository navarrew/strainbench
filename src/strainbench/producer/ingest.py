"""Ingest parsed StrainRecords into the strainbench SQLite database.

The ingestion step is the bridge between the parser layer (which produces
in-memory `StrainRecord` objects) and everything downstream (clustering,
table views, web app). It performs two coordinated writes per strain:

  1. SQL: one row in `strains`, N rows in `cds`, in a single transaction.
  2. Filesystem: per-strain FASTA sidecars at
     `<fasta_dir>/fna/<locus_prefix>.fna` and
     `<fasta_dir>/faa/<locus_prefix>.faa`.

The schema's `UNIQUE (strain_id, locus_tag)` and `UNIQUE (locus_prefix)`
constraints prevent silent duplicates. Three policies for what to do when
a strain is already in the DB are exposed via `on_duplicate`:

  * 'error'   — raise IngestError (default; safest).
  * 'skip'    — leave the existing row untouched, do nothing.
  * 'replace' — delete the existing strain (cascades to its CDS rows
                via ON DELETE CASCADE) and re-insert.

Sidecar FASTA records use the format
    >LOCUSPREFIX|LOCUS_TAG [protein_id=...] [protein=...] description
so that mmseqs cluster output can be unambiguously joined back to a row in
the `cds` table via (locus_prefix, locus_tag).
"""

from __future__ import annotations

import sqlite3
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from strainbench.core.models import CDSRecord, NonCDSRecord, StrainRecord

_FASTA_LINE_WIDTH = 80


class IngestError(Exception):
    """Raised when a strain cannot be ingested (e.g., duplicate with on_duplicate='error')."""


@dataclass
class IngestResult:
    """Per-strain outcome from a batch ingest run."""

    locus_prefix: str
    source_file: str
    status: str  # 'inserted', 'replaced', 'skipped', 'failed'
    strain_id: int | None = None
    cds_count: int = 0
    error: str | None = None


def ingest_strain(
    conn: sqlite3.Connection,
    record: StrainRecord,
    fasta_dir: str | Path,
    *,
    on_duplicate: str = "error",
) -> tuple[int | None, str]:
    """Ingest one StrainRecord into the database and write its sidecar FASTAs.

    Args:
        conn: Open SQLite connection. The caller controls the connection
            lifecycle (open/close). This function manages its own transaction
            boundary via a SAVEPOINT-style `with conn:` block.
        record: A parsed StrainRecord (e.g., from `parse_gbff`).
        fasta_dir: Parent directory for sidecar FASTA files. Subdirectories
            `fna/` and `faa/` are created if needed.
        on_duplicate: One of 'error', 'skip', 'replace'.

    Returns:
        (strain_id, status) where status is one of 'inserted', 'replaced',
        'skipped'. strain_id is None for 'skipped'.

    Raises:
        IngestError: if the strain already exists and on_duplicate='error',
            or if on_duplicate is unrecognized.
    """
    if on_duplicate not in {"error", "skip", "replace"}:
        raise IngestError(f"invalid on_duplicate policy: {on_duplicate!r}")

    fasta_dir = Path(fasta_dir)
    fna_path = fasta_dir / "fna" / f"{record.locus_prefix}.fna"
    faa_path = fasta_dir / "faa" / f"{record.locus_prefix}.faa"
    rna_path = fasta_dir / "rna" / f"{record.locus_prefix}.fna"
    fna_path.parent.mkdir(parents=True, exist_ok=True)
    faa_path.parent.mkdir(parents=True, exist_ok=True)
    rna_path.parent.mkdir(parents=True, exist_ok=True)

    with conn:
        existing = conn.execute(
            "SELECT strain_id FROM strains WHERE locus_prefix = ?",
            (record.locus_prefix,),
        ).fetchone()

        status = "inserted"
        if existing is not None:
            if on_duplicate == "error":
                raise IngestError(
                    f"strain '{record.locus_prefix}' already exists "
                    f"(strain_id={existing['strain_id']}); "
                    f"use on_duplicate='skip' or 'replace'"
                )
            if on_duplicate == "skip":
                return None, "skipped"
            # replace: cascades to cds rows
            conn.execute(
                "DELETE FROM strains WHERE strain_id = ?",
                (existing["strain_id"],),
            )
            status = "replaced"

        strain_id = _insert_strain_row(conn, record)
        _insert_cds_rows(conn, strain_id, record.cds_records)
        _write_fasta(fna_path, record.cds_records, kind="nt", locus_prefix=record.locus_prefix)
        _write_fasta(faa_path, record.cds_records, kind="aa", locus_prefix=record.locus_prefix)
        # Sidecar: non-CDS features (tRNA/rRNA/ncRNA/tmRNA/CRISPR). Only
        # write the file if this strain actually has any — no point creating
        # 124 empty placeholder files for strains with no annotated RNAs.
        if record.non_cds_records:
            _write_rna_fasta(rna_path, record.non_cds_records,
                             locus_prefix=record.locus_prefix)
        elif rna_path.exists():
            # Re-ingest case: previous run had RNAs, this one doesn't (or
            # they got filtered). Drop the stale file.
            rna_path.unlink()

    return strain_id, status


def ingest_records(
    conn: sqlite3.Connection,
    records: Iterable[StrainRecord],
    fasta_dir: str | Path,
    *,
    on_duplicate: str = "error",
    progress: bool = False,
) -> list[IngestResult]:
    """Ingest a sequence of pre-parsed StrainRecords.

    Errors on individual strains are caught and recorded in the result list
    rather than aborting the batch — so one malformed strain doesn't sink
    a 1000-strain run.
    """
    results: list[IngestResult] = []
    for i, record in enumerate(records, start=1):
        result = IngestResult(
            locus_prefix=record.locus_prefix,
            source_file=record.source_file,
            status="failed",
            cds_count=len(record.cds_records),
        )
        try:
            strain_id, status = ingest_strain(
                conn, record, fasta_dir, on_duplicate=on_duplicate
            )
            result.strain_id = strain_id
            result.status = status
        except (IngestError, sqlite3.Error) as exc:
            result.error = str(exc)
        results.append(result)
        if progress:
            print(
                f"  [{i:4d}] {result.status:9s}  {record.locus_prefix:14s}  "
                f"({result.cds_count:5d} CDS)"
                + (f"  ERROR: {result.error}" if result.error else "")
            )
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Internals
# ─────────────────────────────────────────────────────────────────────────────


def _insert_strain_row(conn: sqlite3.Connection, record: StrainRecord) -> int:
    cur = conn.execute(
        """
        INSERT INTO strains (
            locus_prefix, species, strain_name, assembly_id,
            biosample_id, bioproject_id, assembly_level,
            source_format, source_file, cds_count
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        RETURNING strain_id
        """,
        (
            record.locus_prefix,
            record.species,
            record.strain_name,
            record.assembly_id,
            record.biosample_id,
            record.bioproject_id,
            record.assembly_level,
            record.source_format,
            record.source_file,
            len(record.cds_records),
        ),
    )
    return int(cur.fetchone()["strain_id"])


def _insert_cds_rows(
    conn: sqlite3.Connection,
    strain_id: int,
    cds_records: list[CDSRecord],
) -> None:
    rows: list[tuple[Any, ...]] = [
        (
            strain_id,
            c.locus_tag,
            c.protein_id,
            c.nuc_accession,
            c.location,
            c.direction,
            c.nt_length,
            c.aa_length,
            c.gc_pct,
            c.gene_name,
            c.annotation,
            c.notes,
        )
        for c in cds_records
    ]
    conn.executemany(
        """
        INSERT INTO cds (
            strain_id, locus_tag, protein_id, nuc_accession, location, direction,
            nt_length, aa_length, gc_pct, gene_name, annotation, notes
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )


def _write_fasta(
    path: Path,
    cds_records: list[CDSRecord],
    *,
    kind: str,
    locus_prefix: str,
) -> None:
    """Write per-strain FASTA. `kind` is 'nt' or 'aa'."""
    with path.open("w") as f:
        for c in cds_records:
            seq = c.nt_sequence if kind == "nt" else c.aa_sequence
            descriptor_bits = [f"[locus_prefix={locus_prefix}]"]
            if c.protein_id:
                descriptor_bits.append(f"[protein_id={c.protein_id}]")
            if c.gene_name:
                descriptor_bits.append(f"[gene={c.gene_name}]")
            descriptor_bits.append(f"[{c.annotation}]")
            f.write(f">{locus_prefix}|{c.locus_tag} {' '.join(descriptor_bits)}\n")
            for line_start in range(0, len(seq), _FASTA_LINE_WIDTH):
                f.write(seq[line_start : line_start + _FASTA_LINE_WIDTH] + "\n")


def _write_rna_fasta(
    path: Path,
    records: list[NonCDSRecord],
    *,
    locus_prefix: str,
) -> None:
    """Write per-strain non-CDS FASTA (tRNA/rRNA/ncRNA/tmRNA/CRISPR).

    Header format mirrors the CDS sidecars but tags each record with its
    feature_type and product so users can grep for e.g. '16S' or 'CRISPR'
    across all strains' files.

        >LOCUSPREFIX|LOCUS_TAG [tRNA] [locus_prefix=…] [product=tRNA-Ala(GGC)] [location=38420..38493]
        GGGGGCATAGCTC...
    """
    with path.open("w") as f:
        for r in records:
            descriptor_bits = [
                f"[{r.feature_type}]",
                f"[locus_prefix={locus_prefix}]",
            ]
            if r.product:
                descriptor_bits.append(f"[product={r.product}]")
            descriptor_bits.append(f"[location={r.location}]")
            if r.notes:
                descriptor_bits.append(f"[note={r.notes}]")
            f.write(f">{locus_prefix}|{r.locus_tag} {' '.join(descriptor_bits)}\n")
            for line_start in range(0, len(r.nt_sequence), _FASTA_LINE_WIDTH):
                f.write(r.nt_sequence[line_start : line_start + _FASTA_LINE_WIDTH] + "\n")
