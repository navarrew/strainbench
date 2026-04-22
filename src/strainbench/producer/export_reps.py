"""Export cluster representative sequences as a FASTA file.

Two modes:

  * **Default** (fast, DB-only). Reads the stored representative sequence
    for each cluster from the `clusters` table (the one mmseqs picked at
    clustering time). Output has one record per cluster; record IDs are
    cluster names. This is what you feed to external annotators
    (eggNOG-mapper, kofamscan, DeepNOG, etc.).

  * **`full_length_only=True`** (slower, joins the sidecar FASTAs). For each
    cluster, picks any un-truncated member via the `cds.notes` column instead
    of trusting the stored rep (which mmseqs picked without respect to
    truncation). Clusters whose members are ALL truncated are skipped.
    Produces the curated "unique full-length ORF catalog" for downstream
    phylogenetic or comparative analysis.

Record IDs are always the cluster name (e.g. 'INERS_NT_000001') so output
can be joined back to `clusters` without any extra lookup.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path

from strainbench.core.db import resolve_cluster_run_id as _resolve_cluster_run_id

_FASTA_LINE_WIDTH = 80


@dataclass
class ExportRepsSummary:
    output_path: str
    cluster_run_id: int
    n_reps_written: int
    n_clusters_skipped: int              # clusters with no sequence / all truncated
    sequence_type: str
    full_length_only: bool = False
    n_clusters_all_truncated: int = 0    # only populated when full_length_only=True


def export_cluster_representatives(
    db_path: str | Path,
    output_path: str | Path,
    *,
    cluster_run_id: int | None = None,
    sequence_type: str | None = None,
    full_length_only: bool = False,
    fasta_dir: str | Path | None = None,
) -> ExportRepsSummary:
    """Write cluster rep sequences to a FASTA file.

    Args:
        db_path: Path to a populated strainbench DB.
        output_path: Where to write the .faa or .fna file.
        cluster_run_id: Which cluster_run to use (None → most recent active).
        sequence_type: 'protein' or 'nucleotide'. None → match the cluster_run's
            type (protein clustering → aa reps, nucleotide clustering → nt reps).
        full_length_only: If True, bypass the stored rep and instead pick an
            un-truncated cluster member via `cds.notes`. Skips clusters whose
            members are all truncated. Requires `fasta_dir` so sequences can
            be read from the sidecar FASTAs.
        fasta_dir: Path to the sidecar FASTA parent dir (with fna/ and faa/
            subdirs). Required when `full_length_only=True`, ignored otherwise.
    """
    db_path = Path(db_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        # Pass sequence_type as a hint so --sequence-type nucleotide (without
        # an explicit --cluster-run-id) picks the nt run instead of the
        # protein run. Avoids the footgun where mismatched type + run silently
        # produces an empty FASTA.
        cluster_run_id = _resolve_cluster_run_id(
            conn, cluster_run_id, sequence_type_hint=sequence_type,
        )
        run_row = conn.execute(
            "SELECT sequence_type FROM cluster_runs WHERE cluster_run_id = ?",
            (cluster_run_id,),
        ).fetchone()
        if run_row is None:
            raise ValueError(f"cluster_run_id={cluster_run_id} not found")
        sequence_type = sequence_type or run_row["sequence_type"]

        if full_length_only:
            summary = _write_full_length_only(
                conn, cluster_run_id, sequence_type, fasta_dir, output_path,
            )
        else:
            summary = _write_stored_reps(
                conn, cluster_run_id, sequence_type, output_path,
            )
    finally:
        conn.close()
    return summary


# ─────────────────────────────────────────────────────────────────────────────
# Default path: dump stored rep sequences straight from the DB
# ─────────────────────────────────────────────────────────────────────────────


def _write_stored_reps(
    conn: sqlite3.Connection,
    cluster_run_id: int,
    sequence_type: str,
    output_path: Path,
) -> ExportRepsSummary:
    seq_column = (
        "representative_aa_seq" if sequence_type == "protein"
        else "representative_nt_seq"
    )
    rows = conn.execute(
        f"""
        SELECT cluster_name, consensus_annotation, {seq_column} AS seq
        FROM clusters
        WHERE cluster_run_id = ?
        ORDER BY cluster_name
        """,
        (cluster_run_id,),
    ).fetchall()

    n_written = n_skipped = 0
    with output_path.open("w") as f:
        for row in rows:
            seq = row["seq"]
            if not seq:
                n_skipped += 1
                continue
            annotation = row["consensus_annotation"] or ""
            f.write(f">{row['cluster_name']} {annotation}\n")
            for line_start in range(0, len(seq), _FASTA_LINE_WIDTH):
                f.write(seq[line_start:line_start + _FASTA_LINE_WIDTH] + "\n")
            n_written += 1

    return ExportRepsSummary(
        output_path=str(output_path),
        cluster_run_id=cluster_run_id,
        n_reps_written=n_written,
        n_clusters_skipped=n_skipped,
        sequence_type=sequence_type,
        full_length_only=False,
    )


# ─────────────────────────────────────────────────────────────────────────────
# full_length_only path: pick un-truncated members; read from sidecar FASTAs
# ─────────────────────────────────────────────────────────────────────────────


def _write_full_length_only(
    conn: sqlite3.Connection,
    cluster_run_id: int,
    sequence_type: str,
    fasta_dir: str | Path | None,
    output_path: Path,
) -> ExportRepsSummary:
    if fasta_dir is None:
        raise ValueError(
            "fasta_dir is required when full_length_only=True (sequences are "
            "read from the sidecar FASTAs rather than from the DB)."
        )
    fasta_dir = Path(fasta_dir)
    subdir = fasta_dir / ("faa" if sequence_type == "protein" else "fna")
    if not subdir.is_dir():
        raise ValueError(
            f"sidecar dir not found: {subdir} "
            f"(expected <fasta_dir>/{subdir.name}/ with per-strain .{subdir.name} files)"
        )

    # Totals for reporting — how many clusters exist at all, and how many
    # are all-truncated (= skipped in full-length-only mode).
    total_clusters = int(conn.execute(
        "SELECT COUNT(*) AS n FROM clusters WHERE cluster_run_id = ?",
        (cluster_run_id,),
    ).fetchone()["n"])

    # One pick per cluster: the MIN cds_id among its un-truncated members.
    # MIN is arbitrary-but-stable — re-runs produce the same picks.
    picks = conn.execute(
        """
        WITH untruncated AS (
            SELECT m.cluster_id, MIN(cds.cds_id) AS pick_cds_id
            FROM cluster_membership m
            JOIN cds ON cds.cds_id = m.cds_id
            WHERE m.cluster_run_id = ?
              AND (cds.notes IS NULL OR cds.notes NOT LIKE '%truncated%')
            GROUP BY m.cluster_id
        )
        SELECT
            nt.cluster_name,
            nt.consensus_annotation,
            s.locus_prefix,
            cds.locus_tag
        FROM untruncated u
        JOIN cds       ON cds.cds_id = u.pick_cds_id
        JOIN strains s ON s.strain_id = cds.strain_id
        JOIN clusters nt ON nt.cluster_id = u.cluster_id
        ORDER BY nt.cluster_name
        """,
        (cluster_run_id,),
    ).fetchall()

    n_all_truncated = total_clusters - len(picks)

    # Build a BioPython byte-offset index over all strain sidecar files
    # (one-time cost; sub-millisecond lookup after that).
    try:
        from Bio import SeqIO
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "BioPython is required for full_length_only mode; install via "
            "`pip install '.[producer]'`"
        ) from exc

    sidecar_files = sorted(subdir.glob(f"*.{subdir.name}"))
    if not sidecar_files:
        raise ValueError(f"no *.{subdir.name} files found under {subdir}")

    # Index lives in a tempdir for the duration of this call. Cheap to rebuild.
    # Use a tempdir (not NamedTemporaryFile, which creates an empty file that
    # SeqIO.index_db would then try to read as a valid index and fail).
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix="strainbench_fasta_idx_")
    idx_path = str(Path(tmpdir) / "idx.sqlite")
    try:
        index = SeqIO.index_db(
            idx_path, [str(p) for p in sidecar_files], "fasta",
        )

        n_written = n_missing = 0
        with output_path.open("w") as out:
            for row in picks:
                # strainbench's ingest writes record IDs as "<locus_prefix>|<locus_tag>"
                rec_id = f"{row['locus_prefix']}|{row['locus_tag']}"
                if rec_id not in index:
                    n_missing += 1
                    continue
                seq = str(index[rec_id].seq)
                ann = row['consensus_annotation'] or ''
                out.write(f">{row['cluster_name']} {rec_id} {ann}\n".rstrip() + "\n")
                for i in range(0, len(seq), _FASTA_LINE_WIDTH):
                    out.write(seq[i:i + _FASTA_LINE_WIDTH] + "\n")
                n_written += 1
        index.close()
    finally:
        # BioPython's index_db leaves a SQLite file behind; clean up the whole tempdir.
        import shutil as _shutil
        _shutil.rmtree(tmpdir, ignore_errors=True)

    return ExportRepsSummary(
        output_path=str(output_path),
        cluster_run_id=cluster_run_id,
        n_reps_written=n_written,
        n_clusters_skipped=n_missing,
        sequence_type=sequence_type,
        full_length_only=True,
        n_clusters_all_truncated=n_all_truncated,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Internals
# ─────────────────────────────────────────────────────────────────────────────


# resolve_cluster_run_id moved to core.db — imported at the top of this file.
