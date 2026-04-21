"""Cluster CDS sequences with mmseqs2 and load the results into the DB.

The clustering step takes a populated strainbench DB plus its sidecar FASTA
files and produces three coordinated sets of rows:

    cluster_runs:        1 row capturing the parameters (pct_id, coverage,
                         sequence_type, prefix, optional parent for nesting)
    clusters:            N rows, one per output cluster, ordered by member
                         count (largest cluster gets <prefix>_000001)
    cluster_membership:  M rows total, one per CDS-in-a-cluster (M ≈ total CDS)

Per-cluster aggregates (member_count, strain_count, avg_gc_pct, gc_spread,
avg_aa_length) are derived after membership inserts via a single SQL pass.

Representative sequences for each cluster are stored inline on the `clusters`
row (representative_nt_seq, representative_aa_seq) so the web app can render
cluster detail pages without touching the sidecar FASTA files.

mmseqs2 itself is invoked as a subprocess. It must be on PATH; install with
`conda install -c bioconda mmseqs2` or `brew install mmseqs2`.
"""

from __future__ import annotations

import shutil
import sqlite3
import subprocess
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

_DEFAULT_PCT_ID = 80
_DEFAULT_COVERAGE = 90
_DEFAULT_NAME_PREFIX = "CLUSTER"


class MMseqsError(RuntimeError):
    """Raised when an mmseqs2 subprocess returns non-zero."""


class ClusterError(RuntimeError):
    """Raised for clustering problems unrelated to the mmseqs subprocess
    (e.g., empty input, unknown FASTA record IDs in the cluster output)."""


@dataclass
class ClusterRunResult:
    cluster_run_id: int
    label: str
    cluster_count: int
    member_count: int
    sequence_type: str


# ─────────────────────────────────────────────────────────────────────────────
# Public entry points
# ─────────────────────────────────────────────────────────────────────────────


def cluster_strains(
    conn: sqlite3.Connection,
    fasta_dir: str | Path,
    *,
    sequence_type: str = "protein",
    pct_identity: int = _DEFAULT_PCT_ID,
    coverage: int = _DEFAULT_COVERAGE,
    name_prefix: str = _DEFAULT_NAME_PREFIX,
    label: str | None = None,
    parent_run_id: int | None = None,
    parent_cluster_id: int | None = None,
    work_dir: str | Path | None = None,
    keep_work_dir: bool = False,
) -> ClusterRunResult:
    """Cluster all CDSs currently in the DB and store the results.

    Args:
        conn: Open SQLite connection (read-write).
        fasta_dir: Parent of `fna/` and `faa/` sidecar dirs (the ones written
            by ingestion). Clustering reads from `<fasta_dir>/faa/*.faa` for
            sequence_type='protein', and `<fasta_dir>/fna/*.fna` for 'nucleotide'.
        sequence_type: 'protein' or 'nucleotide'. mmseqs2 auto-detects the
            input type, but we pass --search-type 3 explicitly for nucleotide.
        pct_identity: Minimum sequence identity for cluster membership
            (mmseqs --min-seq-id). 0–99.
        coverage: Minimum mutual coverage of aligned region (mmseqs -c).
        name_prefix: Cluster names will be '<prefix>_000001', etc., zero-padded.
        label: Human-readable label for this run. Must be UNIQUE across runs in
            this DB. Defaults to '<sequence_type>-<pct>%-<coverage>%-<timestamp>'.
        parent_run_id, parent_cluster_id: For nested clustering. Currently
            recorded but the input filtering is not yet implemented; will be
            added when the GUI needs sub-clustering.
        work_dir: Where mmseqs writes its scratch files. Default: a tempdir
            that's deleted on success.
        keep_work_dir: If True, do not delete `work_dir` on success
            (useful for debugging).
    """
    if sequence_type not in {"protein", "nucleotide"}:
        raise ValueError(f"sequence_type must be 'protein' or 'nucleotide', got {sequence_type!r}")

    if shutil.which("mmseqs") is None:
        raise MMseqsError(
            "mmseqs2 not found on PATH. Install via "
            "`conda install -c bioconda mmseqs2` or `brew install mmseqs2`."
        )

    fasta_dir = Path(fasta_dir)
    sidecar_glob = "faa/*.faa" if sequence_type == "protein" else "fna/*.fna"
    sidecar_files = sorted(fasta_dir.glob(sidecar_glob))
    if not sidecar_files:
        raise ClusterError(
            f"No sidecar FASTA files found at {fasta_dir}/{sidecar_glob}. "
            "Run `strainbench ingest` first."
        )

    label = label or _default_label(sequence_type, pct_identity, coverage)

    # work_dir setup — keep around if caller asked, else auto-clean.
    if work_dir is None:
        work_dir_obj = tempfile.TemporaryDirectory(prefix="strainbench_mmseqs_")
        work_path = Path(work_dir_obj.name)
        cleanup = work_dir_obj.cleanup
    else:
        work_path = Path(work_dir)
        work_path.mkdir(parents=True, exist_ok=True)
        cleanup = (lambda: None) if keep_work_dir else (lambda: shutil.rmtree(work_path, ignore_errors=True))

    try:
        # 1. Concatenate sidecars into one input file (mmseqs reads one fasta).
        concat_path = work_path / "all_input.fasta"
        _concat_files(sidecar_files, concat_path)

        # 2. Run mmseqs.
        clustered_path = _run_mmseqs(
            input_fasta=concat_path,
            work_dir=work_path,
            pct_identity=pct_identity,
            coverage=coverage,
            sequence_type=sequence_type,
        )

        # 3. Parse the result2flat output into clusters.
        clusters = parse_clustered_fasta(clustered_path)
        if not clusters:
            raise ClusterError(f"mmseqs produced no clusters; check {clustered_path}")

        # 4. Load into DB (one transaction).
        return _load_clusters_into_db(
            conn,
            clusters=clusters,
            sequence_type=sequence_type,
            pct_identity=pct_identity,
            coverage=coverage,
            name_prefix=name_prefix,
            label=label,
            parent_run_id=parent_run_id,
            parent_cluster_id=parent_cluster_id,
        )
    finally:
        cleanup()


# ─────────────────────────────────────────────────────────────────────────────
# Output parser (public — exercised directly by tests)
# ─────────────────────────────────────────────────────────────────────────────


def parse_clustered_fasta(path: str | Path) -> list[list[tuple[str, str]]]:
    """Parse mmseqs2 result2flat output into a list of clusters.

    Each cluster is a list of (header, sequence) tuples. The first tuple is
    the cluster representative.

    The mmseqs2 result2flat format marks cluster boundaries with two
    consecutive '>' lines (the representative's header appearing twice).
    Sequences may span multiple lines; this parser handles both wrapped
    and unwrapped output.
    """
    clusters: list[list[tuple[str, str]]] = []
    current_cluster: list[tuple[str, str]] = []
    pending_header: str | None = None
    pending_seq_lines: list[str] = []
    last_was_header = False

    def flush_record() -> None:
        nonlocal pending_header, pending_seq_lines
        if pending_header is not None and pending_seq_lines:
            current_cluster.append((pending_header, "".join(pending_seq_lines)))
        pending_header = None
        pending_seq_lines = []

    with open(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if last_was_header:
                    # Two '>' in a row → cluster boundary. The first '>' was a
                    # marker (no sequence followed); the current header overrides it.
                    if current_cluster:
                        clusters.append(current_cluster)
                        current_cluster = []
                    pending_header = line[1:]
                    pending_seq_lines = []
                    last_was_header = False
                else:
                    flush_record()
                    pending_header = line[1:]
                    pending_seq_lines = []
                    last_was_header = True
            else:
                pending_seq_lines.append(line)
                last_was_header = False

    flush_record()
    if current_cluster:
        clusters.append(current_cluster)

    # Sort largest-first so cluster names reflect abundance ranking.
    clusters.sort(key=len, reverse=True)
    return clusters


# ─────────────────────────────────────────────────────────────────────────────
# Internals
# ─────────────────────────────────────────────────────────────────────────────


def _default_label(sequence_type: str, pct_identity: int, coverage: int) -> str:
    from datetime import datetime

    return (
        f"{sequence_type}-{pct_identity}id-{coverage}cov-"
        f"{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    )


def _concat_files(files: Iterable[Path], output: Path) -> None:
    with output.open("wb") as out:
        for path in files:
            with path.open("rb") as src:
                shutil.copyfileobj(src, out)


def _run_mmseqs(
    *,
    input_fasta: Path,
    work_dir: Path,
    pct_identity: int,
    coverage: int,
    sequence_type: str,
) -> Path:
    """Run mmseqs createdb / cluster / createseqfiledb / result2flat.

    Returns the path to the final result2flat output (clustered_sequences.fasta).
    """
    db = work_dir / "DB"
    clu_db = work_dir / "clusteredDB"
    clu_seq = work_dir / "clu_seq"
    tmp = work_dir / "tmp"
    tmp.mkdir(exist_ok=True)
    out_fasta = work_dir / "clustered_sequences.fasta"

    min_seq_id = pct_identity / 100.0
    cov_frac = coverage / 100.0

    # mmseqs defaults to using all available cores; no explicit --threads needed.
    # mmseqs 16+: nucleotide vs protein is set via --dbtype on `createdb` (1=aa, 2=nt).
    # The `cluster` step then auto-detects from the DB type — no extra flag needed.
    createdb_cmd = ["mmseqs", "createdb", str(input_fasta), str(db)]
    if sequence_type == "nucleotide":
        createdb_cmd += ["--dbtype", "2"]
    else:
        createdb_cmd += ["--dbtype", "1"]

    cluster_cmd = [
        "mmseqs", "cluster", str(db), str(clu_db), str(tmp),
        "--min-seq-id", str(min_seq_id),
        "-c", str(cov_frac),
        "--cov-mode", "0",
    ]

    pipeline = [
        createdb_cmd,
        cluster_cmd,
        ["mmseqs", "createseqfiledb", str(db), str(clu_db), str(clu_seq)],
        ["mmseqs", "result2flat", str(db), str(db), str(clu_seq), str(out_fasta)],
    ]

    for cmd in pipeline:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise MMseqsError(
                f"mmseqs command failed (exit {result.returncode}):\n"
                f"  cmd: {' '.join(cmd)}\n"
                f"  stderr: {result.stderr.strip()[:1000]}"
            )

    if not out_fasta.exists():
        raise MMseqsError(f"mmseqs finished but produced no output at {out_fasta}")
    return out_fasta


def _load_clusters_into_db(
    conn: sqlite3.Connection,
    *,
    clusters: list[list[tuple[str, str]]],
    sequence_type: str,
    pct_identity: int,
    coverage: int,
    name_prefix: str,
    label: str,
    parent_run_id: int | None,
    parent_cluster_id: int | None,
) -> ClusterRunResult:
    """Insert one cluster_run + its clusters + memberships into the DB.

    All inserts happen in a single transaction. Aggregates (member_count,
    strain_count, avg_gc_pct, gc_spread, avg_aa_length) are computed in SQL
    after membership rows land.
    """
    # Build (locus_prefix, locus_tag) → (cds_id, aa_seq_lookup_key) map.
    # We can't store aa_seq in the dict (we don't have it in DB), but cds_id
    # is what we need for membership inserts.
    lookup: dict[tuple[str, str], int] = {}
    for row in conn.execute(
        "SELECT c.cds_id, c.locus_tag, s.locus_prefix "
        "FROM cds c JOIN strains s ON s.strain_id = c.strain_id"
    ):
        lookup[(row["locus_prefix"], row["locus_tag"])] = row["cds_id"]

    # Map cluster output member IDs → cds_ids. Capture rep separately.
    # Drop clusters whose members we can't resolve (shouldn't happen, but defensive).
    resolved: list[tuple[int, list[int], str, str]] = []
    # tuples: (rep_cds_id, [member_cds_ids], rep_aa_or_nt_seq_from_output, rep_header)
    skipped_unresolved = 0
    for cluster in clusters:
        member_cds_ids: list[int] = []
        rep_cds_id: int | None = None
        rep_seq = ""
        rep_header = cluster[0][0] if cluster else ""
        for header, seq in cluster:
            cds_id = _resolve_member_id(header, lookup)
            if cds_id is None:
                continue
            if rep_cds_id is None:
                rep_cds_id = cds_id
                rep_seq = seq
            member_cds_ids.append(cds_id)
        if rep_cds_id is None or not member_cds_ids:
            skipped_unresolved += 1
            continue
        resolved.append((rep_cds_id, member_cds_ids, rep_seq, rep_header))

    if not resolved:
        raise ClusterError(
            f"None of the {len(clusters)} mmseqs clusters could be resolved to "
            "DB rows. Check that ingest was run on the same fasta_dir."
        )

    with conn:
        # 1. cluster_runs row
        cur = conn.execute(
            """
            INSERT INTO cluster_runs (
                label, sequence_type, pct_identity, coverage, name_prefix,
                parent_run_id, parent_cluster_id, is_active
            ) VALUES (?, ?, ?, ?, ?, ?, ?, 1)
            RETURNING cluster_run_id
            """,
            (label, sequence_type, pct_identity, coverage, name_prefix,
             parent_run_id, parent_cluster_id),
        )
        cluster_run_id = int(cur.fetchone()["cluster_run_id"])

        # 2. clusters rows (one per cluster, ordered largest-first → 000001 etc.)
        cluster_rows: list[tuple] = []
        for i, (rep_cds_id, _members, rep_seq, _rep_header) in enumerate(resolved, start=1):
            cluster_name = f"{name_prefix}_{i:06d}"
            nt_seq = rep_seq if sequence_type == "nucleotide" else None
            aa_seq = rep_seq if sequence_type == "protein" else None
            cluster_rows.append((
                cluster_run_id, cluster_name, rep_cds_id, nt_seq, aa_seq,
            ))
        # Use executemany + RETURNING isn't supported uniformly, so do single inserts
        # and capture cluster_ids.
        cluster_ids: list[int] = []
        for row in cluster_rows:
            cur = conn.execute(
                """
                INSERT INTO clusters (
                    cluster_run_id, cluster_name, representative_cds_id,
                    representative_nt_seq, representative_aa_seq
                ) VALUES (?, ?, ?, ?, ?)
                RETURNING cluster_id
                """,
                row,
            )
            cluster_ids.append(int(cur.fetchone()["cluster_id"]))

        # 3. cluster_membership rows in bulk.
        membership_rows = [
            (cluster_run_id, cluster_id, cds_id)
            for cluster_id, (_rep, members, _seq, _header) in zip(cluster_ids, resolved, strict=True)
            for cds_id in members
        ]
        conn.executemany(
            "INSERT INTO cluster_membership (cluster_run_id, cluster_id, cds_id) "
            "VALUES (?, ?, ?)",
            membership_rows,
        )

        # 4. Aggregate updates. One pass, one query, scoped to this run.
        conn.execute(
            """
            UPDATE clusters
            SET member_count  = agg.member_count,
                strain_count  = agg.strain_count,
                avg_gc_pct    = agg.avg_gc_pct,
                gc_spread     = agg.gc_spread,
                avg_aa_length = agg.avg_aa_length,
                consensus_annotation = agg.modal_annotation
            FROM (
                SELECT
                    m.cluster_id,
                    COUNT(*)                      AS member_count,
                    COUNT(DISTINCT c.strain_id)   AS strain_count,
                    AVG(c.gc_pct)                 AS avg_gc_pct,
                    MAX(c.gc_pct) - MIN(c.gc_pct) AS gc_spread,
                    CAST(AVG(c.aa_length) AS INTEGER) AS avg_aa_length,
                    -- Modal annotation: most common annotation among members.
                    (SELECT c2.annotation
                     FROM cluster_membership m2
                     JOIN cds c2 ON c2.cds_id = m2.cds_id
                     WHERE m2.cluster_id = m.cluster_id
                     GROUP BY c2.annotation
                     ORDER BY COUNT(*) DESC, c2.annotation
                     LIMIT 1) AS modal_annotation
                FROM cluster_membership m
                JOIN cds c ON c.cds_id = m.cds_id
                WHERE m.cluster_run_id = ?
                GROUP BY m.cluster_id
            ) AS agg
            WHERE clusters.cluster_id = agg.cluster_id
            """,
            (cluster_run_id,),
        )

    return ClusterRunResult(
        cluster_run_id=cluster_run_id,
        label=label,
        cluster_count=len(resolved),
        member_count=sum(len(members) for _, members, _, _ in resolved),
        sequence_type=sequence_type,
    )


def _resolve_member_id(
    header: str,
    lookup: dict[tuple[str, str], int],
) -> int | None:
    """Map an mmseqs member header back to a cds_id via the (prefix, locus_tag) key.

    The header format we wrote in ingest.py is:
        LOCUSPREFIX|LOCUS_TAG [locus_prefix=...] [protein_id=...] ...
    so the first whitespace-delimited token, split on '|', gives us the key.
    """
    first_token = header.split(None, 1)[0]
    if "|" not in first_token:
        return None
    locus_prefix, locus_tag = first_token.split("|", 1)
    return lookup.get((locus_prefix, locus_tag))
