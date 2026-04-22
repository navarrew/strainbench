"""Human-readable snapshot of a strainbench database.

Answers "what have I done to this dataset, and what's missing?" — the
command you run as the first thing when returning to a project after a
break. Works on any strainbench DB regardless of schema version (gracefully
degrades when newer columns are absent).
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path


# Annotation sources that try to label every protein. Low coverage on these
# generally means the tool didn't finish or was run on the wrong input.
_GENERAL_PURPOSE_SOURCES = frozenset({
    "COG", "KEGG", "GO", "Pfam", "TIGRFAM", "PANTHER", "SUPERFAMILY", "SMART",
})

# Other sources (AMRFinder, CAZy, DefenseFinder, PADLOC, EC, …) are
# specialized — they label a deliberately narrow class of proteins (AMR
# genes, carb-active enzymes, defense systems, enzymes). Low coverage on
# these is expected biology, not a problem.


@dataclass
class ClusterRunStatus:
    cluster_run_id: int
    label: str
    sequence_type: str
    pct_identity: int
    coverage: int
    is_active: bool
    cluster_count: int
    has_display_order: bool
    annotation_summary: list[dict]   # [{source, n_clusters, coverage_pct, latest_file, latest_at}]


@dataclass
class DatasetStatus:
    db_path: str
    schema_version: int
    strains_count: int
    cds_count: int
    strains_ingested_at: str | None
    cluster_runs: list[ClusterRunStatus]


def collect_status(db_path: str | Path) -> DatasetStatus:
    """Gather counts and annotation coverage from a strainbench DB."""
    db_path = Path(db_path)
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        schema_version = conn.execute(
            "SELECT MAX(version) AS v FROM schema_version"
        ).fetchone()["v"] or 0

        strains_count = conn.execute("SELECT COUNT(*) AS n FROM strains").fetchone()["n"]
        cds_count = conn.execute("SELECT COUNT(*) AS n FROM cds").fetchone()["n"]
        ingested_row = conn.execute(
            "SELECT MIN(added_at) AS first_added FROM strains"
        ).fetchone()
        ingested_at = ingested_row["first_added"]

        runs: list[ClusterRunStatus] = []
        for r in conn.execute(
            "SELECT cluster_run_id, label, sequence_type, pct_identity, coverage, is_active "
            "FROM cluster_runs ORDER BY cluster_run_id"
        ).fetchall():
            cluster_ids = [
                row["cluster_id"]
                for row in conn.execute(
                    "SELECT cluster_id FROM clusters WHERE cluster_run_id = ?",
                    (r["cluster_run_id"],),
                )
            ]
            cluster_count = len(cluster_ids)
            has_order = False
            if cluster_count > 0:
                has_order = conn.execute(
                    "SELECT COUNT(*) AS n FROM clusters "
                    "WHERE cluster_run_id = ? AND display_order IS NOT NULL",
                    (r["cluster_run_id"],),
                ).fetchone()["n"] > 0

            annotation_summary = _annotation_summary(conn, r["cluster_run_id"], cluster_count)
            runs.append(ClusterRunStatus(
                cluster_run_id=int(r["cluster_run_id"]),
                label=r["label"],
                sequence_type=r["sequence_type"],
                pct_identity=int(r["pct_identity"]),
                coverage=int(r["coverage"]),
                is_active=bool(r["is_active"]),
                cluster_count=cluster_count,
                has_display_order=has_order,
                annotation_summary=annotation_summary,
            ))
    finally:
        conn.close()

    return DatasetStatus(
        db_path=str(db_path),
        schema_version=int(schema_version),
        strains_count=strains_count,
        cds_count=cds_count,
        strains_ingested_at=ingested_at,
        cluster_runs=runs,
    )


def _annotation_summary(
    conn: sqlite3.Connection, cluster_run_id: int, cluster_count: int
) -> list[dict]:
    """Return one row per annotation source for this cluster_run."""
    # Older schemas may not have source_file / tool_format — coalesce handles that.
    try:
        cur = conn.execute(
            """
            SELECT
                a.source                       AS source,
                COUNT(DISTINCT a.cluster_id)   AS n_clusters,
                COUNT(*)                       AS n_rows,
                MAX(a.applied_at)              AS latest_at,
                COALESCE(
                    (SELECT a2.source_file FROM cluster_annotations a2
                     WHERE a2.source = a.source AND a2.cluster_id IN (
                        SELECT cluster_id FROM clusters WHERE cluster_run_id = ?
                     )
                     ORDER BY a2.applied_at DESC LIMIT 1),
                    NULL
                )                              AS latest_file
            FROM cluster_annotations a
            JOIN clusters c ON c.cluster_id = a.cluster_id
            WHERE c.cluster_run_id = ?
            GROUP BY a.source
            ORDER BY a.source
            """,
            (cluster_run_id, cluster_run_id),
        )
        rows = cur.fetchall()
    except sqlite3.OperationalError:
        # cluster_annotations table doesn't exist yet on this DB.
        return []

    return [
        {
            "source": row["source"],
            "n_clusters": int(row["n_clusters"]),
            "n_rows": int(row["n_rows"]),
            "coverage_pct": (row["n_clusters"] * 100.0 / cluster_count) if cluster_count else 0.0,
            "latest_at": row["latest_at"],
            "latest_file": row["latest_file"],
        }
        for row in rows
    ]


def format_status(status: DatasetStatus) -> str:
    """Render a DatasetStatus as a human-readable multi-line string."""
    lines = []
    lines.append(f"Database: {status.db_path}")
    lines.append(f"Schema version: v{status.schema_version}")
    lines.append("")
    lines.append(f"Strains:  {status.strains_count:>10,}"
                 + (f"   (first added: {status.strains_ingested_at})" if status.strains_ingested_at else ""))
    lines.append(f"CDSs:     {status.cds_count:>10,}")
    lines.append("")

    if not status.cluster_runs:
        lines.append("No cluster_runs yet. Run `strainbench cluster` next.")
        return "\n".join(lines)

    for run in status.cluster_runs:
        marker = "●" if run.is_active else "○"
        lines.append(
            f"{marker} cluster_run_id={run.cluster_run_id}  '{run.label}'"
        )
        lines.append(
            f"    {run.sequence_type}, id≥{run.pct_identity}%, cov≥{run.coverage}%  "
            f"— {run.cluster_count:,} clusters"
        )
        if run.has_display_order:
            lines.append("    display_order populated (heatmap step ran) ✓")
        else:
            lines.append("    display_order empty (run `strainbench heatmap` to order strains/clusters)")

        if not run.annotation_summary:
            lines.append("    annotations: none")
            lines.append("      → run `strainbench annotate-prep` to start")
        else:
            lines.append("    annotations:")
            for ann in run.annotation_summary:
                file_note = f"  ({ann['latest_file']})" if ann["latest_file"] else ""
                lines.append(
                    f"      {ann['source']:10s} "
                    f"{ann['n_clusters']:>6,} clusters ({ann['coverage_pct']:>5.1f}%)  "
                    f"{ann['n_rows']:>6,} rows   latest {ann['latest_at']}{file_note}"
                )
            # Only warn for general-purpose annotators (COG, KEGG, GO, Pfam, …).
            # Specialized annotators (AMRFinder, CAZy, DefenseFinder, …) legitimately
            # hit a narrow slice of proteins — low coverage there is biology, not a bug.
            for ann in run.annotation_summary:
                if (
                    ann["source"] in _GENERAL_PURPOSE_SOURCES
                    and ann["coverage_pct"] < 10.0
                    and run.cluster_count > 100
                ):
                    lines.append(
                        f"      ❗ {ann['source']} coverage is only {ann['coverage_pct']:.1f}% — "
                        f"expected >60% from a full annotation run."
                    )
        lines.append("")

    return "\n".join(lines)
