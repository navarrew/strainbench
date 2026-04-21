"""Hierarchical clustering of strains and protein clusters from the DB.

This is the SQL-backed equivalent of strain-comp's `5_heatmap.py`. It
reads the presence/absence matrix from `v_cluster_cell` for a given
cluster_run, runs seaborn's clustermap to compute strain and cluster
dendrograms, and writes the resulting orderings back to the DB:

    strains.display_order   ← position along the strain axis
    clusters.display_order  ← position along the cluster axis

Downstream views (`export_xlsx`, the future web app) read those columns to
render strains and clusters in dendrogram order rather than insertion
order.

Optionally writes a heatmap PNG matching strain-comp's visual style so
biologists have something familiar to look at while we figure out the
GUI.

Multi-copy clusters (where a single strain has more than one CDS in the
cluster — typically transposases) are excluded from clustering by default
because copy-count noise dominates the dendrogram. They get a `display_order`
appended after the clustered rows so they're still visible at the bottom of
the spreadsheet.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path


@dataclass
class HeatmapResult:
    cluster_run_id: int
    method: str
    n_clusters_in_clustering: int
    n_clusters_excluded: int
    n_strains: int
    output_png: str | None


def compute_hierarchical_order(
    conn: sqlite3.Connection,
    *,
    cluster_run_id: int | None = None,
    method: str = "average",
    output_png: str | Path | None = None,
    strain_range: tuple[int, int] | None = None,
    color: str = "Blues",
    figsize: tuple[float, float] = (50.0, 50.0),
    write_orderings: bool = True,
) -> HeatmapResult:
    """Run hierarchical clustering on the strain × cluster presence matrix.

    Args:
        conn: Open SQLite connection (read-write if write_orderings=True).
        cluster_run_id: Which cluster_run to operate on. None → most recent
            active run.
        method: Linkage method passed to seaborn.clustermap (default 'average';
            other options: 'ward', 'complete', 'single', 'centroid').
        output_png: If set, write the heatmap figure here.
        strain_range: (min, max) — only include clusters whose strain_count is
            in this range. Default: include all single-copy clusters.
        color: matplotlib colormap name for the heatmap.
        figsize: matplotlib figure size in inches.
        write_orderings: If True, update `strains.display_order` and
            `clusters.display_order` with the dendrogram ordering. Set False
            for dry-run / inspection.
    """
    # Heavy imports localized so non-heatmap CLI commands don't pay the cost.
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    cluster_run_id = _resolve_cluster_run_id(conn, cluster_run_id)

    matrix, excluded_cluster_ids = _build_presence_matrix(conn, cluster_run_id, strain_range)

    if matrix.empty:
        raise ValueError(
            "No clusters left after filtering. Try widening --range or check "
            "that this cluster_run has any single-copy clusters."
        )

    # Run clustermap. yticklabels=False to avoid drawing 3000+ cluster labels.
    sns.set(font_scale=0.05)
    g = sns.clustermap(
        matrix,
        method=method,
        cmap=color,
        vmin=0, vmax=1,
        figsize=figsize,
        cbar=False,
        xticklabels=True,
        yticklabels=False,
    )
    g.ax_heatmap.set(xlabel="STRAINS", ylabel="CLUSTERS")

    if output_png:
        output_png = Path(output_png)
        output_png.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_png, dpi=200, bbox_inches="tight", pad_inches=0)
    plt.close("all")

    cluster_order = matrix.index[g.dendrogram_row.reordered_ind].tolist()
    strain_order = matrix.columns[g.dendrogram_col.reordered_ind].tolist()

    if write_orderings:
        _write_orderings(conn, cluster_run_id, cluster_order, strain_order, excluded_cluster_ids)

    return HeatmapResult(
        cluster_run_id=cluster_run_id,
        method=method,
        n_clusters_in_clustering=len(cluster_order),
        n_clusters_excluded=len(excluded_cluster_ids),
        n_strains=len(strain_order),
        output_png=str(output_png) if output_png else None,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────


def _resolve_cluster_run_id(conn: sqlite3.Connection, requested: int | None) -> int:
    if requested is not None:
        row = conn.execute(
            "SELECT cluster_run_id FROM cluster_runs WHERE cluster_run_id = ?",
            (requested,),
        ).fetchone()
        if row is None:
            raise ValueError(f"cluster_run_id={requested} not found")
        return requested
    row = conn.execute(
        "SELECT cluster_run_id FROM cluster_runs "
        "WHERE is_active = 1 ORDER BY cluster_run_id DESC LIMIT 1"
    ).fetchone()
    if row is None:
        raise ValueError("No active cluster_runs in DB. Run `strainbench cluster` first.")
    return int(row["cluster_run_id"])


def _build_presence_matrix(
    conn: sqlite3.Connection,
    cluster_run_id: int,
    strain_range: tuple[int, int] | None,
):
    """Return (matrix, excluded_cluster_ids).

    Matrix: rows=cluster_id, cols=strain_id, cells = 0/1 presence.
    Excluded cluster_ids: clusters dropped from clustering (multi-copy or
    out-of-range), to be appended to display_order at the end.
    """
    import pandas as pd

    cells = pd.read_sql_query(
        """
        SELECT vc.cluster_id, vc.strain_id, vc.hit_count
        FROM v_cluster_cell vc
        WHERE vc.cluster_run_id = ?
        """,
        conn,
        params=(cluster_run_id,),
    )
    matrix = (
        cells.pivot(index="cluster_id", columns="strain_id", values="hit_count")
        .fillna(0)
        .clip(upper=1)
        .astype(int)
    )

    # Filter: drop multi-copy clusters and apply strain_range if given.
    cluster_meta = pd.read_sql_query(
        """
        SELECT cluster_id, member_count, strain_count
        FROM clusters
        WHERE cluster_run_id = ?
        """,
        conn,
        params=(cluster_run_id,),
    )

    keep_mask = cluster_meta["member_count"] == cluster_meta["strain_count"]
    if strain_range is not None:
        lo, hi = strain_range
        keep_mask &= cluster_meta["strain_count"].between(lo, hi)
    keep_ids = set(cluster_meta.loc[keep_mask, "cluster_id"])
    excluded_ids = [cid for cid in cluster_meta["cluster_id"] if cid not in keep_ids]

    matrix = matrix.loc[matrix.index.isin(keep_ids)]
    return matrix, excluded_ids


def _write_orderings(
    conn: sqlite3.Connection,
    cluster_run_id: int,
    cluster_order: list[int],
    strain_order: list[int],
    excluded_cluster_ids: list[int],
) -> None:
    """Persist dendrogram orderings to display_order columns.

    Strategy:
        - Reset display_order for this run's clusters and for all strains.
        - Assign 1..N to clusters in dendrogram order.
        - Append excluded clusters after, in cluster_id order, continuing
          the numbering. They appear at the bottom of the spreadsheet.
        - Assign 1..M to strains in dendrogram order.
    """
    with conn:
        conn.execute(
            "UPDATE clusters SET display_order = NULL WHERE cluster_run_id = ?",
            (cluster_run_id,),
        )
        conn.execute("UPDATE strains SET display_order = NULL")

        conn.executemany(
            "UPDATE clusters SET display_order = ? WHERE cluster_id = ?",
            [(pos, cid) for pos, cid in enumerate(cluster_order, start=1)],
        )
        conn.executemany(
            "UPDATE strains SET display_order = ? WHERE strain_id = ?",
            [(pos, sid) for pos, sid in enumerate(strain_order, start=1)],
        )
        # Excluded clusters get appended at the end in stable cluster_id order.
        next_pos = len(cluster_order) + 1
        conn.executemany(
            "UPDATE clusters SET display_order = ? WHERE cluster_id = ?",
            [(next_pos + i, cid) for i, cid in enumerate(sorted(excluded_cluster_ids))],
        )
