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
from dataclasses import dataclass, field
from pathlib import Path

from strainbench.core.db import resolve_cluster_run_id as _resolve_cluster_run_id


@dataclass
class HeatmapResult:
    cluster_run_id: int
    method: str
    n_clusters_in_clustering: int
    n_clusters_excluded: int
    n_strains: int
    output_png: str | None
    # Dendrogram orderings — populated whether or not we persist them to the
    # DB. Consumers (e.g. export-strainlist with an explicit --cluster-run-id)
    # can call compute_hierarchical_order with write_orderings=False and read
    # strain_order straight from the result.
    strain_order: list[int] = field(default_factory=list)
    cluster_order: list[int] = field(default_factory=list)


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
    write_strain_order: bool = True,
    strain_label_stride: int | None = None,
    cluster_label_stride: int | None = None,
    label_density_threshold: int = 1250,
) -> HeatmapResult:
    """Run hierarchical clustering on the strain × cluster presence matrix.

    Args:
        conn: Open SQLite connection (read-write if write_orderings=True).
        cluster_run_id: Which cluster_run to operate on. None → most recent
            active protein run (see `_resolve_cluster_run_id`).
        method: Linkage method passed to seaborn.clustermap (default 'average';
            other options: 'ward', 'complete', 'single', 'centroid').
        output_png: If set, write the heatmap figure here.
        strain_range: (min, max) — only include clusters whose strain_count is
            in this range. Default: include all single-copy clusters.
        color: matplotlib colormap name for the heatmap.
        figsize: matplotlib figure size in inches.
        write_orderings: If True, update `clusters.display_order` with the
            dendrogram ordering. Set False for dry-run / inspection.
        write_strain_order: If True (and write_orderings is also True), also
            update `strains.display_order`. Set False when running multiple
            heatmaps (e.g. protein + nucleotide) — only one of them should
            drive the strain axis of the spreadsheet, by convention the
            protein run (canonical pan-genome view).
        strain_label_stride: Show every Nth strain label along the X-axis.
            None → auto (stride chosen so visible labels ≤ threshold).
            0 → no labels. 1 → every label. N → every Nth (after dendrogram
            reordering, so labels are evenly spaced visually).
        cluster_label_stride: Same idea for Y-axis cluster labels. Cluster
            labels show 'CLUSTER_NAME (gene_name)' when a /gene qualifier
            exists in any member, else just 'CLUSTER_NAME'. Auto-stride is
            usually what you want — Gardnerella's 12k clusters get stride 7,
            iners's 3.5k gets stride 2, small datasets get every label.
        label_density_threshold: Auto-stride target. The auto formula picks
            the smallest stride such that visible_labels ≤ threshold. 2000
            is a rough sweet spot on a 50"-tall figure rendered at 200 dpi.
    """
    # Heavy imports localized so non-heatmap CLI commands don't pay the cost.
    import sys

    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    # scipy's hierarchy.dendrogram walks the cluster tree with pure-Python
    # recursion. Python's default limit (1000) blows up on any matrix with
    # more than ~500 items. Bumping the limit here matches what strain-comp's
    # 5_heatmap.py did. For a performance upgrade (avoids the recursion
    # fallback AND is ~10× faster), install `fastcluster` in this env:
    #     pip install fastcluster
    # seaborn auto-detects it.
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 100_000))

    cluster_run_id = _resolve_cluster_run_id(conn, cluster_run_id)

    # Pull metadata for the title bar before running clustermap.
    run_meta = conn.execute(
        "SELECT label, sequence_type, pct_identity, coverage "
        "FROM cluster_runs WHERE cluster_run_id = ?",
        (cluster_run_id,),
    ).fetchone()

    matrix, excluded_cluster_ids = _build_presence_matrix(conn, cluster_run_id, strain_range)

    if matrix.empty:
        raise ValueError(
            "No clusters left after filtering. Try widening --range or check "
            "that this cluster_run has any single-copy clusters."
        )

    n_clusters, n_strains = matrix.shape

    # Auto-pick strides if caller didn't specify. Formula: smallest stride
    # such that the number of visible labels stays below the threshold.
    if strain_label_stride is None:
        strain_label_stride = _auto_stride(n_strains, label_density_threshold)
    if cluster_label_stride is None:
        cluster_label_stride = _auto_stride(n_clusters, label_density_threshold)

    show_strain_labels = strain_label_stride > 0
    show_cluster_labels = cluster_label_stride > 0

    # Capture the ID-keyed axes BEFORE rename — we need these to map
    # dendrogram positions back to strain_id / cluster_id for persistence
    # below (otherwise the post-rename matrix.columns is locus_prefix
    # strings and the UPDATE WHERE strain_id=… would silently match nothing).
    original_strain_ids = list(matrix.columns)
    original_cluster_ids = list(matrix.index)

    # Map matrix index/columns from raw IDs to human-readable labels so
    # seaborn's tick rendering shows strain names + cluster (+gene) names
    # (or blanks where the stride says to suppress).
    strain_label_map = _strain_labels(conn) if show_strain_labels else {}
    cluster_label_map = _cluster_labels(conn, cluster_run_id) if show_cluster_labels else {}
    matrix = matrix.rename(columns=strain_label_map, index=cluster_label_map)

    # Font scales calibrated to figure size and visible label density.
    # The goal: each VISIBLE label gets enough room to be readable when the
    # PNG is viewed at native resolution (200 dpi → ~10,000 px square).
    font_scale = _pick_font_scale(
        figsize,
        n_strains // max(strain_label_stride, 1) if show_strain_labels else 0,
        n_clusters // max(cluster_label_stride, 1) if show_cluster_labels else 0,
    )
    sns.set(font_scale=font_scale)

    g = sns.clustermap(
        matrix,
        method=method,
        cmap=color,
        vmin=0, vmax=1,
        figsize=figsize,
        cbar_pos=None,   # fully remove the colorbar axes — cbar=False alone
                         # can leave a ghost placeholder rectangle in the
                         # top-left corner in some seaborn versions
        xticklabels=show_strain_labels,
        yticklabels=show_cluster_labels,
    )
    g.ax_heatmap.set(xlabel="STRAINS", ylabel="CLUSTERS")

    # Rasterize ONLY the heatmap's cell-grid image — labels, title, axes,
    # and dendrograms stay as vectors. This is the difference between a
    # 100 MB PDF (every cell is its own vector rectangle, e.g. 4.6M cells
    # for Gardnerella) and a ~2 MB PDF with the grid as an embedded raster
    # and everything text-like infinitely zoomable. Only the Axes' images
    # (from pcolormesh/imshow) get rasterized, not its children.
    for child in g.ax_heatmap.get_children():
        # QuadMesh (from pcolormesh) is what seaborn uses for the grid.
        # Rasterize everything that's an image-like artist; leave text alone.
        if hasattr(child, "set_rasterized"):
            name = type(child).__name__
            if name in {"QuadMesh", "AxesImage", "PolyCollection"}:
                child.set_rasterized(True)

    # Apply stride: blank every label that isn't on a stride boundary AFTER
    # dendrogram reordering, so labels are evenly spaced in the displayed
    # image (not in the original matrix order).
    if show_strain_labels and strain_label_stride > 1:
        _apply_stride_to_axis(g.ax_heatmap, axis="x", stride=strain_label_stride)
    if show_cluster_labels and cluster_label_stride > 1:
        _apply_stride_to_axis(g.ax_heatmap, axis="y", stride=cluster_label_stride)

    # Top-of-image title — tells you which cluster_run this PNG is for when
    # you have multiple (protein + nucleotide) side by side.
    if run_meta is not None:
        title = (
            f"{run_meta['label']}   —   "
            f"{run_meta['sequence_type']} clustering "
            f"(id≥{run_meta['pct_identity']}%, cov≥{run_meta['coverage']}%)   —   "
            f"{n_clusters:,} clusters × {n_strains} strains   "
            f"[{method} linkage]"
        )
        # Note in the title when stride is in use, so viewers know labels
        # are subsampled (not random).
        stride_notes = []
        if show_cluster_labels and cluster_label_stride > 1:
            stride_notes.append(f"cluster labels: every {cluster_label_stride}th")
        if show_strain_labels and strain_label_stride > 1:
            stride_notes.append(f"strain labels: every {strain_label_stride}th")
        if not show_cluster_labels:
            stride_notes.append("cluster labels off")
        if stride_notes:
            title += "   (" + ", ".join(stride_notes) + ")"
        # Figure is 50" by default; title sits above the figure at ~1.5× the
        # axis label scale. Big enough to read when the PNG is rendered at
        # 200 dpi but not comically large when zoomed out.
        g.figure.suptitle(title, fontsize=60, y=0.995)

    if output_png:
        output_png = Path(output_png)
        output_png.parent.mkdir(parents=True, exist_ok=True)
        # Format is auto-detected from the path's extension by matplotlib.
        # For vector formats (.pdf, .svg, .eps), dpi is irrelevant — the text
        # and shapes are stored as vectors that scale infinitely. For raster
        # formats (.png, .jpg, .tiff), use 200 dpi for a 10k×10k image.
        # For vector formats (.pdf, .svg) we still pass dpi — it controls
        # the resolution of any RASTERIZED elements inside (the heatmap grid
        # we rasterized above). 200 dpi on a 50" figure = ~10,000 px, plenty
        # for zooming. Text and dendrograms remain vector regardless.
        plt.savefig(output_png, dpi=200, bbox_inches="tight", pad_inches=0.2)
    plt.close("all")

    # Map dendrogram positions back to the original IDs (not the renamed
    # display labels). reordered_ind is positional, so it indexes both views.
    cluster_order = [original_cluster_ids[i] for i in g.dendrogram_row.reordered_ind]
    strain_order = [original_strain_ids[i] for i in g.dendrogram_col.reordered_ind]

    if write_orderings:
        _write_orderings(
            conn, cluster_run_id, cluster_order, strain_order,
            excluded_cluster_ids,
            write_strain_order=write_strain_order,
        )

    return HeatmapResult(
        cluster_run_id=cluster_run_id,
        method=method,
        n_clusters_in_clustering=len(cluster_order),
        n_clusters_excluded=len(excluded_cluster_ids),
        n_strains=len(strain_order),
        output_png=str(output_png) if output_png else None,
        strain_order=strain_order,
        cluster_order=cluster_order,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────


def _strain_labels(conn: sqlite3.Connection) -> dict[int, str]:
    """Map strain_id → locus_prefix for X-axis tick labels.

    Locus prefix is the friendly strain ID (e.g. 'HMPREF0520', 'HXT39') —
    short enough to print on a heatmap column without overlap.
    """
    return {
        int(row["strain_id"]): row["locus_prefix"]
        for row in conn.execute("SELECT strain_id, locus_prefix FROM strains")
    }


def _cluster_labels(conn: sqlite3.Connection, cluster_run_id: int) -> dict[int, str]:
    """Map cluster_id → 'CLUSTER_NAME (gene_name)' or just 'CLUSTER_NAME'.

    The pooled gene-name aggregation matches the xlsx export's column —
    so a biologist seeing 'INERS_000005 (nrdF)' in the heatmap recognizes
    the same identifier they'd see in cluster_table.xlsx.
    """
    rows = conn.execute(
        """
        SELECT
            c.cluster_id,
            c.cluster_name,
            (SELECT GROUP_CONCAT(DISTINCT cds.gene_name)
             FROM cluster_membership m
             JOIN cds ON cds.cds_id = m.cds_id
             WHERE m.cluster_id = c.cluster_id
               AND cds.gene_name IS NOT NULL
               AND cds.gene_name != '')        AS gene_names
        FROM clusters c
        WHERE c.cluster_run_id = ?
        """,
        (cluster_run_id,),
    ).fetchall()
    out: dict[int, str] = {}
    for row in rows:
        cid = int(row["cluster_id"])
        name = row["cluster_name"]
        genes = row["gene_names"]
        if genes:
            # Take just the first gene name to keep the label short — full
            # pooled set is in the xlsx; here we just need a recognizable hint.
            first_gene = genes.split(",")[0].strip()
            out[cid] = f"{name} ({first_gene})"
        else:
            out[cid] = name
    return out


def _auto_stride(n_items: int, threshold: int) -> int:
    """Pick the smallest stride such that visible label count ≤ threshold.

    Stride 1 = every label. Stride 5 = every 5th. The point is to keep the
    LABELED rows/columns evenly spread without overcrowding.
    """
    if n_items <= 0 or threshold <= 0:
        return 1
    if n_items <= threshold:
        return 1
    # ceil(n_items / threshold)
    return -(-n_items // threshold)


def _apply_stride_to_axis(ax, axis: str, stride: int) -> None:
    """Blank every tick label that isn't a multiple of `stride` on `axis`.

    Called AFTER seaborn renders the clustermap, so the dendrogram ordering
    is already applied and the stride is applied to the DISPLAYED rows
    (not the matrix-order rows). That keeps labels evenly spaced in the
    final image.
    """
    if axis == "x":
        labels = [t.get_text() for t in ax.get_xticklabels()]
        kept = [lbl if i % stride == 0 else "" for i, lbl in enumerate(labels)]
        ax.set_xticklabels(kept)
    elif axis == "y":
        labels = [t.get_text() for t in ax.get_yticklabels()]
        kept = [lbl if i % stride == 0 else "" for i, lbl in enumerate(labels)]
        ax.set_yticklabels(kept)


def _pick_font_scale(
    figsize: tuple[float, float],
    n_visible_x_labels: int,
    n_visible_y_labels: int,
) -> float:
    """Compute a seaborn font_scale that sizes labels to fit their per-tick slot.

    seaborn's default base font is ~12pt. font_scale=1.0 → 12pt, 0.5 → 6pt.
    Each VISIBLE tick label needs roughly `figsize_inches / n_visible_labels`
    of axis space. A label at f points needs roughly f/72 inches of slot. So:
        max readable scale = (slot_inches × 72) / 12 = slot_inches × 6
    We take the tighter of the two axes' constraints, then cap at a sensible
    maximum so labels don't dominate the figure on small heatmaps.
    """
    width_in, height_in = figsize
    constraints = [0.6]  # cap so labels never get absurdly huge
    if n_visible_x_labels > 0:
        constraints.append((width_in / n_visible_x_labels) * 6)
    if n_visible_y_labels > 0:
        constraints.append((height_in / n_visible_y_labels) * 6)
    return max(min(constraints), 0.05)  # never go below the historical floor


# resolve_cluster_run_id moved to core.db — imported at the top of this file.
# Protein-first default matches the biology here: nucleotide runs capture
# within-species lineage structure at a finer grain, but the canonical
# pan-genome heatmap is the protein-cluster view.


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
    *,
    write_strain_order: bool = True,
) -> None:
    """Persist dendrogram orderings to display_order columns.

    Strategy:
        - Reset this run's clusters.display_order (scoped to the run — other
          runs' cluster orderings are preserved).
        - Assign 1..N to clusters in dendrogram order; excluded clusters
          append at the bottom in stable cluster_id order.
        - Optionally also reset + reassign strains.display_order. Skip that
          when called from a multi-run driver where only one heatmap
          (the protein one) should drive the strain axis.
    """
    with conn:
        conn.execute(
            "UPDATE clusters SET display_order = NULL WHERE cluster_run_id = ?",
            (cluster_run_id,),
        )
        conn.executemany(
            "UPDATE clusters SET display_order = ? WHERE cluster_id = ?",
            [(pos, cid) for pos, cid in enumerate(cluster_order, start=1)],
        )
        # Excluded clusters get appended at the end in stable cluster_id order.
        next_pos = len(cluster_order) + 1
        conn.executemany(
            "UPDATE clusters SET display_order = ? WHERE cluster_id = ?",
            [(next_pos + i, cid) for i, cid in enumerate(sorted(excluded_cluster_ids))],
        )

        if write_strain_order:
            conn.execute("UPDATE strains SET display_order = NULL")
            conn.executemany(
                "UPDATE strains SET display_order = ? WHERE strain_id = ?",
                [(pos, sid) for pos, sid in enumerate(strain_order, start=1)],
            )
