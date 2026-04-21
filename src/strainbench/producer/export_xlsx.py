"""Export the strain × cluster matrix from a strainbench DB to an Excel file.

This is the SQL-backed equivalent of strain-comp's `4_maketable.py` +
`6_formatxl.py` chain. It reads the populated `clusters`, `cds`, `strains`,
and `v_cluster_cell` tables for a given cluster_run and produces a single
xlsx with:

    * One row per cluster, ordered by member_count DESC.
    * One column per strain, headers rotated 90° and color-shaded by
      assembly_level (Complete=bright green, Scaffold=medium, Contig=light).
    * Per-cluster metadata columns (annotation, GC, length, counts) on the
      left, frozen so they stay visible while scrolling strain columns.
    * Cells filled with '[N]|locus_tag|protein_id, ...' for present hits,
      '*' for absent — matching strain-comp's cell convention so existing
      workflows (sorting, filtering in Excel) feel the same.

The intent is that a biologist can open the output in Excel and have the
exact same browsing experience as a strain-comp output.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import xlsxwriter

# Color hexes lifted from strain-comp's 6_formatxl.py so the visual style
# is unchanged for users moving over.
_COLOR_CLUSTER_BG = "black"
_COLOR_CLUSTER_FG = "white"
_COLOR_LEAD_INFO = "#5DEFF0"     # cyan — GC%, length, etc.
_COLOR_NCBI = "#FCD5B4"          # peach — NCBI annotations column
_COLOR_COUNTS = "#00F900"        # bright green — total count / strain count
_COLOR_ANNOT = "#FFBF00"         # amber — functional annotations (COG/KEGG/GO/Pfam/etc.)
_COLOR_LEVEL = {
    "Complete":   "#00F900",     # bright green (chromosomes too)
    "Chromosome": "#00F900",
    "Scaffold":   "#D4FB79",     # medium green
    "Contig":     "#E3FBD6",     # light green
    "Unknown":    "#EEEEEE",     # gray
}
_COLOR_ANCHOR = "#FFFF00"        # yellow — strain-anchored view: anchor's column


@dataclass
class ExportSummary:
    output_path: str
    cluster_run_id: int
    cluster_count: int
    strain_count: int


def export_cluster_table_xlsx(
    db_path: str | Path,
    output_path: str | Path,
    *,
    cluster_run_id: int | None = None,
) -> ExportSummary:
    """Pivot v_cluster_cell into a wide cluster × strain xlsx.

    Args:
        db_path: Path to a populated strainbench DB.
        output_path: Where to write the .xlsx file. Parent dir created if
            necessary.
        cluster_run_id: Which cluster_run to export. If None, picks the
            most recent active run.
    """
    db_path = Path(db_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        cluster_run_id = _resolve_cluster_run_id(conn, cluster_run_id)
        strains_df = _load_strains(conn)
        clusters_df = _load_clusters(conn, cluster_run_id)
        wide = _build_pivot(conn, cluster_run_id, clusters_df, strains_df)
        annotation_columns = _build_annotation_columns(conn, cluster_run_id, clusters_df)
    finally:
        conn.close()

    _write_xlsx(output_path, clusters_df, strains_df, wide, annotation_columns)

    return ExportSummary(
        output_path=str(output_path),
        cluster_run_id=cluster_run_id,
        cluster_count=len(clusters_df),
        strain_count=len(strains_df),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────


def _resolve_cluster_run_id(
    conn: sqlite3.Connection, requested: int | None
) -> int:
    """Pick a cluster_run for export. Prefers protein runs over nucleotide.

    With both protein and nucleotide runs in the DB, the protein run is the
    canonical "gene families" view that the cluster_table and strain-anchored
    xlsx are organized around. Nucleotide runs are an overlay (see the
    `nt_subcluster` column). User can always force a specific run via
    `cluster_run_id=N`.
    """
    if requested is not None:
        row = conn.execute(
            "SELECT cluster_run_id FROM cluster_runs WHERE cluster_run_id = ?",
            (requested,),
        ).fetchone()
        if row is None:
            raise ValueError(f"cluster_run_id={requested} not found in DB")
        return requested

    row = conn.execute(
        """
        SELECT cluster_run_id FROM cluster_runs
        WHERE is_active = 1
        ORDER BY (sequence_type = 'protein') DESC, cluster_run_id DESC
        LIMIT 1
        """
    ).fetchone()
    if row is None:
        raise ValueError(
            "No active cluster_runs in DB. Run `strainbench cluster` first."
        )
    return int(row["cluster_run_id"])


def _resolve_nt_cluster_run_id(
    conn: sqlite3.Connection, requested: int | None
) -> int | None:
    """Pick a nucleotide cluster_run to overlay (or None if none exists).

    Returns None if the user passed `requested=None` and the DB contains no
    active nucleotide run — caller should treat that as "no nt overlay column".
    """
    if requested is not None:
        return _resolve_cluster_run_id(conn, requested)
    row = conn.execute(
        """
        SELECT cluster_run_id FROM cluster_runs
        WHERE is_active = 1 AND sequence_type = 'nucleotide'
        ORDER BY cluster_run_id DESC LIMIT 1
        """
    ).fetchone()
    return int(row["cluster_run_id"]) if row else None


def _load_strains(conn: sqlite3.Connection) -> pd.DataFrame:
    # COALESCE so strains without a display_order (no heatmap run yet) sort
    # after those that do, in insertion order.
    df = pd.read_sql_query(
        """
        SELECT strain_id, locus_prefix, species, strain_name, assembly_id,
               biosample_id, bioproject_id, assembly_level
        FROM strains
        ORDER BY COALESCE(display_order, 999999), strain_id
        """,
        conn,
    )
    df["header"] = df.apply(_format_strain_header, axis=1)
    return df


def _format_strain_header(row: pd.Series) -> str:
    """Build the spreadsheet column header for one strain.

    Format mirrors strain-comp's strainlist.txt convention so the existing
    `; Complete]` / `; Scaffold]` / `; Contig]` markers can drive header
    coloring downstream.
    """
    return (
        f"{row['locus_prefix']} | {row['species']} {row['strain_name']} "
        f"[{row['assembly_id'] or '?'}; "
        f"{row['biosample_id'] or '?'}; "
        f"{row['bioproject_id'] or '?'}; "
        f"{row['assembly_level'] or 'Unknown'}]"
    )


def _load_clusters(conn: sqlite3.Connection, cluster_run_id: int) -> pd.DataFrame:
    # If display_order is set (from a heatmap run), it dominates.
    # Otherwise sort by abundance, matching the pre-heatmap default.
    return pd.read_sql_query(
        """
        SELECT
            cluster_id,
            cluster_name,
            consensus_annotation,
            member_count   AS total_count,
            strain_count,
            ROUND(avg_gc_pct, 2)  AS gc_pct,
            ROUND(gc_spread, 2)   AS gc_spread,
            avg_aa_length        AS aa_length,
            flags
        FROM clusters
        WHERE cluster_run_id = ?
        ORDER BY
            CASE WHEN display_order IS NOT NULL THEN 0 ELSE 1 END,
            display_order,
            member_count DESC,
            cluster_name
        """,
        conn,
        params=(cluster_run_id,),
    )


def _build_pivot(
    conn: sqlite3.Connection,
    cluster_run_id: int,
    clusters_df: pd.DataFrame,
    strains_df: pd.DataFrame,
) -> pd.DataFrame:
    """Return a cluster_id × strain_id DataFrame of cell strings.

    Cell format: '[N]|hit_summary' for present, '*' for absent.
    """
    cells_df = pd.read_sql_query(
        """
        SELECT cluster_id, strain_id, hit_count, hit_summary
        FROM v_cluster_cell
        WHERE cluster_run_id = ?
        """,
        conn,
        params=(cluster_run_id,),
    )
    # Format cells with strain-comp's [N]|... prefix
    cells_df["cell"] = (
        "[" + cells_df["hit_count"].astype(str) + "]|" + cells_df["hit_summary"]
    )
    pivot = cells_df.pivot(index="cluster_id", columns="strain_id", values="cell")
    # Reindex so missing (cluster, strain) cells become "*", and column/row
    # order matches the metadata DataFrames.
    pivot = pivot.reindex(
        index=clusters_df["cluster_id"],
        columns=strains_df["strain_id"],
    )
    return pivot.fillna("*")


def _build_annotation_columns(
    conn: sqlite3.Connection,
    cluster_run_id: int,
    clusters_df: pd.DataFrame,
) -> list[tuple[str, dict[int, str]]]:
    """Build extra columns for the xlsx from the cluster_annotations table.

    Returns a list of (column_header, {cluster_id: cell_string}) pairs,
    one per (source, category) combination that has any rows. Returns []
    if cluster_annotations is empty for this run.

    Column naming:
        source='GO', category='MF'  → 'GO_MF'
        source='COG', category=None or any → 'COG'
        source='KEGG', category='pathway' → 'KEGG_pathway'
        source='KEGG', category=None (KO) → 'KEGG_KO'
    """
    try:
        ann = pd.read_sql_query(
            """
            SELECT
                a.cluster_id, a.source, a.category, a.code, a.name,
                a.description, a.extra
            FROM cluster_annotations a
            JOIN clusters c ON c.cluster_id = a.cluster_id
            WHERE c.cluster_run_id = ?
            """,
            conn,
            params=(cluster_run_id,),
        )
    except Exception:
        return []
    if ann.empty:
        return []

    import json as _json

    def _entry_text(row: pd.Series) -> str:
        """Build the cell entry for this annotation row.

        Below-threshold KEGG hits get formatted as 'KOxxxxx (e: 1.2e-50)' so
        biologists can see at a glance which calls are guesses vs. confident.
        Other rows get 'code name' (or just 'code' if name missing).
        """
        code = str(row.get("code") or "").strip()
        if not code:
            return ""
        # Below-threshold flagging (currently only kofamscan emits this).
        extra_raw = row.get("extra")
        if extra_raw:
            try:
                extra = _json.loads(extra_raw)
            except (TypeError, ValueError):
                extra = {}
            if extra.get("below_threshold"):
                evalue = extra.get("evalue")
                if evalue is not None:
                    return f"{code} (e: {evalue:.1e})"
                return f"{code} (sub-significant)"
        name = str(row.get("name") or "").strip()
        return f"{code} {name}".strip()

    ann["entry"] = ann.apply(_entry_text, axis=1)

    # Column name derivation
    def _col_name(source: str, category: str | None) -> str:
        if source == "GO":
            return f"GO_{category}" if category in {"MF", "BP", "CC"} else "GO"
        if source == "KEGG":
            if category in {"pathway", "module"}:
                return f"KEGG_{category}"
            return "KEGG_KO"
        if source == "COG":
            # eggNOG-mapper produces single-letter COG categories
            # (category in {'J','K','L',...}); DeepNOG produces COG groups
            # (category='group', code='COG1234'). Split into two columns.
            return "COG_group" if category == "group" else "COG"
        return source

    ann["col"] = ann.apply(lambda r: _col_name(r["source"], r.get("category")), axis=1)

    grouped = (
        ann.groupby(["col", "cluster_id"])["entry"]
        .apply(lambda s: "; ".join(sorted(x for x in s if x)))
        .reset_index()
    )

    # Build a parallel "description" column for KEGG, populated from the
    # `description` field (kofamscan's KO definition text). This gives the
    # user a column showing what each KO actually IS, alongside the bare KO ID.
    kegg_desc_rows = ann[(ann["source"] == "KEGG") & ann["description"].notna() & (ann["description"].str.strip() != "")]
    kegg_desc_map: dict[int, str] = {}
    if not kegg_desc_rows.empty:
        descs_per_cluster = (
            kegg_desc_rows.groupby("cluster_id")["description"]
            .apply(lambda s: "; ".join(sorted({x.strip() for x in s if x.strip()})))
        )
        kegg_desc_map = descs_per_cluster.to_dict()

    # Stable column order: put KEGG_KO first, then KEGG_description right next
    # to it (so users see KO + what it is in adjacent cells).
    preferred = ["KEGG_KO", "KEGG_description", "KEGG_pathway", "KEGG_module",
                 "COG", "COG_group", "GO_MF", "GO_BP", "GO_CC",
                 "Pfam", "TIGRFAM", "PANTHER", "SUPERFAMILY", "SMART",
                 "EC", "CAZy",
                 "AMRFinder", "DefenseFinder", "PADLOC"]
    seen = set(grouped["col"].unique())
    if kegg_desc_map:
        seen.add("KEGG_description")
    ordered = [c for c in preferred if c in seen] + sorted(c for c in seen if c not in preferred)

    columns: list[tuple[str, dict[int, str]]] = []
    for col in ordered:
        if col == "KEGG_description":
            columns.append((col, kegg_desc_map))
            continue
        sub = grouped[grouped["col"] == col]
        per_cluster = dict(zip(sub["cluster_id"], sub["entry"], strict=True))
        columns.append((col, per_cluster))
    return columns


# ─────────────────────────────────────────────────────────────────────────────
# Excel writing
# ─────────────────────────────────────────────────────────────────────────────


def _write_xlsx(
    output_path: Path,
    clusters_df: pd.DataFrame,
    strains_df: pd.DataFrame,
    wide: pd.DataFrame,
    annotation_columns: list[tuple[str, dict[int, str]]] | None = None,
) -> None:
    workbook = xlsxwriter.Workbook(str(output_path), {"strings_to_formulas": False})
    worksheet = workbook.add_worksheet("cluster_table")

    base = {"bold": True, "text_wrap": True, "valign": "bottom", "border": 1}

    fmt_cluster = workbook.add_format({**base, "fg_color": _COLOR_CLUSTER_BG, "font_color": _COLOR_CLUSTER_FG})
    fmt_lead = workbook.add_format({**base, "fg_color": _COLOR_LEAD_INFO})
    fmt_ncbi = workbook.add_format({**base, "fg_color": _COLOR_NCBI})
    fmt_annot = workbook.add_format({**base, "fg_color": _COLOR_ANNOT})
    fmt_counts = workbook.add_format({**base, "fg_color": _COLOR_COUNTS})
    fmt_counts_border = workbook.add_format({**base, "fg_color": _COLOR_COUNTS, "right": 2})
    border_right = workbook.add_format({"right": 2})

    # Strain header formats per assembly level (rotated 90°)
    strain_header_formats = {
        level: workbook.add_format({**base, "rotation": 90, "fg_color": color})
        for level, color in _COLOR_LEVEL.items()
    }

    # ── Header row layout ────────────────────────────────────────────────────
    # Columns: original_order | CLUSTER | annotation | GC% | spread | aa_len |
    #          flags | total | strain_count | [annotation columns…] |
    #          <strain1> | <strain2> | …
    metadata_columns = [
        ("original_order",       fmt_cluster, 6),
        ("CLUSTER",              fmt_cluster, max(20, clusters_df["cluster_name"].str.len().max() + 2 if len(clusters_df) else 20)),
        ("NCBI annotation",      fmt_ncbi,    45),
        ("GC%",                  fmt_lead,    6),
        ("GC spread",            fmt_lead,    6),
        ("aa length",            fmt_lead,    6),
        ("flags",                fmt_lead,    8),
        ("total count",          fmt_counts,  6),
        ("strain count",         fmt_counts_border, 6),
    ]
    # Append annotation columns if present (one per source/aspect pair).
    annotation_columns = annotation_columns or []
    for col_name, _ in annotation_columns:
        metadata_columns.append((col_name, fmt_annot, 24))

    for col_idx, (name, fmt, width) in enumerate(metadata_columns):
        worksheet.write(0, col_idx, name, fmt)
        worksheet.set_column(col_idx, col_idx, width)

    n_meta = len(metadata_columns)

    # Strain headers (rotated, colored by assembly_level)
    for i, row in enumerate(strains_df.itertuples(index=False)):
        col_idx = n_meta + i
        fmt = strain_header_formats.get(row.assembly_level, strain_header_formats["Unknown"])
        worksheet.write(0, col_idx, row.header, fmt)
        worksheet.set_column(col_idx, col_idx, 2)  # narrow strain columns

    # Freeze: 1 header row + (n_meta) metadata columns
    worksheet.freeze_panes(1, n_meta)

    # ── Data rows ────────────────────────────────────────────────────────────
    for r, (cluster_row, wide_row) in enumerate(
        zip(clusters_df.itertuples(index=False), wide.itertuples(index=False), strict=True),
        start=1,
    ):
        worksheet.write(r, 0, r)                                          # original_order
        worksheet.write(r, 1, cluster_row.cluster_name)
        worksheet.write(r, 2, cluster_row.consensus_annotation or "")
        worksheet.write(r, 3, cluster_row.gc_pct if cluster_row.gc_pct is not None else "")
        worksheet.write(r, 4, cluster_row.gc_spread if cluster_row.gc_spread is not None else "")
        worksheet.write(r, 5, cluster_row.aa_length if cluster_row.aa_length is not None else "")
        worksheet.write(r, 6, cluster_row.flags or "")
        worksheet.write(r, 7, cluster_row.total_count)
        worksheet.write(r, 8, cluster_row.strain_count, border_right)
        # Annotation columns, keyed by cluster_id
        col_idx = 9
        for _col_name, per_cluster in annotation_columns:
            worksheet.write(r, col_idx, per_cluster.get(cluster_row.cluster_id, ""))
            col_idx += 1
        # Strain columns start at col_idx (== n_meta after the annotation block)
        for c, value in enumerate(wide_row, start=col_idx):
            worksheet.write(r, c, value)

    workbook.close()


# ─────────────────────────────────────────────────────────────────────────────
# Strain-anchored view
# ─────────────────────────────────────────────────────────────────────────────


def export_strain_anchored_xlsx(
    db_path: str | Path,
    output_path: str | Path,
    strain_locus_prefix: str,
    *,
    cluster_run_id: int | None = None,
    nt_cluster_run_id: int | None = None,
) -> ExportSummary:
    """Export a CDS-anchored matrix for ONE strain.

    Each row is one CDS in the anchor strain, in genomic order
    (replicon-major + position-within-replicon, derived from `cds_id` order).
    Columns are all strains (anchor first, then dendrogram order if a heatmap
    has been run). Cells show what each strain has in the cluster that the
    row's anchor CDS belongs to.

    Multi-copy genes get separate rows — solving the "[N] cells clump together
    when sorting by strain column" problem in the pan-genome view.

    Args:
        db_path: Path to a populated strainbench DB.
        output_path: Where to write the .xlsx.
        strain_locus_prefix: e.g. 'HMPREF0520'. Must exist in the strains table.
        cluster_run_id: Which cluster_run to use. None → most recent active.
    """
    db_path = Path(db_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        cluster_run_id = _resolve_cluster_run_id(conn, cluster_run_id)
        nt_cluster_run_id = _resolve_nt_cluster_run_id(conn, nt_cluster_run_id)
        anchor = conn.execute(
            "SELECT strain_id, locus_prefix, strain_name, species "
            "FROM strains WHERE locus_prefix = ?",
            (strain_locus_prefix,),
        ).fetchone()
        if anchor is None:
            raise ValueError(f"strain {strain_locus_prefix!r} not found in DB")

        anchor_cds = _load_anchor_cds(
            conn, anchor["strain_id"], cluster_run_id, nt_cluster_run_id,
        )
        if anchor_cds.empty:
            raise ValueError(
                f"strain {strain_locus_prefix!r} has no clustered CDSs in "
                f"cluster_run_id={cluster_run_id}. Did the cluster step run?"
            )

        strains_df = _load_strains(conn)
        cell_matrix = _load_cells_for_clusters(
            conn, cluster_run_id, anchor_cds["cluster_id"].unique().tolist(), strains_df,
        )
    finally:
        conn.close()

    _write_strain_anchored_xlsx(
        output_path, anchor_cds, cell_matrix, strains_df, anchor=anchor,
        has_nt_subcluster=(nt_cluster_run_id is not None),
    )

    return ExportSummary(
        output_path=str(output_path),
        cluster_run_id=cluster_run_id,
        cluster_count=len(anchor_cds),  # rows = anchor CDS count
        strain_count=len(strains_df),
    )


def _load_anchor_cds(
    conn: sqlite3.Connection,
    anchor_strain_id: int,
    cluster_run_id: int,
    nt_cluster_run_id: int | None = None,
) -> pd.DataFrame:
    """Load the anchor strain's CDSs joined to their protein cluster.

    If `nt_cluster_run_id` is provided, also LEFT JOINs the matching nucleotide
    cluster name as the column `nt_cluster_name`. CDSs not found in the nt
    run (rare; only if some CDSs were missing from the nt clustering input)
    get NULL.
    """
    if nt_cluster_run_id is None:
        return pd.read_sql_query(
            """
            SELECT
                c.cds_id, c.locus_tag, c.protein_id, c.nuc_accession,
                c.location, c.gene_name, c.aa_length,
                cl.cluster_id, cl.cluster_name, cl.consensus_annotation,
                cl.member_count, cl.strain_count,
                NULL AS nt_cluster_name
            FROM cds c
            JOIN cluster_membership m ON m.cds_id = c.cds_id AND m.cluster_run_id = ?
            JOIN clusters cl          ON cl.cluster_id = m.cluster_id
            WHERE c.strain_id = ?
            ORDER BY c.cds_id
            """,
            conn,
            params=(cluster_run_id, anchor_strain_id),
        )
    return pd.read_sql_query(
        """
        SELECT
            c.cds_id, c.locus_tag, c.protein_id, c.nuc_accession,
            c.location, c.gene_name, c.aa_length,
            cl.cluster_id, cl.cluster_name, cl.consensus_annotation,
            cl.member_count, cl.strain_count,
            nt_cl.cluster_name AS nt_cluster_name
        FROM cds c
        JOIN cluster_membership m  ON m.cds_id = c.cds_id  AND m.cluster_run_id = ?
        JOIN clusters cl           ON cl.cluster_id = m.cluster_id
        LEFT JOIN cluster_membership nm ON nm.cds_id = c.cds_id AND nm.cluster_run_id = ?
        LEFT JOIN clusters nt_cl        ON nt_cl.cluster_id = nm.cluster_id
        WHERE c.strain_id = ?
        ORDER BY c.cds_id
        """,
        conn,
        params=(cluster_run_id, nt_cluster_run_id, anchor_strain_id),
    )


def _load_cells_for_clusters(
    conn: sqlite3.Connection,
    cluster_run_id: int,
    cluster_ids: list[int],
    strains_df: pd.DataFrame,
) -> pd.DataFrame:
    """Return a DataFrame indexed by cluster_id, columns = strain_ids,
    cells = '[N]|hit_summary' or '*'."""
    placeholders = ",".join("?" * len(cluster_ids))
    cells = pd.read_sql_query(
        f"""
        SELECT vc.cluster_id, vc.strain_id, vc.hit_count, vc.hit_summary
        FROM v_cluster_cell vc
        WHERE vc.cluster_run_id = ?
          AND vc.cluster_id IN ({placeholders})
        """,
        conn,
        params=[cluster_run_id, *cluster_ids],
    )
    cells["cell"] = "[" + cells["hit_count"].astype(str) + "]|" + cells["hit_summary"]
    pivot = cells.pivot(index="cluster_id", columns="strain_id", values="cell")
    return pivot.reindex(columns=strains_df["strain_id"]).fillna("*")


def _write_strain_anchored_xlsx(
    output_path: Path,
    anchor_cds: pd.DataFrame,
    cell_matrix: pd.DataFrame,
    strains_df: pd.DataFrame,
    *,
    anchor: sqlite3.Row,
    has_nt_subcluster: bool = False,
) -> None:
    workbook = xlsxwriter.Workbook(str(output_path), {"strings_to_formulas": False})
    sheet_name = f"{anchor['locus_prefix']}_anchored"[:31]  # Excel sheet-name cap
    worksheet = workbook.add_worksheet(sheet_name)

    base = {"bold": True, "text_wrap": True, "valign": "bottom", "border": 1}
    rotated = {**base, "rotation": 90}

    fmt_cluster = workbook.add_format({**base, "fg_color": _COLOR_CLUSTER_BG, "font_color": _COLOR_CLUSTER_FG})
    fmt_lead = workbook.add_format({**base, "fg_color": _COLOR_LEAD_INFO})
    fmt_ncbi = workbook.add_format({**base, "fg_color": _COLOR_NCBI})
    fmt_counts = workbook.add_format({**base, "fg_color": _COLOR_COUNTS})
    fmt_counts_border = workbook.add_format({**base, "fg_color": _COLOR_COUNTS, "right": 2})

    fmt_anchor_header = workbook.add_format({**rotated, "fg_color": _COLOR_ANCHOR})
    fmt_anchor_cell = workbook.add_format({"fg_color": _COLOR_ANCHOR})
    strain_header_formats = {
        level: workbook.add_format({**rotated, "fg_color": color})
        for level, color in _COLOR_LEVEL.items()
    }

    # Metadata column layout — different from the pan-genome view.
    # nt_subcluster (when present) sits right after the protein CLUSTER
    # so a biologist can see at a glance: protein family + nt lineage.
    metadata_columns: list[tuple[str, object, int]] = [
        ("position",          fmt_cluster, 6),
        ("anchor_locus_tag",  fmt_cluster, 22),
        ("nuc_accession",     fmt_lead,    18),
        ("location",          fmt_lead,    18),
        ("gene",              fmt_lead,    8),
        ("CLUSTER",           fmt_cluster, 16),
    ]
    if has_nt_subcluster:
        metadata_columns.append(("nt_subcluster", fmt_cluster, 18))
    metadata_columns.extend([
        ("NCBI annotation",   fmt_ncbi,    45),
        ("aa length",         fmt_lead,    6),
        ("cluster total",     fmt_counts,  7),
        ("cluster strains",   fmt_counts_border, 7),
    ])
    for col_idx, (name, fmt, width) in enumerate(metadata_columns):
        worksheet.write(0, col_idx, name, fmt)
        worksheet.set_column(col_idx, col_idx, width)
    n_meta = len(metadata_columns)

    # Strain columns — anchor first (highlighted), then the rest in dendrogram
    # / display order. Caller's strains_df is already in that order (display_order
    # ASC, falling back to strain_id) thanks to _load_strains().
    anchor_strain_id = anchor["strain_id"]
    other_strains = strains_df[strains_df["strain_id"] != anchor_strain_id]
    ordered_strains = pd.concat(
        [strains_df[strains_df["strain_id"] == anchor_strain_id], other_strains],
        ignore_index=True,
    )

    strain_col_index: dict[int, int] = {}
    for i, row in enumerate(ordered_strains.itertuples(index=False)):
        col_idx = n_meta + i
        if row.strain_id == anchor_strain_id:
            fmt = fmt_anchor_header
            header = "★ " + row.header  # mark the anchor visually too
        else:
            fmt = strain_header_formats.get(row.assembly_level, strain_header_formats["Unknown"])
            header = row.header
        worksheet.write(0, col_idx, header, fmt)
        worksheet.set_column(col_idx, col_idx, 2)
        strain_col_index[row.strain_id] = col_idx

    worksheet.freeze_panes(1, n_meta + 1)  # freeze metadata + anchor strain column

    # Data rows
    for r, anchor_row in enumerate(anchor_cds.itertuples(index=False), start=1):
        cells_for_strains = cell_matrix.loc[anchor_row.cluster_id]

        worksheet.write(r, 0, r)                                       # position (1..N)
        worksheet.write(r, 1, anchor_row.locus_tag)
        worksheet.write(r, 2, anchor_row.nuc_accession or "")
        worksheet.write(r, 3, anchor_row.location or "")
        worksheet.write(r, 4, anchor_row.gene_name or "")
        worksheet.write(r, 5, anchor_row.cluster_name)
        col = 6
        if has_nt_subcluster:
            nt_name = getattr(anchor_row, "nt_cluster_name", None)
            worksheet.write(r, col, nt_name if nt_name else "")
            col += 1
        worksheet.write(r, col,     anchor_row.consensus_annotation or "")
        worksheet.write(r, col + 1, anchor_row.aa_length if anchor_row.aa_length else "")
        worksheet.write(r, col + 2, anchor_row.member_count)
        worksheet.write(r, col + 3, anchor_row.strain_count)

        for strain_id, col_idx in strain_col_index.items():
            if strain_id == anchor_strain_id:
                # Override: show this row's specific CDS (so multi-copy rows
                # are visually distinguishable), highlighted yellow.
                anchor_cell = (
                    f"[1]|{anchor_row.locus_tag}|{anchor_row.protein_id or ''}"
                )
                worksheet.write(r, col_idx, anchor_cell, fmt_anchor_cell)
            else:
                worksheet.write(r, col_idx, cells_for_strains[strain_id])

    workbook.close()
