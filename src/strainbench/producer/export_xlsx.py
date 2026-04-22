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

from strainbench.core.db import (
    resolve_cluster_run_id as _resolve_cluster_run_id,
    resolve_nt_cluster_run_id as _resolve_nt_cluster_run_id,
)

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
        annotation_columns, extra_gene_names = _build_annotation_columns(
            conn, cluster_run_id, clusters_df
        )
    finally:
        conn.close()

    _write_xlsx(
        output_path, clusters_df, strains_df, wide,
        annotation_columns, extra_gene_names,
    )

    return ExportSummary(
        output_path=str(output_path),
        cluster_run_id=cluster_run_id,
        cluster_count=len(clusters_df),
        strain_count=len(strains_df),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────


# The resolver functions used to live here; they've moved to core.db to
# avoid the four-copies problem. The aliased imports at the top of this
# file keep existing call sites working without any further changes.


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

    NCBI frequently bakes the strain identifier into the organism name:
       species  = 'Lactobacillus iners AB-1'    + strain  = 'AB-1'
       species  = 'Gardnerella sp. DNF01159'    + strain  = 'DNF01159'
    Naively concatenating produces 'Lactobacillus iners AB-1 AB-1' which
    is ugly and confusing. We de-dupe here: if the species string already
    ends with the strain name, just use the species verbatim.
    """
    species = (row["species"] or "Unknown").strip()
    strain = (row["strain_name"] or "").strip()
    if strain and (species == strain or species.endswith(" " + strain)):
        display = species
    elif strain:
        display = f"{species} {strain}"
    else:
        display = species
    return (
        f"{row['locus_prefix']} | {display} "
        f"[{row['assembly_id'] or '?'}; "
        f"{row['biosample_id'] or '?'}; "
        f"{row['bioproject_id'] or '?'}; "
        f"{row['assembly_level'] or 'Unknown'}]"
    )


def _load_clusters(conn: sqlite3.Connection, cluster_run_id: int) -> pd.DataFrame:
    """Load cluster metadata plus the pooled NCBI gene names per cluster.

    The correlated subquery aggregates DISTINCT `cds.gene_name` across each
    cluster's members via GROUP_CONCAT. NULL/empty gene names are excluded,
    so clusters where no member had a /gene qualifier get a NULL result
    (rendered as 'none' at write time, matching strain-comp's convention).

    Sort priority:
      1. display_order (populated by hierarchical clustering, if run)
      2. member_count DESC (abundance)
      3. cluster_name (stable tiebreaker)
    """
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
            flags,
            (SELECT GROUP_CONCAT(DISTINCT c.gene_name)
             FROM cluster_membership m
             JOIN cds c ON c.cds_id = m.cds_id
             WHERE m.cluster_id = clusters.cluster_id
               AND c.gene_name IS NOT NULL
               AND c.gene_name != '')           AS pooled_gene_names
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
) -> tuple[list[tuple[str, dict[int, str]]], dict[int, set[str]]]:
    """Build extra columns for the xlsx from the cluster_annotations table.

    Returns:
      columns: list of (column_header, {cluster_id: cell_string}) pairs,
        one per (source, category) combination that has any rows.
      extra_gene_names: {cluster_id: set[str]} — gene symbols extracted from
        KEGG_KO rows' `name` field (eggNOG-mapper's Preferred_name). The
        caller merges these into the pan-genome 'NCBI gene names' column so
        biologists searching for short symbols (lacZ, dnaA) find them in
        the same place regardless of whether NCBI or eggNOG provided them.

    Returns ([], {}) if cluster_annotations is empty for this run.

    Column naming:
        source='GO',  category='MF'                    → 'GO_MF'
        source='COG', category=letter (J, K, L, …)     → 'COG'
        source='COG', category='group' (DeepNOG)       → 'COG_group'
        source='KEGG', category='pathway' / 'module'   → 'KEGG_{cat}'
        source='KEGG', category=NULL (KO entries)      → 'KEGG_KO'

    KEGG_KO special handling (the bit that handles emapper × kofamscan dup):
      - Same K-number from multiple tools → merged to a single entry.
      - When merging, prefer the form carrying an e-value: 'K00001 (e: 1e-50)'
        wins over the bare 'K00001' for the same cluster.
      - Gene names that eggNOG attaches to KEGG KOs (Preferred_name) are
        stripped from the KEGG_KO cell and surfaced in the NCBI gene names
        column instead.
    """
    import json as _json
    from collections import defaultdict

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
        return [], {}
    if ann.empty:
        return [], {}

    extra_gene_names: dict[int, set[str]] = defaultdict(set)

    def _kegg_evalue(extra_raw: object) -> float | None:
        """Pull below-threshold e-value from a KEGG row's `extra` JSON, or None."""
        if not extra_raw:
            return None
        try:
            extra = _json.loads(extra_raw)
        except (TypeError, ValueError):
            return None
        if extra.get("below_threshold") and extra.get("evalue") is not None:
            try:
                return float(extra["evalue"])
            except (TypeError, ValueError):
                return None
        return None

    def _entry_text(row: pd.Series) -> str:
        """Build the cell entry for non-KEGG_KO rows.

        For sources without code+name conflicts (COG categories, GO terms,
        Pfam IDs, etc.) the simple 'code name' form is correct. The KEGG_KO
        column gets dedicated dedup-by-code logic below — entry_text isn't
        used for it.
        """
        code = str(row.get("code") or "").strip()
        if not code:
            return ""
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

    # ── KEGG_KO: special dedup-by-code logic ─────────────────────────────────
    # For each cluster, group rows by K-number (code). When the same K-number
    # comes from multiple tools (emapper + kofamscan, common case), the row
    # whose `extra` carries a below-threshold e-value wins the tiebreak so
    # the e-value annotation 'K00001 (e: 1.2e-50)' is preserved. Gene names
    # (Preferred_name from emapper) are pulled out into extra_gene_names —
    # they belong in the NCBI gene names column, not stuck to the KO code.
    kegg_ko_rows = ann[ann["col"] == "KEGG_KO"]
    kegg_ko_map: dict[int, str] = {}
    if not kegg_ko_rows.empty:
        # Per-cluster: code → best-formatted entry for that code
        per_cluster: dict[int, dict[str, str]] = defaultdict(dict)
        for row in kegg_ko_rows.itertuples(index=False):
            code = (row.code or "").strip()
            if not code:
                continue
            cluster_id = int(row.cluster_id)
            # Pull gene name → extras (later merged into NCBI gene names column)
            name = (row.name or "").strip() if row.name else ""
            if name:
                extra_gene_names[cluster_id].add(name)
            # Build candidate entry: bare or with e-value
            evalue = _kegg_evalue(row.extra)
            candidate = f"{code} (e: {evalue:.1e})" if evalue is not None else code
            existing = per_cluster[cluster_id].get(code)
            if existing is None:
                per_cluster[cluster_id][code] = candidate
            elif "(e:" in candidate and "(e:" not in existing:
                # New row carries e-value, existing doesn't → upgrade to new
                per_cluster[cluster_id][code] = candidate
            # else: existing already has e-value (or both lack it) → keep existing
        kegg_ko_map = {
            cid: "; ".join(code_to_entry[c] for c in sorted(code_to_entry))
            for cid, code_to_entry in per_cluster.items()
        }

    # ── Non-KEGG_KO columns: existing aggregate-then-dedup-by-string logic ──
    grouped = (
        ann[ann["col"] != "KEGG_KO"]
        .groupby(["col", "cluster_id"])["entry"]
        .apply(lambda s: "; ".join(sorted({x for x in s if x})))
        .reset_index()
    )

    # ── KEGG_description column (KO definition text from kofamscan) ─────────
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
    if kegg_ko_map:
        seen.add("KEGG_KO")
    ordered = [c for c in preferred if c in seen] + sorted(c for c in seen if c not in preferred)

    columns: list[tuple[str, dict[int, str]]] = []
    for col in ordered:
        if col == "KEGG_description":
            columns.append((col, kegg_desc_map))
            continue
        if col == "KEGG_KO":
            columns.append((col, kegg_ko_map))
            continue
        sub = grouped[grouped["col"] == col]
        per_cluster = dict(zip(sub["cluster_id"], sub["entry"], strict=True))
        columns.append((col, per_cluster))
    return columns, dict(extra_gene_names)


# ─────────────────────────────────────────────────────────────────────────────
# Excel writing
# ─────────────────────────────────────────────────────────────────────────────


def _write_xlsx(
    output_path: Path,
    clusters_df: pd.DataFrame,
    strains_df: pd.DataFrame,
    wide: pd.DataFrame,
    annotation_columns: list[tuple[str, dict[int, str]]] | None = None,
    extra_gene_names: dict[int, set[str]] | None = None,
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
    # Three blocks left-to-right:
    #   1. Core metadata (10 columns, ending at 'strain count').
    #   2. Strain columns (one per strain, narrow + rotated 90° headers).
    #   3. Annotation columns (KEGG_KO, GO_MF, COG, etc.) at the far right.
    #
    # The frozen pane sits at the end of block 1 so you can scroll horizontally
    # through the strain block AND the annotations while keeping the cluster
    # name + count info pinned. Annotations are at the right edge because
    # they're the "look-at-me-once-I've-decided-which-cluster-matters" data —
    # putting them between metadata and strains made the frozen section
    # uncomfortably wide.
    core_metadata = [
        ("original_order",       fmt_cluster, 6),
        ("CLUSTER",              fmt_cluster, max(20, clusters_df["cluster_name"].str.len().max() + 2 if len(clusters_df) else 20)),
        # Pooled NCBI /gene qualifiers across all cluster members + KEGG-extracted
        # gene symbols. Biologists search spreadsheets for short gene symbols
        # (lacZ, dnaA, recA) more often than long product strings — so this
        # column sits left of NCBI annotation to make sort-then-scroll easy.
        ("NCBI gene names",      fmt_ncbi,    14),
        ("NCBI annotation",      fmt_ncbi,    45),
        ("GC%",                  fmt_lead,    6),
        ("GC spread",            fmt_lead,    6),
        ("aa length",            fmt_lead,    6),
        ("flags",                fmt_lead,    8),
        ("total count",          fmt_counts,  6),
        ("strain count",         fmt_counts_border, 6),
    ]
    n_core = len(core_metadata)

    for col_idx, (name, fmt, width) in enumerate(core_metadata):
        worksheet.write(0, col_idx, name, fmt)
        worksheet.set_column(col_idx, col_idx, width)

    # Block 2: strain headers (rotated, colored by assembly_level)
    for i, row in enumerate(strains_df.itertuples(index=False)):
        col_idx = n_core + i
        fmt = strain_header_formats.get(row.assembly_level, strain_header_formats["Unknown"])
        worksheet.write(0, col_idx, row.header, fmt)
        worksheet.set_column(col_idx, col_idx, 2)  # narrow strain columns

    n_strain = len(strains_df)
    annotation_section_start = n_core + n_strain

    # Block 3: annotation columns at the right edge
    annotation_columns = annotation_columns or []
    for offset, (col_name, _per_cluster) in enumerate(annotation_columns):
        col_idx = annotation_section_start + offset
        worksheet.write(0, col_idx, col_name, fmt_annot)
        worksheet.set_column(col_idx, col_idx, 24)

    # Freeze the cluster-identifier metadata only; strains + annotations scroll.
    worksheet.freeze_panes(1, n_core)

    # ── Data rows ────────────────────────────────────────────────────────────
    for r, (cluster_row, wide_row) in enumerate(
        zip(clusters_df.itertuples(index=False), wide.itertuples(index=False), strict=True),
        start=1,
    ):
        cluster_id = int(cluster_row.cluster_id)
        kegg_extras = extra_gene_names.get(cluster_id, set()) if extra_gene_names else set()
        # Block 1: core metadata (cols 0..9)
        worksheet.write(r, 0, r)                                          # original_order
        worksheet.write(r, 1, cluster_row.cluster_name)
        worksheet.write(r, 2, _format_pooled_gene_names(
            cluster_row.pooled_gene_names, extra=kegg_extras,
        ))
        worksheet.write(r, 3, cluster_row.consensus_annotation or "")
        worksheet.write(r, 4, cluster_row.gc_pct if cluster_row.gc_pct is not None else "")
        worksheet.write(r, 5, cluster_row.gc_spread if cluster_row.gc_spread is not None else "")
        worksheet.write(r, 6, cluster_row.aa_length if cluster_row.aa_length is not None else "")
        worksheet.write(r, 7, cluster_row.flags or "")
        worksheet.write(r, 8, cluster_row.total_count)
        worksheet.write(r, 9, cluster_row.strain_count, border_right)
        # Block 2: strain columns immediately after strain count
        for c, value in enumerate(wide_row, start=n_core):
            worksheet.write(r, c, value)
        # Block 3: annotation columns at the far right
        for offset, (_col_name, per_cluster) in enumerate(annotation_columns):
            worksheet.write(r, annotation_section_start + offset,
                            per_cluster.get(cluster_id, ""))

    # xlsxwriter only flushes to disk on close(); without this the file
    # is silently never created (Workbook's destructor doesn't close).
    workbook.close()


def _format_pooled_gene_names(
    raw: str | float | None,
    *,
    extra: set[str] | None = None,
) -> str:
    """Render the 'NCBI gene names' cell, merging gene names from two sources.

    Inputs:
      raw: GROUP_CONCAT output of distinct `cds.gene_name` values for this
        cluster (from gbff /gene qualifiers across all members in all strains).
        Comma-separated string from SQLite, or None/NaN when no NCBI gene
        names exist.
      extra: gene symbols pulled from KEGG annotations (e.g. eggNOG-mapper's
        Preferred_name). Already a Python set, sourced from
        `_build_annotation_columns`.

    Both inputs are merged into a single sorted set so duplicates collapse
    (e.g., NCBI's 'dnaA' + eggNOG's 'dnaA' → 'dnaA' once). When all sources
    are empty the cell renders as the literal 'none', matching strain-comp's
    sentinel so biologists can still filter for it.
    """
    pooled: set[str] = set()
    if raw is not None and not (isinstance(raw, float) and raw != raw):  # not NaN
        s = str(raw).strip()
        if s:
            pooled.update(part.strip() for part in s.split(",") if part.strip())
    if extra:
        pooled.update(extra)
    if not pooled:
        return "none"
    return ", ".join(sorted(pooled))

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
        # Same gene-symbol extraction as the pan-genome view, so the
        # strain-anchored 'NCBI gene names' column gets the same enrichment.
        _annotation_columns_unused, extra_gene_names = _build_annotation_columns(
            conn, cluster_run_id, clusters_df=anchor_cds[["cluster_id"]].drop_duplicates(),
        )
    finally:
        conn.close()

    _write_strain_anchored_xlsx(
        output_path, anchor_cds, cell_matrix, strains_df, anchor=anchor,
        has_nt_subcluster=(nt_cluster_run_id is not None),
        extra_gene_names=extra_gene_names,
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

    Also pools each cluster's NCBI gene names (across all members in all
    strains) into the `pooled_gene_names` column. The anchor CDS's own
    `gene_name` column is kept separately — so a biologist can see both
    "what NCBI called THIS copy" and "what NCBI called ANY copy of this
    gene family in the whole dataset".

    If `nt_cluster_run_id` is provided, also LEFT JOINs the matching nucleotide
    cluster name as the column `nt_cluster_name`. CDSs not found in the nt
    run (rare; only if some CDSs were missing from the nt clustering input)
    get NULL.
    """
    pooled_subquery = """(
        SELECT GROUP_CONCAT(DISTINCT c2.gene_name)
        FROM cluster_membership m2
        JOIN cds c2 ON c2.cds_id = m2.cds_id
        WHERE m2.cluster_id = cl.cluster_id
          AND c2.gene_name IS NOT NULL
          AND c2.gene_name != ''
    )"""
    if nt_cluster_run_id is None:
        return pd.read_sql_query(
            f"""
            SELECT
                c.cds_id, c.locus_tag, c.protein_id, c.nuc_accession,
                c.location, c.gene_name, c.aa_length,
                cl.cluster_id, cl.cluster_name, cl.consensus_annotation,
                cl.member_count, cl.strain_count,
                {pooled_subquery} AS pooled_gene_names,
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
        f"""
        SELECT
            c.cds_id, c.locus_tag, c.protein_id, c.nuc_accession,
            c.location, c.gene_name, c.aa_length,
            cl.cluster_id, cl.cluster_name, cl.consensus_annotation,
            cl.member_count, cl.strain_count,
            {pooled_subquery} AS pooled_gene_names,
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
    extra_gene_names: dict[int, set[str]] | None = None,
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
    # NCBI gene names is the pooled /gene qualifier column (distinct from
    # the per-row `gene` column which shows the anchor CDS's own value).
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
        ("NCBI gene names",   fmt_ncbi,    14),
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
        # Pooled NCBI gene names — merged with KEGG-extracted gene symbols
        # so a biologist scrolling row-by-row sees a unified pool. May be
        # absent from older anchor-CDS rows if the query wasn't updated;
        # guard with getattr to be safe.
        pooled = getattr(anchor_row, "pooled_gene_names", None)
        cluster_id = int(anchor_row.cluster_id)
        kegg_extras = (
            extra_gene_names.get(cluster_id, set()) if extra_gene_names else set()
        )
        worksheet.write(r, col,     _format_pooled_gene_names(pooled, extra=kegg_extras))
        worksheet.write(r, col + 1, anchor_row.consensus_annotation or "")
        worksheet.write(r, col + 2, anchor_row.aa_length if anchor_row.aa_length else "")
        worksheet.write(r, col + 3, anchor_row.member_count)
        worksheet.write(r, col + 4, anchor_row.strain_count)

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
