"""Import functional annotations onto cluster representatives.

The external annotation workflow is:
    1. `strainbench export-reps` writes cluster rep sequences as FASTA.
    2. User runs an annotator in its own conda env (kofamscan, eggNOG-mapper,
       DeepNOG, etc.) on that FASTA.
    3. `strainbench import-annotation --format <fmt> --input <tsv>` parses the
       output and writes rows into the `cluster_annotations` table.
    3'. OR: `strainbench import-annotation --workdir <dir>` scans a directory
        produced by `annotate-prep`, auto-detects formats, imports all at once.

Each format parser is a function that reads the tool's output file and yields
`ClusterAnnotationRow` tuples ready for bulk insertion.

Formats supported:
    emapper             — eggNOG-mapper v2 .emapper.annotations (TSV).
                          One row yields multiple annotations (COG, KEGG, GO,
                          Pfam, EC, CAZy), each with source= the respective
                          database name.
    kofamscan           — kofamscan detail-tsv format. One cluster → one KEGG KO.
    deepnog             — DeepNOG classify output (TSV). One cluster → one COG group.
    amrfinder           — NCBI AMRFinderPlus TSV. AMR/STRESS/VIRULENCE gene hits.
    defensefinder       — DefenseFinder defense_finder_genes.tsv. Antiphage systems.
    padloc              — PADLOC CSV output. Antiphage + toxin-antitoxin systems.
    interproscan        — InterProScan TSV. Domain-level annotations from many
                          databases (Pfam, TIGRFAM, PANTHER, …) plus InterPro2GO.
"""

from __future__ import annotations

import csv
import sqlite3
import sys
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class ClusterAnnotationRow:
    cluster_name: str
    source: str
    code: str | None = None
    category: str | None = None
    name: str | None = None
    description: str | None = None
    score: float | None = None
    extra: str | None = None


@dataclass
class ImportSummary:
    source_formats: list[str]
    input_path: str
    rows_inserted: int
    sources_written: dict[str, int]   # 'KEGG' → 1520, 'GO' → 18432, etc.
    clusters_touched: int
    unmatched_clusters: int           # cluster_names in the TSV that aren't in the DB


def import_annotation_file(
    db_path: str | Path,
    input_path: str | Path,
    *,
    source_format: str,
    cluster_run_id: int | None = None,
    replace: bool = False,
) -> ImportSummary:
    """Parse an annotation tool's output and insert rows into `cluster_annotations`.

    Args:
        db_path: Path to a populated strainbench DB.
        input_path: Path to the annotator's output file.
        source_format: 'emapper' or 'kofamscan'.
        cluster_run_id: None → most recent active run.
        replace: If True, delete all existing annotations from the same
            source(s) on these clusters before inserting. Use when re-running
            an annotator with different parameters.
    """
    parser = _PARSERS.get(source_format)
    if parser is None:
        raise ValueError(
            f"unknown source_format={source_format!r}; "
            f"available: {sorted(_PARSERS)}"
        )

    db_path = Path(db_path)
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"annotation file not found: {input_path}")

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        cluster_run_id = _resolve_cluster_run_id(conn, cluster_run_id)
        name_to_cluster_id = _load_cluster_name_map(conn, cluster_run_id)

        # Accumulate rows, track stats.
        rows_to_insert: list[tuple] = []
        sources_written: dict[str, int] = {}
        clusters_touched: set[int] = set()
        unmatched: set[str] = set()

        for ann in parser(input_path):
            cluster_id = name_to_cluster_id.get(ann.cluster_name)
            if cluster_id is None:
                unmatched.add(ann.cluster_name)
                continue
            rows_to_insert.append(
                (
                    cluster_id, ann.source, ann.code, ann.category,
                    ann.name, ann.description, ann.score, ann.extra,
                )
            )
            sources_written[ann.source] = sources_written.get(ann.source, 0) + 1
            clusters_touched.add(cluster_id)

        # Tag each row with provenance (source_file, tool_format) so we can
        # later answer "where did this annotation come from?".
        source_file_str = str(input_path.resolve())
        rows_with_provenance = [
            (*row, source_file_str, source_format) for row in rows_to_insert
        ]

        with conn:
            if replace:
                # Scope the delete by tool_format rather than source so that
                # re-importing emapper's COG categories doesn't also wipe
                # DeepNOG's COG groups (both have source='COG' but different
                # tool_format). "Replace this tool's previous contribution,
                # preserve other tools' contributions to the same source."
                conn.execute(
                    """
                    DELETE FROM cluster_annotations
                    WHERE tool_format = ?
                      AND cluster_id IN (
                          SELECT cluster_id FROM clusters WHERE cluster_run_id = ?
                      )
                    """,
                    (source_format, cluster_run_id),
                )

            conn.executemany(
                """
                INSERT INTO cluster_annotations
                    (cluster_id, source, code, category, name, description,
                     score, extra, source_file, tool_format)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                rows_with_provenance,
            )
    finally:
        conn.close()

    return ImportSummary(
        source_formats=[source_format],
        input_path=str(input_path),
        rows_inserted=len(rows_to_insert),
        sources_written=sources_written,
        clusters_touched=len(clusters_touched),
        unmatched_clusters=len(unmatched),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Parsers
# ─────────────────────────────────────────────────────────────────────────────


def parse_emapper(path: Path) -> Iterable[ClusterAnnotationRow]:
    """Parse eggNOG-mapper v2 `.emapper.annotations` output.

    Each row contains multiple annotation fields; we emit one
    ClusterAnnotationRow per non-empty field so each annotation type ends up
    as its own DB row.

    Column layout (eggNOG-mapper v2 default):
        #query  seed_ortholog  evalue  score  eggNOG_OGs  max_annot_lvl
        COG_category  Description  Preferred_name  GOs  EC  KEGG_ko
        KEGG_Pathway  KEGG_Module  KEGG_Reaction  KEGG_rclass  BRITE
        KEGG_TC  CAZy  BiGG_Reaction  PFAMs
    """
    go_lookup = _load_go_aspects()

    with path.open() as f:
        # eggNOG-mapper prefixes comments with '##' and header with '#query...'.
        lines = [line.rstrip("\n") for line in f if not line.startswith("##")]
    if not lines:
        return
    # Header has a leading '#'. strip it.
    header = lines[0].lstrip("#").split("\t")
    reader = csv.DictReader(lines[1:], fieldnames=header, delimiter="\t")

    def _clean(value: str | None) -> str | None:
        if value is None:
            return None
        v = value.strip()
        return None if v in {"", "-"} else v

    for row in reader:
        cluster_name = _clean(row.get("query"))
        if cluster_name is None:
            continue

        try:
            score = float(row.get("score", "") or "nan")
        except ValueError:
            score = None
        if score != score:  # NaN
            score = None

        seed_ortholog = _clean(row.get("seed_ortholog"))
        description = _clean(row.get("Description"))
        preferred_name = _clean(row.get("Preferred_name"))

        # 1. COG category (single letters or compound like 'KLJ')
        cog = _clean(row.get("COG_category"))
        if cog:
            for letter in cog:
                yield ClusterAnnotationRow(
                    cluster_name=cluster_name,
                    source="COG",
                    code=letter,
                    category=letter,
                    name=_COG_CATEGORY_NAMES.get(letter),
                    description=description,
                    score=score,
                    extra=f'{{"seed_ortholog": "{seed_ortholog or ""}"}}',
                )

        # 2. KEGG KOs
        kegg_ko = _clean(row.get("KEGG_ko"))
        if kegg_ko:
            for ko in kegg_ko.split(","):
                ko_clean = ko.strip().removeprefix("ko:")
                if ko_clean:
                    yield ClusterAnnotationRow(
                        cluster_name=cluster_name, source="KEGG", code=ko_clean,
                        category=None, name=preferred_name, description=description,
                        score=score, extra=None,
                    )

        # 3. KEGG Pathway, Module
        for col, cat_label in (("KEGG_Pathway", "pathway"), ("KEGG_Module", "module")):
            v = _clean(row.get(col))
            if v:
                for entry in v.split(","):
                    e = entry.strip()
                    if e:
                        yield ClusterAnnotationRow(
                            cluster_name=cluster_name, source="KEGG", code=e,
                            category=cat_label, name=None, description=None,
                            score=None, extra=None,
                        )

        # 4. GO terms (split by aspect)
        go_ids = _clean(row.get("GOs"))
        if go_ids:
            for go_id in go_ids.split(","):
                g = go_id.strip()
                if not g:
                    continue
                aspect, go_name = go_lookup.get(g, (None, None))
                yield ClusterAnnotationRow(
                    cluster_name=cluster_name, source="GO", code=g,
                    category=aspect, name=go_name, description=None,
                    score=None, extra=None,
                )

        # 5. EC numbers
        ec = _clean(row.get("EC"))
        if ec:
            for e in ec.split(","):
                e_clean = e.strip()
                if e_clean:
                    yield ClusterAnnotationRow(
                        cluster_name=cluster_name, source="EC", code=e_clean,
                        category=None, name=None, description=None,
                        score=None, extra=None,
                    )

        # 6. Pfam domains
        pfams = _clean(row.get("PFAMs"))
        if pfams:
            for pf in pfams.split(","):
                pf_clean = pf.strip()
                if pf_clean:
                    yield ClusterAnnotationRow(
                        cluster_name=cluster_name, source="Pfam", code=pf_clean,
                        category=None, name=pf_clean, description=None,
                        score=None, extra=None,
                    )

        # 7. CAZy families
        cazy = _clean(row.get("CAZy"))
        if cazy:
            for c in cazy.split(","):
                c_clean = c.strip()
                if c_clean:
                    yield ClusterAnnotationRow(
                        cluster_name=cluster_name, source="CAZy", code=c_clean,
                        category=None, name=None, description=None,
                        score=None, extra=None,
                    )


def parse_kofamscan(path: Path) -> Iterable[ClusterAnnotationRow]:
    """Parse kofamscan detail-tsv. Two-pass: assigned hits + best-candidate fallback.

    Format:
        # gene name       KO      thrshld score   E-value "KO definition"
        * INERS_000001    K00001  256.30  312.5   1.2e-90 "alcohol dehydrogenase"
          INERS_000001    K00002  100.00   50.5   1.0e-10 "below-threshold candidate"

    Strategy (matches strain-comp's existing behavior):
      * For clusters with one or more `*` hits (above kofamscan's per-KO
        threshold): emit ALL of them.
      * For clusters with NO `*` hit: emit the single best (highest score)
        below-threshold hit, tagged via extra='{"below_threshold": true}'.
        That way every cluster gets some KEGG signal — strong calls are
        clearly distinguishable, weak guesses are present for hand-checking.

    The naive "all hits including below-threshold" approach swamps the DB
    with ~50 candidates per assigned hit. This strategy keeps the cluster
    table dense without drowning it in noise.
    """
    import json as _json

    @dataclass
    class _Hit:
        cluster: str
        ko: str
        score: float | None
        evalue: float | None
        definition: str
        is_assigned: bool

    # First pass: collect all parsed hits, separated by cluster.
    by_cluster: dict[str, dict[str, list[_Hit]]] = {}
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            is_assigned = line.startswith("*")
            payload = line[1:].strip() if is_assigned else line.strip()
            if '"' in payload:
                left, definition = payload.split('"', 1)
                definition = definition.rstrip('"').strip()
            else:
                left, definition = payload, ""
            parts = left.split()
            if len(parts) < 5:
                continue
            cluster_name, ko, _thrshld, score_str, evalue_str = parts[:5]
            try:
                score = float(score_str)
            except ValueError:
                score = None
            try:
                evalue = float(evalue_str)
            except ValueError:
                evalue = None
            bucket = by_cluster.setdefault(
                cluster_name, {"assigned": [], "candidate": []}
            )
            bucket["assigned" if is_assigned else "candidate"].append(
                _Hit(cluster_name, ko, score, evalue, definition, is_assigned)
            )

    # Second pass: yield assigned hits where present; else the single best candidate.
    for cluster, hits in by_cluster.items():
        if hits["assigned"]:
            for h in hits["assigned"]:
                # Always record evalue in extra so xlsx can show it on demand;
                # NO below_threshold flag, so cells render as bare KO codes.
                extra = _json.dumps({"evalue": h.evalue}) if h.evalue is not None else None
                yield ClusterAnnotationRow(
                    cluster_name=h.cluster, source="KEGG", code=h.ko,
                    category=None, name=None, description=h.definition,
                    score=h.score, extra=extra,
                )
        elif hits["candidate"]:
            best = max(
                hits["candidate"],
                key=lambda h: (h.score if h.score is not None else float("-inf")),
            )
            # Below-threshold flag tells the xlsx writer to format the cell as
            # 'KOxxxxx (e: 1.2e-50)' so biologists can see at a glance that
            # this was a sub-significant guess.
            yield ClusterAnnotationRow(
                cluster_name=best.cluster, source="KEGG", code=best.ko,
                category=None, name=None, description=best.definition,
                score=best.score,
                extra=_json.dumps({
                    "below_threshold": True,
                    "evalue": best.evalue,
                }),
            )


def parse_deepnog(path: Path) -> Iterable[ClusterAnnotationRow]:
    """Parse DeepNOG `infer` output.

    Format: 3-column header `sequence_id,prediction,confidence`. DeepNOG 1.2.x
    actually emits CSV (despite `.tsv` extension); older versions used TSV.
    We sniff the delimiter from the header to handle both.

    DeepNOG predicts orthology group IDs (COG#### for cog2020 db, eggNOG OG
    IDs like ENOG#### or COG#### for eggNOG5 db). Stored under source='COG'
    with category='group' to distinguish from eggNOG-mapper's single-letter
    COG categories.
    """
    with path.open() as f:
        first = f.readline()
        # Sniff delimiter from header
        delimiter = "," if first.count(",") > first.count("\t") else "\t"
        f.seek(0)
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            cluster_name = (row.get("sequence_id") or "").strip()
            prediction = (row.get("prediction") or "").strip()
            if not cluster_name or not prediction or prediction in {"-", "none"}:
                continue
            try:
                score = float(row.get("confidence", "") or "nan")
            except ValueError:
                score = None
            if score != score:  # NaN
                score = None
            yield ClusterAnnotationRow(
                cluster_name=cluster_name,
                source="COG",
                code=prediction,
                category="group",
                name=None,
                description=None,
                score=score,
                extra='{"tool": "deepnog"}',
            )


def parse_amrfinder(path: Path) -> Iterable[ClusterAnnotationRow]:
    """Parse NCBI AMRFinderPlus output (TSV).

    Column names vary noisily across AMRFinder versions. Some examples:
        v3.x:   'Protein identifier' / 'Gene symbol' / 'Sequence name' / 'Element type'
        v4.x:   'Protein id'         / 'Element symbol' / 'Element name' / 'Type'

    We look up by trying several aliases per logical field. The cluster_name
    always comes from the first column, whatever it's called.
    """
    import json as _json

    def _first(row: dict, *keys: str) -> str:
        for k in keys:
            v = row.get(k)
            if v is not None and v.strip():
                return v.strip()
        return ""

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            return
        id_col = reader.fieldnames[0]
        for row in reader:
            cluster_name = (row.get(id_col) or "").strip()
            if not cluster_name:
                continue
            gene = _first(row, "Element symbol", "Gene symbol")
            seq_name = _first(row, "Element name", "Sequence name")
            element_type = _first(row, "Type", "Element type")          # AMR/STRESS/VIRULENCE
            element_sub = _first(row, "Subtype", "Element subtype")
            drug_class = _first(row, "Class")
            drug_subcl = _first(row, "Subclass")
            pct_id_str = _first(row, "% Identity to reference",
                                "% Identity to reference sequence")
            try:
                pct_id = float(pct_id_str) if pct_id_str else None
            except ValueError:
                pct_id = None

            extras = {
                k: row[k].strip() for k in (
                    "Method", "Class", "Subclass",
                    "% Coverage of reference", "% Coverage of reference sequence",
                    "Closest reference accession", "Accession of closest reference",
                ) if row.get(k) and row[k].strip()
            }
            yield ClusterAnnotationRow(
                cluster_name=cluster_name,
                source="AMRFinder",
                code=gene or None,
                category=element_type or None,
                name=gene or None,
                description=seq_name or (f"{drug_class}/{drug_subcl}" if drug_class else None),
                score=pct_id,
                extra=_json.dumps(extras) if extras else None,
            )


def parse_defensefinder(path: Path) -> Iterable[ClusterAnnotationRow]:
    """Parse DefenseFinder defense_finder_genes.tsv.

    One row per gene-in-system. The gene-level granularity is what we want:
    each row maps directly to one cluster rep (hit_id).
    """
    import json as _json
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cluster_name = (row.get("hit_id") or "").strip()
            if not cluster_name:
                continue
            sys_type = (row.get("type") or "").strip()
            sys_sub = (row.get("subtype") or "").strip()
            gene = (row.get("gene_name") or "").strip()
            try:
                score = float(row.get("hit_score") or "nan")
            except ValueError:
                score = None
            if score != score:
                score = None
            extras = {
                k: row[k].strip() for k in ("sys_id", "hit_status", "hit_seq_cov")
                if row.get(k)
            }
            yield ClusterAnnotationRow(
                cluster_name=cluster_name,
                source="DefenseFinder",
                code=sys_type or None,
                category=sys_sub or None,
                name=gene or None,
                description=f"{sys_type} / {sys_sub}" if sys_type and sys_sub else sys_type or None,
                score=score,
                extra=_json.dumps(extras) if extras else None,
            )


def parse_padloc(path: Path) -> Iterable[ClusterAnnotationRow]:
    """Parse PADLOC CSV output (one row per protein-in-system)."""
    import json as _json
    with path.open() as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            cluster_name = (row.get("target.name") or "").strip()
            if not cluster_name:
                continue
            system = (row.get("system") or "").strip()
            protein = (row.get("protein.name") or "").strip()
            hmm_name = (row.get("hmm.name") or "").strip()
            try:
                score = float(row.get("full.seq.score") or "nan")
            except ValueError:
                score = None
            if score != score:
                score = None
            extras = {
                k: row[k].strip() for k in ("hmm.accession", "full.seq.E.value",
                                            "target.coverage", "hmm.coverage")
                if row.get(k)
            }
            yield ClusterAnnotationRow(
                cluster_name=cluster_name,
                source="PADLOC",
                code=system or None,
                category=protein or None,
                name=hmm_name or protein or None,
                description=(row.get("target.description") or "").strip() or None,
                score=score,
                extra=_json.dumps(extras) if extras else None,
            )


# InterProScan TSV column positions (it has no header).
_IPS_COLS = [
    "protein_accession", "md5", "length", "analysis", "signature_accession",
    "signature_description", "start", "stop", "evalue", "status", "date",
    "interpro_accession", "interpro_description", "go_annotations", "pathways",
]
_IPS_ANALYSES = {
    "Pfam", "TIGRFAM", "NCBIfam", "PANTHER", "SUPERFAMILY", "CDD", "SMART",
    "PRINTS", "ProSiteProfiles", "ProSitePatterns", "HAMAP", "Gene3D",
    "SFLD", "FunFam", "PIRSF", "PIRSR", "AntiFam", "MobiDBLite", "Coils",
    "TMHMM", "Phobius", "SignalP_GRAM_POSITIVE", "SignalP_GRAM_NEGATIVE",
    "SignalP_EUK",
}


def parse_interproscan(path: Path) -> Iterable[ClusterAnnotationRow]:
    """Parse InterProScan TSV output.

    TSV has no header. Column count is 11 (no lookup), 13 (+ InterPro), or
    15 (+ GO + pathways). One row per signature hit, one protein → many rows.
    Emits one annotation row per hit. GO annotations (column 14) are split
    out as separate source='GO' rows.
    """
    go_lookup = _load_go_aspects()
    import json as _json

    with path.open() as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            # Pad to 15 so all slots are indexable.
            while len(fields) < 15:
                fields.append("")
            row = dict(zip(_IPS_COLS, fields, strict=True))

            cluster_name = row["protein_accession"].strip()
            analysis = row["analysis"].strip()
            if not cluster_name or not analysis:
                continue
            if analysis not in _IPS_ANALYSES:
                # Unknown analysis column — probably malformed line; skip.
                continue

            sig_acc = row["signature_accession"].strip() or None
            sig_desc = row["signature_description"].strip() or None
            ipr_acc = row["interpro_accession"].strip() or None
            ipr_desc = row["interpro_description"].strip() or None
            try:
                evalue = float(row["evalue"])
            except ValueError:
                evalue = None

            extras: dict = {}
            if ipr_acc:
                extras["interpro_accession"] = ipr_acc
            if ipr_desc:
                extras["interpro_description"] = ipr_desc

            yield ClusterAnnotationRow(
                cluster_name=cluster_name,
                source=analysis,          # e.g. 'Pfam', 'TIGRFAM', 'PANTHER'
                code=sig_acc,
                category=None,
                name=sig_desc,
                description=ipr_desc,
                score=evalue,
                extra=_json.dumps(extras) if extras else None,
            )

            # GO annotations carried by InterProScan (column 14 in the 15-col format)
            go_field = row.get("go_annotations", "").strip()
            if go_field and go_field != "-":
                # Format like 'GO:0005515|...(...|GO:0003824)' — split loosely.
                for token in go_field.replace("|", ",").split(","):
                    tok = token.strip().split("(", 1)[0].strip()
                    if tok.startswith("GO:"):
                        aspect, go_name = go_lookup.get(tok, (None, None))
                        yield ClusterAnnotationRow(
                            cluster_name=cluster_name,
                            source="GO",
                            code=tok,
                            category=aspect,
                            name=go_name,
                            description=None,
                            score=None,
                            extra='{"via": "interproscan"}',
                        )


_PARSERS = {
    "emapper": parse_emapper,
    "kofamscan": parse_kofamscan,
    "deepnog": parse_deepnog,
    "amrfinder": parse_amrfinder,
    "defensefinder": parse_defensefinder,
    "padloc": parse_padloc,
    "interproscan": parse_interproscan,
}


def detect_format(path: Path) -> str | None:
    """Peek at the file and guess which annotator produced it.

    Returns the format key (one of `_PARSERS`), or None if unrecognized.
    """
    try:
        with path.open() as f:
            head = [line.rstrip("\n") for line in f.readlines()[:20]]
    except (OSError, UnicodeDecodeError):
        return None

    # First pass: header-bearing formats (most reliable)
    for line in head:
        if not line or line.startswith("##"):
            continue
        stripped = line.lstrip("#").strip()
        if stripped.startswith("query\t") or line.startswith("#query\t"):
            return "emapper"
        if "gene name" in line and "\tKO\t" in line:
            return "kofamscan"
        # DeepNOG: 1.2.x emits CSV ('sequence_id,prediction,confidence');
        # older versions used TSV. Accept both.
        if line.startswith(("sequence_id,", "sequence_id\t")) and "prediction" in line:
            return "deepnog"
        # AMRFinderPlus header: starts with a protein-id column ("Protein id"
        # in v4+, "Protein identifier" in older versions, sometimes "Name")
        # and contains element/symbol/scope vocabulary.
        if line.startswith(("Protein id\t", "Protein identifier\t",
                            "Protein Identifier\t", "Name\t")):
            if any(k in line for k in
                   ("Element symbol", "Element type", "Gene symbol", "\tType\t")):
                return "amrfinder"
        # DefenseFinder genes.tsv: header columns include 'hit_id' + 'gene_name' + 'sys_id'
        if "\thit_id\t" in line and "gene_name" in line and "sys_id" in line:
            return "defensefinder"
        # PADLOC CSV: comma-separated, first cells are 'system.number,seqid,system,target.name'
        if line.startswith("system.number,") or line.startswith("system.number,seqid"):
            return "padloc"

    # Second pass: InterProScan has no header. Look for a line whose 4th
    # tab-delimited field is one of the known analyses.
    for line in head:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) >= 11 and fields[3] in _IPS_ANALYSES:
            return "interproscan"

    return None


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────


_COG_CATEGORY_NAMES = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "Y": "Nuclear structure",
}


def _load_go_aspects() -> dict[str, tuple[str, str]]:
    """Load GO_ID → (aspect, name) from the shipped go_aspects.tsv.

    The shipped file covers the most common terms; unknowns come through as
    (None, None) and get stored without aspect. Users with a full go-basic.obo
    can regenerate this file and drop it into the package directory.
    """
    try:
        path = files("strainbench.producer") / "data" / "go_aspects.tsv"
        text = path.read_text()
    except (FileNotFoundError, ModuleNotFoundError):
        return {}
    lookup: dict[str, tuple[str, str]] = {}
    for line in text.splitlines()[1:]:  # skip header
        parts = line.split("\t")
        if len(parts) >= 3:
            lookup[parts[0]] = (parts[1], parts[2])
    return lookup


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
        raise ValueError("No active cluster_runs in DB")
    return int(row["cluster_run_id"])


def _load_cluster_name_map(
    conn: sqlite3.Connection, cluster_run_id: int
) -> dict[str, int]:
    return {
        row["cluster_name"]: int(row["cluster_id"])
        for row in conn.execute(
            "SELECT cluster_id, cluster_name FROM clusters WHERE cluster_run_id = ?",
            (cluster_run_id,),
        )
    }


@dataclass
class WorkdirImportSummary:
    """Aggregate outcome of importing all detected files in a workdir."""

    workdir: str
    per_file: list[tuple[str, str, ImportSummary | str]]
    # list of (relative_path, detected_format_or_reason, result_or_error_string)
    total_rows_inserted: int
    total_files_scanned: int
    total_files_imported: int
    total_files_skipped: int


def import_workdir(
    db_path: str | Path,
    workdir: str | Path,
    *,
    cluster_run_id: int | None = None,
    replace: bool = False,
) -> WorkdirImportSummary:
    """Scan a workdir for recognizable annotation outputs and import them all.

    Walks the directory recursively. Every non-tiny TSV / txt file is passed
    through `detect_format`; recognized formats are handed to the matching
    parser. Unrecognized files are reported in the summary but not imported.

    Zero-effect if run twice without `replace=True` — each (source, cluster)
    pair gets duplicate rows otherwise. Typically the caller wants
    `replace=True` for idempotent re-imports.
    """
    workdir = Path(workdir)
    if not workdir.is_dir():
        raise FileNotFoundError(f"workdir not found or not a directory: {workdir}")

    per_file: list[tuple[str, str, ImportSummary | str]] = []
    total_rows = 0
    total_imported = 0
    total_skipped = 0
    total_scanned = 0

    # Candidate files: TSV, CSV (PADLOC), TXT, and the specific .annotations
    # extension that eggNOG-mapper writes.
    candidates: list[Path] = []
    for pattern in ("*.tsv", "*.csv", "*.txt", "*.annotations", "*.emapper.annotations"):
        candidates.extend(workdir.rglob(pattern))
    # De-duplicate while preserving order.
    seen: set[Path] = set()
    ordered_candidates = []
    for p in sorted(candidates):
        if p.resolve() not in seen:
            seen.add(p.resolve())
            ordered_candidates.append(p)

    for path in ordered_candidates:
        total_scanned += 1
        rel = str(path.relative_to(workdir))
        fmt = detect_format(path)
        if fmt is None:
            per_file.append((rel, "unrecognized", "skipped: unknown format"))
            total_skipped += 1
            continue
        try:
            summary = import_annotation_file(
                db_path, path,
                source_format=fmt,
                cluster_run_id=cluster_run_id,
                replace=replace,
            )
        except (ValueError, FileNotFoundError, sqlite3.Error) as exc:
            per_file.append((rel, fmt, f"error: {exc}"))
            total_skipped += 1
            continue
        per_file.append((rel, fmt, summary))
        total_rows += summary.rows_inserted
        total_imported += 1

    return WorkdirImportSummary(
        workdir=str(workdir),
        per_file=per_file,
        total_rows_inserted=total_rows,
        total_files_scanned=total_scanned,
        total_files_imported=total_imported,
        total_files_skipped=total_skipped,
    )
