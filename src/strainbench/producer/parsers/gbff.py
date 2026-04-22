"""Parse NCBI GenBank flat files (.gbff / .gbk) into StrainRecord objects.

The parser is intentionally tolerant: a single .gbff may hold one record
(complete chromosome) or hundreds (a draft assembly's contigs). All records
are scanned for CDS features and merged into one StrainRecord.

Pseudogenes are skipped — defined as a CDS feature that either lacks a
/translation qualifier, has /pseudo or /pseudogene set, or whose translation
is shorter than 30 amino acids. This matches the behavior of strain-comp's
`pseudoremover()`.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation
from Bio.SeqUtils import gc_fraction

from strainbench.core.models import CDSRecord, NonCDSRecord, StrainRecord

_VALID_ASSEMBLY_LEVELS = {"Complete", "Chromosome", "Scaffold", "Contig", "Unknown"}
_MIN_AA_LENGTH = 30

# Non-CDS feature types we extract verbatim into the rna/ sidecar.
_NON_CDS_RNA_TYPES = {"tRNA", "rRNA", "ncRNA", "tmRNA"}


def parse_gbff(
    path: str | Path,
    *,
    metadata: dict[str, Any] | None = None,
) -> StrainRecord:
    """Parse one .gbff file and return a populated StrainRecord.

    Args:
        path: Path to the .gbff file.
        metadata: Optional dict supplied by the caller (typically from
            `ncbi_strain_assembly_table.tab` via load_assembly_table()).
            Recognized keys: assembly_id, species, strain_name, biosample_id,
            bioproject_id, assembly_level, locus_prefix. Any missing key
            falls back to a best-effort value extracted from the gbff
            (organism annotation, source-feature /strain) or 'Unknown'.

    Raises:
        ValueError: if the file holds no GenBank records, or no translatable
            CDS features.
    """
    path = Path(path)
    metadata = dict(metadata or {})

    records = list(SeqIO.parse(path, "genbank"))
    if not records:
        raise ValueError(f"No GenBank records found in {path}")

    cds_records: list[CDSRecord] = []
    non_cds_records: list[NonCDSRecord] = []
    non_cds_index = 0   # for synthesizing locus_tags when /locus_tag is missing
    for record in records:
        for feature in record.features:
            if feature.type == "CDS":
                cds = _build_cds_record(feature, record)
                if cds is not None:
                    cds_records.append(cds)
            else:
                non_cds = _build_non_cds_record(feature, record, non_cds_index)
                if non_cds is not None:
                    non_cds_records.append(non_cds)
                    non_cds_index += 1

    if not cds_records:
        raise ValueError(f"No translatable CDS features found in {path}")

    first = records[0]
    locus_prefix = metadata.get("locus_prefix") or _derive_locus_prefix(cds_records)
    species = metadata.get("species") or first.annotations.get("organism", "Unknown")
    strain_name = (
        metadata.get("strain_name") or _extract_strain_from_source(first) or "Unknown"
    )
    assembly_id = metadata.get("assembly_id") or first.id
    assembly_level = metadata.get("assembly_level") or "Unknown"
    if assembly_level not in _VALID_ASSEMBLY_LEVELS:
        assembly_level = "Unknown"

    return StrainRecord(
        locus_prefix=locus_prefix,
        species=species,
        strain_name=strain_name,
        assembly_id=assembly_id,
        biosample_id=metadata.get("biosample_id"),
        bioproject_id=metadata.get("bioproject_id"),
        assembly_level=assembly_level,
        source_format="gbff",
        source_file=str(path),
        cds_records=cds_records,
        non_cds_records=non_cds_records,
    )


def _build_cds_record(feature, record) -> CDSRecord | None:
    """Build a CDSRecord from one CDS feature, or None to skip it."""
    qualifiers = feature.qualifiers

    if "translation" not in qualifiers:
        return None
    if "pseudo" in qualifiers or "pseudogene" in qualifiers:
        return None

    aa_sequence = qualifiers["translation"][0]
    if len(aa_sequence) < _MIN_AA_LENGTH:
        return None

    nt_sequence = str(feature.extract(record.seq))

    locus_tag = qualifiers.get("locus_tag", [None])[0]
    if not locus_tag:
        # Synthesize a stable identifier so the row can still be ingested
        # without violating the (strain_id, locus_tag) UNIQUE constraint.
        locus_tag = f"{record.id}_{int(feature.location.start)}_{int(feature.location.end)}"

    return CDSRecord(
        locus_tag=locus_tag,
        protein_id=qualifiers.get("protein_id", [None])[0],
        nuc_accession=record.id,
        location=_format_location(feature.location),
        direction="R" if feature.location.strand == -1 else "F",
        nt_sequence=nt_sequence,
        aa_sequence=aa_sequence,
        nt_length=len(nt_sequence),
        aa_length=len(aa_sequence),
        gc_pct=round(100 * gc_fraction(nt_sequence), 2),
        annotation=qualifiers.get("product", ["hypothetical protein"])[0],
        gene_name=qualifiers.get("gene", [None])[0],
        notes=_describe_special_cases(feature),
    )


def _build_non_cds_record(feature, record, index: int) -> NonCDSRecord | None:
    """Build a NonCDSRecord from a tRNA/rRNA/ncRNA/tmRNA or CRISPR repeat_region.

    Returns None for any feature type we don't capture (CDS, gene, source,
    misc_feature, gap, regulatory, …) and for repeat_region features that
    aren't CRISPR-related (tandem repeats, etc.).

    NCBI's CRISPR annotations are unreliable — we capture them anyway as a
    convenience for browsing, but biologists doing real CRISPR analysis
    should use dedicated tools (CRISPRCasFinder, MinCED, etc.) rather than
    trusting these.
    """
    qualifiers = feature.qualifiers
    ftype = feature.type

    if ftype in _NON_CDS_RNA_TYPES:
        product = qualifiers.get("product", [""])[0].strip() or ftype
        normalized_type = ftype
    elif ftype == "repeat_region":
        # Only keep CRISPR-flavored repeat_regions; skip tandem repeats etc.
        markers = " ".join(
            qualifiers.get("rpt_family", []) + qualifiers.get("note", [])
        ).upper()
        if "CRISPR" not in markers:
            return None
        normalized_type = "CRISPR"
        # Use the most informative free-text we have
        product = (qualifiers.get("rpt_family", [""])[0]
                   or qualifiers.get("note", [""])[0]
                   or "CRISPR repeat region").strip()
    else:
        return None

    # Synthesize locus_tag when missing (CRISPR repeats often have none).
    locus_tag = qualifiers.get("locus_tag", [None])[0]
    if not locus_tag:
        locus_tag = f"{record.id}_{normalized_type}_{index:04d}"

    try:
        nt_sequence = str(feature.extract(record.seq))
    except Exception:  # noqa: BLE001 - origin-spanning joins, etc.
        return None

    return NonCDSRecord(
        feature_type=normalized_type,
        locus_tag=locus_tag,
        nuc_accession=record.id,
        location=_format_location(feature.location),
        direction="R" if feature.location.strand == -1 else "F",
        nt_sequence=nt_sequence,
        nt_length=len(nt_sequence),
        product=product,
        notes=_describe_special_cases(feature),
    )


def _format_location(loc) -> str:
    """Render a BioPython location in GenBank-style text.

    BioPython locations are 0-based half-open; GenBank text is 1-based
    inclusive. So start_text = start + 1, end_text = end (the BioPython
    end already equals the GenBank inclusive end).

    Examples:
        SimpleLocation(155, 547, strand=1)            -> '156..547'
        SimpleLocation(155, 547, strand=-1)           -> 'complement(156..547)'
        SimpleLocation(BeforePosition(155), 547, +1)  -> '<156..547'
        SimpleLocation(155, AfterPosition(547), +1)   -> '156..>547'
        Compound (join of two parts)                  -> 'join(156..200,300..547)'
    """
    if isinstance(loc, CompoundLocation):
        return "join(" + ",".join(_format_location(p) for p in loc.parts) + ")"

    start_prefix = "<" if isinstance(loc.start, BeforePosition) else ""
    end_prefix = ">" if isinstance(loc.end, AfterPosition) else ""
    span = f"{start_prefix}{int(loc.start) + 1}..{end_prefix}{int(loc.end)}"
    return f"complement({span})" if loc.strand == -1 else span


def _describe_special_cases(feature) -> str | None:
    """Return a human note for truncated/joined/origin-crossing CDSs, else None."""
    loc = feature.location
    notes: list[str] = []

    fuzzy_start = isinstance(loc.start, BeforePosition)
    fuzzy_end = isinstance(loc.end, AfterPosition)
    is_reverse = loc.strand == -1

    if fuzzy_start and fuzzy_end:
        notes.append("truncated at both ends")
    elif fuzzy_start:
        notes.append("truncated at 3' end" if is_reverse else "truncated at 5' end")
    elif fuzzy_end:
        notes.append("truncated at 5' end" if is_reverse else "truncated at 3' end")

    if isinstance(loc, CompoundLocation):
        # Heuristic: if any non-first part starts at coordinate 0, the join
        # likely wraps the origin of a circular replicon.
        non_first_starts = [int(p.start) for p in loc.parts[1:]]
        if 0 in non_first_starts:
            notes.append("may cross origin of circular DNA")
        else:
            notes.append("compound location (possible programmed frameshift)")

    return "; ".join(notes) if notes else None


def _derive_locus_prefix(cds_records: list[CDSRecord]) -> str:
    """Extract the strain's locus prefix from the first CDS's locus_tag.

    Two formats are common in NCBI annotations:
        'AB0Q63_RS00005' -> 'AB0Q63'   (split on underscore)
        'STM0042'        -> 'STM'       (letters before digits)
    Mirrors strain-comp's get_locus_tag() heuristic for compatibility with
    existing pipelines.
    """
    first_tag = cds_records[0].locus_tag
    if "_" in first_tag:
        return first_tag.split("_")[0]
    match = re.match(r"\D+", first_tag)
    return match.group(0) if match else first_tag


def _extract_strain_from_source(record) -> str | None:
    """Pull a strain identifier from the source feature.

    NCBI submitters are inconsistent about where the strain identifier lives.
    The conventional qualifier is /strain, but a meaningful fraction of
    submissions put it in /isolate instead (especially for clinical isolates
    and metagenome-derived assemblies). We accept either, preferring /strain.
    """
    for feature in record.features:
        if feature.type == "source":
            qualifiers = feature.qualifiers
            for key in ("strain", "isolate"):
                values = qualifiers.get(key)
                if values and values[0].strip():
                    return values[0]
    return None
