"""Load NCBI assembly metadata tables for use as parser inputs.

NCBI's `dataformat tsv genome` command can emit one of two column-header
conventions depending on how the user invoked it:

  SHORT (strain-comp docs' convention, older tools):
    Accession   Species   Strain   BioProject   BioSample   Level

  LONG (newer dataformat default):
    Assembly Accession   ANI Submitted species   Assembly BioSample Strain
    Assembly BioProject Accession   Assembly BioSample Accession   Assembly Level

This parser accepts either — it looks up each logical field by trying the
known aliases in order. Users never have to reformat their metadata file.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

# NCBI uses 'Complete Genome' in some places and 'Complete' in others.
# We canonicalize on the shorter form, which is what our SQLite schema's
# CHECK constraint accepts.
_LEVEL_NORMALIZATION = {
    "Complete Genome": "Complete",
    "Complete": "Complete",
    "Chromosome": "Chromosome",
    "Scaffold": "Scaffold",
    "Contig": "Contig",
}

# Header aliases per logical field. Order matters — first match wins, so put
# the more-common / more-recently-seen variants near the front.
_COL_ALIASES: dict[str, list[str]] = {
    "accession":  ["Accession", "Assembly Accession"],
    "species":    ["Species", "ANI Submitted species",
                   "Organism Name", "Organism Scientific Name"],
    "strain":     ["Strain", "Assembly BioSample Strain",
                   "Organism Infraspecific Names Strain"],
    "bioproject": ["BioProject", "Assembly BioProject Accession"],
    "biosample":  ["BioSample", "Assembly BioSample Accession"],
    "level":      ["Level", "Assembly Level"],
}


def load_assembly_table(path: str | Path) -> dict[str, dict[str, Any]]:
    """Load an NCBI assembly metadata TSV, keyed by assembly accession.

    Returns:
        Dict mapping accession (e.g. 'GCF_022456355.1') → metadata dict
        with keys ready to pass as the `metadata` argument to `parse_gbff`:
        assembly_id, species, strain_name, biosample_id, bioproject_id,
        assembly_level.

        Empty cells become None. The 'Level' column is normalized to one of
        the schema's permitted values, falling back to 'Unknown' for
        unrecognized inputs.

    Raises:
        ValueError: if no recognized accession column can be found in the
            header — almost always means the TSV is in an unexpected format.
    """
    path = Path(path)
    table: dict[str, dict[str, Any]] = {}
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []
        # Map logical field → actual header name present in this file.
        resolved: dict[str, str | None] = {
            field: _pick_header(fieldnames, aliases)
            for field, aliases in _COL_ALIASES.items()
        }
        if resolved["accession"] is None:
            raise ValueError(
                f"Could not find an accession column in {path}. "
                f"Looked for any of {_COL_ALIASES['accession']}. "
                f"Got columns: {fieldnames}"
            )

        for row in reader:
            accession = (row.get(resolved["accession"]) or "").strip()
            if not accession:
                continue
            raw_level = (row.get(resolved["level"]) or "").strip() if resolved["level"] else ""
            table[accession] = {
                "assembly_id": accession,
                "species":       _get_if(row, resolved["species"]),
                "strain_name":   _get_if(row, resolved["strain"]),
                "bioproject_id": _get_if(row, resolved["bioproject"]),
                "biosample_id":  _get_if(row, resolved["biosample"]),
                "assembly_level": _LEVEL_NORMALIZATION.get(raw_level, "Unknown"),
            }
    return table


def _pick_header(fieldnames: list[str], aliases: list[str]) -> str | None:
    """Return the first alias that exists in `fieldnames`, or None if none do."""
    for alias in aliases:
        if alias in fieldnames:
            return alias
    return None


def _get_if(row: dict, col: str | None) -> str | None:
    if col is None:
        return None
    value = row.get(col)
    if value is None:
        return None
    stripped = value.strip()
    return stripped or None


def _strip_or_none(value: str | None) -> str | None:
    """Deprecated — kept for backward compatibility with any external callers."""
    if value is None:
        return None
    stripped = value.strip()
    return stripped or None
