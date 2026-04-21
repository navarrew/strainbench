"""Load NCBI assembly metadata tables for use as parser inputs.

The strain-comp pipeline produces (or consumes) a tab-delimited file of the
form:

    Accession   Species   Strain   BioProject   BioSample   Level
    GCF_…       …         …        …            …           Complete Genome / Scaffold / Contig

The same table can be passed to `parse_gbff(..., metadata=row)` so the
parser doesn't have to guess these values from the gbff itself.
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


def load_assembly_table(path: str | Path) -> dict[str, dict[str, Any]]:
    """Load an `ncbi_strain_assembly_table.tab` file keyed by accession.

    Returns:
        Dict mapping accession (e.g. 'GCF_022456355.1') -> metadata dict
        with keys ready to be passed as the `metadata` argument of
        `parse_gbff`: assembly_id, species, strain_name, biosample_id,
        bioproject_id, assembly_level.

        Empty cells in the source TSV become None. The 'Level' column is
        normalized to one of the schema's permitted values, falling back
        to 'Unknown' for unrecognized inputs.
    """
    path = Path(path)
    table: dict[str, dict[str, Any]] = {}
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            accession = (row.get("Accession") or "").strip()
            if not accession:
                continue
            raw_level = (row.get("Level") or "").strip()
            table[accession] = {
                "assembly_id": accession,
                "species": _strip_or_none(row.get("Species")),
                "strain_name": _strip_or_none(row.get("Strain")),
                "bioproject_id": _strip_or_none(row.get("BioProject")),
                "biosample_id": _strip_or_none(row.get("BioSample")),
                "assembly_level": _LEVEL_NORMALIZATION.get(raw_level, "Unknown"),
            }
    return table


def _strip_or_none(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    return stripped or None
