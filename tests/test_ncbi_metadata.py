"""Tests for the NCBI assembly metadata-table loader."""

from __future__ import annotations

import pytest

from strainbench.producer.parsers.ncbi_metadata import load_assembly_table


def test_loads_rows_keyed_by_accession(tmp_path):
    p = tmp_path / "table.tab"
    p.write_text(
        "Accession\tSpecies\tStrain\tBioProject\tBioSample\tLevel\n"
        "GCF_000001.1\tEscherichia coli\tK-12\tPRJ1\tSAM1\tComplete Genome\n"
        "GCF_000002.1\tEscherichia coli\tBL21\tPRJ2\tSAM2\tContig\n"
    )
    table = load_assembly_table(p)
    assert set(table) == {"GCF_000001.1", "GCF_000002.1"}


def test_complete_genome_normalized_to_complete(tmp_path):
    p = tmp_path / "table.tab"
    p.write_text(
        "Accession\tSpecies\tStrain\tBioProject\tBioSample\tLevel\n"
        "GCF_X.1\tFoo\tBAR\tPRJ\tSAM\tComplete Genome\n"
    )
    assert load_assembly_table(p)["GCF_X.1"]["assembly_level"] == "Complete"


def test_unknown_level_falls_back(tmp_path):
    p = tmp_path / "table.tab"
    p.write_text(
        "Accession\tSpecies\tStrain\tBioProject\tBioSample\tLevel\n"
        "GCF_X.1\tFoo\tBAR\tPRJ\tSAM\tWeirdLevel\n"
    )
    assert load_assembly_table(p)["GCF_X.1"]["assembly_level"] == "Unknown"


def test_long_form_dataformat_headers_also_work(tmp_path):
    """Newer `dataformat tsv genome` emits longer human-readable headers
    by default. The parser must accept both header schemas — silently
    dropping data because of a header-name mismatch was a real bug."""
    p = tmp_path / "table.tab"
    p.write_text(
        "Assembly Accession\tANI Submitted species\tAssembly BioSample Strain\t"
        "Assembly BioProject Accession\tAssembly BioSample Accession\tAssembly Level\n"
        "GCF_049244215.1\tGardnerella sp. DNF01159\tDNF01159\tPRJNA639145\tSAMN15393823\tComplete Genome\n"
    )
    table = load_assembly_table(p)
    assert "GCF_049244215.1" in table
    row = table["GCF_049244215.1"]
    assert row["species"] == "Gardnerella sp. DNF01159"
    assert row["strain_name"] == "DNF01159"
    assert row["bioproject_id"] == "PRJNA639145"
    assert row["biosample_id"] == "SAMN15393823"
    assert row["assembly_level"] == "Complete"


def test_unrecognized_header_format_raises_clearly(tmp_path):
    p = tmp_path / "weird.tab"
    p.write_text(
        "AccNumber\tOrganism\tStrainName\n"
        "GCF_X.1\tFoo\tbar\n"
    )
    with pytest.raises(ValueError, match="Could not find an accession column"):
        load_assembly_table(p)


def test_empty_cells_become_none(tmp_path):
    p = tmp_path / "table.tab"
    p.write_text(
        "Accession\tSpecies\tStrain\tBioProject\tBioSample\tLevel\n"
        "GCF_X.1\tFoo\t\t\t\tContig\n"
    )
    row = load_assembly_table(p)["GCF_X.1"]
    assert row["strain_name"] is None
    assert row["bioproject_id"] is None
    assert row["biosample_id"] is None
