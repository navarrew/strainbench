"""Tests for the NCBI assembly metadata-table loader."""

from __future__ import annotations

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
