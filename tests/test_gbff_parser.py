"""Tests for the gbff parser.

These run against a small programmatically-built fixture rather than a real
genome. The fixture covers: a normal forward CDS, a reverse-strand CDS, a
5'-truncated CDS, and a pseudogene that should be filtered out.
"""

from __future__ import annotations

import pytest

# All tests need BioPython via the parser; skip cleanly if it isn't installed.
pytest.importorskip("Bio")

from strainbench.core.models import StrainRecord
from strainbench.producer.parsers.gbff import parse_gbff


def test_parse_returns_strain_record(synthetic_gbff):
    rec = parse_gbff(synthetic_gbff)
    assert isinstance(rec, StrainRecord)
    assert rec.source_format == "gbff"
    assert rec.source_file == str(synthetic_gbff)


def test_pseudogenes_are_skipped(synthetic_gbff):
    # 4 CDS features in the fixture; one is a pseudogene → 3 expected.
    rec = parse_gbff(synthetic_gbff)
    assert len(rec.cds_records) == 3


def test_locus_prefix_derived_from_first_cds(synthetic_gbff):
    rec = parse_gbff(synthetic_gbff)
    assert rec.locus_prefix == "TESTPFX"
    assert all(c.locus_tag.startswith("TESTPFX") for c in rec.cds_records)


def test_directions_and_locations(synthetic_gbff):
    rec = parse_gbff(synthetic_gbff)
    by_tag = {c.locus_tag: c for c in rec.cds_records}

    fwd = by_tag["TESTPFX_RS00005"]
    assert fwd.direction == "F"
    assert fwd.location == "100..300"
    assert fwd.gene_name == "dnaA"
    assert fwd.protein_id == "WP_000000001.1"

    rev = by_tag["TESTPFX_RS00010"]
    assert rev.direction == "R"
    assert rev.location == "complement(500..700)"


def test_truncation_note_and_fuzzy_position(synthetic_gbff):
    rec = parse_gbff(synthetic_gbff)
    truncated = next(c for c in rec.cds_records if c.locus_tag == "TESTPFX_RS00015")
    assert truncated.notes is not None
    assert "truncated" in truncated.notes
    assert truncated.location.startswith("<")


def test_metadata_overrides_inferred_values(synthetic_gbff):
    rec = parse_gbff(
        synthetic_gbff,
        metadata={
            "assembly_id": "GCF_999999999.1",
            "species": "Override species",
            "strain_name": "OVERRIDE",
            "biosample_id": "SAMN12345678",
            "bioproject_id": "PRJNA123456",
            "assembly_level": "Complete",
        },
    )
    assert rec.assembly_id == "GCF_999999999.1"
    assert rec.species == "Override species"
    assert rec.strain_name == "OVERRIDE"
    assert rec.biosample_id == "SAMN12345678"
    assert rec.bioproject_id == "PRJNA123456"
    assert rec.assembly_level == "Complete"


def test_inferred_values_when_no_metadata(synthetic_gbff):
    rec = parse_gbff(synthetic_gbff)
    assert rec.species == "Test bacterium"
    assert rec.strain_name == "TESTSTRAIN"
    assert rec.assembly_level == "Unknown"


def test_unknown_assembly_level_normalized(synthetic_gbff):
    rec = parse_gbff(synthetic_gbff, metadata={"assembly_level": "Banana"})
    assert rec.assembly_level == "Unknown"


def test_gc_and_lengths_are_populated(synthetic_gbff):
    # The nt-vs-aa length relationship isn't checked here because the synthetic
    # fixture supplies translations independently of the nt sequence. Real gbff
    # files (where /translation decodes the nt) are exercised by the hand
    # validation against ~/Documents/iners/data/input/gbff via the CLI.
    rec = parse_gbff(synthetic_gbff)
    for cds in rec.cds_records:
        assert cds.aa_length > 0
        assert cds.nt_length > 0
        assert 0 <= cds.gc_pct <= 100


def test_isolate_qualifier_used_when_strain_absent(tmp_path):
    """Some NCBI submissions put the strain ID in /isolate instead of /strain.

    The parser must fall back to /isolate so we don't lose the identifier.
    Mirrors the case of GCF_978044465.1 in the iners dataset.
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    seq = Seq("ATG" + "ATCG" * 100 + "TAA")
    features = [
        SeqFeature(
            location=SimpleLocation(0, len(seq), strand=1),
            type="source",
            qualifiers={
                "organism": ["Lactobacillus iners"],
                # Note: NO /strain qualifier; only /isolate.
                "isolate": ["AMBV-2211"],
            },
        ),
        SeqFeature(
            location=SimpleLocation(2, 200, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": ["TESTPFX_RS00005"],
                "product": ["hypothetical protein"],
                "translation": ["MFKAINKGTSL" * 6],
            },
        ),
    ]
    record = SeqRecord(
        seq, id="TEST.1", name="TEST",
        description="Test isolate-only fallback.",
        annotations={"molecule_type": "DNA", "organism": "Lactobacillus iners",
                     "topology": "linear"},
        features=features,
    )
    gbff_path = tmp_path / "isolate_only.gbff"
    SeqIO.write([record], gbff_path, "genbank")

    rec = parse_gbff(gbff_path)
    assert rec.strain_name == "AMBV-2211"


def test_empty_gbff_raises(tmp_path):
    empty = tmp_path / "empty.gbff"
    empty.write_text("")
    with pytest.raises(ValueError, match="No GenBank records"):
        parse_gbff(empty)
