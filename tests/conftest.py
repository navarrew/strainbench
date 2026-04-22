"""Shared pytest fixtures.

`synthetic_gbff` builds a small but format-valid GenBank file containing
assorted CDS shapes (forward, reverse, truncated, pseudogene) so parser
tests don't need to depend on real ~3 MB genome files.
"""

from __future__ import annotations

from pathlib import Path

import pytest


def _build_synthetic_gbff(path: Path) -> None:
    # Imported lazily so tests that don't need BioPython still collect.
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import BeforePosition, SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    seq = Seq("ATGAAA" + "ATCG" * 500 + "TAA")  # ~2 kb of synthetic sequence

    features = [
        SeqFeature(
            location=SimpleLocation(0, len(seq), strand=1),
            type="source",
            qualifiers={
                "organism": ["Test bacterium"],
                "strain": ["TESTSTRAIN"],
            },
        ),
        # Forward CDS, normal
        SeqFeature(
            location=SimpleLocation(99, 300, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": ["TESTPFX_RS00005"],
                "gene": ["dnaA"],
                "product": ["chromosomal replication initiator protein DnaA"],
                "protein_id": ["WP_000000001.1"],
                "translation": ["MFKAINKGTSL" * 6],  # 66 aa, > 30 cutoff
            },
        ),
        # Reverse-strand CDS
        SeqFeature(
            location=SimpleLocation(499, 700, strand=-1),
            type="CDS",
            qualifiers={
                "locus_tag": ["TESTPFX_RS00010"],
                "product": ["hypothetical protein"],
                "protein_id": ["WP_000000002.1"],
                "translation": ["MVHCVPYGSI" * 6],
            },
        ),
        # 5'-truncated CDS (BeforePosition start)
        SeqFeature(
            location=SimpleLocation(BeforePosition(799), 1000, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": ["TESTPFX_RS00015"],
                "product": ["truncated transposase"],
                "protein_id": ["WP_000000003.1"],
                "translation": ["MNQRPAEFWNRALQ" * 4],
            },
        ),
        # Pseudogene — should be skipped (no /translation, /pseudo set)
        SeqFeature(
            location=SimpleLocation(1099, 1300, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": ["TESTPFX_RS00020"],
                "product": ["pseudogene fragment"],
                "pseudo": [""],
            },
        ),
        # tRNA — should land in non_cds_records / rna sidecar
        SeqFeature(
            location=SimpleLocation(1400, 1474, strand=1),
            type="tRNA",
            qualifiers={
                "locus_tag": ["TESTPFX_RS00025"],
                "product": ["tRNA-Ala(GGC)"],
            },
        ),
        # rRNA on the reverse strand — common pattern for 16S in many bacteria
        SeqFeature(
            location=SimpleLocation(1500, 1700, strand=-1),
            type="rRNA",
            qualifiers={
                "locus_tag": ["TESTPFX_RS00030"],
                "product": ["16S ribosomal RNA"],
            },
        ),
        # CRISPR repeat_region — captured even without /locus_tag (synthesized)
        SeqFeature(
            location=SimpleLocation(1750, 1900, strand=1),
            type="repeat_region",
            qualifiers={
                "rpt_family": ["CRISPR"],
                "note": ["CRISPR array consisting of 5 repeat units"],
            },
        ),
        # NON-CRISPR repeat_region — must NOT be captured
        SeqFeature(
            location=SimpleLocation(1920, 1960, strand=1),
            type="repeat_region",
            qualifiers={
                "rpt_type": ["tandem"],
                "note": ["short tandem repeat"],
            },
        ),
    ]

    record = SeqRecord(
        seq,
        id="TEST_ACC.1",
        name="TESTACC",
        description="Test bacterium strain TESTSTRAIN chromosome.",
        annotations={
            "molecule_type": "DNA",
            "organism": "Test bacterium",
            "topology": "linear",
            "taxonomy": ["Bacteria", "Firmicutes"],
        },
        features=features,
    )
    SeqIO.write([record], path, "genbank")


@pytest.fixture(scope="session")
def synthetic_gbff(tmp_path_factory) -> Path:
    """Return the path to a session-scoped synthetic .gbff file."""
    path = tmp_path_factory.mktemp("fixtures") / "synthetic.gbff"
    _build_synthetic_gbff(path)
    return path
