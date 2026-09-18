"""Tests for `strainbench export-cluster-region`.

Builds a real two-strain DB from synthetic gbff files, ingests both, fakes
a cluster grouping their forward-strand CDSs, then verifies that the
export pulls the correct flanking sequence in the right orientation.

Why this test uses real gbff files (not just SQL fixtures): the whole point
of export-cluster-region is to re-extract sequence with strand-aware flanking
from the source gbff, so a SQL-only fixture would mock away the thing we
actually care about.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("Bio")

from strainbench.core import db as core_db
from strainbench.producer.export_cluster_region import (
    ExportRegionSummary,
    export_cluster_region,
)
from strainbench.producer.ingest import ingest_strain
from strainbench.producer.parsers.gbff import parse_gbff


def _build_two_strains_with_cluster(tmp_path: Path):
    """Create db.db + two ingested strains + a fake cluster with their forward CDSs.

    Returns (db_path, gbff_root). gbff_root is the directory the
    strains.source_file relative paths resolve against.
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    gbff_root = tmp_path / "gbff_root"
    gbff_root.mkdir()
    db_path = tmp_path / "test.db"
    fasta_dir = tmp_path / "fasta"
    core_db.initialize(db_path)

    def make_gbff(filename: str, locus_prefix: str, cds_strand: int) -> Path:
        # 1 kb of contig with a CDS at coords 200..400 on the requested strand,
        # a uniquely-identifiable upstream region we can grep for, and a
        # uniquely-identifiable downstream region.
        upstream_marker = "TTTAAACCCGGG" * 8     # 96 bp marker — easy to find
        downstream_marker = "GGGCCCAAATTT" * 8
        cds = "ATG" + ("AAA" * 65) + "TAA"       # 201 bp CDS, MKK...K*
        # We want CDS at coords 200..401. Build:
        # padding to bring upstream marker to coords 104..200
        # then the CDS at 200..401
        # then downstream marker at 401..497
        # then padding to 1000.
        before = "N" * (200 - len(upstream_marker))   # 200 - 96 = 104 N's
        after_marker = "N" * (1000 - 401 - len(downstream_marker))
        seq_str = before + upstream_marker + cds + downstream_marker + after_marker
        assert len(seq_str) == 1000
        if cds_strand == -1:
            # Make this strain's CDS reverse-strand by reverse-complementing the
            # whole construct: the CDS now sits at the 'end' of the contig and
            # reads on the - strand at the new coordinates.
            seq_str = str(Seq(seq_str).reverse_complement())
            # New coords for the CDS: it was 200..401 on +, becomes
            # (1000-401)..(1000-200) = 599..800 on -.
            cds_loc = SimpleLocation(599, 800, strand=-1)
        else:
            cds_loc = SimpleLocation(200, 401, strand=1)

        seq = Seq(seq_str)
        features = [
            SeqFeature(
                location=SimpleLocation(0, 1000, strand=1),
                type="source",
                qualifiers={
                    "organism": ["Test bacterium"],
                    "strain": [locus_prefix],
                },
            ),
            SeqFeature(
                location=cds_loc,
                type="CDS",
                qualifiers={
                    "locus_tag": [f"{locus_prefix}_RS00010"],
                    "product": ["target gene"],
                    "protein_id": [f"WP_{locus_prefix}.1"],
                    "translation": ["MKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK"
                                    "KKKKKKKKKKKKKKKKKKKKKKKKK"],
                },
            ),
        ]
        record = SeqRecord(
            seq, id=f"{locus_prefix}_CTG", name=f"{locus_prefix}CTG",
            description=f"{locus_prefix} synthetic.",
            annotations={"molecule_type": "DNA", "organism": "Test bacterium",
                         "topology": "linear"},
            features=features,
        )
        path = gbff_root / filename
        SeqIO.write([record], path, "genbank")
        return path

    # FWD strain — CDS on the + strand, upstream marker BEFORE start in genome coords
    fwd_path = make_gbff("FWDSTRAIN.gbff", "FWDSTRAIN", cds_strand=+1)
    # REV strain — CDS on the - strand, upstream marker AFTER end in genome coords
    #              (reverse-complemented relative to the natural read direction)
    rev_path = make_gbff("REVSTRAIN.gbff", "REVSTRAIN", cds_strand=-1)

    # Ingest with relative source_file paths so gbff_root resolution mirrors prod.
    for path in (fwd_path, rev_path):
        record = parse_gbff(path)
        record.source_file = path.name   # store JUST the filename — relative
        with core_db.connect(db_path) as conn:
            ingest_strain(conn, record, fasta_dir)

    # Fake a protein cluster_run + cluster grouping the two CDSs.
    with core_db.connect(db_path) as conn:
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, "
            "coverage, name_prefix) "
            "VALUES ('t', 'protein', 80, 80, 'CL')"
        )
        run_id = conn.execute(
            "SELECT MAX(cluster_run_id) FROM cluster_runs"
        ).fetchone()[0]
        conn.execute(
            "INSERT INTO clusters (cluster_run_id, cluster_name, member_count, "
            "strain_count, consensus_annotation) "
            "VALUES (?, 'CL_001', 2, 2, 'target gene')",
            (run_id,),
        )
        cluster_id = conn.execute(
            "SELECT cluster_id FROM clusters WHERE cluster_name = 'CL_001'"
        ).fetchone()[0]
        for cds_id in conn.execute("SELECT cds_id FROM cds").fetchall():
            conn.execute(
                "INSERT INTO cluster_membership (cluster_run_id, cluster_id, cds_id) "
                "VALUES (?, ?, ?)",
                (run_id, cluster_id, cds_id["cds_id"]),
            )
        conn.commit()

    return db_path, gbff_root


def _read_fasta(path: Path) -> dict[str, str]:
    """Tiny FASTA reader for test assertions."""
    records: dict[str, str] = {}
    current = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            records[current] = ""
        elif current is not None:
            records[current] += line
    return records


def test_extracts_upstream_in_natural_orientation_for_both_strands(tmp_path):
    """The 96-bp upstream marker must appear at the 5' end of BOTH strands' output.

    This is the core correctness check: forward-strand and reverse-strand
    CDSs must both have the SAME upstream marker at the start of the
    extracted sequence (because both are written 5'→3' relative to the gene).
    """
    db_path, gbff_root = _build_two_strains_with_cluster(tmp_path)
    out = tmp_path / "region.fna"

    summary = export_cluster_region(
        db_path, out,
        cluster_name="CL_001",
        upstream_bp=96,    # exactly the marker length
        downstream_bp=0,
        gbff_root=gbff_root,
    )

    assert isinstance(summary, ExportRegionSummary)
    assert summary.n_records_written == 2
    assert summary.n_strains_missing_gbff == 0
    assert summary.n_members_missing_locus_tag == 0

    records = _read_fasta(out)
    assert set(records) == {"FWDSTRAIN|FWDSTRAIN_RS00010",
                             "REVSTRAIN|REVSTRAIN_RS00010"}
    upstream_marker = "TTTAAACCCGGG" * 8
    for name, seq in records.items():
        # 96 bp marker + 201 bp CDS = 297 bp total
        assert len(seq) == 297, f"{name}: expected 297 bp, got {len(seq)}"
        # marker first, then ATG start codon
        assert seq.startswith(upstream_marker), \
            f"{name}: upstream marker not at 5' end (got {seq[:30]}...)"
        assert seq[96:99] == "ATG", f"{name}: ATG not where expected"


def test_truncation_at_contig_boundary_is_reported(tmp_path):
    """Asking for more upstream than the contig has should still succeed but flag it."""
    db_path, gbff_root = _build_two_strains_with_cluster(tmp_path)
    out = tmp_path / "region.fna"

    # CDS starts at coord 200 → max possible upstream is 200 bp. Ask for 500.
    summary = export_cluster_region(
        db_path, out,
        cluster_name="CL_001",
        upstream_bp=500,
        gbff_root=gbff_root,
    )
    assert summary.n_records_written == 2
    assert summary.n_members_truncated_window == 2
    # And the headers should show the actual amount
    text = out.read_text()
    assert "[upstream=200/500bp]" in text


def test_gene_name_lookup_resolves_unique_cluster(tmp_path):
    """--gene-name should find the cluster when only one cluster matches."""
    db_path, gbff_root = _build_two_strains_with_cluster(tmp_path)
    # The fixture's CDSs have no gene_name, so set one for the test
    with core_db.connect(db_path) as conn:
        conn.execute("UPDATE cds SET gene_name = 'targ'")
        conn.commit()

    out = tmp_path / "region.fna"
    summary = export_cluster_region(
        db_path, out,
        gene_name="targ",
        upstream_bp=20,
        gbff_root=gbff_root,
    )
    assert summary.cluster_name == "CL_001"
    assert summary.n_records_written == 2


def test_missing_gbff_is_counted_not_fatal(tmp_path):
    """If a strain's source gbff is missing, we should report it but not crash."""
    db_path, gbff_root = _build_two_strains_with_cluster(tmp_path)
    # Delete one of the gbffs to simulate it being missing
    (gbff_root / "FWDSTRAIN.gbff").unlink()

    out = tmp_path / "region.fna"
    summary = export_cluster_region(
        db_path, out,
        cluster_name="CL_001",
        upstream_bp=10,
        gbff_root=gbff_root,
    )
    assert summary.n_strains_missing_gbff == 1
    assert summary.n_records_written == 1   # the surviving one still got written
    assert "FWDSTRAIN.gbff" in summary.missing_gbff_paths[0]


def test_ambiguous_lookup_lists_candidates(tmp_path):
    """--annotation-pattern matching multiple clusters should error with a list."""
    db_path, gbff_root = _build_two_strains_with_cluster(tmp_path)
    # Add a second cluster whose annotation also contains "target"
    with core_db.connect(db_path) as conn:
        run_id = conn.execute(
            "SELECT MAX(cluster_run_id) FROM cluster_runs"
        ).fetchone()[0]
        conn.execute(
            "INSERT INTO clusters (cluster_run_id, cluster_name, member_count, "
            "strain_count, consensus_annotation) "
            "VALUES (?, 'CL_002', 1, 1, 'other target gene variant')",
            (run_id,),
        )
        new_cluster_id = conn.execute(
            "SELECT cluster_id FROM clusters WHERE cluster_name = 'CL_002'"
        ).fetchone()[0]
        # Insert a 2nd cds for one of the strains in the 2nd cluster
        sid = conn.execute("SELECT strain_id FROM strains LIMIT 1").fetchone()[0]
        conn.execute(
            "INSERT INTO cds (strain_id, locus_tag, location, direction, "
            "annotation, aa_length) "
            "VALUES (?, 'X_RS99999', '1..30', 'F', 'target second copy', 10)",
            (sid,),
        )
        cds2 = conn.execute(
            "SELECT cds_id FROM cds WHERE locus_tag = 'X_RS99999'"
        ).fetchone()[0]
        conn.execute(
            "INSERT INTO cluster_membership (cluster_run_id, cluster_id, cds_id) "
            "VALUES (?, ?, ?)",
            (run_id, new_cluster_id, cds2),
        )
        conn.commit()

    out = tmp_path / "region.fna"
    with pytest.raises(ValueError, match="matches 2 clusters"):
        export_cluster_region(
            db_path, out,
            annotation_pattern="target",
            gbff_root=gbff_root,
        )
