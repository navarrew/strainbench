"""Tests for the xlsx export.

We don't try to verify visual formatting (colors, frozen panes, etc.) — that
would require comparing against a snapshot file and is brittle. Instead we
verify the structural shape of the exported workbook by reading it back
with pandas:

    * Correct dimensions (clusters × (metadata cols + strain cols))
    * Metadata columns present in the right places
    * Strain columns named with the locus_prefix-prefixed full header
    * Cells filled with '[N]|...' for present hits and '*' for absent
"""

from __future__ import annotations

import shutil

import pytest

pytest.importorskip("Bio")
pytest.importorskip("pandas")
pytest.importorskip("openpyxl")

import pandas as pd

from strainbench.core import db as core_db
from strainbench.producer.cluster import cluster_strains
from strainbench.producer.export_xlsx import (
    export_cluster_table_xlsx,
    export_strain_anchored_xlsx,
)
from strainbench.producer.ingest import ingest_strain
from strainbench.producer.parsers.gbff import parse_gbff


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs2 not on PATH")
def test_export_produces_readable_xlsx(tmp_path, synthetic_gbff):
    """End-to-end: ingest synthetic gbff → cluster → export → read back."""
    db_path = tmp_path / "test.db"
    fasta_dir = tmp_path / "fasta"
    out_path = tmp_path / "cluster_table.xlsx"

    core_db.initialize(db_path)
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(db_path) as conn:
        ingest_strain(conn, record, fasta_dir)
        cluster_strains(conn, fasta_dir, sequence_type="protein",
                        pct_identity=80, coverage=80,
                        name_prefix="TEST", label="export-test")

    summary = export_cluster_table_xlsx(db_path, out_path)
    assert out_path.exists()
    assert summary.cluster_count >= 1
    assert summary.strain_count == 1

    # Read back. Header is row 0; pandas treats it as column names.
    df = pd.read_excel(out_path, sheet_name="cluster_table")
    expected_meta = [
        "original_order", "CLUSTER", "NCBI annotation", "GC%", "GC spread",
        "aa length", "flags", "total count", "strain count",
    ]
    assert list(df.columns[: len(expected_meta)]) == expected_meta

    # The remaining columns are strain headers — exactly summary.strain_count of them.
    strain_cols = df.columns[len(expected_meta):]
    assert len(strain_cols) == summary.strain_count
    # The synthetic strain's locus_prefix is 'TESTPFX'
    assert any(col.startswith("TESTPFX | ") for col in strain_cols)

    # Cluster rows match the cluster count
    assert len(df) == summary.cluster_count

    # Cell format: every value in strain columns is either '*' or starts with '['
    for col in strain_cols:
        for value in df[col]:
            assert value == "*" or str(value).startswith("[")


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs2 not on PATH")
def test_strain_anchored_export(tmp_path, synthetic_gbff):
    """End-to-end: ingest → cluster → export-strain-table → read back."""
    db_path = tmp_path / "test.db"
    fasta_dir = tmp_path / "fasta"
    out_path = tmp_path / "anchor.xlsx"

    core_db.initialize(db_path)
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(db_path) as conn:
        ingest_strain(conn, record, fasta_dir)
        cluster_strains(conn, fasta_dir, sequence_type="protein",
                        pct_identity=80, coverage=80,
                        name_prefix="T", label="anchor-test")

    summary = export_strain_anchored_xlsx(db_path, out_path, "TESTPFX")
    assert out_path.exists()
    # The synthetic gbff has 3 translatable CDSs → 3 rows.
    assert summary.cluster_count == 3
    assert summary.strain_count == 1

    df = pd.read_excel(out_path, sheet_name="TESTPFX_anchored")
    expected_meta = [
        "position", "anchor_locus_tag", "nuc_accession", "location", "gene",
        "CLUSTER", "NCBI annotation", "aa length", "cluster total", "cluster strains",
    ]
    assert list(df.columns[: len(expected_meta)]) == expected_meta
    # 3 rows, one per anchor CDS, in genomic order
    assert len(df) == 3
    assert list(df["position"]) == [1, 2, 3]
    # locus_tags are the synthetic ones, in cds_id order (= ingest order)
    assert all(tag.startswith("TESTPFX_RS") for tag in df["anchor_locus_tag"])

    # Anchor's column should be the FIRST strain column and start with the ★ marker
    strain_cols = list(df.columns[len(expected_meta):])
    assert strain_cols[0].startswith("★ TESTPFX | ")

    # Anchor cells should be '[1]|<that row's locus_tag>|...' (NOT the cluster summary)
    anchor_col = strain_cols[0]
    for _, row in df.iterrows():
        assert row[anchor_col].startswith(f"[1]|{row['anchor_locus_tag']}|")


def test_strain_anchored_nt_subcluster_column_present_when_nt_run_exists(tmp_path):
    """The nt_subcluster column should appear when a nucleotide cluster_run
    exists in the DB. We exercise the writer directly with hand-built data
    instead of running mmseqs nt clustering — the synthetic fixture has too
    few sequences for mmseqs's nucleotide linclust pre-step to succeed."""
    db_path = tmp_path / "test.db"
    out_path = tmp_path / "anchor.xlsx"
    core_db.initialize(db_path)

    with core_db.connect(db_path) as conn:
        # Minimal: one strain, two CDS, two cluster_runs (protein + nt),
        # one cluster per run, with both CDS as members.
        conn.execute(
            "INSERT INTO strains (locus_prefix, species, strain_name, assembly_id, "
            "biosample_id, bioproject_id, assembly_level, source_format, source_file) "
            "VALUES ('TESTP', 'Test bact', 'TESTSTRAIN', 'GCF_X', '?', '?', 'Complete', 'gbff', '?')"
        )
        sid = conn.execute("SELECT strain_id FROM strains").fetchone()["strain_id"]
        for i, lt in enumerate(("TESTP_001", "TESTP_002"), start=1):
            conn.execute(
                "INSERT INTO cds (strain_id, locus_tag, location, direction, "
                "nt_length, aa_length, gc_pct, annotation) "
                "VALUES (?, ?, ?, 'F', 300, 100, 35.0, 'test protein')",
                (sid, lt, f"{i*100}..{i*100+300}"),
            )
        # Protein run
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('p', 'protein', 80, 80, 'P')"
        )
        prot_run = conn.execute("SELECT MAX(cluster_run_id) FROM cluster_runs").fetchone()[0]
        conn.execute(
            "INSERT INTO clusters (cluster_run_id, cluster_name, member_count, strain_count) "
            "VALUES (?, 'P_000001', 2, 1)",
            (prot_run,),
        )
        prot_cid = conn.execute("SELECT cluster_id FROM clusters WHERE cluster_name='P_000001'").fetchone()[0]
        # Nucleotide run
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('n', 'nucleotide', 95, 90, 'N')"
        )
        nt_run = conn.execute("SELECT MAX(cluster_run_id) FROM cluster_runs").fetchone()[0]
        conn.execute(
            "INSERT INTO clusters (cluster_run_id, cluster_name, member_count, strain_count) "
            "VALUES (?, 'N_000001', 2, 1)",
            (nt_run,),
        )
        nt_cid = conn.execute("SELECT cluster_id FROM clusters WHERE cluster_name='N_000001'").fetchone()[0]
        for cid in conn.execute("SELECT cds_id FROM cds"):
            conn.execute(
                "INSERT INTO cluster_membership (cluster_run_id, cluster_id, cds_id) "
                "VALUES (?, ?, ?)",
                (prot_run, prot_cid, cid["cds_id"]),
            )
            conn.execute(
                "INSERT INTO cluster_membership (cluster_run_id, cluster_id, cds_id) "
                "VALUES (?, ?, ?)",
                (nt_run, nt_cid, cid["cds_id"]),
            )
        conn.commit()

    summary = export_strain_anchored_xlsx(db_path, out_path, "TESTP")
    assert summary.cluster_count == 2

    df = pd.read_excel(out_path, sheet_name="TESTP_anchored")
    expected_meta = [
        "position", "anchor_locus_tag", "nuc_accession", "location", "gene",
        "CLUSTER", "nt_subcluster", "NCBI annotation", "aa length",
        "cluster total", "cluster strains",
    ]
    assert list(df.columns[: len(expected_meta)]) == expected_meta
    assert all(v == "N_000001" for v in df["nt_subcluster"])


def test_strain_anchored_unknown_strain_raises(tmp_path):
    db_path = tmp_path / "empty.db"
    out_path = tmp_path / "out.xlsx"
    core_db.initialize(db_path)
    # Stub a cluster_run so we get past _resolve_cluster_run_id
    with core_db.connect(db_path) as conn:
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('stub', 'protein', 80, 80, 'X')"
        )
        conn.commit()
    with pytest.raises(ValueError, match="not found"):
        export_strain_anchored_xlsx(db_path, out_path, "NOPE")


def test_export_raises_when_no_cluster_run(tmp_path):
    """Trying to export from a DB that has no cluster_runs yet should fail clearly."""
    db_path = tmp_path / "empty.db"
    out_path = tmp_path / "out.xlsx"
    core_db.initialize(db_path)
    with pytest.raises(ValueError, match="No active cluster_runs"):
        export_cluster_table_xlsx(db_path, out_path)
