"""Tests for the clustering layer.

Two flavors:

  * Unit tests for the parser and DB-loader helpers — fast, no mmseqs needed.
  * One integration test that actually invokes mmseqs2; auto-skipped when
    mmseqs is not on PATH.
"""

from __future__ import annotations

import shutil
import sqlite3
import textwrap

import pytest

pytest.importorskip("Bio")

from strainbench.core import db as core_db
from strainbench.producer.cluster import (
    ClusterError,
    cluster_strains,
    parse_clustered_fasta,
)
from strainbench.producer.ingest import ingest_strain
from strainbench.producer.parsers.gbff import parse_gbff


# ─────────────────────────────────────────────────────────────────────────────
# Parser tests — exercise the duplicated-header cluster boundary logic.
# ─────────────────────────────────────────────────────────────────────────────


def test_parser_handles_two_clusters_unwrapped(tmp_path):
    sample = textwrap.dedent("""\
        >ALPHA|ALPHA_001 desc A1
        >ALPHA|ALPHA_001 desc A1
        MFKAINKGTSL
        >ALPHA|ALPHA_002 desc A2
        MVHCVPYGSI
        >BETA|BETA_001 desc B1
        >BETA|BETA_001 desc B1
        MNQRPAEFWNRALQ
    """)
    path = tmp_path / "out.fasta"
    path.write_text(sample)

    clusters = parse_clustered_fasta(path)
    assert len(clusters) == 2

    # Cluster sorted largest first; first one has 2 members, second has 1.
    assert len(clusters[0]) == 2
    assert len(clusters[1]) == 1
    assert clusters[0][0][0].startswith("ALPHA|")
    assert clusters[0][0][1] == "MFKAINKGTSL"
    assert clusters[0][1][1] == "MVHCVPYGSI"
    assert clusters[1][0][0].startswith("BETA|")


def test_parser_handles_wrapped_sequences(tmp_path):
    """If mmseqs ever outputs wrapped sequences, they should still parse."""
    sample = textwrap.dedent("""\
        >X|X_001 desc
        >X|X_001 desc
        MFKAINKGTSL
        MVHCVPYGSI
        MNQRPAEFWN
        >X|X_002 desc
        MAAAA
        MBBBB
    """)
    path = tmp_path / "out.fasta"
    path.write_text(sample)

    clusters = parse_clustered_fasta(path)
    assert len(clusters) == 1
    assert clusters[0][0][1] == "MFKAINKGTSLMVHCVPYGSIMNQRPAEFWN"
    assert clusters[0][1][1] == "MAAAAMBBBB"


def test_parser_clusters_sorted_largest_first(tmp_path):
    sample = textwrap.dedent("""\
        >SMALL|S_001 small
        >SMALL|S_001 small
        SEQ_S
        >BIG|B_001 big
        >BIG|B_001 big
        SEQ_B1
        >BIG|B_002 big
        SEQ_B2
        >BIG|B_003 big
        SEQ_B3
    """)
    path = tmp_path / "out.fasta"
    path.write_text(sample)

    clusters = parse_clustered_fasta(path)
    assert [len(c) for c in clusters] == [3, 1]


# ─────────────────────────────────────────────────────────────────────────────
# Integration test — auto-skipped if mmseqs2 not on PATH.
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def populated_db_with_fasta(tmp_path, synthetic_gbff):
    """Initialize a DB, ingest the synthetic gbff, return (db_path, fasta_dir)."""
    db_path = tmp_path / "test.db"
    fasta_dir = tmp_path / "fasta"
    core_db.initialize(db_path)
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(db_path) as conn:
        ingest_strain(conn, record, fasta_dir)
    return db_path, fasta_dir


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs2 not on PATH")
def test_cluster_strains_end_to_end(populated_db_with_fasta):
    db_path, fasta_dir = populated_db_with_fasta

    with core_db.connect(db_path) as conn:
        result = cluster_strains(
            conn,
            fasta_dir,
            sequence_type="protein",
            pct_identity=80,
            coverage=80,
            name_prefix="TEST",
            label="unit-test-run",
        )

    assert result.cluster_count >= 1
    assert result.member_count >= 1

    with core_db.connect(db_path, read_only=True) as conn:
        runs = conn.execute("SELECT * FROM cluster_runs").fetchall()
        assert len(runs) == 1
        assert runs[0]["label"] == "unit-test-run"
        assert runs[0]["sequence_type"] == "protein"

        clusters = conn.execute("SELECT * FROM clusters").fetchall()
        # Each cluster should have aggregates populated and a name like TEST_000001
        for c in clusters:
            assert c["cluster_name"].startswith("TEST_")
            assert c["member_count"] is not None and c["member_count"] >= 1
            assert c["strain_count"] is not None and c["strain_count"] >= 1
            assert c["representative_aa_seq"] is not None

        memberships = conn.execute("SELECT COUNT(*) AS n FROM cluster_membership").fetchone()
        assert memberships["n"] == result.member_count

        # v_cluster_cell view should now return rows.
        cells = conn.execute("SELECT * FROM v_cluster_cell").fetchall()
        assert len(cells) >= 1


# ─────────────────────────────────────────────────────────────────────────────
# Negative test — empty fasta dir should raise.
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs2 not on PATH")
def test_empty_fasta_dir_raises(tmp_path):
    db_path = tmp_path / "empty.db"
    core_db.initialize(db_path)
    fasta_dir = tmp_path / "nothing"
    fasta_dir.mkdir()
    with core_db.connect(db_path) as conn, pytest.raises(ClusterError, match="No sidecar"):
        cluster_strains(conn, fasta_dir)
