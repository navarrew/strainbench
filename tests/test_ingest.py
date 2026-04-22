"""Tests for the ingestion layer.

These exercise both the SQL side (rows land in the right tables with the
right shape and FK relationships hold) and the filesystem side (sidecar
FASTA files appear with the expected format).

Uses the synthetic_gbff fixture from conftest.py + the real gbff parser
to produce a StrainRecord, so we're testing the end-to-end parser->ingest
flow rather than ingestion in isolation.
"""

from __future__ import annotations

import sqlite3

import pytest

pytest.importorskip("Bio")

from strainbench.core import db as core_db
from strainbench.producer.ingest import (
    IngestError,
    ingest_records,
    ingest_strain,
)
from strainbench.producer.parsers.gbff import parse_gbff


@pytest.fixture
def fresh_db(tmp_path):
    """Yield a freshly-initialized DB path. Closed cleanly between tests."""
    db_path = tmp_path / "test.db"
    core_db.initialize(db_path)
    return db_path


@pytest.fixture
def fasta_dir(tmp_path):
    return tmp_path / "fasta"


def test_inserts_strain_and_cds_rows(fresh_db, fasta_dir, synthetic_gbff):
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(fresh_db) as conn:
        strain_id, status = ingest_strain(conn, record, fasta_dir)
    assert status == "inserted"
    assert isinstance(strain_id, int)

    with core_db.connect(fresh_db, read_only=True) as conn:
        strains = conn.execute("SELECT * FROM strains").fetchall()
        assert len(strains) == 1
        assert strains[0]["locus_prefix"] == "TESTPFX"
        assert strains[0]["cds_count"] == len(record.cds_records)
        assert strains[0]["source_format"] == "gbff"

        cds_rows = conn.execute(
            "SELECT * FROM cds WHERE strain_id = ?", (strain_id,)
        ).fetchall()
        assert len(cds_rows) == len(record.cds_records)
        assert {r["locus_tag"] for r in cds_rows} == {
            c.locus_tag for c in record.cds_records
        }


def test_writes_rna_sidecar_when_non_cds_features_present(fresh_db, fasta_dir, synthetic_gbff):
    """tRNA/rRNA/CRISPR features should land in <fasta_dir>/rna/<prefix>.fna —
    a parallel to the fna/ and faa/ sidecars but for non-CDS data."""
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(fresh_db) as conn:
        ingest_strain(conn, record, fasta_dir)

    rna_path = fasta_dir / "rna" / "TESTPFX.fna"
    assert rna_path.exists(), "rna/ sidecar should be written when non-CDS features exist"

    text = rna_path.read_text()
    # 3 records — tRNA + rRNA + CRISPR — should each have a header line
    assert text.count(">") == 3
    # Headers carry the feature_type tag so users can grep across files
    assert "[tRNA]" in text
    assert "[rRNA]" in text
    assert "[CRISPR]" in text
    # Product info preserved
    assert "16S ribosomal RNA" in text
    # Non-CRISPR repeat_region should NOT appear — tandem-only filter at parser level
    assert "tandem" not in text.lower()


def test_writes_sidecar_fasta_files(fresh_db, fasta_dir, synthetic_gbff):
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(fresh_db) as conn:
        ingest_strain(conn, record, fasta_dir)

    fna = fasta_dir / "fna" / "TESTPFX.fna"
    faa = fasta_dir / "faa" / "TESTPFX.faa"
    assert fna.exists() and faa.exists()

    fna_text = fna.read_text()
    faa_text = faa.read_text()
    # One '>' line per CDS in each file
    assert fna_text.count(">") == len(record.cds_records)
    assert faa_text.count(">") == len(record.cds_records)
    # Header format: >LOCUSPREFIX|LOCUS_TAG ...
    for cds in record.cds_records:
        assert f">TESTPFX|{cds.locus_tag} " in faa_text


def test_duplicate_raises_by_default(fresh_db, fasta_dir, synthetic_gbff):
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(fresh_db) as conn:
        ingest_strain(conn, record, fasta_dir)
        with pytest.raises(IngestError, match="already exists"):
            ingest_strain(conn, record, fasta_dir)


def test_duplicate_skip_returns_skipped(fresh_db, fasta_dir, synthetic_gbff):
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(fresh_db) as conn:
        first_id, first_status = ingest_strain(conn, record, fasta_dir)
        second_id, second_status = ingest_strain(
            conn, record, fasta_dir, on_duplicate="skip"
        )
    assert first_status == "inserted"
    assert second_status == "skipped"
    assert second_id is None


def test_duplicate_replace_drops_old_rows_and_reinserts(
    fresh_db, fasta_dir, synthetic_gbff
):
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(fresh_db) as conn:
        old_id, _ = ingest_strain(conn, record, fasta_dir)
        new_id, status = ingest_strain(
            conn, record, fasta_dir, on_duplicate="replace"
        )

    assert status == "replaced"
    # Note: SQLite may reuse the rowid after delete (no AUTOINCREMENT), so we
    # don't assert new_id != old_id — what matters is that the old row+cds
    # were dropped and the new strain owns all current cds rows.

    with core_db.connect(fresh_db, read_only=True) as conn:
        # Only one strain row exists
        n = conn.execute("SELECT COUNT(*) AS n FROM strains").fetchone()["n"]
        assert n == 1
        # All cds rows belong to the new strain_id (cascade dropped the old)
        rows = conn.execute(
            "SELECT DISTINCT strain_id FROM cds"
        ).fetchall()
        assert {r["strain_id"] for r in rows} == {new_id}


def test_failed_ingest_does_not_leave_partial_db_rows(fresh_db, fasta_dir, synthetic_gbff):
    """Transaction must roll back on error so we never see half-ingested strains."""
    record = parse_gbff(synthetic_gbff)
    # Inject a bogus locus_tag duplicate within the same record to violate
    # UNIQUE(strain_id, locus_tag) on the second CDS insert.
    record.cds_records[1].locus_tag = record.cds_records[0].locus_tag

    with core_db.connect(fresh_db) as conn:
        with pytest.raises(sqlite3.IntegrityError):
            ingest_strain(conn, record, fasta_dir)

        # Strain row must NOT exist — rollback should have undone it.
        n = conn.execute("SELECT COUNT(*) AS n FROM strains").fetchone()["n"]
        assert n == 0
        n = conn.execute("SELECT COUNT(*) AS n FROM cds").fetchone()["n"]
        assert n == 0


def test_ingest_records_continues_on_per_strain_failure(
    fresh_db, fasta_dir, synthetic_gbff
):
    """Batch ingest collects failures rather than aborting."""
    good = parse_gbff(synthetic_gbff)
    # Make a "bad" record that will fail on the second CDS insert
    bad = parse_gbff(synthetic_gbff)
    bad.locus_prefix = "BADPFX"
    bad.cds_records[1].locus_tag = bad.cds_records[0].locus_tag

    with core_db.connect(fresh_db) as conn:
        results = ingest_records(conn, [good, bad], fasta_dir)

    assert [r.status for r in results] == ["inserted", "failed"]
    assert results[1].error is not None
    # Good record made it in; bad one didn't pollute the DB
    with core_db.connect(fresh_db, read_only=True) as conn:
        names = {
            r["locus_prefix"]
            for r in conn.execute("SELECT locus_prefix FROM strains").fetchall()
        }
    assert names == {"TESTPFX"}


def test_cds_metadata_round_trips(fresh_db, fasta_dir, synthetic_gbff):
    """Verify a representative CDS preserves its full metadata through ingest."""
    record = parse_gbff(synthetic_gbff)
    expected = next(c for c in record.cds_records if c.locus_tag == "TESTPFX_RS00005")
    with core_db.connect(fresh_db) as conn:
        ingest_strain(conn, record, fasta_dir)

    with core_db.connect(fresh_db, read_only=True) as conn:
        row = conn.execute(
            "SELECT * FROM cds WHERE locus_tag = ?", (expected.locus_tag,)
        ).fetchone()
    assert row["protein_id"] == expected.protein_id
    assert row["nuc_accession"] == expected.nuc_accession
    assert row["direction"] == expected.direction
    assert row["gc_pct"] == expected.gc_pct
    assert row["aa_length"] == expected.aa_length
    assert row["nt_length"] == expected.nt_length
    assert row["gene_name"] == expected.gene_name
    assert row["annotation"] == expected.annotation
