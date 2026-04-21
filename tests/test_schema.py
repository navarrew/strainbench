"""Smoke tests for the core schema and DB initialization.

These tests don't exercise business logic — they verify that the schema
file is parseable by SQLite, that every expected table/view is created,
that foreign-key constraints are actually enforced at runtime, and that
the schema version is recorded.
"""

from __future__ import annotations

import sqlite3

import pytest

from strainbench.core import db as core_db


EXPECTED_TABLES = {
    "schema_version",
    "strains",
    "cds",
    "cluster_runs",
    "clusters",
    "cluster_membership",
    "annotations",
}
EXPECTED_VIEWS = {"v_cluster_cell"}


def test_initialize_creates_expected_tables(tmp_path):
    db_path = tmp_path / "test.db"
    core_db.initialize(db_path)

    with core_db.connect(db_path, read_only=True) as conn:
        tables = {
            row["name"]
            for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()
        }
        views = {
            row["name"]
            for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='view'"
            ).fetchall()
        }

    assert EXPECTED_TABLES.issubset(tables), f"Missing tables: {EXPECTED_TABLES - tables}"
    assert EXPECTED_VIEWS.issubset(views), f"Missing views: {EXPECTED_VIEWS - views}"


def test_schema_version_is_recorded(tmp_path):
    db_path = tmp_path / "test.db"
    core_db.initialize(db_path)
    # Latest version recorded equals the constant the module advertises.
    assert core_db.schema_version(db_path) == core_db.CURRENT_SCHEMA_VERSION


def test_initialize_is_idempotent(tmp_path):
    db_path = tmp_path / "test.db"
    core_db.initialize(db_path)
    rows_first = list(core_db.connect(db_path, read_only=True)
                      .execute("SELECT version FROM schema_version ORDER BY version"))
    core_db.initialize(db_path)  # should not duplicate any version rows
    rows_second = list(core_db.connect(db_path, read_only=True)
                       .execute("SELECT version FROM schema_version ORDER BY version"))
    assert rows_first == rows_second


def test_foreign_keys_are_enforced(tmp_path):
    db_path = tmp_path / "test.db"
    core_db.initialize(db_path)

    # Inserting a CDS row with a strain_id that doesn't exist must fail.
    with core_db.connect(db_path) as conn, pytest.raises(sqlite3.IntegrityError):
        conn.execute(
            "INSERT INTO cds (strain_id, locus_tag) VALUES (?, ?)",
            (999, "FAKE_0001"),
        )


def test_cascade_delete_removes_cds(tmp_path):
    db_path = tmp_path / "test.db"
    core_db.initialize(db_path)

    with core_db.connect(db_path) as conn:
        conn.execute(
            "INSERT INTO strains (locus_prefix, species, assembly_level, source_format) "
            "VALUES (?, ?, ?, ?)",
            ("TESTPFX", "Escherichia coli", "Complete", "gbff"),
        )
        strain_id = conn.execute(
            "SELECT strain_id FROM strains WHERE locus_prefix = ?", ("TESTPFX",)
        ).fetchone()["strain_id"]
        conn.execute(
            "INSERT INTO cds (strain_id, locus_tag) VALUES (?, ?)",
            (strain_id, "TESTPFX_0001"),
        )
        conn.commit()

        # Deleting the strain should cascade to its CDS rows.
        conn.execute("DELETE FROM strains WHERE strain_id = ?", (strain_id,))
        conn.commit()

        remaining = conn.execute(
            "SELECT COUNT(*) AS n FROM cds WHERE strain_id = ?", (strain_id,)
        ).fetchone()["n"]
        assert remaining == 0


def test_unique_locus_prefix_is_enforced(tmp_path):
    db_path = tmp_path / "test.db"
    core_db.initialize(db_path)

    with core_db.connect(db_path) as conn:
        conn.execute(
            "INSERT INTO strains (locus_prefix, source_format, assembly_level) "
            "VALUES (?, ?, ?)",
            ("DUPEPFX", "gbff", "Complete"),
        )
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute(
                "INSERT INTO strains (locus_prefix, source_format, assembly_level) "
                "VALUES (?, ?, ?)",
                ("DUPEPFX", "gbff", "Contig"),
            )
