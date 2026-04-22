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


def test_resolve_cluster_run_id_prefers_protein(tmp_path):
    """The canonical resolver in core.db picks protein runs by default when
    multiple active cluster_runs exist. Used by every export subcommand."""
    from strainbench.core.db import resolve_cluster_run_id

    db_path = tmp_path / "t.db"
    core_db.initialize(db_path)
    with core_db.connect(db_path) as conn:
        # Insert runs in reverse preference order to exercise the ORDER BY:
        # a newer nt run (id=1) and an older protein run (id=2 we'll fake
        # by setting it last). To make the protein run NEWER we add it second.
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('nt', 'nucleotide', 95, 90, 'NT')"
        )
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('prot', 'protein', 80, 90, 'P')"
        )
        conn.commit()

        # Default (no hint): should pick the protein run even though nt was created first
        assert resolve_cluster_run_id(conn, None) == 2

        # sequence_type hint 'nucleotide' → finds nt run
        assert resolve_cluster_run_id(conn, None, sequence_type_hint="nucleotide") == 1

        # Explicit requested run is used verbatim
        assert resolve_cluster_run_id(conn, 1) == 1
        assert resolve_cluster_run_id(conn, 2) == 2

        # Unknown run raises
        import pytest as _pytest
        with _pytest.raises(ValueError, match="not found"):
            resolve_cluster_run_id(conn, 99)


def test_resolve_nt_cluster_run_id_returns_none_when_no_nt_run(tmp_path):
    """Unlike the main resolver, the nt-specific one returns None rather than
    raising — it's used as an optional overlay for the strain-anchored xlsx."""
    from strainbench.core.db import resolve_nt_cluster_run_id

    db_path = tmp_path / "t.db"
    core_db.initialize(db_path)
    with core_db.connect(db_path) as conn:
        # No runs at all → None
        assert resolve_nt_cluster_run_id(conn, None) is None

        # Add a protein-only run — still no nt, still None
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('p', 'protein', 80, 90, 'P')"
        )
        conn.commit()
        assert resolve_nt_cluster_run_id(conn, None) is None

        # Add an nt run — now returned
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('n', 'nucleotide', 95, 90, 'N')"
        )
        conn.commit()
        assert resolve_nt_cluster_run_id(conn, None) == 2


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
