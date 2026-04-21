"""Tests for `strainbench status` and `strainbench annotate-prep`."""

from __future__ import annotations

import json

import pytest

pytest.importorskip("Bio")

from strainbench.core import db as core_db
from strainbench.producer.annotate_prep import prepare_annotation_workdir
from strainbench.producer.ingest import ingest_strain
from strainbench.producer.parsers.gbff import parse_gbff
from strainbench.producer.status import collect_status, format_status


def _seed_db_with_synthetic_strain(db_path, fasta_dir, synthetic_gbff):
    core_db.initialize(db_path)
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(db_path) as conn:
        ingest_strain(conn, record, fasta_dir)
    # Add a stub cluster_run and one cluster so we can test prep + status
    with core_db.connect(db_path) as conn:
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('stub', 'protein', 80, 80, 'X')"
        )
        run_id = conn.execute("SELECT MAX(cluster_run_id) FROM cluster_runs").fetchone()[0]
        conn.execute(
            "INSERT INTO clusters (cluster_run_id, cluster_name, representative_aa_seq) "
            "VALUES (?, ?, ?)",
            (run_id, "X_000001", "MKFNS"),
        )
        conn.commit()


def test_status_fresh_db_no_strains(tmp_path):
    db_path = tmp_path / "empty.db"
    core_db.initialize(db_path)
    status = collect_status(db_path)
    assert status.schema_version >= 4
    assert status.strains_count == 0
    assert status.cds_count == 0
    assert status.cluster_runs == []
    rendered = format_status(status)
    assert "No cluster_runs yet" in rendered


def test_status_with_clusters_and_no_annotations(tmp_path, synthetic_gbff):
    db_path = tmp_path / "t.db"
    fasta_dir = tmp_path / "fasta"
    _seed_db_with_synthetic_strain(db_path, fasta_dir, synthetic_gbff)

    status = collect_status(db_path)
    assert status.strains_count == 1
    assert status.cds_count >= 1
    assert len(status.cluster_runs) == 1
    run = status.cluster_runs[0]
    assert run.is_active is True
    assert run.annotation_summary == []
    rendered = format_status(status)
    assert "annotations: none" in rendered


def test_annotate_prep_creates_workdir_with_scripts(tmp_path, synthetic_gbff):
    db_path = tmp_path / "t.db"
    fasta_dir = tmp_path / "fasta"
    _seed_db_with_synthetic_strain(db_path, fasta_dir, synthetic_gbff)

    workdir = tmp_path / "prep"
    summary = prepare_annotation_workdir(db_path, workdir)

    # All expected files present
    expected_scripts = (
        "eggnog.sh", "kofamscan.sh", "deepnog.sh",
        "amrfinder.sh", "defensefinder.sh", "padloc.sh", "interproscan.sh",
    )
    assert (workdir / "reps.faa").exists()
    for name in expected_scripts:
        assert (workdir / name).exists()
    assert (workdir / "README.md").exists()
    assert (workdir / ".strainbench_workdir.json").exists()

    # Each script is executable, env-var-driven, and has no hard-coded user paths.
    expected_envvars = {
        "eggnog.sh": "EMAPPER_DATA_DIR",
        "kofamscan.sh": "KOFAM_PROFILES",
        "deepnog.sh": "DEEPNOG_DB",
        "amrfinder.sh": "AMRFINDER_THREADS",
        "defensefinder.sh": "DEFENSEFINDER_THREADS",
        "padloc.sh": "PADLOC_CPUS",
        "interproscan.sh": "INTERPROSCAN_BIN",
    }
    for name in expected_scripts:
        content = (workdir / name).read_text()
        assert content.startswith("#!/usr/bin/env bash")
        assert expected_envvars[name] in content
        assert "/Users/williamnavarre" not in content

    # Manifest has usable structure
    manifest = json.loads((workdir / ".strainbench_workdir.json").read_text())
    assert manifest["db_path"].endswith("t.db")
    assert manifest["cluster_run_id"] == summary.cluster_run_id
    assert manifest["reps_file"] == "reps.faa"
    assert len(manifest["scripts"]) == 7
    assert {s["import_format"] for s in manifest["scripts"]} == {
        "emapper", "kofamscan", "deepnog",
        "amrfinder", "defensefinder", "padloc", "interproscan",
    }


def test_annotate_prep_keeps_existing_reps_without_overwrite(tmp_path, synthetic_gbff):
    db_path = tmp_path / "t.db"
    fasta_dir = tmp_path / "fasta"
    _seed_db_with_synthetic_strain(db_path, fasta_dir, synthetic_gbff)
    workdir = tmp_path / "prep"
    prepare_annotation_workdir(db_path, workdir)
    reps = workdir / "reps.faa"
    reps.write_text("manually edited\n")
    # Re-running without overwrite should keep the manually-edited file.
    prepare_annotation_workdir(db_path, workdir)
    assert reps.read_text() == "manually edited\n"
    # But with overwrite=True the file gets regenerated.
    prepare_annotation_workdir(db_path, workdir, overwrite=True)
    assert ">X_000001" in reps.read_text()
