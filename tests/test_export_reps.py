"""Tests for `strainbench export-reps`.

The default stored-rep path is exercised indirectly by other tests. These
tests focus on `full_length_only=True`, which:
  - must read from sidecar FASTAs (needs fasta_dir)
  - must skip clusters where every member is truncated
  - must pick an un-truncated member where one exists, even if the stored
    representative would have been truncated
"""

from __future__ import annotations

import sqlite3

import pytest

pytest.importorskip("Bio")

from strainbench.core import db as core_db
from strainbench.producer.export_reps import export_cluster_representatives


def _build_db_and_sidecar(tmp_path):
    """Build a minimal DB with three clusters and matching sidecar FASTAs.

    Cluster A: rep un-truncated + second member un-truncated → use either.
    Cluster B: rep is truncated, another member un-truncated → must swap to the
               un-truncated one when full_length_only=True.
    Cluster C: every member truncated → must be skipped entirely.
    """
    db_path = tmp_path / "t.db"
    fasta_dir = tmp_path / "fasta"
    (fasta_dir / "faa").mkdir(parents=True)
    (fasta_dir / "fna").mkdir(parents=True)
    core_db.initialize(db_path)

    # Sidecar FASTAs. Record IDs follow strainbench's convention: <prefix>|<locus_tag>
    aa_records = {
        "A_001": "MAAAAAAA",   "A_002": "MAAAAAAB",
        "B_001": "MBBBBBBB",   "B_002": "MBBBBBBC",
        "C_001": "MCCCCCCC",   "C_002": "MCCCCCCD",
    }
    with (fasta_dir / "faa" / "PFX.faa").open("w") as f:
        for tag, seq in aa_records.items():
            f.write(f">PFX|{tag} fixture\n{seq}\n")

    with core_db.connect(db_path) as conn:
        conn.execute(
            "INSERT INTO strains (locus_prefix, source_format, assembly_level) "
            "VALUES ('PFX', 'gbff', 'Complete')"
        )
        sid = conn.execute("SELECT strain_id FROM strains").fetchone()["strain_id"]

        # CDSs — notes are the key: NULL = full-length, '...truncated...' = cut off.
        cds_specs = [
            ("A_001", None),
            ("A_002", None),
            ("B_001", "truncated at 5' end"),   # rep for cluster B — will be "wrong" pick
            ("B_002", None),
            ("C_001", "truncated at both ends"),
            ("C_002", "truncated at 3' end"),
        ]
        for tag, notes in cds_specs:
            conn.execute(
                "INSERT INTO cds (strain_id, locus_tag, location, direction, "
                "aa_length, notes) VALUES (?, ?, '1..24', 'F', 8, ?)",
                (sid, tag, notes),
            )

        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('t', 'protein', 80, 80, 'CL')"
        )
        run_id = conn.execute("SELECT MAX(cluster_run_id) FROM cluster_runs").fetchone()[0]

        cluster_map = {
            "CL_A": ["A_001", "A_002"],
            "CL_B": ["B_001", "B_002"],
            "CL_C": ["C_001", "C_002"],
        }
        for name, members in cluster_map.items():
            # Stored rep is always the FIRST member listed — which for CL_B
            # is a truncated one (hence the point of this test).
            rep_tag = members[0]
            rep_aa = aa_records[rep_tag]
            rep_cds_id = conn.execute(
                "SELECT cds_id FROM cds WHERE locus_tag = ?", (rep_tag,)
            ).fetchone()[0]
            conn.execute(
                "INSERT INTO clusters (cluster_run_id, cluster_name, "
                "representative_cds_id, representative_aa_seq, member_count, strain_count) "
                "VALUES (?, ?, ?, ?, ?, 1)",
                (run_id, name, rep_cds_id, rep_aa, len(members)),
            )
            cluster_id = conn.execute(
                "SELECT cluster_id FROM clusters WHERE cluster_name = ?", (name,)
            ).fetchone()[0]
            for tag in members:
                cds_id = conn.execute(
                    "SELECT cds_id FROM cds WHERE locus_tag = ?", (tag,)
                ).fetchone()[0]
                conn.execute(
                    "INSERT INTO cluster_membership (cluster_run_id, cluster_id, cds_id) "
                    "VALUES (?, ?, ?)",
                    (run_id, cluster_id, cds_id),
                )
        conn.commit()

    return db_path, fasta_dir


def _read_fasta(path):
    records = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = ""
            elif current is not None:
                records[current] += line
    return records


def test_default_writes_stored_rep_for_every_cluster(tmp_path):
    db_path, _ = _build_db_and_sidecar(tmp_path)
    out = tmp_path / "reps.faa"

    summary = export_cluster_representatives(db_path, out)
    records = _read_fasta(out)

    # All 3 clusters present — default path doesn't filter.
    assert summary.n_reps_written == 3
    assert set(records) == {"CL_A", "CL_B", "CL_C"}
    # Cluster B's rep was the truncated B_001 (stored as-is in default mode).
    assert records["CL_B"] == "MBBBBBBB"


def test_full_length_only_swaps_truncated_rep_and_skips_all_truncated(tmp_path):
    db_path, fasta_dir = _build_db_and_sidecar(tmp_path)
    out = tmp_path / "reps.faa"

    summary = export_cluster_representatives(
        db_path, out,
        full_length_only=True,
        fasta_dir=fasta_dir,
    )
    records = _read_fasta(out)

    # Only CL_A and CL_B present; CL_C fully-truncated → skipped.
    assert summary.n_reps_written == 2
    assert summary.n_clusters_all_truncated == 1
    assert set(records) == {"CL_A", "CL_B"}

    # For CL_B, the rep was truncated B_001 — full_length_only must have
    # swapped in B_002 (the only un-truncated member).
    assert records["CL_B"] == "MBBBBBBC"


def test_full_length_only_requires_fasta_dir(tmp_path):
    db_path, _ = _build_db_and_sidecar(tmp_path)
    out = tmp_path / "reps.faa"
    with pytest.raises(ValueError, match="fasta_dir is required"):
        export_cluster_representatives(db_path, out, full_length_only=True)


def test_full_length_only_raises_when_sidecar_subdir_missing(tmp_path):
    db_path, _ = _build_db_and_sidecar(tmp_path)
    empty_dir = tmp_path / "no_sidecar"
    empty_dir.mkdir()
    out = tmp_path / "reps.faa"
    with pytest.raises(ValueError, match="sidecar dir not found"):
        export_cluster_representatives(
            db_path, out, full_length_only=True, fasta_dir=empty_dir,
        )
