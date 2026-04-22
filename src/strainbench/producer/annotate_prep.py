"""Generate a self-contained workdir for external annotation.

The workdir holds:
    reps.faa                      — cluster rep sequences (input for every tool)
    eggnog.sh / kofamscan.sh / deepnog.sh
                                  — portable shell scripts (env-var-driven)
    README_workdir.md             — human instructions for the workflow
    .strainbench_workdir.json     — machine-readable manifest (db path, date,
                                    cluster_run used, expected outputs)

The scripts are *templates but runnable as-is* provided your databases are in
the default locations (~/databases/<tool>/) or you've set the relevant
env vars. This means: one `annotate-prep` call, then just `bash eggnog.sh`
(with the right conda env active) — no editing.
"""

from __future__ import annotations

import json
import os
import stat
from dataclasses import dataclass
from datetime import datetime
from importlib.resources import files
from pathlib import Path

from strainbench.producer.export_reps import export_cluster_representatives


@dataclass
class AnnotatePrepSummary:
    workdir: str
    cluster_run_id: int
    reps_written: int
    scripts: list[str]


SCRIPTS = (
    ("eggnog.sh",       "eggnog_mapper", "emapper_out",       "emapper"),
    ("kofamscan.sh",    "kofamscan",     "kofamscan_out",     "kofamscan"),
    ("deepnog.sh",      "deepnog",       "deepnog_out",       "deepnog"),
    ("amrfinder.sh",    "amrfinder",     "amrfinder_out",     "amrfinder"),
    ("defensefinder.sh","defensefinder", "defensefinder_out", "defensefinder"),
    ("padloc.sh",       "padloc",        "padloc_out",        "padloc"),
    ("interproscan.sh", "interproscan",  "interproscan_out",  "interproscan"),
)


def prepare_annotation_workdir(
    db_path: str | Path,
    workdir: str | Path,
    *,
    cluster_run_id: int | None = None,
    overwrite: bool = False,
) -> AnnotatePrepSummary:
    """Build an annotation-ready workdir for the given DB's most recent cluster run."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    reps_path = workdir / "reps.faa"
    if reps_path.exists() and not overwrite:
        # Keep the existing reps.faa — user may have partial outputs referencing it.
        # Still refresh the shell templates + manifest (cheap, no data loss).
        pass
    else:
        export_cluster_representatives(
            db_path, reps_path,
            cluster_run_id=cluster_run_id,
            sequence_type="protein",
        )

    # Re-resolve the cluster_run_id for the manifest. Uses the canonical
    # protein-preferring resolver (before this consolidation, this was an
    # inline query with the old "most-recent-any" semantics — a footgun
    # waiting to happen once users had protein + nucleotide runs).
    import sqlite3
    from strainbench.core.db import resolve_cluster_run_id
    with sqlite3.connect(db_path) as conn:
        conn.row_factory = sqlite3.Row
        cluster_run_id = resolve_cluster_run_id(conn, cluster_run_id)
        reps_written = conn.execute(
            "SELECT COUNT(*) AS n FROM clusters "
            "WHERE cluster_run_id = ? AND representative_aa_seq IS NOT NULL",
            (cluster_run_id,),
        ).fetchone()["n"]

    # Copy shell templates from the package data dir.
    templates_dir = files("strainbench.producer") / "templates"
    written_scripts: list[str] = []
    for script_name, _tool, _output_dir, _fmt in SCRIPTS:
        src_text = (templates_dir / script_name).read_text()
        dest = workdir / script_name
        dest.write_text(src_text)
        dest.chmod(dest.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        written_scripts.append(script_name)

    readme_src = (templates_dir / "README_workdir.md").read_text()
    (workdir / "README.md").write_text(readme_src)

    # Manifest for machine-readable lookups ('what DB does this belong to?').
    manifest = {
        "strainbench_workdir_version": 1,
        "db_path": str(Path(db_path).resolve()),
        "cluster_run_id": cluster_run_id,
        "reps_written": reps_written,
        "reps_file": "reps.faa",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "scripts": [
            {
                "script": s,
                "tool": tool,
                "output_dir": out_dir,
                "import_format": fmt,
            }
            for s, tool, out_dir, fmt in SCRIPTS
        ],
    }
    (workdir / ".strainbench_workdir.json").write_text(json.dumps(manifest, indent=2) + "\n")

    return AnnotatePrepSummary(
        workdir=str(workdir),
        cluster_run_id=cluster_run_id,
        reps_written=reps_written,
        scripts=written_scripts,
    )
