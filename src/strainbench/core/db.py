"""Database connection and initialization helpers.

Both the producer and the web consumer use these — the producer opens
read-write connections, the web consumer opens read-only ones.

`initialize()` is idempotent: it creates the DB if missing AND brings an
existing older DB up to the current schema by applying ALTER TABLE
migrations. Operations that touch a DB should call `initialize()` defensively
to avoid surprising the user with stale-schema errors.
"""

from __future__ import annotations

import sqlite3
from importlib.resources import files
from pathlib import Path

CURRENT_SCHEMA_VERSION = 4

# One entry per schema bump beyond v1. Each statement runs in order; each must
# be idempotent (we swallow "duplicate column" / "already exists" errors so
# re-running against a newer DB is safe).
_MIGRATIONS: dict[int, list[str]] = {
    2: [
        "ALTER TABLE strains  ADD COLUMN display_order INTEGER",
        "ALTER TABLE clusters ADD COLUMN display_order INTEGER",
        "INSERT OR IGNORE INTO schema_version (version, description) "
        "VALUES (2, 'add display_order columns for hierarchical clustering')",
    ],
    3: [
        # cluster_annotations table is declared in schema.sql with CREATE TABLE
        # IF NOT EXISTS, so running schema.sql against an old DB is enough;
        # this migration just records the version bump.
        "INSERT OR IGNORE INTO schema_version (version, description) "
        "VALUES (3, 'add cluster_annotations table for cluster-level functional annotations')",
    ],
    4: [
        "ALTER TABLE cluster_annotations ADD COLUMN source_file TEXT",
        "ALTER TABLE cluster_annotations ADD COLUMN tool_format TEXT",
        "INSERT OR IGNORE INTO schema_version (version, description) "
        "VALUES (4, 'provenance columns on cluster_annotations (source_file, tool_format)')",
    ],
}


def connect(db_path: str | Path, *, read_only: bool = False) -> sqlite3.Connection:
    """Open a SQLite connection with foreign keys enabled and row access by name.

    Read-only mode uses SQLite's URI syntax so the DB file is opened with
    mode=ro — safe for concurrent web serving even while the producer is
    writing to a different copy of the file.
    """
    db_path = Path(db_path)
    if read_only:
        uri = f"file:{db_path.resolve()}?mode=ro"
        conn = sqlite3.connect(uri, uri=True)
    else:
        conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def initialize(db_path: str | Path) -> None:
    """Create or upgrade a strainbench database at the given path.

    Idempotent on both fronts: a fresh DB gets the full schema; an existing
    older DB is brought up to CURRENT_SCHEMA_VERSION via the migrations
    listed in `_MIGRATIONS` (idempotent ALTER TABLE statements).
    """
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    schema_sql = (files("strainbench.core") / "schema.sql").read_text()
    with connect(db_path) as conn:
        conn.executescript(schema_sql)
        _apply_migrations(conn)
        conn.commit()


def _apply_migrations(conn: sqlite3.Connection) -> None:
    """Apply any DDL needed to bring an existing DB up to CURRENT_SCHEMA_VERSION.

    Each statement is run individually with "duplicate column" / "already
    exists" errors swallowed, so this is safe to run on a freshly-created DB
    (where schema.sql already declared the new columns) or an old one.
    """
    for version in sorted(_MIGRATIONS):
        if version > CURRENT_SCHEMA_VERSION:
            continue
        for stmt in _MIGRATIONS[version]:
            try:
                conn.execute(stmt)
            except sqlite3.OperationalError as exc:
                msg = str(exc).lower()
                if "duplicate column" in msg or "already exists" in msg:
                    continue
                raise


def schema_version(db_path: str | Path) -> int | None:
    """Return the highest applied schema version in this database, or None if empty."""
    with connect(db_path, read_only=True) as conn:
        row = conn.execute("SELECT MAX(version) AS v FROM schema_version").fetchone()
        return row["v"] if row and row["v"] is not None else None


# ─────────────────────────────────────────────────────────────────────────────
# Cluster-run resolution — single canonical implementation used everywhere.
# Lives here (not in a producer module) because every producer subcommand
# needs to pick a cluster_run, and previously this logic was copy-pasted in
# four separate files. Each divergence cost at least one user-visible bug.
# ─────────────────────────────────────────────────────────────────────────────


def resolve_cluster_run_id(
    conn: sqlite3.Connection,
    requested: int | None,
    *,
    sequence_type_hint: str | None = None,
) -> int:
    """Pick a cluster_run_id for any operation that targets one.

    Resolution order:
      1. Explicit `requested` (validated to exist) → used as-is.
      2. `sequence_type_hint` ('protein' or 'nucleotide') → most recent
         active run of that type. Raises if no such run exists.
      3. Default → most recent active protein run. Falls back to most
         recent active run of any type if no protein run exists.

    The protein default keeps the cluster_table xlsx and strain-anchored
    views pointed at the canonical "gene families" view; nucleotide runs
    are an overlay (see the `nt_subcluster` column). The hint mode lets
    `--sequence-type nucleotide` find an nt run automatically without
    needing the user to look up its cluster_run_id.
    """
    if requested is not None:
        row = conn.execute(
            "SELECT cluster_run_id FROM cluster_runs WHERE cluster_run_id = ?",
            (requested,),
        ).fetchone()
        if row is None:
            raise ValueError(f"cluster_run_id={requested} not found in DB")
        return int(requested)

    if sequence_type_hint in {"protein", "nucleotide"}:
        row = conn.execute(
            "SELECT cluster_run_id FROM cluster_runs "
            "WHERE is_active = 1 AND sequence_type = ? "
            "ORDER BY cluster_run_id DESC LIMIT 1",
            (sequence_type_hint,),
        ).fetchone()
        if row is not None:
            return int(row["cluster_run_id"])
        raise ValueError(
            f"no active cluster_run of sequence_type={sequence_type_hint!r}. "
            f"Run `strainbench cluster --sequence-type {sequence_type_hint} ...` first."
        )

    row = conn.execute(
        """
        SELECT cluster_run_id FROM cluster_runs
        WHERE is_active = 1
        ORDER BY (sequence_type = 'protein') DESC, cluster_run_id DESC
        LIMIT 1
        """
    ).fetchone()
    if row is None:
        raise ValueError(
            "No active cluster_runs in DB. Run `strainbench cluster` first."
        )
    return int(row["cluster_run_id"])


def resolve_nt_cluster_run_id(
    conn: sqlite3.Connection,
    requested: int | None,
) -> int | None:
    """Pick a nucleotide cluster_run to overlay (or None if none exists).

    Used by the strain-anchored xlsx export to find the nt sub-cluster
    information to overlay alongside the canonical protein view. Returns
    None instead of raising when no nt run exists, because the overlay
    is optional — the strain-anchored view still works without it.

    If `requested` is given, validates and returns it (any cluster_run_id,
    not just nt). If unspecified, returns the most recent active nt run,
    or None if there are none.
    """
    if requested is not None:
        return resolve_cluster_run_id(conn, requested)
    row = conn.execute(
        """
        SELECT cluster_run_id FROM cluster_runs
        WHERE is_active = 1 AND sequence_type = 'nucleotide'
        ORDER BY cluster_run_id DESC LIMIT 1
        """
    ).fetchone()
    return int(row["cluster_run_id"]) if row else None
