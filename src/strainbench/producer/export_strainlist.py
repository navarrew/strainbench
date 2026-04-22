"""Export a strainlist.txt file in dendrogram order.

Format mirrors strain-comp's `strainlist.txt` exactly so any existing
downstream scripts that parse it work unchanged:

    AC0MRG | Lactobacillus iners P022T4 [GCF_055700695.1; SAMN48590876; PRJNA822121; Contig]
    AC02ZQ | Lactobacillus iners In_Kul_1991H [GCF_056138725.1; SAMN49888361; …]
    …

Default ordering: whatever's currently in `strains.display_order` — i.e. the
ordering the most recent `heatmap` step wrote. By convention that's the
protein dendrogram order. To get a different ordering (e.g. the nucleotide
dendrogram), pass `--cluster-run-id N` or `--sequence-type nucleotide`;
strainbench will recompute the dendrogram on the fly without disturbing the
canonical `strains.display_order` value.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from strainbench.core.db import resolve_cluster_run_id as _resolve_cluster_run_id
from strainbench.producer.export_xlsx import _format_strain_header


@dataclass
class ExportStrainlistSummary:
    output_path: str
    n_strains: int
    cluster_run_id: int | None      # None when using strains.display_order directly
    sequence_type: str | None       # which run drove the order, if recomputed
    recomputed: bool                # True if dendrogram was re-run for this export


def export_strainlist(
    db_path: str | Path,
    output_path: str | Path,
    *,
    cluster_run_id: int | None = None,
    sequence_type: str | None = None,
) -> ExportStrainlistSummary:
    """Write strain headers (one per line) in dendrogram order.

    Args:
        db_path: Path to the strainbench DB.
        output_path: Where to write the strainlist.txt.
        cluster_run_id: If specified, recompute the strain dendrogram from
            this run (without persisting). Otherwise use strains.display_order.
        sequence_type: 'protein' or 'nucleotide'. Same effect as
            cluster_run_id but selects the most recent active run of that
            type.
    """
    db_path = Path(db_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        # Two paths: fast (use strains.display_order) vs recompute.
        if cluster_run_id is not None or sequence_type is not None:
            run_id = _resolve_run_for_strainlist(conn, cluster_run_id, sequence_type)
            run_meta = conn.execute(
                "SELECT sequence_type FROM cluster_runs WHERE cluster_run_id = ?",
                (run_id,),
            ).fetchone()
            run_type = run_meta["sequence_type"] if run_meta else None

            # Compute the dendrogram ordering without persisting it. This
            # imports seaborn etc. — slow but only when the user explicitly
            # asks for a non-canonical ordering.
            from strainbench.producer.heatmap import compute_hierarchical_order
            result = compute_hierarchical_order(
                conn, cluster_run_id=run_id, write_orderings=False,
            )
            strain_order_ids = result.strain_order

            strains_df = pd.read_sql_query(
                "SELECT strain_id, locus_prefix, species, strain_name, "
                "       assembly_id, biosample_id, bioproject_id, assembly_level "
                "FROM strains",
                conn,
            ).set_index("strain_id")
            strains_df = strains_df.loc[strain_order_ids].reset_index()
            recomputed = True
        else:
            run_id = None
            run_type = None
            strains_df = pd.read_sql_query(
                "SELECT strain_id, locus_prefix, species, strain_name, "
                "       assembly_id, biosample_id, bioproject_id, assembly_level "
                "FROM strains "
                "ORDER BY COALESCE(display_order, 999999), strain_id",
                conn,
            )
            recomputed = False
    finally:
        conn.close()

    with output_path.open("w") as f:
        for _, row in strains_df.iterrows():
            f.write(_format_strain_header(row) + "\n")

    return ExportStrainlistSummary(
        output_path=str(output_path),
        n_strains=len(strains_df),
        cluster_run_id=run_id,
        sequence_type=run_type,
        recomputed=recomputed,
    )


def _resolve_run_for_strainlist(
    conn: sqlite3.Connection,
    cluster_run_id: int | None,
    sequence_type: str | None,
) -> int:
    """Pick a cluster_run_id for recomputation, accepting either id or type."""
    if cluster_run_id is not None and sequence_type is not None:
        # Both given — verify they agree
        row = conn.execute(
            "SELECT sequence_type FROM cluster_runs WHERE cluster_run_id = ?",
            (cluster_run_id,),
        ).fetchone()
        if row is None:
            raise ValueError(f"cluster_run_id={cluster_run_id} not found")
        if row["sequence_type"] != sequence_type:
            raise ValueError(
                f"cluster_run_id={cluster_run_id} is sequence_type="
                f"{row['sequence_type']!r}, but --sequence-type was "
                f"{sequence_type!r}. Drop one of the flags."
            )
        return cluster_run_id
    return _resolve_cluster_run_id(
        conn, cluster_run_id, sequence_type_hint=sequence_type
    )
