"""Export every member of a cluster as `gene + N bp upstream` (and/or downstream).

Use case: "give me 200 bp of promoter region for all 343 vaginolysin orthologs
across my Gardnerella pan-genome, ready for an MSA." The cluster_id is the
strainbench-native identifier — it survives re-annotation and is unaffected by
which strains happened to have a given gene labeled in their gbff.

Why this can't reuse the per-strain `fna/` sidecars: ingest stores each CDS
sliced exactly to its bounds. Flanking sequence was thrown away. So we go back
to the source gbff for every member and re-extract with a wider window.

Output is multi-FASTA, one record per cluster member, in the gene's natural
5'→3' orientation. Each header carries enough metadata to map back to the
source row in the strainbench DB without any extra lookup.

The cluster can be selected three ways:
    --cluster-name GARD_000458       (exact, recommended for scripts)
    --gene-name vly                  (looks at cds.gene_name)
    --annotation-pattern "vaginolysin"  (LIKE search on cds.annotation,
                                         then resolved to the cluster)

For gene-name / annotation-pattern, we error out if matches resolve to more
than one cluster — strainbench can't guess which one you wanted. The error
message lists the candidates so you can pick one and re-run with
--cluster-name.
"""

from __future__ import annotations

import sqlite3
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from strainbench.core.db import resolve_cluster_run_id as _resolve_cluster_run_id

_FASTA_LINE_WIDTH = 80


@dataclass
class ExportRegionSummary:
    output_path: str
    cluster_run_id: int
    cluster_name: str
    cluster_id: int
    upstream_bp: int
    downstream_bp: int
    n_members_total: int
    n_records_written: int
    n_strains_missing_gbff: int = 0
    n_members_missing_locus_tag: int = 0
    n_members_truncated_window: int = 0       # got fewer flanking bp than requested
    missing_gbff_paths: list[str] = field(default_factory=list)
    missing_locus_tags: list[str] = field(default_factory=list)


def export_cluster_region(
    db_path: str | Path,
    output_path: str | Path,
    *,
    cluster_name: str | None = None,
    gene_name: str | None = None,
    annotation_pattern: str | None = None,
    cluster_run_id: int | None = None,
    upstream_bp: int = 100,
    downstream_bp: int = 0,
    gbff_root: str | Path | None = None,
) -> ExportRegionSummary:
    """Write a multi-FASTA of every cluster member with N bp flanking.

    Args:
        db_path: Path to a strainbench DB.
        output_path: Where to write the FASTA.
        cluster_name: Exact cluster name (e.g. "GARD_000458").
        gene_name: Look up cluster via cds.gene_name (e.g. "vly").
        annotation_pattern: Look up via cds.annotation LIKE '%pattern%'.
            Exactly one of cluster_name / gene_name / annotation_pattern
            must be set.
        cluster_run_id: Constrain cluster lookup to this run. None →
            most-recent active protein run (the natural "ortholog group" axis).
        upstream_bp: Bases of 5' UTR to prepend (strand-aware).
        downstream_bp: Bases of 3' UTR to append (strand-aware).
        gbff_root: Directory the relative paths in `strains.source_file` are
            anchored to. Default: the directory containing the .db file.
            (strainbench's own `ingest` writes relative paths, so as long as
            you keep the gbff layout next to the .db this never has to be set.)

    Raises:
        ValueError: ambiguous lookup, no match, missing BioPython, etc.
    """
    selectors = sum(x is not None for x in (cluster_name, gene_name, annotation_pattern))
    if selectors != 1:
        raise ValueError(
            "Exactly one of --cluster-name, --gene-name, or --annotation-pattern "
            "must be supplied"
        )

    db_path = Path(db_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gbff_root = Path(gbff_root) if gbff_root else db_path.parent

    try:
        from Bio import SeqIO  # noqa: F401  (used in _extract_members)
    except ImportError as exc:
        raise RuntimeError(
            "BioPython is required for export-cluster-region. "
            "Install with `pip install '.[producer]'`."
        ) from exc

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        run_id = _resolve_cluster_run_id(
            conn, cluster_run_id, sequence_type_hint="protein",
        )

        cluster_row = _resolve_cluster(
            conn, run_id,
            cluster_name=cluster_name,
            gene_name=gene_name,
            annotation_pattern=annotation_pattern,
        )

        # All CDS members of the chosen cluster, grouped by strain so we open
        # each gbff at most once even for multi-copy clusters.
        members = conn.execute(
            """
            SELECT
                cds.cds_id,
                cds.locus_tag,
                cds.gene_name,
                cds.annotation,
                cds.location,
                cds.direction,
                cds.notes,
                s.strain_id,
                s.locus_prefix,
                s.source_file
            FROM cluster_membership m
            JOIN cds      ON cds.cds_id     = m.cds_id
            JOIN strains  s ON s.strain_id  = cds.strain_id
            WHERE m.cluster_run_id = ? AND m.cluster_id = ?
            ORDER BY s.locus_prefix, cds.locus_tag
            """,
            (run_id, cluster_row["cluster_id"]),
        ).fetchall()
    finally:
        conn.close()

    if not members:
        raise ValueError(
            f"cluster {cluster_row['cluster_name']!r} (cluster_id={cluster_row['cluster_id']}) "
            f"has no members in cluster_run_id={run_id}"
        )

    summary = ExportRegionSummary(
        output_path=str(output_path),
        cluster_run_id=run_id,
        cluster_name=cluster_row["cluster_name"],
        cluster_id=cluster_row["cluster_id"],
        upstream_bp=upstream_bp,
        downstream_bp=downstream_bp,
        n_members_total=len(members),
        n_records_written=0,
    )

    # Group by strain → one gbff parse per strain, all members extracted in one pass.
    by_strain: dict[int, list[sqlite3.Row]] = defaultdict(list)
    for m in members:
        by_strain[m["strain_id"]].append(m)

    with output_path.open("w") as out:
        for strain_id, strain_members in by_strain.items():
            source_file = strain_members[0]["source_file"]
            gbff_path = (gbff_root / source_file).resolve()
            if not gbff_path.is_file():
                # Try the absolute path as recorded (for back-compat with old DBs
                # that might have written full paths).
                gbff_path = Path(source_file)
                if not gbff_path.is_file():
                    summary.n_strains_missing_gbff += 1
                    summary.missing_gbff_paths.append(str(strain_members[0]["source_file"]))
                    continue

            written, missing, truncated = _extract_and_write(
                gbff_path, strain_members, upstream_bp, downstream_bp, out,
            )
            summary.n_records_written += written
            summary.n_members_missing_locus_tag += missing
            summary.n_members_truncated_window += truncated
            summary.missing_locus_tags.extend(missing for _ in range(0))  # placeholder

    return summary


# ─────────────────────────────────────────────────────────────────────────────
# Internals
# ─────────────────────────────────────────────────────────────────────────────


def _resolve_cluster(
    conn: sqlite3.Connection,
    run_id: int,
    *,
    cluster_name: str | None,
    gene_name: str | None,
    annotation_pattern: str | None,
) -> sqlite3.Row:
    """Return exactly one matching cluster row, else raise ValueError with help text."""
    if cluster_name is not None:
        rows = conn.execute(
            """
            SELECT cluster_id, cluster_name, consensus_annotation, member_count, strain_count
            FROM clusters
            WHERE cluster_run_id = ? AND cluster_name = ?
            """,
            (run_id, cluster_name),
        ).fetchall()
        if not rows:
            raise ValueError(
                f"no cluster named {cluster_name!r} in cluster_run_id={run_id}"
            )
        return rows[0]

    # gene_name / annotation_pattern: find candidate clusters via their members
    if gene_name is not None:
        where = "LOWER(cds.gene_name) = LOWER(?)"
        params: tuple = (gene_name,)
        descr = f"gene_name == {gene_name!r}"
    else:
        where = "LOWER(cds.annotation) LIKE LOWER(?)"
        params = (f"%{annotation_pattern}%",)
        descr = f"annotation LIKE %{annotation_pattern}%"

    rows = conn.execute(
        f"""
        SELECT
            c.cluster_id,
            c.cluster_name,
            c.consensus_annotation,
            c.member_count,
            c.strain_count,
            COUNT(*) AS matching_member_count
        FROM cluster_membership m
        JOIN cds      ON cds.cds_id     = m.cds_id
        JOIN clusters c ON c.cluster_id = m.cluster_id
        WHERE m.cluster_run_id = ?
          AND {where}
        GROUP BY c.cluster_id
        ORDER BY matching_member_count DESC, c.strain_count DESC
        """,
        (run_id, *params),
    ).fetchall()

    if not rows:
        raise ValueError(f"no clusters in cluster_run_id={run_id} have {descr}")
    if len(rows) > 1:
        # Multiple cluster matches — show them so the user can disambiguate.
        lines = [
            f"  {r['cluster_name']:<20s}  members={r['member_count']:<5d}  "
            f"strains={r['strain_count']:<5d}  annotation={(r['consensus_annotation'] or '')!r}"
            for r in rows[:10]
        ]
        more = f"\n  ... and {len(rows) - 10} more" if len(rows) > 10 else ""
        raise ValueError(
            f"{descr} matches {len(rows)} clusters in cluster_run_id={run_id}; "
            f"use --cluster-name to pick one:\n" + "\n".join(lines) + more
        )
    return rows[0]


def _extract_and_write(
    gbff_path: Path,
    members: list[sqlite3.Row],
    upstream_bp: int,
    downstream_bp: int,
    out,
) -> tuple[int, int, int]:
    """Parse one gbff, extract every wanted CDS + flanking, write FASTA records.

    Returns (n_written, n_missing_locus_tag, n_truncated_window).
    """
    from Bio import SeqIO

    wanted_tags = {m["locus_tag"]: m for m in members}
    found_tags: set[str] = set()
    n_written = 0
    n_truncated = 0

    for record in SeqIO.parse(gbff_path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            tag = (feature.qualifiers.get("locus_tag") or [None])[0]
            if tag is None or tag not in wanted_tags or tag in found_tags:
                continue
            found_tags.add(tag)

            member = wanted_tags[tag]
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand
            contig_len = len(record.seq)

            if strand == 1:
                ext_start = max(0, start - upstream_bp)
                ext_end = min(contig_len, end + downstream_bp)
                got_upstream = start - ext_start
                got_downstream = ext_end - end
                seq = record.seq[ext_start:ext_end]
            else:
                ext_start = max(0, start - downstream_bp)  # 3' end (gene-relative)
                ext_end = min(contig_len, end + upstream_bp)  # 5' end (gene-relative)
                got_upstream = ext_end - end
                got_downstream = start - ext_start
                seq = record.seq[ext_start:ext_end].reverse_complement()

            if got_upstream < upstream_bp or got_downstream < downstream_bp:
                n_truncated += 1

            product = (feature.qualifiers.get("product") or ["?"])[0]
            header_bits = [
                f">{member['locus_prefix']}|{tag}",
                f"[upstream={got_upstream}/{upstream_bp}bp]",
                f"[downstream={got_downstream}/{downstream_bp}bp]",
                f"[direction={'F' if strand == 1 else 'R'}]",
                f"[contig={record.id}]",
                f"[product={product}]",
            ]
            out.write(" ".join(header_bits) + "\n")
            seq_str = str(seq)
            for i in range(0, len(seq_str), _FASTA_LINE_WIDTH):
                out.write(seq_str[i:i + _FASTA_LINE_WIDTH] + "\n")
            n_written += 1

            if len(found_tags) == len(wanted_tags):
                # Found everything we need from this gbff — stop early.
                return n_written, 0, n_truncated

    n_missing = len(wanted_tags) - len(found_tags)
    return n_written, n_missing, n_truncated
