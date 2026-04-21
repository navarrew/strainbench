"""`strainbench` console script.

Subcommands so far:
    init        — create a new empty SQLite database
    parse-gbff  — parse one .gbff and print a summary (sanity check)
    ingest      — parse + load gbff(s) into a DB and write sidecar FASTA files
    cluster     — run mmseqs2 over the ingested CDSs and store the result
    heatmap     — hierarchical clustering: writes display_order columns + PNG
    export-xlsx — export the strain × cluster matrix as a formatted .xlsx
    export-strain-table — strain-anchored view: rows = one strain's CDSs in genomic order
    export-reps — export cluster representative sequences as FASTA
    annotate-prep — generate a workdir with shell templates for external annotators
    import-annotation — import annotator output (single file or whole workdir)
    status — human-readable snapshot of a strainbench DB

Producer-side imports are deferred until their subcommand actually runs,
so users with only the web-side install can still run `init`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from strainbench import __version__
from strainbench.core import db as core_db


def cmd_init(args: argparse.Namespace) -> int:
    target = Path(args.path)
    existed = target.exists()
    if existed and args.force:
        target.unlink()
        existed = False
    # initialize() is idempotent: creates fresh if missing, applies migrations
    # if an existing strainbench DB is older than CURRENT_SCHEMA_VERSION.
    core_db.initialize(target)
    version = core_db.schema_version(target)
    action = "upgraded" if existed else "initialized"
    print(f"strainbench database {action} (schema v{version}): {target}")
    return 0


def cmd_parse_gbff(args: argparse.Namespace) -> int:
    # Producer-only deps (BioPython); imported lazily.
    from strainbench.producer.parsers.gbff import parse_gbff
    from strainbench.producer.parsers.ncbi_metadata import load_assembly_table

    metadata = None
    if args.metadata:
        table = load_assembly_table(args.metadata)
        accession = Path(args.path).stem  # e.g. "GCF_000160875.1"
        metadata = table.get(accession)
        if metadata is None:
            print(
                f"Warning: no metadata row found for accession '{accession}' "
                f"in {args.metadata}",
                file=sys.stderr,
            )

    record = parse_gbff(args.path, metadata=metadata)

    print(f"locus_prefix:    {record.locus_prefix}")
    print(f"species:         {record.species}")
    print(f"strain_name:     {record.strain_name}")
    print(f"assembly_id:     {record.assembly_id}")
    print(f"assembly_level:  {record.assembly_level}")
    print(f"biosample_id:    {record.biosample_id}")
    print(f"bioproject_id:   {record.bioproject_id}")
    print(f"source_format:   {record.source_format}")
    print(f"source_file:     {record.source_file}")
    print(f"cds_count:       {len(record.cds_records)}")

    n = max(0, args.show_cds)
    if n:
        print(f"\nFirst {min(n, len(record.cds_records))} CDS records:")
        for cds in record.cds_records[:n]:
            note = f"  [{cds.notes}]" if cds.notes else ""
            print(
                f"  {cds.locus_tag:20}  {cds.protein_id or '?':18}  "
                f"{cds.direction}  {cds.location:32}  "
                f"GC={cds.gc_pct:5.2f}%  aa={cds.aa_length:4d}  "
                f'"{cds.annotation[:50]}"{note}'
            )
    return 0


def cmd_ingest(args: argparse.Namespace) -> int:
    # Producer-only deps; lazy import.
    import glob
    import time

    from strainbench.producer.ingest import ingest_records
    from strainbench.producer.parsers.gbff import parse_gbff
    from strainbench.producer.parsers.ncbi_metadata import load_assembly_table

    db_path = Path(args.db)
    if not db_path.exists():
        print(
            f"Database not found: {db_path}. Run `strainbench init {db_path}` first.",
            file=sys.stderr,
        )
        return 1

    fasta_dir = Path(args.fasta_dir)

    # Resolve the input list: either one --gbff file, or all *.gbff in --gbff-dir.
    if args.gbff and args.gbff_dir:
        print("--gbff and --gbff-dir are mutually exclusive", file=sys.stderr)
        return 2
    if args.gbff:
        gbff_paths = [Path(args.gbff)]
    elif args.gbff_dir:
        gbff_paths = sorted(Path(p) for p in glob.glob(f"{args.gbff_dir}/*.gbff"))
    else:
        print("must supply either --gbff or --gbff-dir", file=sys.stderr)
        return 2

    if args.limit:
        gbff_paths = gbff_paths[: args.limit]
    if not gbff_paths:
        print("no .gbff files matched", file=sys.stderr)
        return 1

    metadata_table = load_assembly_table(args.metadata) if args.metadata else {}

    print(f"Parsing {len(gbff_paths)} gbff file(s)...")
    t0 = time.time()
    records = []
    parse_failures: list[tuple[Path, str]] = []
    for path in gbff_paths:
        accession = path.stem
        try:
            records.append(parse_gbff(path, metadata=metadata_table.get(accession)))
        except Exception as exc:  # noqa: BLE001
            parse_failures.append((path, str(exc)))
    parse_dt = time.time() - t0
    print(f"  parsed {len(records)} ok, {len(parse_failures)} failed ({parse_dt:.1f}s)")
    for path, msg in parse_failures[:5]:
        print(f"  PARSE FAIL  {path.name}: {msg}", file=sys.stderr)

    # Open one connection for the whole batch — much faster than per-strain.
    from strainbench.core import db as core_db

    print(f"Ingesting into {db_path} (sidecar FASTAs in {fasta_dir})...")
    t0 = time.time()
    conn = core_db.connect(db_path)
    try:
        results = ingest_records(
            conn,
            records,
            fasta_dir,
            on_duplicate=args.on_duplicate,
            progress=args.verbose,
        )
    finally:
        conn.close()
    ingest_dt = time.time() - t0

    by_status = {"inserted": 0, "replaced": 0, "skipped": 0, "failed": 0}
    for r in results:
        by_status[r.status] = by_status.get(r.status, 0) + 1
    print(
        f"\nIngest complete in {ingest_dt:.1f}s:  "
        + ", ".join(f"{k}={v}" for k, v in by_status.items())
    )
    if by_status["failed"]:
        print("\nFailures:", file=sys.stderr)
        for r in results:
            if r.status == "failed":
                print(f"  {r.locus_prefix or '?':14s}  {r.source_file}: {r.error}", file=sys.stderr)

    # Non-zero exit if anything failed (parse OR ingest).
    return 1 if (parse_failures or by_status["failed"]) else 0


def cmd_cluster(args: argparse.Namespace) -> int:
    import time

    from strainbench.core import db as core_db
    from strainbench.producer.cluster import (
        ClusterError,
        MMseqsError,
        cluster_strains,
    )

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Database not found: {db_path}", file=sys.stderr)
        return 1

    print(f"Clustering CDSs in {db_path} ({args.sequence_type}, "
          f"id≥{args.pct_identity}%, cov≥{args.coverage}%)")
    t0 = time.time()
    conn = core_db.connect(db_path)
    try:
        result = cluster_strains(
            conn,
            args.fasta_dir,
            sequence_type=args.sequence_type,
            pct_identity=args.pct_identity,
            coverage=args.coverage,
            name_prefix=args.name_prefix,
            label=args.label,
            work_dir=args.work_dir,
            keep_work_dir=args.keep_work_dir,
        )
    except (MMseqsError, ClusterError) as exc:
        print(f"\nClustering failed: {exc}", file=sys.stderr)
        return 1
    finally:
        conn.close()
    dt = time.time() - t0

    print(
        f"\nDone in {dt:.1f}s.\n"
        f"  cluster_run_id: {result.cluster_run_id}\n"
        f"  label:          {result.label}\n"
        f"  clusters:       {result.cluster_count:,}\n"
        f"  total members:  {result.member_count:,}\n"
        f"  sequence_type:  {result.sequence_type}"
    )
    return 0


def cmd_heatmap(args: argparse.Namespace) -> int:
    import time

    from strainbench.core import db as core_db
    from strainbench.producer.heatmap import compute_hierarchical_order

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Database not found: {db_path}", file=sys.stderr)
        return 1
    # Auto-migrate older DBs so display_order columns exist.
    core_db.initialize(db_path)

    strain_range = None
    if args.range:
        try:
            lo, hi = args.range.split(":", 1)
            strain_range = (int(lo), int(hi))
        except ValueError:
            print(f"--range must look like 'min:max' (got {args.range!r})", file=sys.stderr)
            return 2

    output_png = Path(args.output_png) if args.output_png else (db_path.parent / "heatmap.png")

    print(
        f"Hierarchical clustering on {db_path}\n"
        f"  method:        {args.method}\n"
        f"  strain range:  {strain_range or 'all single-copy clusters'}\n"
        f"  output PNG:    {output_png}"
    )
    t0 = time.time()
    conn = core_db.connect(db_path)
    try:
        result = compute_hierarchical_order(
            conn,
            method=args.method,
            output_png=output_png,
            strain_range=strain_range,
            color=args.color,
        )
    except ValueError as exc:
        print(f"\nHeatmap failed: {exc}", file=sys.stderr)
        return 1
    finally:
        conn.close()
    dt = time.time() - t0

    print(
        f"\nDone in {dt:.1f}s.\n"
        f"  cluster_run_id:           {result.cluster_run_id}\n"
        f"  clusters in clustering:   {result.n_clusters_in_clustering:,}\n"
        f"  clusters excluded:        {result.n_clusters_excluded:,} (multi-copy or out of range)\n"
        f"  strains:                  {result.n_strains}\n"
        f"  display_order written to clusters and strains tables.\n"
        f"  Re-run `strainbench export-xlsx` to get a dendrogram-ordered spreadsheet."
    )
    return 0


def cmd_export_xlsx(args: argparse.Namespace) -> int:
    import time

    from strainbench.producer.export_xlsx import export_cluster_table_xlsx

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Database not found: {db_path}", file=sys.stderr)
        return 1

    print(f"Exporting cluster table from {db_path} → {args.output}")
    t0 = time.time()
    try:
        summary = export_cluster_table_xlsx(
            db_path, args.output, cluster_run_id=args.cluster_run_id
        )
    except ValueError as exc:
        print(f"\nExport failed: {exc}", file=sys.stderr)
        return 1
    dt = time.time() - t0

    print(
        f"\nDone in {dt:.1f}s.\n"
        f"  cluster_run_id: {summary.cluster_run_id}\n"
        f"  rows (clusters): {summary.cluster_count:,}\n"
        f"  cols (strains):  {summary.strain_count:,}\n"
        f"  output:         {summary.output_path}"
    )
    return 0


def cmd_export_strain_table(args: argparse.Namespace) -> int:
    import time

    from strainbench.producer.export_xlsx import export_strain_anchored_xlsx

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Database not found: {db_path}", file=sys.stderr)
        return 1

    print(f"Building strain-anchored table for {args.strain} from {db_path}")
    t0 = time.time()
    try:
        summary = export_strain_anchored_xlsx(
            db_path,
            args.output,
            args.strain,
            cluster_run_id=args.cluster_run_id,
            nt_cluster_run_id=args.nt_cluster_run_id,
        )
    except ValueError as exc:
        print(f"\nExport failed: {exc}", file=sys.stderr)
        return 1
    dt = time.time() - t0

    print(
        f"\nDone in {dt:.1f}s.\n"
        f"  cluster_run_id: {summary.cluster_run_id}\n"
        f"  rows (CDSs in {args.strain}): {summary.cluster_count:,}\n"
        f"  cols (strains):                {summary.strain_count}\n"
        f"  output: {summary.output_path}"
    )
    return 0


def cmd_export_reps(args: argparse.Namespace) -> int:
    import time

    from strainbench.producer.export_reps import export_cluster_representatives

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Database not found: {db_path}", file=sys.stderr)
        return 1

    if args.full_length_only and not args.fasta_dir:
        print("--full-length-only requires --fasta-dir", file=sys.stderr)
        return 2

    mode = "full-length-only" if args.full_length_only else "stored reps"
    print(f"Exporting cluster representatives from {db_path} ({mode}) → {args.output}")
    t0 = time.time()
    try:
        summary = export_cluster_representatives(
            db_path, args.output,
            cluster_run_id=args.cluster_run_id,
            sequence_type=args.sequence_type,
            full_length_only=args.full_length_only,
            fasta_dir=args.fasta_dir,
        )
    except ValueError as exc:
        print(f"\nExport failed: {exc}", file=sys.stderr)
        return 1
    dt = time.time() - t0
    print(f"\nDone in {dt:.1f}s.")
    print(f"  cluster_run_id: {summary.cluster_run_id}")
    print(f"  sequence_type:  {summary.sequence_type}")
    print(f"  reps written:   {summary.n_reps_written:,}")
    if summary.full_length_only:
        print(f"  all-truncated:  {summary.n_clusters_all_truncated:,} clusters skipped "
              f"(no untruncated member)")
        if summary.n_clusters_skipped:
            print(f"  sidecar misses: {summary.n_clusters_skipped} clusters (record not in FASTA)")
    else:
        print(f"  reps skipped:   {summary.n_clusters_skipped} (no sequence stored)")
    print(f"  output:         {summary.output_path}")
    return 0


def cmd_import_annotation(args: argparse.Namespace) -> int:
    import time

    from strainbench.core import db as core_db
    from strainbench.producer.annotations import (
        import_annotation_file,
        import_workdir,
    )

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Database not found: {db_path}", file=sys.stderr)
        return 1
    # Auto-migrate older DBs so cluster_annotations + provenance columns exist.
    core_db.initialize(db_path)

    if args.workdir:
        print(f"Scanning {args.workdir} for annotation outputs...")
        t0 = time.time()
        try:
            summary = import_workdir(
                db_path, args.workdir,
                cluster_run_id=args.cluster_run_id,
                replace=args.replace,
            )
        except (ValueError, FileNotFoundError) as exc:
            print(f"\nImport failed: {exc}", file=sys.stderr)
            return 1
        dt = time.time() - t0

        print(f"\nScanned in {dt:.1f}s:")
        print(f"  files scanned:  {summary.total_files_scanned}")
        print(f"  files imported: {summary.total_files_imported}")
        print(f"  files skipped:  {summary.total_files_skipped}")
        print(f"  rows inserted:  {summary.total_rows_inserted:,}")
        if summary.per_file:
            print("\nPer-file detail:")
            for rel, fmt, result in summary.per_file:
                if hasattr(result, "rows_inserted"):
                    sources = ", ".join(f"{k}={v}" for k, v in result.sources_written.items())
                    print(f"  [{fmt:9s}] {rel}: {result.rows_inserted:,} rows  ({sources})")
                else:
                    print(f"  [{fmt:9s}] {rel}: {result}")
        return 0 if summary.total_files_imported > 0 else 1

    # Single-file mode
    if not args.input or not args.format:
        print("--input and --format are required (or use --workdir)", file=sys.stderr)
        return 2

    print(f"Importing annotations: format={args.format}, input={args.input}")
    t0 = time.time()
    try:
        summary = import_annotation_file(
            db_path, args.input,
            source_format=args.format,
            cluster_run_id=args.cluster_run_id,
            replace=args.replace,
        )
    except (ValueError, FileNotFoundError) as exc:
        print(f"\nImport failed: {exc}", file=sys.stderr)
        return 1
    dt = time.time() - t0
    print(f"\nDone in {dt:.1f}s.")
    print(f"  rows inserted:      {summary.rows_inserted:,}")
    print(f"  clusters touched:   {summary.clusters_touched:,}")
    print(f"  unmatched clusters: {summary.unmatched_clusters} (not found in DB)")
    print("  per source:")
    for source, count in sorted(summary.sources_written.items(), key=lambda kv: -kv[1]):
        print(f"    {source:10s} {count:>7,} rows")
    return 0


def cmd_annotate_prep(args: argparse.Namespace) -> int:
    import time

    from strainbench.core import db as core_db
    from strainbench.producer.annotate_prep import prepare_annotation_workdir

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Database not found: {db_path}", file=sys.stderr)
        return 1
    core_db.initialize(db_path)  # upgrades old DBs so we have cluster_annotations

    print(f"Preparing annotation workdir at {args.workdir}")
    t0 = time.time()
    summary = prepare_annotation_workdir(
        db_path, args.workdir,
        cluster_run_id=args.cluster_run_id,
        overwrite=args.overwrite,
    )
    dt = time.time() - t0
    print(
        f"\nDone in {dt:.1f}s.\n"
        f"  workdir:         {summary.workdir}\n"
        f"  cluster_run_id:  {summary.cluster_run_id}\n"
        f"  reps written:    {summary.reps_written:,}\n"
        f"  scripts:         {', '.join(summary.scripts)}\n\n"
        f"Next:\n"
        f"  1. cd {summary.workdir}\n"
        f"  2. Read README.md to configure env vars for your databases.\n"
        f"  3. Activate each tool's conda env and run the matching script:\n"
        f"       conda activate emapper && bash eggnog.sh\n"
        f"       conda activate kofam    && bash kofamscan.sh\n"
        f"       conda activate deepnog  && bash deepnog.sh\n"
        f"  4. Back in the strainbench env:\n"
        f"       strainbench import-annotation --db {db_path} --workdir {summary.workdir}\n"
    )
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    from strainbench.producer.status import collect_status, format_status

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Database not found: {db_path}", file=sys.stderr)
        return 1
    status = collect_status(db_path)
    print(format_status(status))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="strainbench",
        description="Pan-genome analysis workbench.",
    )
    parser.add_argument("--version", action="version", version=f"strainbench {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    p_init = subparsers.add_parser(
        "init", help="Initialize a new empty strainbench database."
    )
    p_init.add_argument("path", help="Path where the new .db file should be created.")
    p_init.add_argument(
        "--force", action="store_true", help="Overwrite the target file if it exists."
    )
    p_init.set_defaults(func=cmd_init)

    p_parse = subparsers.add_parser(
        "parse-gbff",
        help="Parse one .gbff file and print a summary (sanity check; no DB write).",
    )
    p_parse.add_argument("path", help="Path to the .gbff file.")
    p_parse.add_argument(
        "--metadata",
        help="Optional path to an ncbi_strain_assembly_table.tab file. "
             "The .gbff filename's stem (e.g. 'GCF_000160875.1') is used "
             "as the lookup key.",
    )
    p_parse.add_argument(
        "--show-cds",
        type=int,
        default=5,
        help="Number of CDS records to display (default: 5; pass 0 to suppress).",
    )
    p_parse.set_defaults(func=cmd_parse_gbff)

    p_ingest = subparsers.add_parser(
        "ingest",
        help="Parse and load gbff(s) into a strainbench DB plus sidecar FASTA files.",
    )
    p_ingest.add_argument("--db", required=True, help="Path to an existing strainbench .db file.")
    p_ingest.add_argument(
        "--fasta-dir",
        required=True,
        help="Parent directory for sidecar FASTA files. "
             "Per-strain files land at <fasta-dir>/fna/<prefix>.fna and "
             "<fasta-dir>/faa/<prefix>.faa.",
    )
    src = p_ingest.add_mutually_exclusive_group(required=True)
    src.add_argument("--gbff", help="Path to a single .gbff file.")
    src.add_argument("--gbff-dir", help="Directory of .gbff files (all *.gbff are ingested).")
    p_ingest.add_argument(
        "--metadata",
        help="Optional ncbi_strain_assembly_table.tab to enrich strain metadata. "
             "Lookup key is the gbff filename stem (e.g. 'GCF_000160875.1').",
    )
    p_ingest.add_argument(
        "--on-duplicate",
        choices=["error", "skip", "replace"],
        default="error",
        help="What to do when a strain's locus_prefix is already in the DB. Default: error.",
    )
    p_ingest.add_argument("--limit", type=int, help="Process only the first N gbff files.")
    p_ingest.add_argument(
        "-v", "--verbose", action="store_true", help="Print one line per strain."
    )
    p_ingest.set_defaults(func=cmd_ingest)

    p_cluster = subparsers.add_parser(
        "cluster",
        help="Cluster CDSs in the DB with mmseqs2 and store the results.",
    )
    p_cluster.add_argument("--db", required=True, help="Path to a populated strainbench .db file.")
    p_cluster.add_argument(
        "--fasta-dir",
        required=True,
        help="Sidecar FASTA parent dir (containing fna/ and faa/ subdirs).",
    )
    p_cluster.add_argument(
        "--sequence-type",
        choices=["protein", "nucleotide"],
        default="protein",
        help="Cluster protein sequences (default) or nucleotide sequences.",
    )
    p_cluster.add_argument(
        "--pct-identity", type=int, default=80,
        help="mmseqs --min-seq-id as a percentage (default: 80).",
    )
    p_cluster.add_argument(
        "--coverage", type=int, default=90,
        help="mmseqs -c (mutual coverage) as a percentage (default: 90).",
    )
    p_cluster.add_argument(
        "--name-prefix", default="CLUSTER",
        help="Prefix for cluster names; e.g. 'INERS' → 'INERS_000001' (default: CLUSTER).",
    )
    p_cluster.add_argument(
        "--label", default=None,
        help="Human-readable label for this run; must be unique. "
             "Defaults to '<seqtype>-<pct>id-<cov>cov-<timestamp>'.",
    )
    p_cluster.add_argument(
        "--work-dir", default=None,
        help="Directory for mmseqs scratch files. Default: an auto-cleaned tempdir.",
    )
    p_cluster.add_argument(
        "--keep-work-dir", action="store_true",
        help="Keep the mmseqs scratch directory after a successful run "
             "(useful for debugging).",
    )
    p_cluster.set_defaults(func=cmd_cluster)

    p_heatmap = subparsers.add_parser(
        "heatmap",
        help="Hierarchical clustering of strains and clusters; "
             "writes display_order columns and an optional PNG.",
    )
    p_heatmap.add_argument("--db", required=True, help="Path to a populated strainbench .db file.")
    p_heatmap.add_argument(
        "--method", default="average",
        help="Linkage method: average (default), ward, complete, single, centroid.",
    )
    p_heatmap.add_argument(
        "--color", default="Blues",
        help="matplotlib colormap (default: Blues).",
    )
    p_heatmap.add_argument(
        "--range", default=None,
        help="Restrict clustering to clusters whose strain_count is in MIN:MAX. "
             "Default: all single-copy clusters.",
    )
    p_heatmap.add_argument(
        "--output-png", default=None,
        help="Path for the heatmap PNG (default: <db_dir>/heatmap.png).",
    )
    p_heatmap.set_defaults(func=cmd_heatmap)

    p_export = subparsers.add_parser(
        "export-xlsx",
        help="Export the strain × cluster matrix as a formatted .xlsx file.",
    )
    p_export.add_argument("--db", required=True, help="Path to a populated strainbench .db file.")
    p_export.add_argument(
        "-o", "--output", required=True,
        help="Output .xlsx file path. Parent directory created if needed.",
    )
    p_export.add_argument(
        "--cluster-run-id", type=int, default=None,
        help="Which cluster_run to export (default: most recent active run).",
    )
    p_export.set_defaults(func=cmd_export_xlsx)

    p_strain = subparsers.add_parser(
        "export-strain-table",
        help="Strain-anchored xlsx: rows = one strain's CDSs in genomic order, "
             "columns = all strains. Multi-copy genes get separate rows.",
    )
    p_strain.add_argument("--db", required=True, help="Path to a populated strainbench .db file.")
    p_strain.add_argument(
        "--strain", required=True,
        help="Locus prefix of the anchor strain, e.g. 'HMPREF0520'.",
    )
    p_strain.add_argument(
        "-o", "--output", required=True,
        help="Output .xlsx file path.",
    )
    p_strain.add_argument(
        "--cluster-run-id", type=int, default=None,
        help="Which protein cluster_run to use (default: most recent active "
             "protein run).",
    )
    p_strain.add_argument(
        "--nt-cluster-run-id", type=int, default=None,
        help="Which nucleotide cluster_run to overlay as the nt_subcluster "
             "column (default: most recent active nucleotide run; column "
             "is omitted if no nucleotide run exists).",
    )
    p_strain.set_defaults(func=cmd_export_strain_table)

    p_reps = subparsers.add_parser(
        "export-reps",
        help="Export cluster representative sequences as FASTA (for external annotators).",
    )
    p_reps.add_argument("--db", required=True)
    p_reps.add_argument("-o", "--output", required=True,
                         help="Output .faa or .fna file path.")
    p_reps.add_argument("--cluster-run-id", type=int, default=None)
    p_reps.add_argument(
        "--sequence-type", choices=["protein", "nucleotide"], default=None,
        help="Which sequence type to export (default: match the cluster_run).",
    )
    p_reps.add_argument(
        "--full-length-only", action="store_true",
        help="Pick any un-truncated cluster member instead of the stored rep, "
             "and skip clusters whose members are all truncated. "
             "Produces the curated 'unique full-length ORF catalog'. "
             "Requires --fasta-dir.",
    )
    p_reps.add_argument(
        "--fasta-dir", default=None,
        help="Sidecar FASTA parent dir (with fna/ and faa/ subdirs). "
             "Required when --full-length-only is set.",
    )
    p_reps.set_defaults(func=cmd_export_reps)

    p_import = subparsers.add_parser(
        "import-annotation",
        help="Parse annotator output(s) and write rows into cluster_annotations. "
             "Pass --input/--format for a single file, or --workdir for a directory "
             "produced by `strainbench annotate-prep`.",
    )
    p_import.add_argument("--db", required=True)
    p_import.add_argument(
        "--workdir", default=None,
        help="Directory to scan for annotation outputs (formats auto-detected).",
    )
    p_import.add_argument(
        "--input", default=None,
        help="Single annotation file (use with --format).",
    )
    p_import.add_argument(
        "--format", default=None,
        choices=["emapper", "kofamscan", "deepnog", "amrfinder",
                 "defensefinder", "padloc", "interproscan"],
        help="Source format of --input. Only needed when importing a single file.",
    )
    p_import.add_argument("--cluster-run-id", type=int, default=None)
    p_import.add_argument(
        "--replace", action="store_true",
        help="Delete existing rows from the same source(s) before inserting.",
    )
    p_import.set_defaults(func=cmd_import_annotation)

    p_prep = subparsers.add_parser(
        "annotate-prep",
        help="Generate a workdir with cluster reps + ready-to-run annotation shell scripts "
             "for eggNOG-mapper, kofamscan, and DeepNOG.",
    )
    p_prep.add_argument("--db", required=True)
    p_prep.add_argument(
        "--workdir", required=True,
        help="Directory to create (or reuse). Will contain reps.faa, shell scripts, README.",
    )
    p_prep.add_argument("--cluster-run-id", type=int, default=None)
    p_prep.add_argument(
        "--overwrite", action="store_true",
        help="Re-export reps.faa even if the file already exists.",
    )
    p_prep.set_defaults(func=cmd_annotate_prep)

    p_status = subparsers.add_parser(
        "status",
        help="Print a human-readable snapshot of what's in a strainbench DB.",
    )
    p_status.add_argument("--db", required=True)
    p_status.set_defaults(func=cmd_status)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
