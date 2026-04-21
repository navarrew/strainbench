# strainbench

A pan-genome analysis workbench for bacterial genomes. Successor to
[strain-comp](../strain-comp) with a SQLite-backed data model that supports
iterative strain addition/removal, nested clustering, and browsable
server-side publication of datasets.

> **Status: early scaffolding.** The schema and contract layer are in place.
> Parsers, clustering, ingestion, export, and the web app are next.

## Architecture at a glance

strainbench is two cooperating halves that share a SQLite database:

```
┌────────────────────────┐          ┌──────────────────────────┐
│       PRODUCER         │          │       WEB CONSUMER       │
│   (local, CPU-heavy)   │          │   (server, read-only)    │
├────────────────────────┤          ├──────────────────────────┤
│ • parse gbff / NCBI    │          │ • browse cluster tables  │
│ • mmseqs2 clustering   │          │ • nested expand/collapse │
│ • write SQLite DB      │──┐  ┌───▶│ • search & filter        │
│ • bundle FASTA sidecars│  │  │    │ • download cluster FASTA │
└────────────────────────┘  │  │    └──────────────────────────┘
                            ▼  │
                  ┌──────────────────────┐
                  │  strainbench.db      │   ← portable artifact
                  │  + sequences.tar.gz  │     (scp this and you're done)
                  └──────────────────────┘
```

The portable artifact between them is a single SQLite file plus optional
FASTA bundles. No server-side compute, no background job queue, no auth
required beyond optional access gating.

## Data model (overview)

- **strains** — one row per input genome.
- **cds** — the atomic unit; one row per coding sequence. Sequences
  themselves live in sidecar FASTA files (`data/fna/`, `data/faa/`),
  indexed by `(strain_id, locus_tag)`.
- **cluster_runs** — a single clustering invocation with its parameters.
  Runs can nest via `parent_run_id` / `parent_cluster_id` to support
  nucleotide sub-clustering within protein clusters.
- **clusters** — a named group of CDSs produced by one run.
- **cluster_membership** — many-to-many between CDSs and clusters, scoped
  by run. One CDS can belong to many clusters across runs.
- **annotations** — pluggable table for COG, KEGG, deepnog, and future
  annotation sources.

See [`src/strainbench/core/schema.sql`](src/strainbench/core/schema.sql)
for the authoritative schema.

## Install

This package uses optional dependency groups so producer-side heavy deps
don't bloat the web server.

```bash
# Producer machine (your Mac or a workstation):
pip install -e ".[producer]"

# Web server:
pip install -e ".[web]"

# Development:
pip install -e ".[producer,web,dev]"
```

Producer-side also requires `mmseqs2` installed on `PATH`.

## Quick smoke test

```bash
strainbench init /tmp/test.db
# → "Initialized empty strainbench database: /tmp/test.db"
```

## Relationship to strain-comp

strain-comp is preserved as the working reference implementation. strainbench
is a significant rewrite, not a refactor — the data model moves from
tab-delimited flat files to SQLite, and the input layer decouples from NCBI's
cds_from_genomic format so gbff (and eventually Prokka, Bakta, etc.) work
as first-class inputs.

Where possible, strainbench reuses battle-tested logic from strain-comp —
but through a canonical `StrainRecord` / `CDSRecord` contract rather than
file-format-specific parsing scattered across scripts.
