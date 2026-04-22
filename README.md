# strainbench

A pan-genome analysis workbench for bacterial genomes. Successor to
[strain-comp](https://github.com/navarrew/strain-comp) with a SQLite-backed
data model that supports iterative strain addition/removal, nested clustering
(protein × nucleotide), pluggable functional annotation from many tools, and
browsable Excel + future web outputs.

> **Status: producer pipeline complete.** The full backend — from raw GenBank
> files to a populated SQLite database with annotated, hierarchically-clustered,
> Excel-exported pan-genomes — is built and tested on real data. The
> student-facing web consumer is the next phase.

---

## Architecture at a glance

strainbench is two cooperating halves that share a SQLite database:

```
┌────────────────────────┐          ┌──────────────────────────┐
│       PRODUCER         │          │      WEB CONSUMER        │
│   (local, CPU-heavy)   │          │  (server, read-only)     │
├────────────────────────┤          ├──────────────────────────┤
│ • parse gbff           │          │ • browse cluster tables  │
│ • mmseqs2 clustering   │          │ • search by gene name    │
│ • run annotators       │──┐  ┌───▶│ • download cluster FASTA │
│ • write SQLite DB      │  │  │    │ • per-species dropdown   │
│ • write Excel exports  │  │  │    │                          │
└────────────────────────┘  │  │    └──────────────────────────┘
                            ▼  │           ↑
              ┌─────────────────────────┐  │ COMING NEXT
              │  <species>.db           │  │
              │  + sidecar FASTAs/      │──┘ (you ship this folder
              │  + .xlsx exports        │     to a webserver)
              └─────────────────────────┘
```

The producer is fully working today (this README documents it). The web
consumer is aspirational — it'll read the same `.db` files the producer
writes, with no schema changes needed.

---

## Pipeline overview

The producer is **one command per stage**, run in order. Each stage reads
the artifacts of earlier stages and adds new ones. You can re-run later
stages without redoing earlier ones.

```
                                 ┌──────────────────────────────────┐
  raw NCBI gbff files            │  external annotators             │
   (one per genome)              │  (each in its own conda env)     │
        │                        │  • DeepNOG                       │
        ▼                        │  • kofamscan                     │
  ┌──────────────┐               │  • AMRFinderPlus                 │
  │ ingest       │  →  strains   │  • eggNOG-mapper (web/local)     │
  │              │     cds       │  • InterProScan (web/local)      │
  │              │     sidecar   │                                  │
  │              │     FASTAs    │                                  │
  └──────────────┘               └──────────────────────────────────┘
        │                                ▲              │
        ▼                                │ runs on      │ produce TSVs
  ┌──────────────┐                       │              │
  │ cluster      │  →  cluster_runs ─────┴──────────┐   │
  │ (mmseqs2)    │     clusters    ┌───────────┐    │   │
  │              │     cluster_    │ export-   │    │   │
  │              │     membership  │ reps      │────┘   │
  └──────────────┘                 └───────────┘        │
        │                                ▲              │
        ▼                                │              │
  ┌──────────────┐               ┌─────────────┐        │
  │ heatmap      │  → display_   │ annotate-   │        │
  │ (seaborn)    │    order on   │ prep        │        │
  │              │    strains +  │ (writes     │        │
  │              │    clusters   │  shell      │        │
  └──────────────┘               │  templates) │        │
        │                        └─────────────┘        │
        │                                ▼              │
        │                         ┌─────────────┐       │
        │                         │ import-     │ ◀─────┘
        │                         │ annotation  │
        │                         │             │ → cluster_annotations
        │                         └─────────────┘
        │                                │
        ▼                                ▼
  ┌──────────────────────────────────────────────────┐
  │ export-xlsx           export-strain-table        │
  │ (pan-genome view)     (per-strain genome walk)   │
  └──────────────────────────────────────────────────┘
        │                                │
        ▼                                ▼
  cluster_table.xlsx          strain_anchored_<name>.xlsx
  (3,661 × 124 grid           (1,175 × 124 grid for L. iners
   for L. iners pan-genome)    HMPREF0520, in genomic order)
```

At every stage `strainbench status --db <file>` tells you what's been done
and what's missing.

---

## Quickstart: build a database from scratch

This walks through building an L. iners pan-genome end-to-end. Adapt the
paths and species name for your own data.

### Prerequisites (one-time setup)

You'll want at least **two conda environments**:

```bash
# Strainbench itself + mmseqs2:
conda create -n strainbench -c conda-forge -c bioconda \
    python=3.11 biopython pandas numpy scipy seaborn matplotlib \
    xlsxwriter mmseqs2
conda activate strainbench
pip install -e /path/to/strainbench   # editable install of this repo
pip install '.[dev]'                  # pulls in openpyxl + pytest if you want tests

# Each annotator wants its own env (they fight each other for dependencies).
# Set these up only as you need them — see "Annotation workflow" below.
```

### 1. Initialize an empty database

```bash
strainbench init iners.db
```

Creates an empty SQLite file with the strainbench schema. Idempotent —
re-running on an existing DB applies any pending migrations without
losing data.

### 2. Get gbff files

Use the NCBI `datasets` CLI tool (your existing strain-comp workflow works
fine here):

```bash
datasets download genome taxon "Lactobacillus iners" \
    --include cds --assembly-source RefSeq --exclude-atypical \
    --annotated --exclude-multi-isolate --mag exclude
unzip ncbi_dataset.zip
# you should now have ncbi_dataset/data/GCF_*/genomic.gbff files
```

Also generate the assembly metadata table:

```bash
dataformat tsv genome \
    --inputfile ncbi_dataset/data/assembly_data_report.jsonl \
    --fields accession,ani-submitted-species,assminfo-biosample-strain,\
assminfo-bioproject,assminfo-biosample-accession,assminfo-level \
    > assembly_metadata.tab
```

### 3. Parse + ingest into the DB

```bash
strainbench ingest \
    --db iners.db \
    --fasta-dir iners_fasta/ \
    --gbff-dir ncbi_dataset/data/   \   # all */genomic.gbff get picked up
    --metadata assembly_metadata.tab
```

What this does:
- For each gbff: extract every CDS feature with a `/translation`, capture
  metadata (locus_tag, protein_id, location, GC%, NCBI annotation, gene_name).
- Write per-strain FASTA files at `iners_fasta/fna/<prefix>.fna` and
  `iners_fasta/faa/<prefix>.faa`.
- Insert one row in `strains` and N rows in `cds` per genome.
- All wrapped in transactions — if one strain fails, the others still load.

Typical timing: ~70 ms per strain. For 124 iners strains: ~10 s.

### 4. Cluster the proteins (and optionally nucleotides)

```bash
# Protein clustering at 80% identity, 90% mutual coverage:
strainbench cluster \
    --db iners.db \
    --fasta-dir iners_fasta/ \
    --sequence-type protein \
    --pct-identity 80 \
    --coverage 90 \
    --name-prefix INERS \
    --label "iners-protein-80id-90cov"

# Optional: also cluster nucleotides at 95% for within-species lineage detection
strainbench cluster \
    --db iners.db \
    --fasta-dir iners_fasta/ \
    --sequence-type nucleotide \
    --pct-identity 95 --coverage 90 \
    --name-prefix INERS_NT --label "iners-nt-95id-90cov"
```

The two cluster runs coexist — every CDS gets a row in `cluster_membership`
for each run. Queries can join across both runs to ask things like "for
each protein cluster, how many distinct nucleotide sub-clusters does it
contain?"

Typical timing: ~10 s per run for ~150k CDSs.

### 5. Hierarchical clustering (heatmap)

```bash
strainbench heatmap --db iners.db --output-png iners_heatmap.png
```

Runs seaborn clustermap on the strain × cluster presence/absence matrix.
Writes the resulting strain ordering and cluster ordering back to the DB
as `strains.display_order` and `clusters.display_order`. Subsequent xlsx
exports use this ordering automatically.

### 6. Functional annotation (the big section — see below)

```bash
# Export protein representatives (one per cluster — far fewer than total CDSs)
strainbench export-reps --db iners.db -o iners_reps.faa

# Generate a workdir with shell scripts for every annotator strainbench knows:
strainbench annotate-prep --db iners.db --workdir iners_annot/

# Then for each annotator (in its own conda env):
cd iners_annot/
conda activate kofam        && bash kofamscan.sh
conda activate deepnog      && bash deepnog.sh
conda activate amrfinder    && bash amrfinder.sh
# eggNOG-mapper: upload reps.faa to https://usegalaxy.ca, download result here
# InterProScan:  upload reps.faa to https://www.ebi.ac.uk/interpro/, download here

# Pull all annotation outputs back into the DB at once:
strainbench import-annotation --db iners.db --workdir iners_annot/ --replace
```

See **Annotation workflow** below for details.

### 7. Export the spreadsheets

```bash
# Pan-genome view: strain × cluster matrix with annotation columns
strainbench export-xlsx --db iners.db -o iners_cluster_table.xlsx

# Strain-anchored view: walk one strain's genome, see each gene's cluster + sub-cluster
strainbench export-strain-table --db iners.db --strain HMPREF0520 \
    -o iners_HMPREF0520_anchored.xlsx
```

### 8. Sanity-check what you've got

```bash
strainbench status --db iners.db
```

Prints a summary of all populated tables, annotation coverage per source,
and warnings (e.g., low coverage that suggests an annotator didn't finish).

---

## Annotation workflow

This is the most failure-prone part of the pipeline because it depends on
external tools. strainbench's job is to **make running them painless and
to import whatever output they produce** — it doesn't run the tools itself
(too brittle across machines).

### Recommended tools

| Tool | What it adds | Where it runs |
|---|---|---|
| **DeepNOG** | COG group (e.g., COG0208) | Local — `deepnog` env, ~10 min |
| **kofamscan** | KEGG KO numbers | Local — `kofam` env, ~30 min |
| **AMRFinderPlus** | AMR / virulence / stress genes | Local — `amrfinder` env, ~15 s |
| **eggNOG-mapper** | COG + KEGG + GO + Pfam + EC + CAZy in one pass | **Web** (usegalaxy.ca or eggnog-mapper.embl.de) |
| **InterProScan** | Pfam/TIGRFAM/PANTHER + InterPro2GO | **Web** (https://www.ebi.ac.uk/interpro/) |

The web options are recommended for eggNOG-mapper and InterProScan because
their local installs require ~50 GB databases each.

### How it works

1. `strainbench export-reps --db <db> -o reps.faa` writes one protein per
   cluster (~3,000 sequences for a typical species — vastly faster to
   annotate than all 150,000 CDSs).
2. `strainbench annotate-prep --db <db> --workdir <dir>/` writes a
   directory containing `reps.faa` plus seven ready-to-run shell scripts
   (`eggnog.sh`, `kofamscan.sh`, `deepnog.sh`, `amrfinder.sh`, etc.). Each
   script reads database paths from environment variables so it works
   without editing across machines.
3. You run whichever annotators you want. Each writes its output to its
   own subdir of the workdir (`emapper_out/`, `kofamscan_out/`, etc.).
4. `strainbench import-annotation --db <db> --workdir <dir>/ --replace`
   auto-detects each output file's format and imports them all at once
   into the `cluster_annotations` table. The `--replace` flag scopes
   deletion by `tool_format`, so re-importing one tool preserves
   contributions from others.

### Updating annotations later

A year from now, you just re-run the same scripts (with refreshed databases)
and re-import. `strainbench status` will tell you what was last imported
and from which file.

---

## Outputs

After a complete pipeline run, you have:

| Artifact | Format | Purpose |
|---|---|---|
| `<species>.db` | SQLite | The source of truth. Ship this to a webserver. |
| `<species>_fasta/fna/*.fna` | FASTA (per strain) | Nucleotide sequences sidecar |
| `<species>_fasta/faa/*.faa` | FASTA (per strain) | Protein sequences sidecar |
| `<species>_reps.faa` | FASTA | One protein per cluster (annotator input) |
| `<species>_cluster_table.xlsx` | Excel | Pan-genome view: clusters × strains, plus annotation columns |
| `<species>_<strain>_anchored.xlsx` | Excel | Walk one strain's genome in order, see each gene's cluster |
| `<species>_heatmap.png` | PNG | Hierarchical clustering visualization |
| `<species>_unique_full_length_nt.fna` | FASTA | One full-length nt rep per nt cluster (for downstream comparative analysis) |

The Excel files have a three-block layout: **frozen core metadata** (cluster
name, gene names, NCBI annotation, counts) | **strain block** (one column
per strain, presence/absence) | **annotation block** (KEGG_KO, COG, GO_MF/BP/CC,
Pfam, EC, CAZy, AMRFinder).

---

## Data model (SQLite schema)

The full DDL is in [`src/strainbench/core/schema.sql`](src/strainbench/core/schema.sql).
Summary:

| Table | Holds | Notes |
|---|---|---|
| `strains` | One row per genome | locus_prefix is the friendly strain ID (e.g., `HMPREF0520`) |
| `cds` | One row per coding sequence — the **atomic unit** | UNIQUE (strain_id, locus_tag) |
| `cluster_runs` | One row per clustering invocation | Protein and nucleotide runs coexist; nested via `parent_run_id` |
| `clusters` | One row per cluster | Includes representative sequences + aggregates (member_count, etc.) |
| `cluster_membership` | Junction: which CDS belongs to which cluster, scoped by run | Lets one CDS belong to multiple clusters across runs |
| `cluster_annotations` | Functional annotations from external tools | Tagged by `source`, `tool_format`, `source_file` for full provenance |
| `annotations` | Per-CDS annotations (currently unused) | Reserved for per-protein facts that differ within a cluster |

A useful view, `v_cluster_cell`, computes the strain × cluster cell content
on the fly — this is what `export-xlsx` pivots into the Excel grid.

---

## Subcommand reference

| Command | What it does |
|---|---|
| `strainbench init <db>` | Create or upgrade a database |
| `strainbench parse-gbff <file>` | Parse one gbff and print a summary (sanity check, no DB write) |
| `strainbench ingest --db --gbff-dir --fasta-dir [--metadata]` | Parse + load gbffs into DB + write FASTA sidecars |
| `strainbench cluster --db --fasta-dir --sequence-type {protein,nucleotide} [--pct-identity N --coverage N --name-prefix X --label Y]` | mmseqs2 clustering |
| `strainbench heatmap --db [--output-png]` | Hierarchical clustering, populates display_order |
| `strainbench export-reps --db -o [--full-length-only --fasta-dir]` | Cluster representative sequences as FASTA |
| `strainbench annotate-prep --db --workdir` | Generate a workdir with `reps.faa` + shell scripts for all annotators |
| `strainbench import-annotation --db {--workdir,--input --format} [--replace]` | Parse annotator output(s) into cluster_annotations |
| `strainbench export-xlsx --db -o` | Pan-genome view (strain × cluster matrix, annotated) |
| `strainbench export-strain-table --db --strain -o` | Strain-anchored view (one strain's genes in genomic order) |
| `strainbench status --db` | Snapshot of what's in the DB + coverage warnings |

Run `strainbench <subcommand> --help` for details on flags.

---

## What's verified working

The full pipeline has been built end-to-end on a 124-strain *Lactobacillus
iners* dataset:

| Stage | Time | Output |
|---|---|---|
| Parse + ingest | 14 s | 124 strains, 148,673 CDSs |
| Protein cluster | 11 s | 3,661 clusters at 80% id |
| Nucleotide cluster | 9 s | 4,064 clusters at 95% id |
| Heatmap | 4 s | 595 KB PNG, display_order populated |
| Annotation: DeepNOG | 10 min | COG group for 100% of clusters |
| Annotation: kofamscan | ~30 min | KEGG KO for 82% of clusters |
| Annotation: AMRFinderPlus | 15 s | 9 AMR genes detected |
| Annotation: eggNOG-mapper (web) | varies | COG/KEGG/GO/Pfam/EC/CAZy for 88% of clusters |
| Export xlsx | 2 s | 2.7 MB cluster_table |

Total wall-clock from raw gbffs to fully-annotated spreadsheet: **~45 minutes**,
of which 35-40 minutes is annotators running.

---

## Roadmap

| Phase | Status |
|---|---|
| Producer pipeline (this README) | ✓ done |
| Per-genome system-detection workflow (DefenseFinder, PADLOC) | planned |
| Centralize duplicate `_resolve_cluster_run_id` helpers | planned (small cleanup) |
| Web consumer: dropdown + read-only browse | aspirational — start when 3+ species datasets exist |
| Multi-species cross-DB queries | requires web phase first |

---

## Relationship to strain-comp

[strain-comp](https://github.com/navarrew/strain-comp) is preserved as the
working reference implementation. strainbench is a significant rewrite, not
a refactor — the data model moves from tab-delimited flat files to SQLite,
and the input layer decouples from NCBI's `cds_from_genomic` format so gbff
(and eventually Prokka, Bakta, etc.) work as first-class inputs.

Where possible, strainbench reuses logic from strain-comp — but through a
canonical `StrainRecord` / `CDSRecord` contract rather than file-format-specific
parsing scattered across scripts.

---

## Install

```bash
# Producer machine (your Mac or a workstation):
git clone https://github.com/navarrew/strainbench
cd strainbench
pip install -e ".[producer]"

# Verify:
strainbench init /tmp/test.db
# → "strainbench database initialized (schema v4): /tmp/test.db"
```

You'll also need `mmseqs2` on your PATH (use the `strainbench` conda env
recipe in the Quickstart section above). Each annotator goes in its own
conda env — see the `annotate-prep` workdir's README for setup commands.
