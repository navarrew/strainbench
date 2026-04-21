-- strainbench schema v1
--
-- Canonical SQLite schema. The producer writes this DB; the web consumer
-- reads it. This file is the single source of truth for the data model.
--
-- Design notes:
--   * CDS is the atomic unit. Strains and clusters are aggregations of CDSs.
--   * cluster_membership is a junction so one CDS can belong to many clusters
--     across different cluster_runs (protein vs nucleotide, different params).
--   * cluster_runs can nest via parent_run_id + parent_cluster_id to express
--     nucleotide sub-clustering within a protein cluster.
--   * Full sequences are NOT stored inline; only representative sequences on
--     the `clusters` table for fast web display. Per-CDS sequences live in
--     sidecar FASTA files under data/fna/ and data/faa/.

PRAGMA foreign_keys = ON;

-- ─────────────────────────────────────────────────────────────────────────────
-- Schema versioning
-- ─────────────────────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS schema_version (
    version      INTEGER PRIMARY KEY,
    applied_at   TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    description  TEXT
);
INSERT OR IGNORE INTO schema_version (version, description)
VALUES (1, 'initial: strains, cds, cluster_runs, clusters, cluster_membership, annotations');
-- v2 and v3 version rows are inserted by migrations in db.py, after their
-- ALTER/CREATE statements succeed.

-- ─────────────────────────────────────────────────────────────────────────────
-- strains: one row per input genome
-- ─────────────────────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS strains (
    strain_id       INTEGER PRIMARY KEY,
    locus_prefix    TEXT    NOT NULL UNIQUE,
    species         TEXT,
    strain_name     TEXT,
    assembly_id     TEXT,
    biosample_id    TEXT,
    bioproject_id   TEXT,
    assembly_level  TEXT    CHECK (assembly_level IN
                        ('Complete', 'Chromosome', 'Scaffold', 'Contig', 'Unknown')),
    source_format   TEXT    CHECK (source_format IN ('gbff', 'ncbi_cds', 'other')),
    source_file     TEXT,
    cds_count       INTEGER,
    display_order   INTEGER,                  -- populated by hierarchical clustering
    added_at        TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX IF NOT EXISTS idx_strains_species ON strains(species);

-- ─────────────────────────────────────────────────────────────────────────────
-- cds: the atomic unit
-- ─────────────────────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS cds (
    cds_id          INTEGER PRIMARY KEY,
    strain_id       INTEGER NOT NULL REFERENCES strains(strain_id) ON DELETE CASCADE,
    locus_tag       TEXT    NOT NULL,
    protein_id      TEXT,
    nuc_accession   TEXT,
    location        TEXT,
    direction       CHAR(1) CHECK (direction IN ('F', 'R')),
    nt_length       INTEGER,
    aa_length       INTEGER,
    gc_pct          REAL,
    gene_name       TEXT,
    annotation      TEXT,
    notes           TEXT,
    UNIQUE (strain_id, locus_tag)
);
CREATE INDEX IF NOT EXISTS idx_cds_strain     ON cds(strain_id);
CREATE INDEX IF NOT EXISTS idx_cds_locus_tag  ON cds(locus_tag);
CREATE INDEX IF NOT EXISTS idx_cds_protein_id ON cds(protein_id);
CREATE INDEX IF NOT EXISTS idx_cds_annotation ON cds(annotation);

-- ─────────────────────────────────────────────────────────────────────────────
-- cluster_runs: one row per clustering invocation.
--
-- Nested runs use parent_run_id (the run whose output we sub-clustered) and
-- parent_cluster_id (the specific cluster within that run we operated on).
-- Top-level runs have both NULL.
--
-- Note: this table has a forward reference to clusters(cluster_id).
-- SQLite does not enforce FKs at DDL time, so creation order is fine; at
-- INSERT time, insert cluster_runs first (with parent_cluster_id NULL),
-- then clusters, then UPDATE cluster_runs.parent_cluster_id for nested runs.
-- ─────────────────────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS cluster_runs (
    cluster_run_id     INTEGER PRIMARY KEY,
    label              TEXT    NOT NULL UNIQUE,
    sequence_type      TEXT    NOT NULL CHECK (sequence_type IN ('protein', 'nucleotide')),
    pct_identity       INTEGER NOT NULL,
    coverage           INTEGER NOT NULL,
    name_prefix        TEXT    NOT NULL,
    parent_run_id      INTEGER REFERENCES cluster_runs(cluster_run_id),
    parent_cluster_id  INTEGER REFERENCES clusters(cluster_id),
    is_active          INTEGER NOT NULL DEFAULT 1,
    created_at         TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    notes              TEXT
);
CREATE INDEX IF NOT EXISTS idx_runs_parent_run     ON cluster_runs(parent_run_id);
CREATE INDEX IF NOT EXISTS idx_runs_parent_cluster ON cluster_runs(parent_cluster_id);

-- ─────────────────────────────────────────────────────────────────────────────
-- clusters: a named group of CDSs, produced by one cluster_run.
--
-- Representative sequences stored inline so cluster detail pages can render
-- without hitting sidecar files. ~65 MB for 50k clusters — negligible.
-- ─────────────────────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS clusters (
    cluster_id              INTEGER PRIMARY KEY,
    cluster_run_id          INTEGER NOT NULL REFERENCES cluster_runs(cluster_run_id)
                                    ON DELETE CASCADE,
    cluster_name            TEXT    NOT NULL,
    representative_cds_id   INTEGER REFERENCES cds(cds_id),
    representative_nt_seq   TEXT,
    representative_aa_seq   TEXT,
    consensus_annotation    TEXT,
    member_count            INTEGER,
    strain_count            INTEGER,
    avg_gc_pct              REAL,
    gc_spread               REAL,
    avg_aa_length           INTEGER,
    flags                   TEXT,
    display_order           INTEGER,                -- populated by hierarchical clustering
    UNIQUE (cluster_run_id, cluster_name)
);
CREATE INDEX IF NOT EXISTS idx_clusters_run        ON clusters(cluster_run_id);
CREATE INDEX IF NOT EXISTS idx_clusters_annotation ON clusters(consensus_annotation);

-- ─────────────────────────────────────────────────────────────────────────────
-- cluster_membership: many-to-many between CDSs and clusters, scoped by run.
-- One CDS has at most one cluster per run, but may belong to clusters in
-- multiple runs (e.g., a protein run AND a nucleotide run).
-- ─────────────────────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS cluster_membership (
    cluster_run_id  INTEGER NOT NULL REFERENCES cluster_runs(cluster_run_id) ON DELETE CASCADE,
    cluster_id      INTEGER NOT NULL REFERENCES clusters(cluster_id)         ON DELETE CASCADE,
    cds_id          INTEGER NOT NULL REFERENCES cds(cds_id)                  ON DELETE CASCADE,
    PRIMARY KEY (cluster_run_id, cds_id)
);
CREATE INDEX IF NOT EXISTS idx_membership_cluster ON cluster_membership(cluster_id);
CREATE INDEX IF NOT EXISTS idx_membership_cds     ON cluster_membership(cds_id);

-- ─────────────────────────────────────────────────────────────────────────────
-- annotations: pluggable table for COG, KEGG, deepnog, and future sources.
-- One CDS can have many annotations (e.g., both a COG code and a KEGG KO).
-- ─────────────────────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS annotations (
    annotation_id  INTEGER PRIMARY KEY,
    cds_id         INTEGER NOT NULL REFERENCES cds(cds_id) ON DELETE CASCADE,
    source         TEXT    NOT NULL,
    code           TEXT,
    description    TEXT,
    score          REAL
);
CREATE INDEX IF NOT EXISTS idx_annotations_cds    ON annotations(cds_id);
CREATE INDEX IF NOT EXISTS idx_annotations_source ON annotations(source);
CREATE INDEX IF NOT EXISTS idx_annotations_code   ON annotations(code);

-- ─────────────────────────────────────────────────────────────────────────────
-- cluster_annotations: functional annotations that apply to a whole CLUSTER
-- rather than to individual CDSs.
--
-- This is the natural home for annotations derived from running a tool on
-- the cluster representative sequences — KEGG/KO via kofamscan, COG, Pfam,
-- GO via eggNOG-mapper, AMRFinderPlus hits, DefenseFinder systems, etc.
-- One cluster typically has many annotations (esp. for GO, where each
-- cluster may have dozens of MF/BP/CC terms) — each lands as its own row.
--
-- `category` is a source-dependent sub-classification:
--    GO:    'MF' / 'BP' / 'CC'
--    COG:   the one-letter category (J, K, L, …)
--    KEGG:  pathway or module ID (when available)
-- ─────────────────────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS cluster_annotations (
    cluster_ann_id  INTEGER PRIMARY KEY,
    cluster_id      INTEGER NOT NULL REFERENCES clusters(cluster_id) ON DELETE CASCADE,
    source          TEXT    NOT NULL,   -- 'KEGG', 'COG', 'GO', 'Pfam', 'eggNOG', 'AMRFinder', …
    code            TEXT,                -- KO number, COG code, GO:XXXXXXX, Pfam ID, …
    category        TEXT,                -- see header note
    name            TEXT,                -- short human-readable name
    description     TEXT,                -- longer description when available
    score           REAL,                -- bit score / e-value / confidence
    extra           TEXT,                -- JSON blob for source-specific fields
    source_file     TEXT,                -- annotator's output file path (for provenance)
    tool_format     TEXT,                -- importer format: 'emapper' / 'kofamscan' / 'deepnog' / …
    applied_at      TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX IF NOT EXISTS idx_cluster_ann_cluster  ON cluster_annotations(cluster_id);
CREATE INDEX IF NOT EXISTS idx_cluster_ann_source   ON cluster_annotations(source);
CREATE INDEX IF NOT EXISTS idx_cluster_ann_code     ON cluster_annotations(code);
CREATE INDEX IF NOT EXISTS idx_cluster_ann_category ON cluster_annotations(category);

-- ─────────────────────────────────────────────────────────────────────────────
-- Views: convenience queries for the consumer side.
-- ─────────────────────────────────────────────────────────────────────────────

-- One row per (cluster_run, cluster, strain) cell, packing CDS hits into a
-- human-readable summary. This is what `4_maketable.py` used to produce as
-- tab-delimited rows; now it's a live query.
CREATE VIEW IF NOT EXISTS v_cluster_cell AS
SELECT
    m.cluster_run_id                                         AS cluster_run_id,
    c.cluster_id                                             AS cluster_id,
    c.cluster_name                                           AS cluster_name,
    s.strain_id                                              AS strain_id,
    s.locus_prefix                                           AS locus_prefix,
    COUNT(m.cds_id)                                          AS hit_count,
    GROUP_CONCAT(cds.locus_tag || '|' || COALESCE(cds.protein_id, ''), ', ')
                                                             AS hit_summary
FROM cluster_membership m
JOIN cds      ON cds.cds_id    = m.cds_id
JOIN strains  s ON s.strain_id = cds.strain_id
JOIN clusters c ON c.cluster_id = m.cluster_id
GROUP BY m.cluster_run_id, c.cluster_id, s.strain_id;
