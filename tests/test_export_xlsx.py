"""Tests for the xlsx export.

We don't try to verify visual formatting (colors, frozen panes, etc.) — that
would require comparing against a snapshot file and is brittle. Instead we
verify the structural shape of the exported workbook by reading it back
with pandas:

    * Correct dimensions (clusters × (metadata cols + strain cols))
    * Metadata columns present in the right places
    * Strain columns named with the locus_prefix-prefixed full header
    * Cells filled with '[N]|...' for present hits and '*' for absent
"""

from __future__ import annotations

import shutil

import pytest

pytest.importorskip("Bio")
pytest.importorskip("pandas")
pytest.importorskip("openpyxl")

import pandas as pd

from strainbench.core import db as core_db
from strainbench.producer.cluster import cluster_strains
from strainbench.producer.export_xlsx import (
    export_cluster_table_xlsx,
    export_strain_anchored_xlsx,
)
from strainbench.producer.ingest import ingest_strain
from strainbench.producer.parsers.gbff import parse_gbff


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs2 not on PATH")
def test_export_produces_readable_xlsx(tmp_path, synthetic_gbff):
    """End-to-end: ingest synthetic gbff → cluster → export → read back."""
    db_path = tmp_path / "test.db"
    fasta_dir = tmp_path / "fasta"
    out_path = tmp_path / "cluster_table.xlsx"

    core_db.initialize(db_path)
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(db_path) as conn:
        ingest_strain(conn, record, fasta_dir)
        cluster_strains(conn, fasta_dir, sequence_type="protein",
                        pct_identity=80, coverage=80,
                        name_prefix="TEST", label="export-test")

    summary = export_cluster_table_xlsx(db_path, out_path)
    assert out_path.exists()
    assert summary.cluster_count >= 1
    assert summary.strain_count == 1

    # Read back. Header is row 0; pandas treats it as column names.
    df = pd.read_excel(out_path, sheet_name="cluster_table")
    # Block 1: core metadata (always 10 columns, ending at strain count)
    expected_core = [
        "original_order", "CLUSTER", "NCBI gene names", "NCBI annotation",
        "GC%", "GC spread", "aa length", "flags", "total count", "strain count",
    ]
    assert list(df.columns[: len(expected_core)]) == expected_core

    # Block 2: exactly summary.strain_count strain header columns immediately after
    strain_cols = list(df.columns[len(expected_core): len(expected_core) + summary.strain_count])
    assert len(strain_cols) == summary.strain_count
    assert any(col.startswith("TESTPFX | ") for col in strain_cols)

    # Block 3: any annotation columns (KEGG_KO etc.) live at the far right.
    # The synthetic gbff has no annotations imported, so this block is empty.
    annotation_cols = list(df.columns[len(expected_core) + summary.strain_count:])
    assert annotation_cols == []

    # Cluster rows match the cluster count
    assert len(df) == summary.cluster_count

    # Cell format: every value in strain columns is either '*' or starts with '['
    for col in strain_cols:
        for value in df[col]:
            assert value == "*" or str(value).startswith("[")


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs2 not on PATH")
def test_strain_anchored_export(tmp_path, synthetic_gbff):
    """End-to-end: ingest → cluster → export-strain-table → read back."""
    db_path = tmp_path / "test.db"
    fasta_dir = tmp_path / "fasta"
    out_path = tmp_path / "anchor.xlsx"

    core_db.initialize(db_path)
    record = parse_gbff(synthetic_gbff)
    with core_db.connect(db_path) as conn:
        ingest_strain(conn, record, fasta_dir)
        cluster_strains(conn, fasta_dir, sequence_type="protein",
                        pct_identity=80, coverage=80,
                        name_prefix="T", label="anchor-test")

    summary = export_strain_anchored_xlsx(db_path, out_path, "TESTPFX")
    assert out_path.exists()
    # The synthetic gbff has 3 translatable CDSs → 3 rows.
    assert summary.cluster_count == 3
    assert summary.strain_count == 1

    df = pd.read_excel(out_path, sheet_name="TESTPFX_anchored")
    expected_meta = [
        "position", "anchor_locus_tag", "nuc_accession", "location", "gene",
        "CLUSTER", "NCBI gene names", "NCBI annotation", "aa length",
        "cluster total", "cluster strains",
    ]
    assert list(df.columns[: len(expected_meta)]) == expected_meta
    # 3 rows, one per anchor CDS, in genomic order
    assert len(df) == 3
    assert list(df["position"]) == [1, 2, 3]
    # locus_tags are the synthetic ones, in cds_id order (= ingest order)
    assert all(tag.startswith("TESTPFX_RS") for tag in df["anchor_locus_tag"])

    # Anchor's column should be the FIRST strain column and start with the ★ marker
    strain_cols = list(df.columns[len(expected_meta):])
    assert strain_cols[0].startswith("★ TESTPFX | ")

    # Anchor cells should be '[1]|<that row's locus_tag>|...' (NOT the cluster summary)
    anchor_col = strain_cols[0]
    for _, row in df.iterrows():
        assert row[anchor_col].startswith(f"[1]|{row['anchor_locus_tag']}|")


def test_strain_anchored_nt_subcluster_column_present_when_nt_run_exists(tmp_path):
    """The nt_subcluster column should appear when a nucleotide cluster_run
    exists in the DB. We exercise the writer directly with hand-built data
    instead of running mmseqs nt clustering — the synthetic fixture has too
    few sequences for mmseqs's nucleotide linclust pre-step to succeed."""
    db_path = tmp_path / "test.db"
    out_path = tmp_path / "anchor.xlsx"
    core_db.initialize(db_path)

    with core_db.connect(db_path) as conn:
        # Minimal: one strain, two CDS, two cluster_runs (protein + nt),
        # one cluster per run, with both CDS as members.
        conn.execute(
            "INSERT INTO strains (locus_prefix, species, strain_name, assembly_id, "
            "biosample_id, bioproject_id, assembly_level, source_format, source_file) "
            "VALUES ('TESTP', 'Test bact', 'TESTSTRAIN', 'GCF_X', '?', '?', 'Complete', 'gbff', '?')"
        )
        sid = conn.execute("SELECT strain_id FROM strains").fetchone()["strain_id"]
        for i, lt in enumerate(("TESTP_001", "TESTP_002"), start=1):
            conn.execute(
                "INSERT INTO cds (strain_id, locus_tag, location, direction, "
                "nt_length, aa_length, gc_pct, annotation) "
                "VALUES (?, ?, ?, 'F', 300, 100, 35.0, 'test protein')",
                (sid, lt, f"{i*100}..{i*100+300}"),
            )
        # Protein run
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('p', 'protein', 80, 80, 'P')"
        )
        prot_run = conn.execute("SELECT MAX(cluster_run_id) FROM cluster_runs").fetchone()[0]
        conn.execute(
            "INSERT INTO clusters (cluster_run_id, cluster_name, member_count, strain_count) "
            "VALUES (?, 'P_000001', 2, 1)",
            (prot_run,),
        )
        prot_cid = conn.execute("SELECT cluster_id FROM clusters WHERE cluster_name='P_000001'").fetchone()[0]
        # Nucleotide run
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('n', 'nucleotide', 95, 90, 'N')"
        )
        nt_run = conn.execute("SELECT MAX(cluster_run_id) FROM cluster_runs").fetchone()[0]
        conn.execute(
            "INSERT INTO clusters (cluster_run_id, cluster_name, member_count, strain_count) "
            "VALUES (?, 'N_000001', 2, 1)",
            (nt_run,),
        )
        nt_cid = conn.execute("SELECT cluster_id FROM clusters WHERE cluster_name='N_000001'").fetchone()[0]
        for cid in conn.execute("SELECT cds_id FROM cds"):
            conn.execute(
                "INSERT INTO cluster_membership (cluster_run_id, cluster_id, cds_id) "
                "VALUES (?, ?, ?)",
                (prot_run, prot_cid, cid["cds_id"]),
            )
            conn.execute(
                "INSERT INTO cluster_membership (cluster_run_id, cluster_id, cds_id) "
                "VALUES (?, ?, ?)",
                (nt_run, nt_cid, cid["cds_id"]),
            )
        conn.commit()

    summary = export_strain_anchored_xlsx(db_path, out_path, "TESTP")
    assert summary.cluster_count == 2

    df = pd.read_excel(out_path, sheet_name="TESTP_anchored")
    expected_meta = [
        "position", "anchor_locus_tag", "nuc_accession", "location", "gene",
        "CLUSTER", "nt_subcluster", "NCBI gene names", "NCBI annotation",
        "aa length", "cluster total", "cluster strains",
    ]
    assert list(df.columns[: len(expected_meta)]) == expected_meta
    assert all(v == "N_000001" for v in df["nt_subcluster"])


def test_strain_anchored_unknown_strain_raises(tmp_path):
    db_path = tmp_path / "empty.db"
    out_path = tmp_path / "out.xlsx"
    core_db.initialize(db_path)
    # Stub a cluster_run so we get past _resolve_cluster_run_id
    with core_db.connect(db_path) as conn:
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('stub', 'protein', 80, 80, 'X')"
        )
        conn.commit()
    with pytest.raises(ValueError, match="not found"):
        export_strain_anchored_xlsx(db_path, out_path, "NOPE")


def test_cluster_table_pools_ncbi_gene_names(tmp_path):
    """Pool all distinct /gene qualifiers from a cluster's members into a single
    cell. Matches strain-comp's get_gene_prot_names() behavior:
      - 1 member has 'lacZ', others NULL  → cell shows 'lacZ'
      - 2 members have different names    → 'lacZ, lacZ2' (sorted)
      - no member has any /gene           → cell shows 'none'
    """
    db_path = tmp_path / "t.db"
    out_path = tmp_path / "ct.xlsx"
    core_db.initialize(db_path)

    with core_db.connect(db_path) as conn:
        conn.execute(
            "INSERT INTO strains (locus_prefix, species, strain_name, assembly_id, "
            "source_format, assembly_level) VALUES ('PFX', 'Test', 'X', 'GCF_X', 'gbff', 'Complete')"
        )
        sid = conn.execute("SELECT strain_id FROM strains").fetchone()["strain_id"]

        # Three clusters demonstrating the three cases:
        #   CL_001 → one member with 'lacZ', one NULL       (single name wins)
        #   CL_002 → two members, 'dnaA' and 'dnaA2'        (concat, sorted)
        #   CL_003 → both members NULL                      (rendered as 'none')
        plan = {
            "CL_001": [("m_1a", "lacZ"), ("m_1b", None)],
            "CL_002": [("m_2a", "dnaA"), ("m_2b", "dnaA2")],
            "CL_003": [("m_3a", None), ("m_3b", None)],
        }
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('t', 'protein', 80, 80, 'CL')"
        )
        run = conn.execute("SELECT MAX(cluster_run_id) FROM cluster_runs").fetchone()[0]
        for cname, members in plan.items():
            conn.execute(
                "INSERT INTO clusters (cluster_run_id, cluster_name, member_count, strain_count) "
                "VALUES (?, ?, ?, 1)",
                (run, cname, len(members)),
            )
            cid = conn.execute("SELECT cluster_id FROM clusters WHERE cluster_name=?", (cname,)).fetchone()[0]
            for locus, gene in members:
                conn.execute(
                    "INSERT INTO cds (strain_id, locus_tag, direction, aa_length, gene_name) "
                    "VALUES (?, ?, 'F', 100, ?)",
                    (sid, locus, gene),
                )
                cds_id = conn.execute("SELECT cds_id FROM cds WHERE locus_tag=?", (locus,)).fetchone()[0]
                conn.execute(
                    "INSERT INTO cluster_membership (cluster_run_id, cluster_id, cds_id) "
                    "VALUES (?, ?, ?)",
                    (run, cid, cds_id),
                )
        conn.commit()

    export_cluster_table_xlsx(db_path, out_path)
    df = pd.read_excel(out_path, sheet_name="cluster_table")

    by_cluster = dict(zip(df["CLUSTER"], df["NCBI gene names"]))
    assert by_cluster["CL_001"] == "lacZ"
    # Multi-value: order is whatever GROUP_CONCAT returned; accept either order.
    assert set(by_cluster["CL_002"].split(", ")) == {"dnaA", "dnaA2"}
    assert by_cluster["CL_003"] == "none"


def test_kegg_ko_dedup_and_gene_name_extraction(tmp_path):
    """The three transformations on the KEGG_KO column:
       1. Same K-number from multiple rows → single entry per cluster
       2. K-number with e-value annotation wins the merge
       3. Gene symbols (eggNOG's Preferred_name) move to NCBI gene names
    """
    db_path = tmp_path / "t.db"
    out_path = tmp_path / "ct.xlsx"
    core_db.initialize(db_path)
    with core_db.connect(db_path) as conn:
        conn.execute(
            "INSERT INTO strains (locus_prefix, source_format, assembly_level) "
            "VALUES ('PFX', 'gbff', 'Complete')"
        )
        sid = conn.execute("SELECT strain_id FROM strains").fetchone()["strain_id"]
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('t', 'protein', 80, 80, 'CL')"
        )
        run = conn.execute("SELECT MAX(cluster_run_id) FROM cluster_runs").fetchone()[0]
        conn.execute(
            "INSERT INTO clusters (cluster_run_id, cluster_name, member_count, strain_count) "
            "VALUES (?, 'CL_001', 1, 1)", (run,)
        )
        cid = conn.execute("SELECT cluster_id FROM clusters").fetchone()[0]
        conn.execute(
            "INSERT INTO cds (strain_id, locus_tag, direction, aa_length) "
            "VALUES (?, 'm_a', 'F', 100)", (sid,)
        )
        cds_id = conn.execute("SELECT cds_id FROM cds").fetchone()[0]
        conn.execute(
            "INSERT INTO cluster_membership (cluster_run_id, cluster_id, cds_id) "
            "VALUES (?, ?, ?)", (run, cid, cds_id),
        )
        # Same KO from emapper (with gene name 'dnaE') and kofamscan (no name)
        conn.execute(
            "INSERT INTO cluster_annotations (cluster_id, source, code, name, tool_format) "
            "VALUES (?, 'KEGG', 'K02337', 'dnaE', 'emapper')", (cid,)
        )
        conn.execute(
            "INSERT INTO cluster_annotations (cluster_id, source, code, name, tool_format) "
            "VALUES (?, 'KEGG', 'K02337', NULL, 'kofamscan')", (cid,)
        )
        # A second KO appearing twice — once bare, once with e-value annotation.
        # The e-value form should win the merge.
        conn.execute(
            "INSERT INTO cluster_annotations (cluster_id, source, code, name, tool_format) "
            "VALUES (?, 'KEGG', 'K00001', NULL, 'emapper')", (cid,)
        )
        conn.execute(
            "INSERT INTO cluster_annotations (cluster_id, source, code, name, extra, tool_format) "
            "VALUES (?, 'KEGG', 'K00001', NULL, '{\"below_threshold\": true, \"evalue\": 1e-50}', 'kofamscan')",
            (cid,),
        )
        conn.commit()

    export_cluster_table_xlsx(db_path, out_path)
    df = pd.read_excel(out_path, sheet_name="cluster_table")
    row = df[df["CLUSTER"] == "CL_001"].iloc[0]

    # (1) + (2): K02337 deduped to one entry (no gene-name suffix), K00001
    # merged with e-value form winning. Order is by K-number sort.
    kegg_cell = row["KEGG_KO"]
    assert "K02337" in kegg_cell
    # Bare K02337 — no 'dnaE' suffix anymore (moved to gene names column).
    assert "dnaE" not in kegg_cell
    # K00001 must carry the e-value
    assert "K00001 (e:" in kegg_cell
    # Each KO appears exactly once in the cell
    assert kegg_cell.count("K02337") == 1
    assert kegg_cell.count("K00001") == 1

    # (3): dnaE shows up in NCBI gene names, set-deduplicated
    gene_cell = row["NCBI gene names"]
    assert "dnaE" in gene_cell


def test_pfam_no_longer_duplicates_id_in_cell(tmp_path):
    """Regression: parser used to set name=code for Pfam → 'PF00001 PF00001'
    cells. Now name=None so the cell shows just 'PF00001'."""
    from strainbench.producer.annotations import parse_emapper
    import textwrap as _tw
    src = _tw.dedent("""\
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        CL_001\t.\t1e-50\t300\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tRibonuc_red_sm
    """)
    p = tmp_path / "e.tsv"
    p.write_text(src)
    pfam_rows = [r for r in parse_emapper(p) if r.source == "Pfam"]
    assert len(pfam_rows) == 1
    assert pfam_rows[0].code == "Ribonuc_red_sm"
    assert pfam_rows[0].name is None   # was the bug — used to equal the code


def test_annotation_columns_appear_after_strain_block(tmp_path):
    """Layout: core metadata → strain columns → annotation columns (far right).
    Frozen pane sits at end of core metadata (col index 10) so the strain
    block + annotations scroll freely."""
    db_path = tmp_path / "t.db"
    out_path = tmp_path / "ct.xlsx"
    core_db.initialize(db_path)
    with core_db.connect(db_path) as conn:
        for prefix in ("ALPHA", "BETA"):
            conn.execute(
                "INSERT INTO strains (locus_prefix, source_format, assembly_level) "
                "VALUES (?, 'gbff', 'Complete')", (prefix,)
            )
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('t', 'protein', 80, 80, 'CL')"
        )
        run = conn.execute("SELECT MAX(cluster_run_id) FROM cluster_runs").fetchone()[0]
        conn.execute(
            "INSERT INTO clusters (cluster_run_id, cluster_name, member_count, strain_count) "
            "VALUES (?, 'CL_001', 1, 1)", (run,)
        )
        cid = conn.execute("SELECT cluster_id FROM clusters").fetchone()[0]
        # Give the cluster a KEGG KO → this triggers a KEGG_KO annotation column
        conn.execute(
            "INSERT INTO cluster_annotations (cluster_id, source, code, tool_format) "
            "VALUES (?, 'KEGG', 'K00001', 'kofamscan')", (cid,),
        )
        # One member CDS so the cluster appears in v_cluster_cell
        sid = conn.execute("SELECT strain_id FROM strains WHERE locus_prefix='ALPHA'").fetchone()[0]
        conn.execute(
            "INSERT INTO cds (strain_id, locus_tag, direction, aa_length) "
            "VALUES (?, 'A_01', 'F', 100)", (sid,)
        )
        cds_id = conn.execute("SELECT cds_id FROM cds").fetchone()[0]
        conn.execute(
            "INSERT INTO cluster_membership (cluster_run_id, cluster_id, cds_id) "
            "VALUES (?, ?, ?)", (run, cid, cds_id)
        )
        conn.commit()

    export_cluster_table_xlsx(db_path, out_path)
    df = pd.read_excel(out_path, sheet_name="cluster_table")

    # Layout: 10 core metadata → 2 strain cols (ALPHA, BETA) → 1 annotation col (KEGG_KO)
    assert list(df.columns[:10]) == [
        "original_order", "CLUSTER", "NCBI gene names", "NCBI annotation",
        "GC%", "GC spread", "aa length", "flags", "total count", "strain count",
    ]
    # Strain columns at positions 10 and 11
    assert df.columns[10].startswith("ALPHA | ")
    assert df.columns[11].startswith("BETA | ")
    # Annotation column at the far right
    assert df.columns[12] == "KEGG_KO"
    assert len(df.columns) == 13  # no more, no less


def test_strain_header_dedupes_when_organism_contains_strain():
    """NCBI bakes the strain into the organism name for unspeciated isolates
    (e.g. 'Gardnerella sp. DNF01159') AND also reports it in /strain. We
    must not produce 'Gardnerella sp. DNF01159 DNF01159' in the header."""
    from strainbench.producer.export_xlsx import _format_strain_header
    import pandas as _pd

    # 1. Strain baked into organism name → no duplication
    row = _pd.Series({
        "locus_prefix": "HXT39",
        "species": "Gardnerella sp. DNF01159",
        "strain_name": "DNF01159",
        "assembly_id": "GCF_049244215.1",
        "biosample_id": "SAMN15393823",
        "bioproject_id": "PRJNA639145",
        "assembly_level": "Complete",
    })
    h = _format_strain_header(row)
    # Should appear once, not twice
    assert h.count("DNF01159") == 1
    # And should still include the assembly metadata
    assert "GCF_049244215.1" in h
    assert "Complete" in h

    # 2. Same identity case
    row_same = _pd.Series({
        "locus_prefix": "X",
        "species": "DNF01159",
        "strain_name": "DNF01159",
        "assembly_id": "GCF_X",
        "biosample_id": "S",
        "bioproject_id": "P",
        "assembly_level": "Contig",
    })
    h_same = _format_strain_header(row_same)
    assert h_same.count("DNF01159") == 1

    # 3. Distinct strain (e.g., classic format) → both appear once each
    row_normal = _pd.Series({
        "locus_prefix": "Y",
        "species": "Lactobacillus iners",
        "strain_name": "AB-1",
        "assembly_id": "GCF_000177755.1",
        "biosample_id": "SAMN02470230",
        "bioproject_id": "PRJNA43549",
        "assembly_level": "Contig",
    })
    h_normal = _format_strain_header(row_normal)
    assert "Lactobacillus iners AB-1" in h_normal


def test_export_raises_when_no_cluster_run(tmp_path):
    """Trying to export from a DB that has no cluster_runs yet should fail clearly."""
    db_path = tmp_path / "empty.db"
    out_path = tmp_path / "out.xlsx"
    core_db.initialize(db_path)
    with pytest.raises(ValueError, match="No active cluster_runs"):
        export_cluster_table_xlsx(db_path, out_path)
