"""Tests for the annotation import layer."""

from __future__ import annotations

import sqlite3
import textwrap

import pytest

from strainbench.core import db as core_db
from strainbench.producer.annotations import (
    detect_format,
    import_annotation_file,
    import_workdir,
    parse_amrfinder,
    parse_deepnog,
    parse_defensefinder,
    parse_emapper,
    parse_interproscan,
    parse_kofamscan,
    parse_padloc,
)


# ─────────────────────────────────────────────────────────────────────────────
# Parser unit tests (no DB needed)
# ─────────────────────────────────────────────────────────────────────────────


def test_parse_emapper_yields_multiple_sources_per_row(tmp_path):
    content = textwrap.dedent("""\
        ##
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        INERS_000001\t123.foo\t1e-50\t300.5\tCOG0001@2|Bacteria\tBacteria\tJ\tAminotransferase\trecA\tGO:0003677,GO:0005524\t2.6.1.1\tko:K00811\tmap00260\t-\t-\t-\t-\t-\t-\t-\tPF00155
        """)
    path = tmp_path / "emapper.tsv"
    path.write_text(content)
    rows = list(parse_emapper(path))
    sources = {r.source for r in rows}
    # We expect all the non-empty fields to produce rows.
    assert "COG" in sources
    assert "GO" in sources
    assert "KEGG" in sources
    assert "EC" in sources
    assert "Pfam" in sources

    # GO terms split into separate rows
    go_rows = [r for r in rows if r.source == "GO"]
    assert len(go_rows) == 2
    assert {r.code for r in go_rows} == {"GO:0003677", "GO:0005524"}

    # KEGG: one KO row + one pathway row
    kegg = [r for r in rows if r.source == "KEGG"]
    assert {r.code for r in kegg} == {"K00811", "map00260"}

    # COG: one row per letter
    cog = [r for r in rows if r.source == "COG"]
    assert len(cog) == 1
    assert cog[0].code == "J"
    assert cog[0].category == "J"
    assert cog[0].name and "Translation" in cog[0].name

    # GO rows have aspect populated for terms in our shipped lookup
    dna_binding = next(r for r in go_rows if r.code == "GO:0003677")
    assert dna_binding.category == "MF"


def test_parse_kofamscan_keeps_assigned_plus_best_fallback(tmp_path):
    """kofamscan parser strategy:
       - Clusters with assigned ('*') hits → emit ALL assigned hits.
       - Clusters with no '*' hit → emit the single best candidate (highest score),
         tagged below_threshold in extra so consumers can flag it visually.
    """
    content = textwrap.dedent("""\
        # gene name\tKO\tthrshld\tscore\tE-value\t"KO definition"
        * INERS_000001\tK00001\t256.30\t312.5\t1.2e-90\t"alcohol dehydrogenase"
        * INERS_000001\tK00121\t260.00\t300.0\t1.0e-85\t"second assigned hit"
          INERS_000002\tK99999\t180.00\t100.0\t1.0e-30\t"weaker candidate"
          INERS_000002\tK88888\t180.00\t150.0\t1.0e-50\t"better candidate"
          INERS_000003\tK11111\t100.00\t70.0\t1.0e-20\t"single candidate, no '*'"
        """)
    path = tmp_path / "kofam.tsv"
    path.write_text(content)
    rows = list(parse_kofamscan(path))

    by_cluster = {r.cluster_name: [] for r in rows}
    for r in rows:
        by_cluster[r.cluster_name].append(r)

    # Cluster 1 has 2 assigned hits → both kept, no below_threshold flag.
    # (extra holds the evalue but should NOT contain below_threshold)
    assert len(by_cluster["INERS_000001"]) == 2
    assert {r.code for r in by_cluster["INERS_000001"]} == {"K00001", "K00121"}
    import json
    for r in by_cluster["INERS_000001"]:
        extra = json.loads(r.extra) if r.extra else {}
        assert "below_threshold" not in extra
        assert "evalue" in extra

    # Cluster 2 has 2 candidates, no assigned → best (K88888) kept with flag.
    assert len(by_cluster["INERS_000002"]) == 1
    fallback2 = by_cluster["INERS_000002"][0]
    assert fallback2.code == "K88888"
    import json
    extra2 = json.loads(fallback2.extra)
    assert extra2["below_threshold"] is True
    # E-value of the fallback hit should be captured for display
    assert extra2["evalue"] == 1e-50

    # Cluster 3 has only 1 candidate → kept with flag and its e-value
    assert len(by_cluster["INERS_000003"]) == 1
    extra3 = json.loads(by_cluster["INERS_000003"][0].extra)
    assert extra3["below_threshold"] is True
    assert extra3["evalue"] == 1e-20


# ─────────────────────────────────────────────────────────────────────────────
# Integration: parse + insert + query
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def db_with_clusters(tmp_path):
    """Create a DB with two fake clusters so we can test import flow without mmseqs."""
    db_path = tmp_path / "test.db"
    core_db.initialize(db_path)
    with core_db.connect(db_path) as conn:
        conn.execute(
            "INSERT INTO cluster_runs (label, sequence_type, pct_identity, coverage, name_prefix) "
            "VALUES ('t', 'protein', 80, 80, 'INERS')"
        )
        run_id = conn.execute("SELECT MAX(cluster_run_id) FROM cluster_runs").fetchone()[0]
        for name in ("INERS_000001", "INERS_000002"):
            conn.execute(
                "INSERT INTO clusters (cluster_run_id, cluster_name) VALUES (?, ?)",
                (run_id, name),
            )
        conn.commit()
    return db_path


def test_import_emapper_writes_rows(tmp_path, db_with_clusters):
    content = textwrap.dedent("""\
        ##
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        INERS_000001\t.\t1e-50\t300.5\t-\t-\tJ\ttest\trecA\tGO:0003677\t-\tko:K00001\t-\t-\t-\t-\t-\t-\t-\t-\t-
        """)
    tsv = tmp_path / "e.tsv"
    tsv.write_text(content)
    summary = import_annotation_file(db_with_clusters, tsv, source_format="emapper")
    assert summary.rows_inserted >= 3   # COG + GO + KEGG
    assert "COG" in summary.sources_written
    assert "GO" in summary.sources_written
    assert "KEGG" in summary.sources_written
    assert summary.clusters_touched == 1
    assert summary.unmatched_clusters == 0

    with sqlite3.connect(db_with_clusters) as conn:
        conn.row_factory = sqlite3.Row
        rows = list(conn.execute(
            "SELECT source, code, category FROM cluster_annotations ORDER BY source, code"
        ))
    sources = [r["source"] for r in rows]
    assert "COG" in sources and "GO" in sources and "KEGG" in sources


def test_import_unmatched_clusters_reported(tmp_path, db_with_clusters):
    content = textwrap.dedent("""\
        # gene name\tKO\tthrshld\tscore\tE-value\t"KO definition"
        * INERS_999999\tK00001\t256.30\t312.5\t1.2e-90\t"fake cluster not in DB"
        * INERS_000001\tK00002\t256.30\t312.5\t1.2e-90\t"real one"
        """)
    tsv = tmp_path / "k.tsv"
    tsv.write_text(content)
    summary = import_annotation_file(db_with_clusters, tsv, source_format="kofamscan")
    assert summary.rows_inserted == 1
    assert summary.unmatched_clusters == 1


def test_parse_deepnog_yields_cog_groups(tmp_path):
    content = textwrap.dedent("""\
        sequence_id\tprediction\tconfidence
        INERS_000001\tCOG1234\t0.89
        INERS_000002\tCOG0150\t0.95
        INERS_000003\t-\t0.0
    """)
    path = tmp_path / "deepnog.tsv"
    path.write_text(content)
    rows = list(parse_deepnog(path))
    # '-' predictions are skipped.
    assert len(rows) == 2
    for row in rows:
        assert row.source == "COG"
        assert row.code.startswith("COG")
        assert row.category == "group"
        assert row.extra and "deepnog" in row.extra


def test_detect_format_identifies_all_three(tmp_path):
    emapper = tmp_path / "e.tsv"
    emapper.write_text("## header comment\n#query\tseed\tscore\nX\t1\t2\n")
    kofam = tmp_path / "k.tsv"
    kofam.write_text('# gene name\tKO\tthrshld\tscore\tE-value\t"KO definition"\n')
    deepnog = tmp_path / "d.tsv"
    deepnog.write_text("sequence_id\tprediction\tconfidence\nFOO\tCOG0001\t0.5\n")
    unknown = tmp_path / "u.tsv"
    unknown.write_text("no annotator header here\nrandom data\n")

    assert detect_format(emapper) == "emapper"
    assert detect_format(kofam) == "kofamscan"
    assert detect_format(deepnog) == "deepnog"
    assert detect_format(unknown) is None


def test_detect_format_identifies_new_four(tmp_path):
    amr = tmp_path / "a.tsv"
    amr.write_text(
        "Protein identifier\tContig id\tStart\tStop\tStrand\tGene symbol\t"
        "Sequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\t"
        "Method\tTarget length\tReference sequence length\t"
        "% Coverage of reference sequence\t% Identity to reference sequence\n"
    )
    df = tmp_path / "df.tsv"
    df.write_text(
        "replicon\thit_id\tgene_name\thit_pos\thit_status\thit_score\tsys_id\ttype\tsubtype\n"
    )
    pl = tmp_path / "pl.csv"
    pl.write_text(
        "system.number,seqid,system,target.name,hmm.accession,hmm.name,protein.name,"
        "full.seq.E.value,full.seq.score,target.coverage\n"
    )
    ips = tmp_path / "ips.tsv"
    # IPS has no header — first data row's column 4 is the analysis name
    ips.write_text(
        "INERS_001\tdeadbeefmd5\t300\tPfam\tPF00001\t7TM receptor\t1\t300\t1e-50\tT\t14-04-2026\n"
    )

    assert detect_format(amr) == "amrfinder"
    assert detect_format(df) == "defensefinder"
    assert detect_format(pl) == "padloc"
    assert detect_format(ips) == "interproscan"


def test_parse_amrfinder(tmp_path):
    p = tmp_path / "a.tsv"
    p.write_text(
        "Protein identifier\tContig id\tStart\tStop\tStrand\tGene symbol\t"
        "Sequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\t"
        "Method\tTarget length\tReference sequence length\t"
        "% Coverage of reference sequence\t% Identity to reference sequence\n"
        "INERS_001234\tNA\tNA\tNA\tNA\ttetM\ttetracycline resistance protein TetM\t"
        "core\tAMR\tAMR\tTETRACYCLINE\tTETRACYCLINE\tHMM\t639\t639\t100.00\t99.84\n"
    )
    rows = list(parse_amrfinder(p))
    assert len(rows) == 1
    r = rows[0]
    assert r.cluster_name == "INERS_001234"
    assert r.source == "AMRFinder"
    assert r.code == "tetM"
    assert r.category == "AMR"
    assert r.score == 99.84
    assert "tetracycline" in r.description.lower()


def test_parse_defensefinder(tmp_path):
    p = tmp_path / "df.tsv"
    p.write_text(
        "replicon\thit_id\tgene_name\thit_pos\thit_status\thit_score\tsys_id\ttype\tsubtype\n"
        "rep1\tINERS_001234\tCas9\t1\tmandatory\t500.5\tSYS_1\tCRISPR-Cas\tType_II-C\n"
        "rep1\tINERS_001235\tCas1\t2\tmandatory\t300.0\tSYS_1\tCRISPR-Cas\tType_II-C\n"
    )
    rows = list(parse_defensefinder(p))
    assert len(rows) == 2
    assert {r.cluster_name for r in rows} == {"INERS_001234", "INERS_001235"}
    assert all(r.source == "DefenseFinder" for r in rows)
    assert all(r.code == "CRISPR-Cas" for r in rows)
    assert all(r.category == "Type_II-C" for r in rows)
    cas9 = next(r for r in rows if r.cluster_name == "INERS_001234")
    assert cas9.name == "Cas9"


def test_parse_padloc(tmp_path):
    p = tmp_path / "pl.csv"
    p.write_text(
        "system.number,seqid,system,target.name,hmm.accession,hmm.name,protein.name,"
        "full.seq.E.value,full.seq.score,target.coverage,target.description\n"
        "1,contig1,CRISPR-Cas_TypeI-B,INERS_002000,PADLOC_HMM_001,Cas3,Cas3,"
        "1e-100,500.5,0.95,Cas3 helicase-nuclease\n"
    )
    rows = list(parse_padloc(p))
    assert len(rows) == 1
    r = rows[0]
    assert r.cluster_name == "INERS_002000"
    assert r.source == "PADLOC"
    assert r.code == "CRISPR-Cas_TypeI-B"
    assert r.category == "Cas3"
    assert r.score == 500.5


def test_parse_interproscan_with_go(tmp_path):
    p = tmp_path / "ips.tsv"
    # 15 columns including GO annotations
    p.write_text(
        "INERS_001\tmd5\t300\tPfam\tPF00001\t7TM receptor\t1\t300\t1e-50\tT\t"
        "14-04-2026\tIPR000001\t7-transmembrane\tGO:0005524|GO:0007165\t-\n"
        "INERS_001\tmd5\t300\tTIGRFAM\tTIGR00010\tDNA polymerase III\t1\t300\t"
        "1e-30\tT\t14-04-2026\t-\t-\t-\t-\n"
    )
    rows = list(parse_interproscan(p))
    # Two signature hits + two GO terms = 4 rows expected
    sources = [r.source for r in rows]
    assert "Pfam" in sources
    assert "TIGRFAM" in sources
    assert sources.count("GO") == 2

    pfam = next(r for r in rows if r.source == "Pfam")
    assert pfam.code == "PF00001"
    go_codes = {r.code for r in rows if r.source == "GO"}
    assert go_codes == {"GO:0005524", "GO:0007165"}
    # GO:0005524 ATP binding is in our shipped go_aspects.tsv → MF
    atp = next(r for r in rows if r.code == "GO:0005524")
    assert atp.category == "MF"


def test_import_workdir_scans_and_imports_all(tmp_path, db_with_clusters):
    workdir = tmp_path / "wd"
    # Drop one of each annotation format into subdirs, as annotate-prep would.
    (workdir / "emapper_out").mkdir(parents=True)
    (workdir / "emapper_out" / "test.emapper.annotations").write_text(textwrap.dedent("""\
        ##
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        INERS_000001\t.\t1e-50\t300.5\t-\t-\tJ\ttest\trecA\tGO:0003677\t-\tko:K00001\t-\t-\t-\t-\t-\t-\t-\t-\t-
    """))
    (workdir / "kofamscan_out").mkdir()
    (workdir / "kofamscan_out" / "detail.tsv").write_text(
        '# gene name\tKO\tthrshld\tscore\tE-value\t"KO definition"\n'
        '* INERS_000002\tK99999\t100\t200\t1e-50\t"test ko"\n'
    )
    (workdir / "deepnog_out").mkdir()
    (workdir / "deepnog_out" / "classify.tsv").write_text(
        "sequence_id\tprediction\tconfidence\n"
        "INERS_000001\tCOG1234\t0.9\n"
    )
    # Sprinkle an unrelated TXT file to confirm it's skipped gracefully
    (workdir / "readme.txt").write_text("just some notes\n")

    summary = import_workdir(db_with_clusters, workdir)
    assert summary.total_files_scanned >= 3
    assert summary.total_files_imported == 3
    assert summary.total_rows_inserted >= 5   # ≥1 from each of 3 files
    formats_seen = {fmt for _path, fmt, _result in summary.per_file if fmt != "unrecognized"}
    assert {"emapper", "kofamscan", "deepnog"}.issubset(formats_seen)


def test_import_workdir_records_provenance(tmp_path, db_with_clusters):
    workdir = tmp_path / "wd"
    workdir.mkdir()
    (workdir / "kofamscan_out").mkdir()
    kofam_path = workdir / "kofamscan_out" / "detail.tsv"
    kofam_path.write_text(
        '# gene name\tKO\tthrshld\tscore\tE-value\t"KO definition"\n'
        '* INERS_000001\tK11111\t100\t200\t1e-50\t"test"\n'
    )
    import_workdir(db_with_clusters, workdir)

    with sqlite3.connect(db_with_clusters) as conn:
        conn.row_factory = sqlite3.Row
        row = conn.execute(
            "SELECT source_file, tool_format FROM cluster_annotations LIMIT 1"
        ).fetchone()
    assert row["source_file"].endswith("detail.tsv")
    assert row["tool_format"] == "kofamscan"


def test_import_replace_deletes_old_rows(tmp_path, db_with_clusters):
    tsv = tmp_path / "k.tsv"
    tsv.write_text(textwrap.dedent("""\
        # header
        * INERS_000001\tK00001\t256.30\t312.5\t1.2e-90\t"first pass"
    """))
    import_annotation_file(db_with_clusters, tsv, source_format="kofamscan")

    tsv.write_text(textwrap.dedent("""\
        # header
        * INERS_000001\tK99999\t256.30\t312.5\t1.2e-90\t"second pass"
    """))
    summary = import_annotation_file(
        db_with_clusters, tsv, source_format="kofamscan", replace=True,
    )
    assert summary.rows_inserted == 1

    with sqlite3.connect(db_with_clusters) as conn:
        codes = [r[0] for r in conn.execute(
            "SELECT code FROM cluster_annotations WHERE source='KEGG'"
        )]
    assert codes == ["K99999"]   # old K00001 is gone


def test_replace_preserves_other_tools_rows_at_same_source(tmp_path, db_with_clusters):
    """Re-importing emapper with --replace shouldn't wipe DeepNOG's COG rows.

    Both tools produce rows with source='COG' but they're different granularities
    (letters vs groups). Scoping --replace by tool_format preserves cross-tool
    contributions to the same source.
    """
    # DeepNOG first: writes COG group rows
    dn = tmp_path / "d.tsv"
    dn.write_text("sequence_id\tprediction\tconfidence\nINERS_000001\tCOG1234\t0.9\n")
    import_annotation_file(db_with_clusters, dn, source_format="deepnog")

    # emapper with --replace: writes COG category (letter) rows
    em = tmp_path / "e.tsv"
    em.write_text(textwrap.dedent("""\
        ##
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        INERS_000001\t.\t1e-50\t300\t-\t-\tJ\ttest\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-
    """))
    import_annotation_file(db_with_clusters, em, source_format="emapper", replace=True)

    with sqlite3.connect(db_with_clusters) as conn:
        rows = list(conn.execute(
            "SELECT source, code, category, tool_format FROM cluster_annotations "
            "WHERE source='COG' ORDER BY tool_format"
        ))
    # Expect BOTH tools' rows to be present: DeepNOG's COG1234 AND emapper's J.
    tool_formats = {r[3] for r in rows}
    assert tool_formats == {"deepnog", "emapper"}
    deepnog_codes = [r[1] for r in rows if r[3] == "deepnog"]
    emapper_codes = [r[1] for r in rows if r[3] == "emapper"]
    assert deepnog_codes == ["COG1234"]
    assert emapper_codes == ["J"]
