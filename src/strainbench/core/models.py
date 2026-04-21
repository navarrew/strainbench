"""The contract layer: canonical in-memory representations of strains and CDSs.

Every input parser produces instances of these dataclasses, and every
downstream step (ingestion into the DB, clustering, export) consumes them.
If a parser honors this contract and a consumer honors this contract, they
can be swapped independently.

Sequences (`nt_sequence`, `aa_sequence`) live on `CDSRecord` during in-memory
handoff from parsers to ingestion, but they are NOT stored inline in the
database — they are written to sidecar FASTA files. The DB row only stores
metadata and length.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class CDSRecord:
    """One coding sequence. The atomic unit of the strainbench data model."""

    locus_tag: str
    """Per-strain unique identifier, e.g. 'JMT79_RS10105'."""

    nuc_accession: str
    """NCBI accession of the replicon (chromosome or plasmid) this CDS sits on."""

    location: str
    """Feature location string, e.g. 'complement(2043708..2044166)'."""

    direction: str
    """Strand: 'F' (forward) or 'R' (reverse/complement)."""

    nt_sequence: str
    """The nucleotide sequence of this CDS. Stored in sidecar FASTA, not the DB."""

    aa_sequence: str
    """The translated protein sequence. Stored in sidecar FASTA, not the DB."""

    nt_length: int
    aa_length: int
    gc_pct: float
    annotation: str
    """The /product qualifier, or equivalent free-text description."""

    protein_id: str | None = None
    """NCBI protein accession, e.g. 'WP_000502119.1'. May be absent for custom genomes."""

    gene_name: str | None = None
    """The /gene qualifier if present (e.g. 'dnaA')."""

    notes: str | None = None
    """Flags such as 'truncated_5_prime', 'frameshift', 'crosses_origin'."""


@dataclass
class StrainRecord:
    """One genome. Produced by a parser, consumed by the ingestion step."""

    locus_prefix: str
    """The shared prefix of this strain's locus tags, e.g. 'JMT79'. Unique per DB."""

    species: str
    strain_name: str
    assembly_id: str
    """Assembly accession, e.g. 'GCF_022662295.1'."""

    assembly_level: str
    """One of: 'Complete', 'Chromosome', 'Scaffold', 'Contig', 'Unknown'."""

    source_format: str
    """One of: 'gbff', 'ncbi_cds', 'other'."""

    source_file: str
    """Absolute or relative path to the original input file."""

    cds_records: list[CDSRecord] = field(default_factory=list)

    biosample_id: str | None = None
    bioproject_id: str | None = None
