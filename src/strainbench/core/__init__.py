"""Core contract and database layer — shared by the producer and web halves."""

from strainbench.core.models import CDSRecord, StrainRecord

__all__ = ["CDSRecord", "StrainRecord"]
