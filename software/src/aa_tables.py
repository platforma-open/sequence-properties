"""Amino-acid lookup tables.

Property functions delegate hydropathy, mass, extinction-coefficient, and
DIWV tables to BioPython's `Bio.SeqUtils.ProtParam` per the spec M1
strategy. The 20 standard single-letter codes stay here because the
sequence-cleanup layer (`properties._prepare`) needs to filter ambiguity
codes before any BioPython call — BioPython rejects non-standard residues.
"""

from __future__ import annotations

# 20 standard single-letter codes used everywhere in this codebase.
STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"
STANDARD_AA_SET = frozenset(STANDARD_AAS)
