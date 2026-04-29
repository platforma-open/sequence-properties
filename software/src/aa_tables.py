"""Amino-acid lookup tables.

Standard 20 single-letter codes, residue properties, Kyte-Doolittle hydropathy,
average residue masses (condensed form), aromatic and aliphatic sets.

References:
- Kyte J, Doolittle RF. *J Mol Biol* 157:105-132 (1982).
- NIST / UniMod average residue masses (free amino acid mass minus 18.0153 for
  the water molecule lost during peptide bond formation).
"""

from __future__ import annotations

# 20 standard single-letter codes used everywhere in this codebase.
STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"
STANDARD_AA_SET = frozenset(STANDARD_AAS)

# Kyte-Doolittle hydropathy scale.
KD_SCALE: dict[str, float] = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

# Average residue masses (Da), condensed form: free amino acid mass − 18.0153.
# Source: NIST / UniMod amino acid average masses.
AVG_RESIDUE_MASS: dict[str, float] = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
}

# Average mass of one water molecule, Da. Added to the residue-mass sum to
# restore the terminal H₂O (N-terminal H + C-terminal OH).
H2O_AVG_MASS = 18.0153

# Aromatic residues — used for aromaticity fraction.
AROMATIC_AAS = frozenset("FWY")

# Pace et al. extinction coefficient parameters (M⁻¹·cm⁻¹ at 280 nm).
EC_TYR = 1490
EC_TRP = 5500
EC_DISULFIDE = 125  # per disulfide bond, i.e. per pair of Cys
