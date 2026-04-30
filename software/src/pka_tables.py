"""IPC 2.0 pKa sets (Kozlowski 2021).

Two contexts:

- **Peptide set**: short, fully solvent-exposed residues. Used for peptide-mode
  inputs and for CDR3 sequences regardless of full-chain availability.
- **Protein set**: residues in a folded globular-domain context. Used for
  reconstructed full VH / VL chains.

Source: Kozlowski LP. *IPC 2.0: prediction of isoelectric point and pKa
dissociation constants.* Nucleic Acids Research 49(W1):W285-W292 (2021).

NOTE FOR IMPLEMENTORS: these constants were transcribed from the published
IPC 2.0 supplementary data. **Verify each value against the paper before
shipping a release** — wrong pKa values will quietly produce systematically
biased pI / charge values.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class PKaSet:
    """A complete pKa set: ionizable side chains plus terminal pKa values."""

    name: str
    side_chain: dict[str, float]  # AA → pKa for ionizable side chains
    n_terminus: float  # pKa of the free α-amino terminus
    c_terminus: float  # pKa of the free α-carboxyl terminus


# Acidic residues — deprotonate to negative charge as pH rises.
ACIDIC_AAS = frozenset("DECY")

# Basic residues — protonate to positive charge as pH falls.
BASIC_AAS = frozenset("HKR")


# IPC 2.0 peptide pKa set.
IPC2_PEPTIDE = PKaSet(
    name="IPC2_peptide",
    side_chain={
        "C": 7.555,  # Cys
        "D": 3.872,  # Asp
        "E": 4.412,  # Glu
        "H": 5.637,  # His
        "K": 9.052,  # Lys
        "R": 11.84,  # Arg
        "Y": 10.85,  # Tyr
    },
    n_terminus=9.094,
    c_terminus=2.869,
)


# IPC 2.0 protein pKa set.
IPC2_PROTEIN = PKaSet(
    name="IPC2_protein",
    side_chain={
        "C": 8.578,
        "D": 3.887,
        "E": 4.317,
        "H": 6.018,
        "K": 10.517,
        "R": 12.503,
        "Y": 10.071,
    },
    n_terminus=9.564,
    c_terminus=2.383,
)
