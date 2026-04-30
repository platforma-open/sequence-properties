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


# IPC 2.0 peptide pKa set.
# Source: Kozlowski 2021, Supplementary Table S1; cross-checked against
# http://ipc2-isoelectric-point.org/theory.html (paper-linked DOI 10.1093/nar/gkab295).
IPC2_PEPTIDE = PKaSet(
    name="IPC2_peptide",
    side_chain={
        "C": 9.439,  # Cys
        "D": 3.969,  # Asp
        "E": 4.507,  # Glu
        "H": 6.439,  # His
        "K": 8.165,  # Lys
        "R": 11.493,  # Arg
        "Y": 9.153,  # Tyr
    },
    n_terminus=7.947,
    c_terminus=2.977,
)


# IPC 2.0 protein pKa set.
# Source: Kozlowski 2021, Supplementary Table S1; cross-checked against
# http://ipc2-isoelectric-point.org/theory.html (paper-linked DOI 10.1093/nar/gkab295).
IPC2_PROTEIN = PKaSet(
    name="IPC2_protein",
    side_chain={
        "C": 7.890,
        "D": 3.766,
        "E": 4.497,
        "H": 5.492,
        "K": 9.247,
        "R": 10.223,
        "Y": 11.491,
    },
    n_terminus=5.779,
    c_terminus=6.065,
)
