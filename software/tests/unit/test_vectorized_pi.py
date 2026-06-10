"""Parity + determinism tests: vectorized isoelectric point vs the scalar oracle.

The vectorized pI (vectorized.isoelectric_point) is a lockstep, fixed-iteration
bisection over the row-wise charge function `_charge_raw`, mirroring the scalar
`properties._bisect_charge_zero` exactly:

* SAME bracket [0, 14], SAME tol 1e-3 -> SAME data-independent iteration count
  (ceil(log2(14/0.001)) == 14 halvings for every sequence);
* SAME branch test `(f_mid > 0) == (f_lo > 0)`;
* `_charge_raw` is bit-identical to the scalar charge (T7: 0.0 residual).

Therefore the vectorized pI is bit-identical to the scalar pI for every valid
sequence, modulo the scalar's astronomically-rare `if f_mid == 0.0: return mid`
early-out (a midpoint landing on EXACT zero charge), which the vectorized loop
ignores and instead converges to the same value within tol.

Two tests:

* PARITY (Hypothesis): vectorized within 1.5e-3 of scalar over generated
  sequences (valid + invalid), both pKa sets, both include_cys; NaN where scalar
  None; float64 dtype. (Expect ~0.0 residual; 1.5e-3 is headroom.)
* STRADDLE (deterministic, seeded): over the committed peptide corpus PLUS >=2000
  seeded random standard-AA sequences, under BOTH (IPC2_PEPTIDE, include_cys=True)
  and (IPC2_PROTEIN, include_cys=False): round(vectorized, 3) == round(scalar, 3)
  for EVERY cell where the scalar is not None. The straddle count must be exactly 0.
"""

from __future__ import annotations

import csv
import math
import random
from pathlib import Path

import numpy as np
from hypothesis import given
from hypothesis import strategies as st

import properties as scalar
from aa_tables import STANDARD_AAS
from pka_tables import IPC2_PEPTIDE, IPC2_PROTEIN
from vectorized import build_counts, isoelectric_point

# Mixed valid (standard-AA-only text) + invalid (empty, stop codon,
# non-standard-only) sequences, so every row exercises one of the branches.
seqs = st.lists(
    st.one_of(
        st.text(alphabet=STANDARD_AAS, min_size=1, max_size=50),
        st.sampled_from(["", "*", "BZXJ"]),
    ),
    min_size=1,
    max_size=25,
)

# pipeline configs: peptides use IPC2_PEPTIDE/include_cys=True;
# full chains use IPC2_PROTEIN/include_cys=False.
PKA = [IPC2_PEPTIDE, IPC2_PROTEIN]
CONFIGS = [(IPC2_PEPTIDE, True), (IPC2_PROTEIN, False)]

_CORPUS = Path(__file__).resolve().parents[1] / "data" / "characterization" / "peptide_input.tsv"


@given(seqs, st.sampled_from(PKA), st.booleans())
def test_pi_parity(seqs, pka, inc):
    sub = build_counts(seqs)
    vec = isoelectric_point(sub, pka, include_cys=inc)
    assert vec.dtype == np.float64
    assert vec.shape == (len(seqs),)
    for i, s in enumerate(seqs):
        exp = scalar.isoelectric_point(s, pka, include_cys=inc)
        if exp is None:
            assert math.isnan(vec[i])
        else:
            assert abs(vec[i] - exp) <= 1.5e-3


def _read_corpus_sequences() -> list[str]:
    with _CORPUS.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [row.get("sequence", "") or "" for row in reader]


def _seeded_random_sequences(rng: random.Random, n: int) -> list[str]:
    out: list[str] = []
    for _ in range(n):
        length = rng.randint(1, 60)
        out.append("".join(rng.choice(STANDARD_AAS) for _ in range(length)))
    return out


def test_pi_no_straddles():
    """The determinism gate: round(vectorized, 3) == round(scalar, 3) for every
    cell where the scalar is not None, across the committed corpus + >=2000
    seeded random sequences, under both pipeline pKa configs. Zero straddles.
    """
    rng = random.Random(0)
    sample = _read_corpus_sequences() + _seeded_random_sequences(rng, 2500)
    assert len(sample) >= 2000

    sub = build_counts(sample)

    straddles: list[tuple[str, str, float, float]] = []
    compared = 0
    for pka, inc in CONFIGS:
        vec = isoelectric_point(sub, pka, include_cys=inc)
        assert vec.dtype == np.float64
        for i, s in enumerate(sample):
            exp = scalar.isoelectric_point(s, pka, include_cys=inc)
            if exp is None:
                # scalar None must correspond to vectorized NaN.
                assert math.isnan(vec[i])
                continue
            compared += 1
            rv = round(float(vec[i]), 3)
            re = round(float(exp), 3)
            if rv != re:
                straddles.append((pka.name, s, float(vec[i]), float(exp)))

    assert compared >= 2000
    assert straddles == [], f"{len(straddles)} straddle(s) at 3 dp; first: {straddles[0]}"
