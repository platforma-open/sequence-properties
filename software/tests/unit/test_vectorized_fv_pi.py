"""Parity + determinism tests: vectorized Fv isoelectric point vs the scalar oracle.

The vectorized Fv pI (vectorized.fv_isoelectric_point) bisects the per-chain
charge SUM, mirroring the scalar `properties.fv_isoelectric_point` (and the
Fv-row path in `pipeline._compute_fv_row_from_ctx`) exactly:

* SAME bracket [0, 14], SAME tol 1e-3 -> SAME data-independent iteration count;
* SAME branch test `(f_mid > 0) == (f_lo > 0)`;
* bisected function is `_charge_raw(vh, ph) + _charge_raw(vl, ph)`, each
  bit-identical to the scalar per-chain charge.

The Fv path always uses the protein pKa set with Cys excluded (IPC2_PROTEIN,
include_cys=False), matching the full-chain rule.

Two tests:

* PARITY (Hypothesis): vectorized within 1.5e-3 of scalar over generated VH/VL
  pairs (valid + invalid); NaN where scalar None (either chain invalid); float64.
* STRADDLE (deterministic, seeded): over >=1000 seeded random VH/VL pairs,
  round(vectorized, 3) == round(scalar, 3) for every pair where the scalar is
  not None. Zero straddles.
"""

from __future__ import annotations

import math
import random

import numpy as np
from hypothesis import given
from hypothesis import strategies as st

import properties as scalar
from aa_tables import STANDARD_AAS
from pka_tables import IPC2_PROTEIN
from vectorized import build_counts, fv_isoelectric_point

# Mixed valid + invalid sequences so every row exercises one of the branches.
_seq = st.one_of(
    st.text(alphabet=STANDARD_AAS, min_size=1, max_size=50),
    st.sampled_from(["", "*", "BZXJ"]),
)
pairs = st.lists(st.tuples(_seq, _seq), min_size=1, max_size=25)


@given(pairs, st.booleans())
def test_fv_pi_parity(pairs, inc):
    vhs = [vh for vh, _ in pairs]
    vls = [vl for _, vl in pairs]
    sub_vh = build_counts(vhs)
    sub_vl = build_counts(vls)
    vec = fv_isoelectric_point(sub_vh, sub_vl, IPC2_PROTEIN, include_cys=inc)
    assert vec.dtype == np.float64
    assert vec.shape == (len(pairs),)
    for i, (vh, vl) in enumerate(pairs):
        # Scalar fv_isoelectric_point hard-codes include_cys=False; mirror its
        # cleaning here and bisect the SUM with the requested include_cys so the
        # parity holds for both settings.
        vh_clean = scalar._prepare(vh)
        vl_clean = scalar._prepare(vl)
        if vh_clean is None or vl_clean is None:
            assert math.isnan(vec[i])
            continue
        ip_vh = scalar._ipc2_isoelectric_point(vh_clean, IPC2_PROTEIN, include_cys=inc)
        ip_vl = scalar._ipc2_isoelectric_point(vl_clean, IPC2_PROTEIN, include_cys=inc)
        exp = scalar._bisect_charge_zero(lambda ph: ip_vh.charge_at_pH(ph) + ip_vl.charge_at_pH(ph))
        if exp is None:
            assert math.isnan(vec[i])
        else:
            assert abs(vec[i] - exp) <= 1.5e-3


def _seeded_random_sequences(rng: random.Random, n: int) -> list[str]:
    out: list[str] = []
    for _ in range(n):
        length = rng.randint(20, 130)  # realistic VH/VL length band
        out.append("".join(rng.choice(STANDARD_AAS) for _ in range(length)))
    return out


def test_fv_pi_no_straddles():
    """Determinism gate: round(vectorized, 3) == round(scalar, 3) for every pair
    where the scalar is not None, across >=1000 seeded random VH/VL pairs under
    the Fv config (IPC2_PROTEIN, include_cys=False). Zero straddles.

    Uses the public scalar `properties.fv_isoelectric_point` (the exact oracle
    the pipeline replaces) as the reference.
    """
    rng = random.Random(0)
    vhs = _seeded_random_sequences(rng, 1500)
    vls = _seeded_random_sequences(rng, 1500)

    sub_vh = build_counts(vhs)
    sub_vl = build_counts(vls)
    vec = fv_isoelectric_point(sub_vh, sub_vl, IPC2_PROTEIN, include_cys=False)
    assert vec.dtype == np.float64

    straddles: list[tuple[str, str, float, float]] = []
    compared = 0
    for i, (vh, vl) in enumerate(zip(vhs, vls)):
        exp = scalar.fv_isoelectric_point(vh, vl, IPC2_PROTEIN)
        if exp is None:
            assert math.isnan(vec[i])
            continue
        compared += 1
        if round(float(vec[i]), 3) != round(float(exp), 3):
            straddles.append((vh, vl, float(vec[i]), float(exp)))

    assert compared >= 1000
    assert straddles == [], f"{len(straddles)} straddle(s) at 3 dp; first: {straddles[0]}"
