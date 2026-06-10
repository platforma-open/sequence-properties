"""Parity tests: vectorized charge / ΔCharge vs the scalar properties.py oracle.

Charge follows BioPython's IsoelectricPoint.charge_at_pH (Henderson-Hasselbalch)
with IPC 2.0 pKa overrides. For every VALID generated sequence each vectorized
value must equal the scalar value within abs=1e-6, across BOTH pKa sets
(IPC2_PEPTIDE, IPC2_PROTEIN) and BOTH include_cys settings, at several pH
points. Where the scalar returns None (invalid: empty, stop codon,
non-standard-only), the vectorized row must be NaN. Output arrays are float64.
"""

from __future__ import annotations

import math

import numpy as np
from hypothesis import given
from hypothesis import strategies as st

import properties as scalar
from aa_tables import STANDARD_AAS
from pka_tables import IPC2_PEPTIDE, IPC2_PROTEIN
from vectorized import build_counts, charge_at_ph, charge_shift

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

PKA = [IPC2_PEPTIDE, IPC2_PROTEIN]


@given(seqs, st.sampled_from(PKA), st.booleans(), st.sampled_from([6.0, 7.0, 7.4]))
def test_charge_parity(seqs, pka, inc, ph):
    sub = build_counts(seqs)
    vec = charge_at_ph(sub, ph, pka, include_cys=inc)
    assert vec.dtype == np.float64
    assert vec.shape == (len(seqs),)
    for i, s in enumerate(seqs):
        exp = scalar.charge_at_ph(s, ph, pka, include_cys=inc)
        if exp is None:
            assert math.isnan(vec[i])
        else:
            assert abs(vec[i] - exp) <= 1e-6


@given(seqs, st.sampled_from(PKA), st.booleans())
def test_charge_shift_parity(seqs, pka, inc):
    sub = build_counts(seqs)
    vec = charge_shift(sub, pka, include_cys=inc)
    assert vec.dtype == np.float64
    assert vec.shape == (len(seqs),)
    for i, s in enumerate(seqs):
        exp = scalar.charge_shift(s, pka, include_cys=inc)
        if exp is None:
            assert math.isnan(vec[i])
        else:
            assert abs(vec[i] - exp) <= 1e-6
