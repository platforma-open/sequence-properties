"""Parity tests: vectorized linear properties vs the scalar properties.py oracle.

For every VALID generated sequence each vectorized property must equal the
scalar value within abs=1e-6. Where the scalar returns None (invalid: empty,
stop codon, non-standard-only), the vectorized row must be NaN (and (NaN, NaN)
for extinction). Output arrays must be float64 so the downstream pipeline's
Float64 quantization rounding still applies.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

import properties as scalar
from aa_tables import STANDARD_AAS
from vectorized import (
    aa_fractions,
    aliphatic_index,
    aromaticity,
    build_counts,
    extinction,
    gravy,
    molecular_weight,
)

_TOL = 1e-6

# Mixed list: valid standard-AA seqs of varied length AND known invalids
# (empty, stop codon, non-standard-only, mixed-with-stop).
seqs_strategy = st.lists(
    st.one_of(
        st.text(alphabet=STANDARD_AAS, min_size=1, max_size=50),
        st.sampled_from(["", "*", "BZXJ", "ACDEF*GH", "X", "---", "acdik"]),
    ),
    min_size=1,
    max_size=25,
)


def _assert_scalar_parity(vec, seqs, scalar_fn):
    """vec[i] == scalar_fn(seqs[i]) within tol; NaN where scalar returns None."""
    for i, s in enumerate(seqs):
        exp = scalar_fn(s)
        if exp is None:
            assert math.isnan(vec[i]), f"seq={s!r}: expected NaN, got {vec[i]}"
        else:
            assert not math.isnan(vec[i]), f"seq={s!r}: expected {exp}, got NaN"
            assert abs(vec[i] - exp) <= _TOL, f"seq={s!r}: vec={vec[i]} scalar={exp} diff={abs(vec[i] - exp)}"


@given(seqs_strategy)
def test_gravy_parity(seqs):
    vec = gravy(build_counts(seqs))
    assert vec.dtype == np.float64
    _assert_scalar_parity(vec, seqs, scalar.gravy)


@given(seqs_strategy)
def test_molecular_weight_parity(seqs):
    vec = molecular_weight(build_counts(seqs))
    assert vec.dtype == np.float64
    _assert_scalar_parity(vec, seqs, scalar.molecular_weight)


@given(seqs_strategy)
def test_aromaticity_parity(seqs):
    vec = aromaticity(build_counts(seqs))
    assert vec.dtype == np.float64
    _assert_scalar_parity(vec, seqs, scalar.aromaticity)


@given(seqs_strategy)
def test_aliphatic_index_parity(seqs):
    vec = aliphatic_index(build_counts(seqs))
    assert vec.dtype == np.float64
    _assert_scalar_parity(vec, seqs, scalar.aliphatic_index)


@given(seqs_strategy)
def test_extinction_parity(seqs):
    ox, red = extinction(build_counts(seqs))
    assert ox.dtype == np.float64
    assert red.dtype == np.float64
    for i, s in enumerate(seqs):
        exp_ox, exp_red = scalar.extinction_coefficients(s)
        if exp_ox is None:
            assert math.isnan(ox[i]), f"seq={s!r}: expected NaN ox, got {ox[i]}"
            assert math.isnan(red[i]), f"seq={s!r}: expected NaN red, got {red[i]}"
        else:
            assert abs(ox[i] - exp_ox) <= _TOL, f"seq={s!r}: ox vec={ox[i]} scalar={exp_ox}"
            assert abs(red[i] - exp_red) <= _TOL, f"seq={s!r}: red vec={red[i]} scalar={exp_red}"


@given(seqs_strategy)
def test_aa_fractions_parity(seqs):
    frac = aa_fractions(build_counts(seqs))
    assert frac.dtype == np.float64
    assert frac.shape == (len(seqs), len(STANDARD_AAS))
    for i, s in enumerate(seqs):
        exp = scalar.aa_fractions(s)
        if exp is None:
            assert np.all(np.isnan(frac[i])), f"seq={s!r}: expected all-NaN row"
        else:
            for j, aa in enumerate(STANDARD_AAS):
                assert abs(frac[i, j] - exp[aa]) <= _TOL, f"seq={s!r} aa={aa}: vec={frac[i, j]} scalar={exp[aa]}"


def test_invalid_rows_are_nan():
    """Explicit invalid-row coverage, decoupled from hypothesis sampling."""
    seqs = ["", "*", "BZXJ", "ACDEF*GH", "ACDEFGHIKLMNPQRSTVWY"]
    sub = build_counts(seqs)
    for fn in (gravy, molecular_weight, aromaticity, aliphatic_index):
        vec = fn(sub)
        assert np.all(np.isnan(vec[:4]))  # invalids
        assert not math.isnan(vec[4])  # the valid seq
    ox, red = extinction(sub)
    assert np.all(np.isnan(ox[:4])) and np.all(np.isnan(red[:4]))
    assert not math.isnan(ox[4]) and not math.isnan(red[4])
    frac = aa_fractions(sub)
    assert np.all(np.isnan(frac[:4]))
    assert pytest.approx(frac[4].sum(), abs=_TOL) == 1.0
