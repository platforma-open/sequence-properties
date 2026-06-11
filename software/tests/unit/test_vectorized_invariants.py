"""Invariant (property-based) tests for the VECTORIZED engine.

`test_invariants.py` asserts these same invariants against the scalar oracle
(properties.py). The vectorized engine is the code that actually ships, so it
must hold the invariants independently — a vectorized charge that is
non-monotonic in pH, an aa-fraction row that does not sum to 1, or a pI outside
[0, 14] is a bug even if every oracle-parity test happened to miss it.

The engine is columnar: one `build_counts(seqs)` call covers a whole list, so
each invariant is asserted across all rows of a generated list at once.

Run from blocks/sequence-properties/software/:
    uv sync
    uv run pytest tests/unit/test_vectorized_invariants.py
"""

from __future__ import annotations

import numpy as np
from hypothesis import given
from hypothesis import strategies as st

from aa_tables import STANDARD_AAS
from pka_tables import IPC2_PEPTIDE
from properties import INSTABILITY_MIN_LENGTH, effective_length
from vectorized import (
    aa_fractions,
    build_counts,
    charge_at_ph,
    instability_index,
    isoelectric_point,
)

# Lists of valid standard-AA sequences. max_size kept modest so Hypothesis
# explores list shapes cheaply.
valid_seqs = st.lists(st.text(alphabet=STANDARD_AAS, min_size=1, max_size=60), min_size=1, max_size=20)

# Same, salted with known invalids (empty / stop codon / non-standard-only /
# mixed-with-stop) to exercise the NaN branches.
mixed_seqs = st.lists(
    st.one_of(
        st.text(alphabet=STANDARD_AAS, min_size=1, max_size=60),
        st.sampled_from(["", "*", "BZXJ", "ACDE*FG"]),
    ),
    min_size=1,
    max_size=20,
)


@given(valid_seqs)
def test_aa_fractions_sum_to_one(seqs):
    # Every row is a valid standard-AA sequence -> its mole fractions sum to 1.
    frac = aa_fractions(build_counts(seqs))
    assert np.allclose(frac.sum(axis=1), 1.0, atol=1e-9)


@given(mixed_seqs)
def test_aa_fractions_valid_sum_to_one_invalid_all_nan(seqs):
    sub = build_counts(seqs)
    frac = aa_fractions(sub)
    for i in range(len(seqs)):
        if sub.valid[i]:
            assert np.isclose(frac[i].sum(), 1.0, atol=1e-9)
        else:
            assert np.all(np.isnan(frac[i]))


@given(mixed_seqs)
def test_pi_in_range_or_nan(seqs):
    # Every defined pI sits in the bisection bracket [0, 14]; invalid -> NaN.
    pi = isoelectric_point(build_counts(seqs), IPC2_PEPTIDE, include_cys=True)
    finite = pi[~np.isnan(pi)]
    assert np.all((finite >= 0.0) & (finite <= 14.0))


@given(
    valid_seqs,
    st.floats(min_value=0.5, max_value=13.5),
    st.floats(min_value=0.5, max_value=13.5),
)
def test_charge_monotonic_in_ph(seqs, ph_a, ph_b):
    lo, hi = sorted((ph_a, ph_b))
    sub = build_counts(seqs)
    c_lo = charge_at_ph(sub, lo, IPC2_PEPTIDE, include_cys=True)
    c_hi = charge_at_ph(sub, hi, IPC2_PEPTIDE, include_cys=True)
    # Net charge is non-increasing in pH for every row.
    assert np.all(c_lo >= c_hi - 1e-9)


@given(mixed_seqs)
def test_instability_floor(seqs):
    # NaN below the effective-length floor (spec R9), finite at/above it.
    ii = instability_index(seqs)
    for i, s in enumerate(seqs):
        if effective_length(s) < INSTABILITY_MIN_LENGTH:
            assert np.isnan(ii[i])
        else:
            assert not np.isnan(ii[i])
