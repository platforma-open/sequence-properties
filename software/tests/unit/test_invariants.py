from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from aa_tables import STANDARD_AAS
from pka_tables import IPC2_PEPTIDE, IPC2_PROTEIN
from properties import (
    INSTABILITY_MIN_LENGTH,
    aa_fractions,
    charge_at_ph,
    effective_length,
    fv_charge,
    instability_index,
    isoelectric_point,
)

aa_seq = st.text(alphabet=STANDARD_AAS, min_size=1, max_size=60)


@given(aa_seq)
def test_aa_fractions_sum_to_one(seq):
    fr = aa_fractions(seq)
    assert fr is not None
    assert sum(fr.values()) == pytest.approx(1.0, abs=1e-9)


@given(aa_seq)
def test_pi_in_range_or_none(seq):
    pi = isoelectric_point(seq, IPC2_PEPTIDE, include_cys=True)
    assert pi is None or (0.0 <= pi <= 14.0)


@given(
    aa_seq,
    st.floats(min_value=0.5, max_value=13.5),
    st.floats(min_value=0.5, max_value=13.5),
)
def test_charge_monotonic_in_ph(seq, ph_a, ph_b):
    lo, hi = sorted((ph_a, ph_b))
    c_lo = charge_at_ph(seq, lo, IPC2_PEPTIDE, include_cys=True)
    c_hi = charge_at_ph(seq, hi, IPC2_PEPTIDE, include_cys=True)
    assert c_lo >= c_hi - 1e-9  # charge is non-increasing in pH


@given(aa_seq, aa_seq)
def test_fv_charge_additive(vh, vl):
    fv = fv_charge(vh, vl, 7.0, IPC2_PROTEIN)
    a = charge_at_ph(vh, 7.0, IPC2_PROTEIN, include_cys=False)
    b = charge_at_ph(vl, 7.0, IPC2_PROTEIN, include_cys=False)
    assert fv == pytest.approx(a + b, abs=1e-9)


@given(aa_seq)
def test_instability_floor(seq):
    ii = instability_index(seq)
    if effective_length(seq) < INSTABILITY_MIN_LENGTH:
        assert ii is None
    else:
        assert ii is not None
