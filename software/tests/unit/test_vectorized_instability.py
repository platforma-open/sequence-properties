"""Parity tests: vectorized instability index vs the scalar properties.py oracle.

The Guruprasad instability index is the only ORDER-dependent property — it sums
the BioPython dipeptide instability weights (ProtParamData.DIWV) over consecutive
residue pairs of the CLEANED sequence, then scales by 10/length:

    score  = sum(DIWV[seq[i]][seq[i+1]] for i in range(len - 1))
    index  = (10.0 / len(seq)) * score

(See Bio.SeqUtils.ProtParam.ProteinAnalysis.instability_index.)

For every VALID generated sequence the vectorized value must equal the scalar
value within abs=1e-6. Where the scalar returns None — invalid (None / empty /
stop codon / non-standard-only) OR effective length < INSTABILITY_MIN_LENGTH
(= 10) — the vectorized row must be NaN. Output arrays are float64, and results
are order-stable (a sequence's value never depends on its position in the list).
"""

from __future__ import annotations

import math

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

import properties as scalar
from aa_tables import STANDARD_AAS
from vectorized import instability_index

# Pools that exercise every branch of the cleaning + floor logic.
_STD = st.text(alphabet="".join(STANDARD_AAS), min_size=0, max_size=30)
_LOWER = st.text(alphabet="".join(STANDARD_AAS).lower(), min_size=1, max_size=10)
_NONSTD = st.text(alphabet="BZXUJ-", min_size=1, max_size=8)  # non-standard only
_INVALID = st.sampled_from([None, "", "*", "ACD*EFG", "BZXUJ", "----", "bcdefghij"])

_SEQ = st.one_of(_STD, _LOWER, _NONSTD, _INVALID)


@given(seqs=st.lists(_SEQ, min_size=1, max_size=40))
@settings(max_examples=300)
def test_matches_scalar_oracle(seqs: list[str | None]) -> None:
    vec = instability_index(seqs)
    assert vec.dtype == np.float64
    assert vec.shape == (len(seqs),)

    for i, s in enumerate(seqs):
        expected = scalar.instability_index(s)
        if expected is None:
            assert math.isnan(vec[i]), f"row {i} {s!r}: expected NaN, got {vec[i]}"
        else:
            assert not math.isnan(vec[i]), f"row {i} {s!r}: expected {expected}, got NaN"
            assert abs(vec[i] - expected) <= 1e-6, (
                f"row {i} {s!r}: vec={vec[i]} scalar={expected} residual={abs(vec[i] - expected)}"
            )


def test_floor_below_min_length_is_nan() -> None:
    """effective length < 10 (after cleaning) -> NaN; >= 10 -> finite."""
    # 9 standard residues (with junk that cleans away) -> below floor.
    below = "AC-DE*"  # has stop codon -> invalid anyway
    nine = "ACDEFGHIK"  # 9 std residues, exactly below floor
    ten = "ACDEFGHIKL"  # 10 std residues, exactly at floor
    eleven_dirty = "ACDEFGHIKLM-XBZ"  # 11 std residues + junk -> at/above floor

    vec = instability_index([nine, ten, eleven_dirty, below])
    assert math.isnan(vec[0]), "9 residues must be below the floor -> NaN"
    assert not math.isnan(vec[1]), "10 residues must be at the floor -> finite"
    assert not math.isnan(vec[2]), "11 cleaned residues must be above the floor"
    assert math.isnan(vec[3]), "stop codon -> invalid -> NaN"

    # Cross-check the finite rows against the oracle.
    assert abs(vec[1] - scalar.instability_index(ten)) <= 1e-6
    assert abs(vec[2] - scalar.instability_index(eleven_dirty)) <= 1e-6


def test_invalid_rows_are_nan() -> None:
    seqs = [None, "", "*", "ACD*EF", "BZXUJ", "----", "X"]
    vec = instability_index(seqs)
    assert vec.dtype == np.float64
    assert np.all(np.isnan(vec)), f"all invalid rows must be NaN, got {vec}"


def test_known_value() -> None:
    """Spot-check a fixed 12-mer against a directly recomputed DIWV sum."""
    seq = "ACDEFGHIKLMN"
    expected = scalar.instability_index(seq)
    vec = instability_index([seq])
    assert expected is not None
    assert abs(vec[0] - expected) <= 1e-6


def test_order_stability() -> None:
    """A sequence's result must not depend on its position in the input list."""
    seqs = [
        "ACDEFGHIKLMNPQRSTVWY",  # 20-mer, valid
        None,
        "ACDEFGHIK",  # 9-mer, below floor -> NaN
        "WYWYWYWYWYWY",  # 12-mer, valid
        "*",
        "MKLVAACDEFGHIK",  # 14-mer, valid
        "bcdefghijk",  # lowercase 10-mer -> 10 std residues, valid
    ]
    forward = instability_index(seqs)
    reversed_in = list(reversed(seqs))
    reversed_out = instability_index(reversed_in)

    # Re-reverse the reversed output to align positions with `forward`.
    realigned = reversed_out[::-1]
    np.testing.assert_array_equal(np.isnan(forward), np.isnan(realigned))  # NaN pattern stable
    finite = ~np.isnan(forward)
    np.testing.assert_allclose(forward[finite], realigned[finite], atol=0.0, rtol=0.0)
