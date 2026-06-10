from __future__ import annotations

from hypothesis import given
from hypothesis import strategies as st

from aa_tables import STANDARD_AA_SET, STANDARD_AAS
from properties import aa_counts, clean_sequence, effective_length, is_invalid_sequence
from vectorized import _BYTE_TO_AA, build_counts

# Include non-standard residues, stop codons, lowercase, gaps, and empties.
raw = st.lists(st.text(alphabet=STANDARD_AAS + "*BZXJUabc-", max_size=40), min_size=1, max_size=30)


@given(raw)
def test_counts_match_scalar_oracle(seqs):
    sub = build_counts(seqs)
    assert sub.counts.shape == (len(seqs), 20)
    for i, s in enumerate(seqs):
        invalid = is_invalid_sequence(s) or not clean_sequence(s)
        assert bool(sub.valid[i]) == (not invalid)
        if not invalid:
            expected = [aa_counts(s)[aa] for aa in STANDARD_AAS]
            assert sub.counts[i].tolist() == expected
            assert int(sub.length[i]) == effective_length(s)


def test_handles_none_and_empty():
    sub = build_counts([None, "", "*", "BZXJ", "A"])
    assert sub.valid.tolist() == [False, False, False, False, True]
    assert int(sub.length[4]) == 1


def test_byte_table_kept_set_matches_oracle():
    """Structural sync guard: the vectorized byte->AA table keeps a byte iff the
    oracle's cleaning rule keeps that character. A future change to STANDARD_AAS
    or the byte table that desyncs them from `properties.clean_sequence`'s
    `c.upper() in STANDARD_AA_SET` rule fails here. Covers upper AND lower case.
    """
    for b in range(256):
        kept = _BYTE_TO_AA[b] >= 0
        oracle_kept = chr(b).upper() in STANDARD_AA_SET
        assert kept == oracle_kept, f"byte {b} ({chr(b)!r}): table kept={kept}, oracle kept={oracle_kept}"
