from __future__ import annotations

from hypothesis import given
from hypothesis import strategies as st

from aa_tables import STANDARD_AAS
from properties import aa_counts, clean_sequence, effective_length, is_invalid_sequence
from vectorized import build_counts

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
