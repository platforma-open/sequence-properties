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


def test_non_ascii_chars_dropped_not_expanded():
    """build_counts cleans via a single latin-1 byte per character, applied
    BEFORE upper-casing — so each character is either one standard-AA byte or
    dropped, and the standard-residue count can never exceed the input length.

    This DIVERGES by construction from the scalar oracle's str.upper()-then-
    filter: Python's str.upper() case-folds some non-ASCII chars to MULTIPLE
    ASCII letters ('ß' -> 'SS', the 'fi' ligature -> 'FI'), which the oracle
    would then count, whereas the single-byte gather drops them. Real input is
    ASCII single-letter AA codes, so this never fires in production. The test
    pins the intentional behavior so a future "match str.upper exactly" change
    is a deliberate decision, not an accident — and `test_counts_match_scalar_
    oracle` above only generates ASCII, so it never exercises this seam.
    """
    # 'ß' (U+00DF) is one latin-1 byte (0xDF); the 'fi' ligature (U+FB01) is not
    # latin-1, so `errors="replace"` maps it to '?' (0x3F). Both map to -1 in the
    # byte table -> dropped. Neither is expanded to the SS / FI str.upper() would.
    sub = build_counts(["ACßDE", "ACﬁDE"])
    assert int(sub.length[0]) == 4  # ß dropped, NOT expanded to SS (would be 6)
    assert int(sub.length[1]) == 4  # fi-ligature dropped, NOT expanded to FI
    for i in range(2):
        counts = dict(zip(STANDARD_AAS, sub.counts[i].tolist()))
        assert counts["A"] == counts["C"] == counts["D"] == counts["E"] == 1
        assert sum(sub.counts[i].tolist()) == 4  # exactly A, C, D, E — no S/F/I


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
