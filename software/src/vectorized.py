"""Vectorized compute substrate. A column of sequences -> per-residue count
matrix + length + validity, matching properties.py's scalar cleaning exactly.
Single-threaded by construction (pure numpy elementwise + per-seq counting)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aa_tables import STANDARD_AAS

_AA_INDEX = {aa: i for i, aa in enumerate(STANDARD_AAS)}


@dataclass(frozen=True)
class Substrate:
    counts: np.ndarray  # (N, 20) int64, STANDARD_AAS order
    length: np.ndarray  # (N,) int64 effective length
    valid: np.ndarray  # (N,) bool


def build_counts(seqs: list[str | None]) -> Substrate:
    n = len(seqs)
    counts = np.zeros((n, 20), dtype=np.int64)
    valid = np.zeros(n, dtype=bool)
    for i, s in enumerate(seqs):
        if s is None or s == "" or "*" in s:
            continue
        row = counts[i]
        any_std = False
        for c in s.upper():
            j = _AA_INDEX.get(c)
            if j is not None:
                row[j] += 1
                any_std = True
        valid[i] = any_std
    length = counts.sum(axis=1)
    return Substrate(counts=counts, length=length, valid=valid)
