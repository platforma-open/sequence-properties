---
'@platforma-open/MiLaboratories.sequence-properties.software': patch
---

Switch python runenv from `3.12.10-scientific-slim` to the full `3.12.10`, and bump `biopython` from 1.84 to 1.87.

The slim variant bundled only polars/numpy/scipy/pyarrow wheels; biopython had to come from elsewhere and was not in any local cache, breaking fresh venv builds (`Could not find a version that satisfies the requirement biopython==1.84`). The full runenv ships biopython==1.87 wheels for every supported platform, matching the workspace majority (23 blocks use the full runenv, only this block + `titeseq-analysis` were on slim). pip still installs only what `requirements.txt` lists, so the venv size is unchanged.

Verified: full python suite (165/165), block build clean, live run on TinyTrees produced `coverageTier: full_chain` with no info errors.
