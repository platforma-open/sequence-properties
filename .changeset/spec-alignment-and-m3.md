---
'@platforma-open/MiLaboratories.sequence-properties.workflow': patch
'@platforma-open/MiLaboratories.sequence-properties.software': patch
---

Align two spec deviations and discharge the M3 external-validation acceptance criteria.

- **SD-004 reconciled.** Peptide-mode TSV column renamed from `peptide_seq` to `sequence` per spec L457. Internal Tengo→Python contract; no PColumn output names changed and no downstream consumer reads the TSV directly.
- **SD-006 reconciled.** Modality detection now scans `axesSpec` for the first axis with a recognized name + domain (per spec R1a) instead of unconditionally picking `axes[len-1]`. Equivalent on every observed `[sampleId, key]` input; first-matching is correct for hypothetical multi-axis layouts.
- **M3 validation tests added.** New `tests/unit/test_m3_validation.py` (38 cases) cross-checks pI / charge / Fv / aliphatic against the IPC 2.0 webserver and an independent textbook Henderson-Hasselbalch reference. Discharges the M3 acceptance criteria for >=5 VH pI, >=2 VL pI, Fv on >=2 paired chains, >=10 CDR-H3 charge, >=3 CDR-L3 charge, and >=3 VH aliphatic. Webserver values pinned (spec L518 says webserver may be unavailable in CI).
