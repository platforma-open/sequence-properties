---
"@platforma-open/milaboratories.sequence-properties.model": minor
"@platforma-open/milaboratories.sequence-properties.ui": minor
"@platforma-open/milaboratories.sequence-properties": minor
---

Add Scatterplot and Histogram tabs to the Sequence Properties block. Both
panels read the existing `propertiesPf` p-frame and pick modality-aware
defaults: peptide charge / hydrophobicity in peptide mode, chain "A" CDR3
charge / hydrophobicity in antibody/TCR mode. Axis pickers list every
numeric scalar PColumn emitted by the run, excluding the 2-axis AA fraction
column. R21 / R21a reference line at GRAVY = 0 deferred — see
docs/spec-deviations.md SD-009.
