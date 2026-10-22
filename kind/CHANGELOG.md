# @platforma-open/milaboratories.sequence-properties.kind

## 1.1.0

### Minor Changes

- d32a5ac: Update the SDK and accept a subset filter (e.g. a Repertoire Labeling label or Lead Selection picks) with the input dataset. The dataset list is built with buildDatasetOptions (requireEnrichments: false). The subset restricts the input before the analysis, so properties are computed only for its entities, and dataset-level summaries (e.g. median CDR3 length) describe the subset. Outputs of a subset run carry pl7.app/inputSubset in their domain and a trace rooted in the subset, so their labels name it; the default subtitle names the subset too. Unfiltered runs are unchanged.
