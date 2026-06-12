---
"@platforma-open/milaboratories.sequence-properties.workflow": patch
---

Namespace exported property columns by a content hash of the source table (pl7.app/contentHash, derived from the backend CanonicalID) instead of the per-block blockId, and drop the blockId-keyed trace step id. Identical results across blocks/projects now produce content-identical columns that dedupe downstream instead of being made unique per block; different results stay distinct. The heavy Python step and Parquet import already deduped; this extends dedup to the exported pframe identity.
