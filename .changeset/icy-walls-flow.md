---
'@platforma-open/MiLaboratories.sequence-properties.software': minor
'@platforma-open/MiLaboratories.sequence-properties.workflow': minor
'@platforma-open/MiLaboratories.sequence-properties': minor
'@platforma-open/MiLaboratories.sequence-properties.model': minor
'@platforma-open/MiLaboratories.sequence-properties.test': minor
'@platforma-open/MiLaboratories.sequence-properties.ui': minor
---

Migrate properties table to `createPlDataTableV3` and bump `@platforma-sdk/*` to 1.69.0. Wire content-addressed dedup: sorted Tengo map iteration plus canonical JSON resources (`plan.json`, `params`, `infoBlob`); Python TSV output sorted on `entity_key` so resource bytes hash identically across runs.
