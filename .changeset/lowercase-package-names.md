---
'@platforma-open/milaboratories.sequence-properties': patch
'@platforma-open/milaboratories.sequence-properties.workflow': patch
'@platforma-open/milaboratories.sequence-properties.model': patch
'@platforma-open/milaboratories.sequence-properties.ui': patch
---

Rename all package scopes from `MiLaboratories.sequence-properties` to `milaboratories.sequence-properties`. npm registry rejects new package names with uppercase letters, which blocked the first publish. Lowercase form aligns with the existing `@platforma-open/milaboratories.*` convention used by sibling blocks. Also corrects the GitHub URL in the block manifest to point at the actual repo (`platforma-open/sequence-properties`).
