---
'@platforma-open/milaboratories.sequence-properties.software': minor
'@platforma-open/milaboratories.sequence-properties.workflow': minor
'@platforma-open/milaboratories.sequence-properties': minor
---

Add ΔCharge (pH 7.4 → 6.0) metric — `pl7.app/chargeShift` — emitted at peptide, CDR3 (per chain), and Fv scopes. Captures pH-switching capacity (FcRn recycling, endosomal release); negative values mean the molecule gains positive charge on acidification, the productive direction for histidine-driven pH switching. Histidine dominates the metric (~−0.46 per His; pKa ~6.0 sits in the window). Domain carries the pH endpoints (`pl7.app/pH/from`, `pl7.app/pH/to`) so additional pH pairs can land later without breaking the v1 column identity. Default-visible alongside the static charge column at each scope; not marked `isScore` (interpretive, not a Lead Selection ranking criterion).
