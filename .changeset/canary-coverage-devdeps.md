---
"@platforma-open/milaboratories.sequence-properties.ui": patch
---

Canary (throwaway): touches only the UI while the catalog bumps build-tool
devDependencies (block-tools, tengo-builder). Verifies the changeset-coverage
check no longer flags the model/workflow packages that merely consume those
build tools via `catalog:`.
