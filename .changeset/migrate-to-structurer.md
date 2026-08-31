---
'@platforma-open/milaboratories.sequence-properties': minor
'@platforma-open/milaboratories.sequence-properties.model': minor
'@platforma-open/milaboratories.sequence-properties.ui': patch
'@platforma-open/milaboratories.sequence-properties.workflow': patch
'@platforma-open/milaboratories.sequence-properties.software': patch
---

Migrate onto the structurer and upgrade the SDK

The block now follows the canonical `block-tools structure` layout on
block-tools 2.14.3, with the SDK bumped to model/ui-vue 1.83.x, workflow-tengo
6.8.3, tengo-builder 4.0.23 and test 1.83.2.

Two changes are visible outside the block:

- **A new `kind` component** declares the block's identity and its init-params
  contract. `BlockParams` is `{ inputAnchor?: PlRef }`, so a project template
  can seed a new instance with the dataset to analyse. Every other field in the
  block's data is view state and still defaults.
- **`block` is now a slim facade.** It publishes with no dependencies, bundles
  the whole block into `dist/` plus `block-pack/`, and exports
  `SequencePropertiesBlockPointer` for consumers that add the block from code.

`@platforma-sdk/ui-vue` is on 1.83.3, which publishes the component
declaration files its own type entry re-exports again — 1.83.1 had dropped
them and broke the facade build.
