# Dedup Investigation — sequence-properties

Switching the input anchor away and back triggers a fresh Python execution.
The platform's content-addressed dedup should recover the previous result for
identical inputs but does not. Our top-level workflow code is canonical;
something deeper in the SDK is not.

**Date:** 2026-04-29
**Project used for testing:** BulkTrees (`NG:0x58857b`)
**Block:** sequence-properties (`f1df7393-3be8-4d9a-bc15-7bdeda726c36`)

---

## Symptom

Round-trip timeline on a single block, MiXCR Clonotyping upstream unchanged:

| Step | Input anchor | Python log timestamp |
|---|---|---|
| 1. Initial run | IG | `11:05:40` |
| 2. Switch + run | TCR Alpha/Beta | `11:09:35` |
| 3. Switch back + run | IG (same as step 1) | `11:09:54` ← new |

If dedup hit on step 3, the log would still read `11:05:40` (the original
blob's content). It read `11:09:54` — Python ran again.

The harness states dedup is the design intent:

> Each template body executes once for a given (template, inputs) pair. The
> same inputs do not re-run the body — the platform recovers the previous
> result. **This is the point.**
> — `mictx-helper/harnesses/block-dev/workflow.md`

## What we ruled out

| Hypothesis | Probe | Result |
|---|---|---|
| `vdjCols` iteration order non-deterministic | Captured `s.key` for all 18 columns across runs | Byte-identical alphabetical order both runs |
| `canonical.encode(plan)` produces different bytes | Captured the encoded string both runs | Byte-identical: `{"chains":["A","B"],"fullChains":["A","B"],"hasFv":true,"mode":"antibody_tcr_legacy_sc","receptor":"IG"}` |
| `chainsFound` / `fullChain` projections non-canonical | Already use `maps.getKeys()` (sorted) | Code path canonical |

## What we did not pin down

The dedup miss happens despite canonical top-level inputs. Probable layers
to investigate next, ranked by likelihood:

1. **`bundle.getColumn(s.key)`** may return per-render resource handles
   even though `s.key` is stable. The exec's `addFile("input.tsv", seqTable)`
   then sees a new handle each time. Need to inspect the bundle accessor in
   `pframes.lib.tengo` and `bundle.lib.tengo`.
2. **`pframes.tsvFileBuilder().build()`** wraps a sub-template. If that
   sub-template's resource construction is non-canonical, every build
   produces a new seqTable resource even when content is byte-equal.
   Need to read `xsv-builder.lib.tengo` for the sub-template definition.
3. **`wf.prepare`'s `bb.build()`** runs each render and may allocate the
   bundle freshly. If the bundle resource is part of the workflow body's
   input hash, every render of the body sees a new bundle handle, breaking
   exec dedup transitively.

`seqTable.id` returned `{ResourceID: 4611686018427387926, Name: "result"}` on
the first IG run and `{ResourceID: 4611686018427387927, Name: "result"}` on
the round-trip. The numeric ID is a sequential resource counter, not a
content hash, so it does not disambiguate (1)/(2)/(3).

## Diagnostic instrumentation currently in the block

Tracking what to revert before merging:

- `workflow/src/main.tpl.tengo`:
  - `debugIterKeys := []` declared near `seqTb` build site.
  - The `for s in vdjCols` loop appends `string(s.key) + " => col.id=" + string(col.id)` to `debugIterKeys`.
  - `debugBlob := canonicalJsonResource({ vdjOrder, planEncoded, seqTableId })` and `debugVdjOrder` exposed in `outputs`.
- `model/src/index.ts`:
  - `.output("debugVdjOrder", ...)` — JSON of the diagnostic blob.
  - `.output("processingLogText", ...)` — `getLastLogs(50)` snapshot of the Python stderr stream.

All four removable in one pass. No production code reads them.

## Other information collected during the session

- The block is correct end-to-end after applying SD-001 / SD-002 / SD-003.
- BulkTrees IG dataset has 65 clones; TCRAB has 0 clones.
- The pl MCP `get_block_logs` tool does not surface our `processingLog`
  output's content. `getLastLogs` exposed as a model output worked as a
  workaround for inspection.
- The `processingLog` URI changes on every render even when the underlying
  exec hits cache (cannot be confirmed without a content read). For this
  block the URI changes correlated with content changes, so URI inequality
  is sufficient evidence of dedup miss here, but is not a general probe.

## Next steps

1. **Read the SDK source.** Start with `bundle.lib.tengo` (column accessor)
   and `pframes.xsv-builder.lib.tengo` (`build()` path). Look for
   per-render allocations or non-canonical map iteration in the sub-template
   construction.
2. **Probe the seqTable's actual content hash.** The Tengo `.id` field is
   sequential. Need a content-derived field (CID or equivalent). If the SDK
   does not expose one, write the `input.tsv` to disk and hash externally
   across the two IG runs — definitive evidence whether the seqTable bytes
   are identical.
3. **If seqTable bytes are byte-equal but exec resource hash differs**, the
   issue is platform-side: exec dedup uses handle identity, not content
   identity. File a reflection in
   `mictx-helper/harnesses/block-dev/_meta/reflections/`.
4. **If seqTable bytes differ**, the issue is in our exec input chain. Most
   likely culprit is the tsvFileBuilder's sub-template producing
   non-canonical bytes. Fix at that layer.
5. **Either way, write a minimal reproducer** — a tiny test block that does
   nothing but build a seqTable from a stable input and feeds an exec.
   Confirms the issue without sequence-properties' real complexity in the
   way.
6. **Clean up the diagnostic instrumentation** listed above before any PR.
   It compiles and runs, but it pollutes outputs and adds work to every
   render.

## References

- Harness: `mictx-helper/harnesses/block-dev/workflow.md` (Canonicality —
  the load-bearing rule).
- Spec deviations applied during this session: `docs/spec-deviations.md`
  (SD-001, SD-002, SD-003).
- Test commands used: `mcp__pl__set_block_data`, `mcp__pl__run_block`,
  `mcp__pl__get_block_state` with transforms.
