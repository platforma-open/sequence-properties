# Sequence-Properties — Dedup Wiring Plan

> **Goal.** Make this block a first-class citizen of the platform's content-addressed dedup machinery, so a duplicate project (or a second project pointing at the same upstream peptide-extraction / MiXCR-clonotyping data) lands on `Done` immediately, with no Python re-execution and shared storage on `exports.properties`.
>
> **Status.** Plan only — no code changes yet. Operator review required before implementation.

## 1. Spec status

The project spec at `docs/text/work/projects/sequence-properties/{README,pcolumn-spec}.md` covers the science (formulas, edge cases, AA tables, quantization) but does **not** discuss dedup, CIDs, canonicality, or `metaInputs`. This plan is the missing piece. It does not modify the project spec.

## 2. What the block does today

`workflow/src/main.tpl.tengo`:

- `wf.prepare` — builds a PBundle with anchored peptide / VDJ region columns. Pure, fine.
- `wf.body` — assembles a `seqTable` (input.tsv) from the bundle, builds a `plan` JSON, runs Python `compute-properties` via `exec.builder()`, hands `propertiesTsv` + `aaFractionTsv` + `params` to `processTpl` via `render.create` (pure ✓).
- `processTpl` constructs PColumn specs and emits `propertiesPf` (model output) and `exportPframe` (downstream `exports.properties`).
- Python software is deterministic by construction (closed-form formulas; pI by bisection; quantization planned).

## 3. Dedup-breaking issues found

All references are to `workflow/src/main.tpl.tengo` unless noted.

### 3.1 Unsorted map iteration affecting workflow inputs (load-bearing)

The harness rule (`harnesses/block-dev/workflow.md` § "Canonicality"): Tengo map iteration order is unspecified; if it shapes a downstream resource's input, dedup silently breaks.

| Site | Symptom |
|---|---|
| `for chain, feats in chainsFound { ... }` (~L174 — chain receptor / infoMessages logic) | `receptor` and `infoMessages` order depend on iteration order |
| `for c, v in chainsFound { chainList += [c] }` (~L223) | `chainList` order varies → goes into `plan.json` → Python input bytes vary → CID mismatch |
| `for c, v in fullChain { fullChainList += [c] }` (~L225) | same — `fullChainList` feeds `plan` and `processTpl.params` |

These are the highest-priority fixes. `chainList` and `fullChainList` are the actual byte payload that breaks Python-step CID stability.

### 3.2 Non-canonical `params` resource for `processTpl`

`processTpl` receives `smart.createJsonResource({ ..., chains: chainList, fullChains: fullChainList, ... })`. Same root cause as 3.1 — the lists must be sorted before being placed into the JSON resource.

### 3.3 `infoMessages` ordering (low priority, model-output only)

`infoMessages` is wrapped in `infoBlob` and emitted as `outputs.info`. It does not feed `processTpl`, so it does not affect Python-step CID, but a non-deterministic `infoBlob` CID still defeats output-level dedup for the model. Build messages from sorted iteration too — cheap to fix while we're there.

### 3.4 Resource-allocation inputs — verify, don't assume

`exec.builder().mem("4GiB").cpu(1)` and `seqTb.mem("4GiB"); seqTb.cpu(1)` should be SERVICE fields excluded from CID per the SDK contract. Read `core/platforma/sdk/workflow-tengo/src/exec/...` and `pframes/...` to confirm both helpers route memory/CPU through the meta-input path. If they don't, follow the documented `metaInputs` pattern from `render.lib.tengo:26-29, 175-178` (third positional arg of `render.create`).

### 3.5 Python output byte-stability

CID merge requires byte-identical `properties.tsv` / `aa_fraction.tsv` between runs of the same input.

- **Row order.** Python must emit rows in a stable order (sort by entity-axis key). Polars default order ≠ guaranteed; verify in `io_layer.py::write_output_tsv` and add a `.sort([key_col])` pass.
- **Float formatting.** Apply the planned quantization (`docs/columnar-and-quantization-plan.md`) before write — fixed-precision string formatting per column.
- **NaN/null tokens.** `null_value=""` is already the policy (per `test_io_layer.py`) — keep.
- **Column order.** Write columns in a fixed declared order, not dict-insertion order.

A non-deterministic Python output is the platform-level `CIDConflictError` failure mode. Catch it in tests (§5).

### 3.6 `processTpl` audit

Repeat the same map-iteration sweep inside `process.tpl.tengo`. Out of scope for this plan until I read it; add as a follow-up task.

## 4. Proposed changes (Tengo)

### 4.1 The SDK surface — `canonical.*` and `maps.*`

Read `core/platforma/sdk/workflow-tengo/src/canonical.lib.tengo` (44 lines) and `maps.lib.tengo` before changing anything — both libraries are small enough to read end-to-end.

**`canonical` exports exactly one function: `canonical.encode(obj) -> string`.** There is no `canonical.decode`, `canonical.hash`, or `canonical.sort`. The full export is literally `ll.toStrict({ encode: encode })`.

What `canonical.encode` does — and doesn't do:

| Behavior | Status |
|---|---|
| Recursive sorted-key JSON for maps (incl. nested) | ✓ |
| Handles `ll.toStrict()` wrapped maps via `ll.fromStrict` | ✓ |
| Preserves array element order | ✓ (this is correct per JSON semantics — but see caveat) |
| **Sorts array contents** | **✗ — `canonical.encode([3,2,1]) == "[3,2,1]"`** |
| **Normalizes numbers per RFC8785** | **✗ — header comment: "partial implementation focusing only on the key sorting aspect"** |
| Strips/normalizes whitespace | ✓ (no whitespace ever emitted) |

So `canonical.encode` solves *map-key-ordering* non-determinism. It does **not** solve *array-content-ordering* non-determinism, and it does **not** solve *float-formatting* non-determinism. Both of those still need to be fixed at the source — i.e. by sorting the arrays before they're encoded, and by quantizing/formatting floats before they're written.

**`maps` exports the helpers we want for sorted iteration:**

- `maps.getKeys(m)` — uses `slices.quickSortInPlace` internally; **sorted output is guaranteed** (verified against source).
- `maps.getValues(m)` — values in key-sorted order.
- `maps.forEach(m, callback)` — callback fires in key-sorted order. Cleanest call site for our use case.

### 4.2 Concrete edits in `main.tpl.tengo`

```tengo
canonical := import("@platforma-sdk/workflow-tengo:canonical")
maps := import("@platforma-sdk/workflow-tengo:maps")

// Build chainList / fullChainList in sorted order.
chainList := []
maps.forEach(chainsFound, func(c, _) { chainList += [c] })

fullChainList := []
maps.forEach(fullChain, func(c, _) { fullChainList += [c] })

// Replace the receptor / infoMessages map loop:
maps.forEach(chainsFound, func(chain, feats) {
    // ... existing logic that consumes chain + feats ...
})

// Sort messages for stable output bytes (cheap, model-output only).
infoMessages = slices.sort(infoMessages)
```

For the `plan.json` file written into the Python step — switch to canonical encoding. This is the highest-value single change in the file:

```tengo
// Before:
//   writeFile("plan.json", json.encode(plan)).
// After:
    writeFile("plan.json", canonical.encode(plan)).
```

This guarantees that the bytes of `plan.json` in the exec workdir are byte-identical across runs with the same logical `plan` — even if Tengo map iteration ever reordered the top-level keys. Real-world precedent: `blocks/repertoire-distance/workflow/src/run-distance.tpl.tengo:22` uses this exact pattern (`writeFile("metrics.json", canonical.encode(metrics))`).

Caveat: this only fixes the **bytes** for `plan.json`. The arrays inside `plan` (`chains`, `fullChains`) must already be sorted before encoding — `canonical.encode` will not reorder them. So the `maps.forEach` fix above is still required.

### 4.3 `smart.createJsonResource` is dedup-broken for maps — verified

> **Correction to an earlier note in this doc:** `createJsonResource` is **not** bound via `plapi`. It's pure Tengo at `core/platforma/sdk/workflow-tengo/src/smart.lib.tengo:1499`:
>
> ```tengo
> createJsonResource = func(value) {
>     ll.assert(value != undefined, "can't create Json value from undefined")
>     ll.assert(!ll.isStrict(value), "can't encode strict map: ", value)
>     encoded := json.encode(value)
>     return createValueResource(constants.RTYPE_JSON, encoded)
> }
> ```
>
> So the bytes of the resulting JSON resource are exactly whatever Tengo's `json.encode` produces.

**Tengo's `json.encode` does not sort map keys.** The pl backend uses a MiLab fork of Tengo (`github.com/milaboratory/tengo/v2 v2.0.0-20250410181927-9a19ce5ac955`, pinned via `replace` in `core/pl/go.mod:437`). Its `Encode` function for `*tengo.Map` at `stdlib/json/encode.go` literally does:

```go
case *tengo.Map:
    b = append(b, '{')
    len1 := len(o.Value) - 1
    idx := 0
    for key, value := range o.Value {       // ← Go's randomized map iteration. No sort.
        b = encodeString(b, key)
        ...
    }
```

`tengo.Map.Value` is a plain `map[string]Object` (`objects.go:1187`), and Go intentionally randomizes `range` iteration order over maps. The fork's own `json_test.go` only round-trips (encode → decode → compare to original), which doesn't catch byte non-determinism.

**The pl backend hashes the raw resource data bytes directly to compute a value resource's CID.** From `core/pl/platform/core/transaction/post_change.go:266-267`:

```go
mihash.MustWriteToHasher(cidHash, rTypeBytes)
mihash.MustWriteToHasher(cidHash, rChange.new.Data())   // ← raw bytes, no normalization
```

There is no JSON re-parse or canonicalization on the way in. `Data()` is whatever `json.encode` wrote.

**Conclusion:** `smart.createJsonResource(someMap)` produces non-deterministic bytes between runs → non-deterministic CID → no dedup. The `canonical` library exists precisely because of this trap. The harness's `workflow.md` § "Canonicality" warns about map iteration in template *bodies* but does not call out the JSON-resource trap explicitly — worth a reflection back to the harness.

### 4.4 The fix — `canonical.encode` + `smart.createValueResource`

Verified that `smart.createValueResource(resourceType, dataBytes)` is exported (`smart.lib.tengo:1551`) and is the lower-level primitive `createJsonResource` itself uses. So the dedup-safe pattern is:

```tengo
smart := import("@platforma-sdk/workflow-tengo:smart")
canonical := import("@platforma-sdk/workflow-tengo:canonical")
constants := import("@platforma-sdk/workflow-tengo:constants")

// Helper, define once at the top of main.tpl.tengo (or in a small lib if used in multiple places).
canonicalJsonResource := func(value) {
    return smart.createValueResource(constants.RTYPE_JSON, canonical.encode(value))
}
```

Then, at the `processTpl` call site:

```tengo
processResult := render.create(processTpl, {
    blockId: blockId,
    propertiesTsv: propertiesTsv,
    aaFractionTsv: aaFractionTsv,
    params: canonicalJsonResource({
        datasetSpec: datasetSpec,
        keyAxisIdx: keyAxisIdx,
        mode: mode,
        receptor: receptor,
        chains: chainList,        // already sorted per 4.2
        fullChains: fullChainList, // already sorted per 4.2
        hasFv: hasFv
    })
})
```

The downstream `processTpl` reads the resource via `getDataAsJson()` exactly as before — `canonical.encode` produces valid JSON that any JSON decoder will parse correctly. No changes required in `process.tpl.tengo` beyond consuming an `inputCanonicalize` audit (see 3.6).

### 4.5 `infoBlob` (model output)

Same pattern, same fix:

```tengo
infoBlob := canonicalJsonResource({
    mode: mode,
    receptor: receptor,
    coverageTier: coverageTier,
    messages: infoMessages   // sort upstream per 4.2
})
```

### 4.6 If 3.4 finds CID-poisoning mem/cpu

Switch the offending `render.create` to use `metaInputs` (third positional arg) — pattern verbatim from `core/platforma/sdk/workflow-tengo/src/render.lib.tengo:26-29, 175-178`. For `exec.builder()` and `seqTb`, if the SDK doesn't already meta-route memory/CPU, file an SDK ticket — that's a platform fix, not a block fix.

## 5. Validation

Per `feedback_verify_before_commit.md` — rebuild + reload via pl MCP + verify symptom is gone before commit.

1. **Unit-level**: Tengo test in `test/src/wf.test.ts` runs the workflow twice on identical inputs; assert `exports.properties` (the PColumn data resource handle) has the same resource-pool entry on second run. The current test runs once — extend it.
2. **Python determinism**: a pytest in `software/tests/integration/` runs `compute-properties` twice on the same `input.tsv` + `plan.json` and asserts `sha256(properties.tsv)` is identical. This catches 3.5 ahead of Tengo work.
3. **Cross-project manual**: via pl MCP — create project P1 with the block + an upstream peptide-extraction. Run to Done. Create project P2 with the same upstream graph. Add the block. Hit Run. Expected: lands on Done within seconds; no `compute-properties` workdir is spawned in P2. Inspect with `get_block_state` (`outputs.propertiesPf?.value`).
4. **Negative test**: change one byte of `plan.json` between runs (e.g., flip a chain order) and confirm a fresh Python run is triggered — proves the CID is actually distinguishing inputs and we're not just always cache-hitting.

## 6. Out of scope

- **Anonymization** (`harnesses/block-dev/anonymization.md`). The block consumes amino-acid sequences and chain labels — no sample names or patient IDs flow into pure-template inputs. Not needed.
- **Hash override** (`workflow.md` § Hash override). Reserved for semantic-stable code refactors; not relevant for new dedup wiring.
- Spec edits under `docs/text/work/projects/sequence-properties/` — left for a separate operator-driven pass.

## 7. Effort

| Step | Estimate |
|---|---|
| Sorted-map fixes in `main.tpl.tengo` (3.1, 3.2, 3.3) | 0.5 d |
| Audit + fix `process.tpl.tengo` (3.6) | 0.5 d |
| Verify mem/cpu meta-input routing (3.4) | 0.25 d |
| Python sort + quantization at write (3.5) | 0.5 – 1 d (pairs with the existing quantization plan) |
| Tests (§5 items 1, 2, 4) | 0.5 d |
| Manual cross-project validation (§5 item 3) | 0.25 d |
| **Total** | **~2.5 – 3 d** |

## 8. References

- `mictx-helper/harnesses/block-dev/workflow.md` — Pure templates and dedup, Canonicality, Hash override, Anti-patterns.
- `mictx-helper/harnesses/block-dev/anonymization.md` — confirms N/A here.
- `core/platforma/sdk/workflow-tengo/src/render.lib.tengo:26-29, 175-178` — `metaInputs` semantics.
- `core/platforma/sdk/workflow-tengo/src/maps.lib.tengo` — `getKeys()` sorted iteration.
- `core/platforma/sdk/workflow-tengo/src/canonical.lib.tengo` — canonical JSON if needed.
- `core/pl/platform/migrations/model_v000/canonical_id.go` — CID + recovery mechanic (the "land on Done immediately" path).
- `core/pl/platform/core/transaction/deduplication.go` — duplicate-resource transition.
- `docs/text/spec/resource.md` § Content Addressing & Deduplication.
- Pattern reference: `blocks/mixcr-clonotyping/workflow/src/main.tpl.tengo` (uses `render.create` for pure preset calc + `render.createEphemeral` for the heavy MiXCR step — read it before deciding pure vs ephemeral on any new sub-templates we add).
