<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnSpec } from "@platforma-sdk/model";
import { computed } from "vue";
import { useApp } from "../app";
import { defaultScatterAxes, isNumericScalar, numericScalarsInOrder } from "../utils/scalarColumns";

const app = useApp();

// R19 default + R19a fallback. Returns null while inputs are loading; null
// keeps GraphMaker in a "no defaults yet" state without applying stale picks.
const defaultOptions = computed((): PredefinedGraphOption<"scatterplot-umap">[] | null => {
  const cols = app.model.outputs.propertiesPfCols;
  const info = app.model.outputs.info;
  if (!cols || !info) return null;

  const { x: defX, y: defY } = defaultScatterAxes(cols, info);
  const scalars = numericScalarsInOrder(cols);
  // Fall back to emission order when the modality default is absent (e.g.
  // light-chain-only input, or full-chain-only mode without CDR3).
  const x = defX ?? scalars[0];
  const y = defY ?? (scalars[1] === defX ? scalars[0] : scalars[1]);
  if (!x || !y || x === y) return null;

  // No default grouping. Grouping by the entity-key axis assigns one colour
  // per entity — thousands of categorical groups explode the palette,
  // legend, and render path on real datasets. The picker still exposes
  // sample-axis and entity-axis meta cols; we leave the choice to the user.
  return [
    { inputName: "x", selectedSource: x },
    { inputName: "y", selectedSource: y },
  ];
});

// Data = own scalar properties. propertiesPfHandle holds our trace-injected
// pCols plus single-axis upstream metadata; the trace match isolates ours
// to drive default X/Y picks.
const dataColumnPredicate = (spec: PColumnSpec) =>
  isNumericScalar(spec) &&
  spec.annotations?.["pl7.app/trace"]?.includes("sequence-properties") === true;

// Meta = every column in the pframe. The model layer already curates the
// set (own scalars + single-axis upstream metadata; aaFraction excluded for
// the cell-count guard), so anything that reaches us is a valid dimension
// for Filter / Grouping-Color / Highlight / Size / Tab / Tooltip / Label /
// Additional-curves. Own scalars deliberately appear in both data and meta
// roles — users color the scatter by, e.g., Aromaticity while plotting
// Charge vs Hydrophobicity. If multi-axis columns later enter the pframe,
// add an `isSingleAxis` guard here.
const metaColumnPredicate = (_spec: PColumnSpec) => true;
</script>

<template>
  <GraphMaker
    v-model="app.model.data.graphStateScatter"
    chart-type="scatterplot-umap"
    :p-frame="app.model.outputs.propertiesPfHandle"
    :default-options="defaultOptions"
    :default-palette="{ categorical: 'bright' }"
    :data-column-predicate="dataColumnPredicate"
    :meta-column-predicate="metaColumnPredicate"
    :status-text="{
      noPframe: { title: 'Select an input dataset on the Main tab to plot.' },
    }"
  />
</template>
