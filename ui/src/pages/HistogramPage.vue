<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnSpec } from "@platforma-sdk/model";
import { computed } from "vue";
import { useApp } from "../app";
import {
  defaultHistogramMetric,
  isNumericScalar,
  numericScalarsInOrder,
} from "../utils/scalarColumns";

const app = useApp();

// R20 default + R20a fallback.
const defaultOptions = computed((): PredefinedGraphOption<"histogram">[] | null => {
  const cols = app.model.outputs.propertiesPfCols;
  const info = app.model.outputs.info;
  if (!cols || !info) return null;

  const metric = defaultHistogramMetric(cols, info) ?? numericScalarsInOrder(cols)[0];
  if (!metric) return null;

  return [{ inputName: "value", selectedSource: metric }];
});

// Data = own scalar properties. propertiesPfHandle holds our trace-injected
// pCols plus single-axis upstream metadata; the trace match isolates ours
// to drive the default value-axis pick.
const dataColumnPredicate = (spec: PColumnSpec) =>
  isNumericScalar(spec) &&
  spec.annotations?.["pl7.app/trace"]?.includes("sequence-properties") === true;

// Meta = every column in the pframe. The model layer already curates the
// set (own scalars + single-axis upstream metadata; aaFraction excluded for
// the cell-count guard), so anything that reaches us is a valid dimension
// for Filter / Grouping-Color / Highlight / Size / Tab / Tooltip / Label /
// Additional-curves. Own scalars deliberately appear in both data and meta
// roles — users bin the histogram by, e.g., chain while plotting Aromaticity
// values. If multi-axis columns later enter the pframe, add an
// `isSingleAxis` guard here.
const metaColumnPredicate = (_spec: PColumnSpec) => true;
</script>

<template>
  <GraphMaker
    v-model="app.model.data.graphStateHistogram"
    chart-type="histogram"
    :p-frame="app.model.outputs.propertiesPfHandle"
    :default-options="defaultOptions"
    :data-column-predicate="dataColumnPredicate"
    :meta-column-predicate="metaColumnPredicate"
    :status-text="{
      noPframe: { title: 'Select an input dataset on the Main tab to plot.' },
    }"
  />
</template>
