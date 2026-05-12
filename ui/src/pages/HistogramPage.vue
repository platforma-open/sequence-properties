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

const dataColumnPredicate = (spec: PColumnSpec) =>
  isNumericScalar(spec) &&
  spec.annotations?.["pl7.app/isOutput"] === "true" &&
  !spec.annotations?.["pl7.app/trace"]?.includes("sequence-properties");

const metaColumnPredicate = (spec: PColumnSpec) =>
  !spec.annotations?.["pl7.app/trace"]?.includes("sequence-properties");
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
