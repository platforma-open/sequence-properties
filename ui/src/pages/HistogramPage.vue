<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnSpec } from "@platforma-sdk/model";
import { PlBlockPage } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import {
  defaultHistogramMetric,
  isNumericScalar,
  numericScalarsInOrder,
} from "../utils/scalarColumns";

const app = useApp();

const dataColumnPredicate = (spec: PColumnSpec) => isNumericScalar(spec);

// R20 default + R20a fallback.
const defaultOptions = computed((): PredefinedGraphOption<"histogram">[] | null => {
  const cols = app.model.outputs.propertiesPfCols;
  const info = app.model.outputs.info;
  if (!cols || !info) return null;

  const metric = defaultHistogramMetric(cols, info) ?? numericScalarsInOrder(cols)[0];
  if (!metric) return null;

  return [{ inputName: "value", selectedSource: metric }];
});
</script>

<template>
  <PlBlockPage>
    <template #title>Histogram</template>
    <GraphMaker
      v-model="app.model.data.graphStateHistogram"
      chart-type="histogram"
      :p-frame="app.model.outputs.propertiesPfHandle"
      :default-options="defaultOptions"
      :data-column-predicate="dataColumnPredicate"
      :status-text="{
        noPframe: { title: 'Select an input dataset on the Main tab to plot.' },
      }"
    />
  </PlBlockPage>
</template>
