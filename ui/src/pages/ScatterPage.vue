<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnSpec } from "@platforma-sdk/model";
import { PlBlockPage } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import { defaultScatterAxes, isNumericScalar, numericScalarsInOrder } from "../utils/scalarColumns";

const app = useApp();

const dataColumnPredicate = (spec: PColumnSpec) => isNumericScalar(spec);

// R19 default + R19a fallback. Returns null while inputs are loading; null
// keeps GraphMaker in a "no defaults yet" state without applying stale picks.
const defaultOptions = computed((): PredefinedGraphOption<"scatterplot">[] | null => {
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

  // Default grouping (point colour) = the entity-key axis. Property columns
  // are 1-axis, so axesSpec[0] is the entity key — peptide variantKey for
  // peptide mode, cloneId / clonotypeKey / scClonotypeKey for the various
  // antibody/TCR modes. Each entity gets its own colour out of the box.
  const groupingAxis = x.axesSpec[0];

  return [
    { inputName: "x", selectedSource: x },
    { inputName: "y", selectedSource: y },
    { inputName: "grouping", selectedSource: groupingAxis },
  ];
});
</script>

<template>
  <PlBlockPage>
    <template #title>Scatterplot</template>
    <GraphMaker
      v-model="app.model.data.graphStateScatter"
      chart-type="scatterplot"
      :p-frame="app.model.outputs.propertiesPfHandle"
      :default-options="defaultOptions"
      :default-palette="{ categorical: 'bright' }"
      :data-column-predicate="dataColumnPredicate"
      :status-text="{
        noPframe: { title: 'Select an input dataset on the Main tab to plot.' },
      }"
    />
  </PlBlockPage>
</template>
