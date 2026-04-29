<script setup lang="ts">
import {
  PlAgDataTableV2,
  PlAlert,
  PlBlockPage,
  PlDropdownRef,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { useApp } from "../app";

const app = useApp();

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.propertiesTable,
});
</script>

<template>
  <PlBlockPage>
    <template #title>Sequence Properties</template>

    <PlDropdownRef
      v-model="app.model.data.inputAnchor"
      :options="app.model.outputs.inputOptions"
      label="Input dataset"
    />

    <template v-if="app.model.outputs.info">
      <PlAlert
        v-for="(message, idx) in app.model.outputs.info?.messages ?? []"
        :key="idx"
        type="info"
      >
        {{ message }}
      </PlAlert>
    </template>

    <PlAgDataTableV2 v-model="app.model.data.tableState" :settings="tableSettings" />
  </PlBlockPage>
</template>
