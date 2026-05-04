<script setup lang="ts">
import {
  PlAgDataTableV2,
  PlAlert,
  PlBlockPage,
  PlBtnGhost,
  PlDropdownRef,
  PlLogView,
  PlMaskIcon24,
  PlSlideModal,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { ref, watch } from "vue";
import { useApp } from "../app";

const app = useApp();

const logOpen = ref(false);
const settingsOpen = ref(app.model.data.inputAnchor === undefined);

watch(
  () => app.model.outputs.isRunning,
  (isRunning) => {
    if (isRunning) settingsOpen.value = false;
  },
);

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.propertiesTable,
});
</script>

<template>
  <PlBlockPage>
    <template #title>Sequence Properties</template>
    <template #append>
      <PlBtnGhost @click.stop="() => (logOpen = true)">
        Logs
        <template #append>
          <PlMaskIcon24 name="file-logs" />
        </template>
      </PlBtnGhost>
      <PlBtnGhost @click.stop="() => (settingsOpen = true)">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>

    <PlAlert
      v-for="(message, idx) in app.model.outputs.info?.messages ?? []"
      :key="idx"
      type="info"
    >
      {{ message }}
    </PlAlert>

    <PlAgDataTableV2
      v-model="app.model.data.tableState"
      :settings="tableSettings"
      not-ready-text="Select an input dataset to compute sequence properties."
      show-export-button
    />
  </PlBlockPage>

  <PlSlideModal v-model="settingsOpen" close-on-outside-click shadow>
    <template #title>Settings</template>
    <PlDropdownRef
      v-model="app.model.data.inputAnchor"
      :options="app.model.outputs.inputOptions"
      label="Input dataset"
      clearable
      required
    >
      <template #tooltip>
        Peptide extraction or MiXCR clonotyping output. Modality is auto-detected.
      </template>
    </PlDropdownRef>
  </PlSlideModal>

  <PlSlideModal v-model="logOpen" width="80%">
    <template #title>Processing Log</template>
    <PlLogView :log-handle="app.model.outputs.processingLog" />
  </PlSlideModal>
</template>
