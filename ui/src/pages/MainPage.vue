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
import { ref } from "vue";
import { useApp } from "../app";

const app = useApp();

const logOpen = ref(false);

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
    </template>

    <PlDropdownRef
      v-model="app.model.data.inputAnchor"
      :options="app.model.outputs.inputOptions"
      label="Input dataset"
    >
      <template #tooltip>
        <div>
          <strong>Sequence dataset</strong><br />
          Peptide extraction or MiXCR clonotyping output. The block detects modality (peptide vs
          antibody/TCR) from the dataset's axes and computes physico-chemical properties
          accordingly.
        </div>
      </template>
    </PlDropdownRef>

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

  <PlSlideModal v-model="logOpen" width="80%">
    <template #title>Processing Log</template>
    <PlLogView :log-handle="app.model.outputs.processingLog" />
  </PlSlideModal>
</template>
