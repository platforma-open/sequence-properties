<script setup lang="ts">
import type { DatasetSelection } from "@platforma-sdk/model";
import {
  createDatasetSelection,
  createPlDataTableStateV2,
  createPrimaryRef,
} from "@platforma-sdk/model";
import {
  PlAgDataTableV2,
  PlAlert,
  PlBlockPage,
  PlBtnGhost,
  PlDatasetSelector,
  PlLogView,
  PlMaskIcon24,
  PlSlideModal,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { computed, ref, watch } from "vue";
import { dismissedInfoMessages, useApp } from "../app";

const app = useApp();

const logOpen = ref(false);
const settingsOpen = ref(app.model.data.inputAnchor === undefined);

watch(
  () => app.model.outputs.isRunning,
  (isRunning) => {
    if (isRunning) settingsOpen.value = false;
  },
);

// Setter, not a watch on `inputAnchor`. The SDK replaces the `data`
// object on server patches (other-client edits, app reopen); a watcher
// would fire on that replacement and reset tableState the user did
// not touch. The setter runs only on the user's gesture.
const inputSelection = computed<DatasetSelection | undefined>({
  get: () => {
    const { inputAnchor, filterRef } = app.model.data;
    if (inputAnchor === undefined) return undefined;
    return createDatasetSelection(createPrimaryRef(inputAnchor, filterRef));
  },
  set: (selection) => {
    app.model.data.inputAnchor = selection?.primary.column;
    app.model.data.filterRef = selection?.primary.filter;
    app.model.data.tableState = createPlDataTableStateV2();
  },
});

// Source-of-truth for session-only dismissal lives at module scope in
// `app.ts` so it survives in-block route changes (Main ↔ Property
// Relationships ↔ Property Distribution). See app.ts for the lifecycle
// notes.
const visibleInfoMessages = computed(() =>
  (app.model.outputs.info?.messages ?? []).filter((m) => !dismissedInfoMessages.value.has(m)),
);
function dismissInfoMessage(message: string) {
  dismissedInfoMessages.value = new Set(dismissedInfoMessages.value).add(message);
}

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.propertiesTable,
});
</script>

<template>
  <PlBlockPage
    v-model:subtitle="app.model.data.customBlockLabel"
    :subtitle-placeholder="app.model.data.defaultBlockLabel"
  >
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
      v-for="message in visibleInfoMessages"
      :key="message"
      type="info"
      closeable
      @update:model-value="() => dismissInfoMessage(message)"
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
    <PlDatasetSelector
      v-model="inputSelection"
      :options="app.model.outputs.inputOptions"
      label="Input dataset"
      clearable
      required
    >
      <template #tooltip>
        Peptide extraction or MiXCR clonotyping output. Modality is auto-detected.
      </template>
    </PlDatasetSelector>
  </PlSlideModal>

  <PlSlideModal v-model="logOpen" width="80%">
    <template #title>Processing Log</template>
    <PlLogView :log-handle="app.model.outputs.processingLog" />
  </PlSlideModal>
</template>
