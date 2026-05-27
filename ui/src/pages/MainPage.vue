<script setup lang="ts">
import { createPlDataTableStateV2, type PlRef } from "@platforma-sdk/model";
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
import { computed, ref, watch } from "vue";
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

// Setter, not a watch on `inputAnchor`. The SDK replaces the `data`
// object on server patches (other-client edits, app reopen); a watcher
// would fire on that replacement and reset tableState the user did
// not touch. The setter runs only on the user's gesture.
function setInput(ref?: PlRef) {
  app.model.data.inputAnchor = ref;
  app.model.data.tableState = createPlDataTableStateV2();
}

// Session-only dismissal of info-message alerts. Local ref → not persisted;
// resets when the block UI unmounts (project close, app reload). The write
// only happens on PlAlert's close-button user gesture, so this is not a
// hairpin.
const dismissedInfoMessages = ref(new Set<string>());
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
      :model-value="true"
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
    <PlDropdownRef
      :model-value="app.model.data.inputAnchor"
      :options="app.model.outputs.inputOptions"
      label="Input dataset"
      clearable
      required
      @update:model-value="setInput"
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
