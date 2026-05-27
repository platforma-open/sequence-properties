import { platforma } from "@platforma-open/milaboratories.sequence-properties.model";
import { defineAppV3 } from "@platforma-sdk/ui-vue";
import { ref, watchEffect } from "vue";
import HistogramPage from "./pages/HistogramPage.vue";
import MainPage from "./pages/MainPage.vue";
import ScatterPage from "./pages/ScatterPage.vue";

// Module-level singleton — session-only dismissal of info-message alerts.
// Survives in-block route changes (Main ↔ Property Relationships ↔
// Property Distribution) and switching away to other blocks (desktop
// LRU-caches recently-used block UIs, limit 4). Resets on project close,
// block reload, LRU eviction after opening several other blocks, or app
// restart. Hairpin-safe: writes only happen on user gesture from
// PlAlert's close-button emit.
export const dismissedInfoMessages = ref(new Set<string>());

export const sdkPlugin = defineAppV3(platforma, (app) => {
  app.model.data.customBlockLabel ??= "";

  watchEffect(() => {
    const anchor = app.model.data.inputAnchor;
    const opts = app.model.outputs.inputOptions ?? [];
    const match = anchor
      ? opts.find((o) => o.ref?.blockId === anchor.blockId && o.ref?.name === anchor.name)
      : undefined;
    app.model.data.defaultBlockLabel = match?.label ?? "";
  });

  return {
    routes: {
      "/": () => MainPage,
      "/scatter": () => ScatterPage,
      "/histogram": () => HistogramPage,
    },
  };
});

export const useApp = sdkPlugin.useApp;
