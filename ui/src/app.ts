import { platforma } from "@platforma-open/milaboratories.sequence-properties.model";
import { defineAppV3 } from "@platforma-sdk/ui-vue";
import { watchEffect } from "vue";
import HistogramPage from "./pages/HistogramPage.vue";
import MainPage from "./pages/MainPage.vue";
import ScatterPage from "./pages/ScatterPage.vue";

export const sdkPlugin = defineAppV3(platforma, (app) => {
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
