import { platforma } from "@platforma-open/milaboratories.sequence-properties.model";
import { defineAppV3 } from "@platforma-sdk/ui-vue";
import { ref, watch, watchEffect } from "vue";
import HistogramPage from "./pages/HistogramPage.vue";
import MainPage from "./pages/MainPage.vue";
import ScatterPage from "./pages/ScatterPage.vue";

// Module-level singleton — session-only dismissal of info-message alerts.
// Tied to currently-firing advisories: an entry stays here only while its
// string is in `outputs.info.messages`. The `watch` inside `sdkPlugin`
// prunes entries that the workflow no longer emits, so a message that
// "naturally" disappears doesn't carry over its dismissal to a future
// re-fire. Survives in-block route changes and switching away to other
// blocks (desktop LRU-caches block UIs, limit 4). Resets on project close,
// block reload, LRU eviction, or app restart. Hairpin-safe: writes happen
// only via user gesture (PlAlert close emit) or output → local-ref prune,
// both sanctioned by `harnesses/block-dev/hairpin.md`.
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

  // Auto-prune the dismissed set to the currently-firing messages. When
  // the workflow stops emitting a previously-dismissed string (re-run
  // with a different input, advisory naturally goes away), drop it from
  // the set so a future re-fire shows fresh. Output → local Vue ref is a
  // sanctioned pattern per `hairpin.md` (state lives in the JS context,
  // not BlockData; no multi-client interleave).
  watch(
    () => app.model.outputs.info?.messages,
    (messages) => {
      const fired = new Set(messages ?? []);
      const current = dismissedInfoMessages.value;
      const filtered = [...current].filter((m) => fired.has(m));
      if (filtered.length !== current.size) {
        dismissedInfoMessages.value = new Set(filtered);
      }
    },
    { immediate: true },
  );

  return {
    routes: {
      "/": () => MainPage,
      "/scatter": () => ScatterPage,
      "/histogram": () => HistogramPage,
    },
  };
});

export const useApp = sdkPlugin.useApp;
