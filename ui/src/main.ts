import { BlockLayout } from "@platforma-sdk/ui-vue";
import { createApp } from "vue";
import { sdkPlugin } from "./app";

// canary: changeset-coverage devDependency-bump fix (throwaway, do not merge)
createApp(BlockLayout).use(sdkPlugin).mount("#app");
