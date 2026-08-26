import { assertParamsObject, defineBlockKind } from "@platforma-sdk/block-kind";
import type { PlRef } from "@platforma-sdk/model";
import { isPlRef } from "@platforma-sdk/model";
import { name, version } from "../package.json" with { type: "json" };

/**
 * This block's init-params contract — the upstream dataset a new instance
 * computes sequence properties on.
 *
 * That is the whole contract. Everything else in the model's `BlockData` is
 * view state the UI owns and always defaults: the table state, the two
 * graph-maker states, and the two label fields (`defaultBlockLabel` is written
 * by the UI from the chosen input's option label, `customBlockLabel` is the
 * subtitle the user types).
 *
 * `inputAnchor` is optional because the projection hands live state back
 * untouched, and a block whose input is not picked yet holds `undefined` there.
 * Requiring it would make the block export a file its own kind refuses to
 * apply, so export and apply would stop being inverses.
 */
export type BlockParams = {
  inputAnchor?: PlRef;
};

/**
 * The same contract at runtime, for params that arrive from a template file
 * rather than from typed code.
 *
 * A readable `{ block, name }` reference is expanded to a full `PlRef` before
 * this runs, so `isPlRef` is the only shape to accept here.
 */
function parseInitializationParams(value: unknown): BlockParams {
  assertParamsObject(value);

  const { inputAnchor } = value;
  if (inputAnchor !== undefined && !isPlRef(inputAnchor)) {
    throw new Error(
      "'inputAnchor' must be a reference to an upstream column, written as { block, name }.",
    );
  }

  return { inputAnchor };
}

// Identity (`name`/`version`) comes from this package's own `package.json`, so
// the on-wire `{name}@{version}` reference can never drift from what npm
// publishes; the bundler inlines the JSON import.
export const kind = defineBlockKind<BlockParams>({
  name,
  version,
  parseInitializationParams,
});
