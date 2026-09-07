---
name: spatial-alignment-compose
description: Compose confirmed spatial alignment recipes into corrected coordinate tables with complete edit history and remote/local consistency checks. Use when Codex needs to apply user-confirmed pairwise, spatial-only, ROI/drop, or manual alignment corrections to downstream slice blocks, write manual_x/manual_y/manual_z, record every operation, generate recipe_chain.json/correction_history.csv/apply_summary.json/input_file_manifest.json, replay locally or remotely, or compare sha256/numeric hashes for reproducibility.
---

# Spatial Alignment Compose

## Package paths

Run the relative commands below from this skill directory. Install this skill
beside `spateo-2d-alignment` (or keep the complete repository checkout).
The shared support tools are the fixed pre-zebrafish snapshot, located in
`../spateo-2d-alignment/pipelines/pairwise-rigid`. Verify entrypoints with its `skill.lock.yaml`.
This packaging restoration does not change alignment algorithms or defaults.


## Overview

Use this skill only after the user confirms one or more alignment candidates.
Composition is non-destructive: preserve original coordinate columns in an
internal audit table and export a separate clean user-facing coordinate CSV.

Recipe creation and remote reference replay happen on the remote server.
Local replay is allowed and expected when the local machine has an immutable
copy of the original coordinate cache. The local and remote outputs must be
generated from the same `recipe_chain.json`, the same transform code version,
the same input coordinate columns, and the same output formatting rules.

Never overwrite the original local coordinate cache. Local compose must write
to a new output directory, preserve input hashes, and compare against the remote
reference output before claiming consistency.

## Invariants

- `recipe_chain.json` is the replay source of truth.
- Source h5ad identity is the source of truth for final cell ids. Preserve or
  restore ids from original `obs_names`, `obs["cell_id"]`, `obs["CellID"]`, and
  slice annotations; never invent clean `cell_id` values during compose.
- Remote execution is mandatory for Spateo/ROI recipe generation and for the
  reference replay.
- Local deterministic replay is allowed from a read-only local coordinate cache
  when the same recipe chain and transform code are used.
- Do not claim reproducibility unless replay check succeeds.
- Record every operation in order, including normal pairwise, spatial-only
  rescue, ROI/drop refine, and manual recipes.
- Apply edge corrections to the intended downstream block, not only the moving
  slice, when the user requests whole-embryo consistency.
- Do not generate HTML here unless the user asks; use `spatial-pointcloud-viewer`
  for viewers.
- Do not overwrite source coordinate files. Use new local and remote output
  directories for every compose/replay.
- Do not treat local and remote outputs as equivalent until checksum or numeric
  replay comparison passes.
- Do not hand off `audit_points.csv` or `corrected_points.csv` as the latest
  user coordinate table. Share `exports/*.clean_coordinates.csv`.
- If input rows use internal ids such as `slice_id:0`, keep those as
  `workflow_cell_id` in audit outputs and pass an h5ad ID map during clean
  export so user-facing `cell_id` values are h5ad-derived ids such as
  `SL38_CELL.18`.

## Compose Command

Use the shared compose script with recipes listed in execution order. First run
it on the remote server to create the reference output:

```bash
python ../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/core/compose_alignment_recipes.py \
  --input SOURCE_POINTS_CSV \
  --operation RECIPE_001.json \
  --operation RECIPE_002.json \
  --x-col X_COL --y-col Y_COL --z-col Z_COL \
  --slice-col sl_number --celltype-col celltype --cell-id-col cell_id \
  --out-x-col manual_x --out-y-col manual_y --out-z-col manual_z \
  --output-dir COMPOSE_OUTPUT_DIR \
  --replay-check
```

Replay an existing chain on the remote server:

```bash
python ../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/core/compose_alignment_recipes.py \
  --input SOURCE_POINTS_CSV \
  --recipe-chain COMPOSE_OUTPUT_DIR/recipe_chain.json \
  --x-col X_COL --y-col Y_COL --z-col Z_COL \
  --slice-col sl_number --celltype-col celltype --cell-id-col cell_id \
  --output-dir REPLAY_OUTPUT_DIR \
  --replay-check
```

Use the current final coordinate columns as the next input columns only when
intentionally composing on a previously corrected remote table. Prefer
replaying from the original remote source CSV plus the full recipe chain when
producing final deliverables.

Then run the same command locally only if the local input cache is a preserved
copy of the same source coordinate table. Use a new local output directory and
do not modify the original cache.

## Required Outputs

Every compose run must produce:

```text
audit_points.csv
corrected_points.csv        # compatibility hardlink/copy of audit_points.csv
clean_coordinates.csv       # or an explicit exports/*.clean_coordinates.csv
clean_export_manifest.json
changed_cells.csv
recipe_chain.json
correction_history.csv
apply_summary.json
input_file_manifest.json
```

The clean CSV schema is fixed:

```text
cell_id,slice_id,stage,chip_id,sl_number,celltype,x,y,z
```

`cell_id` in this schema must be derived from the original h5ad annotations,
typically `<slice>_<CellID>`. It must not be a generated zero-based row index.

Intermediate columns, edit labels, displacement, and component ranks stay in
`audit_points.csv`, `changed_cells.csv`, or QC sidecars.

`recipe_chain.json` must include, for each operation:

- operation id and type: `pairwise_transform`, `spatial_only_rescue`,
  `roi_drop_refine`, or `manual_recipe`;
- source recipe path and hash;
- fixed slice, moving slice, and affected slices;
- input and output coordinate columns;
- transform steps or transform hash;
- expression mode and sigma settings;
- ROI/drop condition when present;
- timestamp and user confirmation note when available.

`correction_history.csv` is the human-readable audit table. It must include:

- step order;
- edge id;
- operation type;
- affected slices;
- points affected;
- coordinate columns written;
- source run directory;
- before/after metrics when available.

`apply_summary.json` must include total points, modified points per step, final
coordinate columns, final slice list, source input hash, recipe chain hash, and
replay status.

`input_file_manifest.json` must include source coordinate file metadata and all
recipe/provenance/audit inputs with path, size, mtime, and sha256 when feasible.

## Consistency Check

After local and remote replay both finish, compare results:

- Prefer full-file `sha256` when the same code version and CSV formatting are
  used on both sides.
- If file hashes differ only because of row endings or float formatting, compare
  stable hashes of `manual_x/manual_y/manual_z` after fixed rounding and report
  `max_abs_delta`.
- Record the comparison in `local_remote_consistency.json` with local path,
  remote path, recipe chain hash, transform code hash, input hashes, file
  hashes, coordinate hashes, max absolute delta, tolerance, and pass/fail.
- Do not claim "consistent" unless this check passes.

## Current Recipe Shape

The compose script accepts existing recipe files with `steps` such as:

```json
{
  "steps": [
    {
      "operation": "rotate_z",
      "scope": {"slice_labels": ["193", "196", "199"]},
      "parameters": {"degrees": 14.55},
      "center": {"x": 0, "y": 0, "z": 0}
    },
    {
      "operation": "translate",
      "scope": {"slice_labels": ["193", "196", "199"]},
      "parameters": {"dx": 434.8, "dy": -840.7, "dz": 0}
    }
  ]
}
```

Supported operations are `translate`, `rotate_z`, `scale`, and 2D affine
matrix operations. Add new operation types to the script before using them in a
final recipe chain.

## Viewer Handoff

After compose succeeds, generate local inspection HTML with
`spatial-pointcloud-viewer` from either the local replay output that passed the
remote consistency check or a copied remote final cache. Verify:

- `sampled=false`;
- total points match the intended full table;
- slices are populated;
- the viewer uses the final coordinate columns.

For manual issue finding or acceptance, render a full-points viewer from the
clean coordinates or component-labeled full table and name it
`review_full_points` or `candidate_full_points`. For post-edit structure checks
and repeated candidate iteration, render a balanced300k CSV created from the
clean coordinates.
