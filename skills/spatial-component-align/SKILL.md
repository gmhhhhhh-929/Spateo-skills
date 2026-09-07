---
name: spatial-component-align
description: Build non-destructive component-wise and component-group spatial alignment candidates from connected-component labels and coordinate tables. Use when Codex needs to match tissue islands across adjacent spatial transcriptomics slices, estimate per-component or per-group rigid transforms, use small anchor components to move larger blocks, preview piecewise alignment, or prepare a recipe for later full replay without overwriting source data.
---

# Spatial Component Align

## Package paths

Run the relative commands below from this skill directory. Install this skill
beside `spateo-2d-alignment` (or keep the complete repository checkout).
The shared support tools are the fixed pre-zebrafish snapshot, located in
`../spateo-2d-alignment/pipelines/pairwise-rigid`. Verify entrypoints with its `skill.lock.yaml`.
This packaging restoration does not change alignment algorithms or defaults.


## Scope

Use this skill after `detect-spatial-components` has produced per-slice
component labels. It creates candidate component/group transforms and review
artifacts only. It must not overwrite source h5ad files, clean coordinate CSVs,
or accepted/final states.

This is the middle layer:

- component detection happens elsewhere;
- component groups and anchors are defined here;
- preview transforms and recipes are written here;
- accepted recipes are replayed to full data with `spatial-alignment-compose`.

## Locked Entrypoint

Use the active release through `../spateo-2d-alignment/pipelines/pairwise-rigid/skill.lock.yaml`:

```text
../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/components/component_group_align_preview.py
```

Record `component_group_align` entrypoint sha256 in
`provenance/skill_run.json`. Stop on hash mismatch.

## Group Plan

Treat components as explicit groups. A group can estimate a transform from one
anchor component and apply that transform to a larger moving block.

```json
{
  "fixed_slice": 59,
  "moving_slice": 62,
  "groups": [
    {
      "group_id": "left_body",
      "fixed_components": [1],
      "moving_components": [1],
      "estimate_from_fixed": [1],
      "estimate_from_moving": [1],
      "apply_to_moving": [1],
      "component_key": "component_rank",
      "method": "icp"
    },
    {
      "group_id": "right_block_with_island_anchor",
      "fixed_components": [2, 7],
      "moving_components": [2, 5],
      "estimate_from_fixed": [7],
      "estimate_from_moving": [5],
      "apply_to_moving": [2, 5],
      "component_key": "component_id",
      "method": "icp"
    }
  ]
}
```

Human notes must name components as `<slice> rank<rank> id<component_id>`, for
example `SL59 rank2 id7`. Component ids are per-slice, not global.

## Workflow

1. Confirm the candidate scope with the user: fixed slice, moving slice,
   component rank/id, and whether an anchor should also move a larger block.
2. Inspect the full-points component viewer first. Use balanced300k only for
   post-edit speed checks and repeated tuning.
3. Build a group plan JSON with explicit `estimate_from_*` and `apply_to_*`
   lists.
4. Run `component_group_align_preview.py` on the chosen coordinate table.
   Balanced300k is acceptable for interactive tuning; full replay is required
   before acceptance.
5. Write preview points, metrics, recipe, and review HTML/summary inside the
   step directory.
6. Mark the result as `candidate` until the user accepts it.

## Required Outputs

```text
recipes/component_group_recipe.json
qc/component_group_metrics.json
light/component_group_preview_points.csv
review/candidate_balanced300k.viewer.html
decision.md
record.json
provenance/skill_run.json
```

For a review milestone, replay the recipe to full data and export:

```text
full/<sample>.<step>.audit_points.csv
exports/<sample>.<step>.clean_coordinates.csv
qc/<sample>.<step>.full_replay_validation.json
```

## Guardrails

- Do not apply a small-anchor transform to an entire slice unless the group
  plan explicitly says so.
- Do not infer component identity from rank alone when ids are available.
- Do not propagate component transforms across long slice ranges without
  explicit user acceptance.
- Do not use component-level candidates as final data until full replay,
  changed-cell validation, and celltype preservation checks pass.
- Do not renumber or invent `cell_id` values.
