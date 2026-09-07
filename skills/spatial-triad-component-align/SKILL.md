---
name: spatial-triad-component-align
description: Scan three-slice spatial transcriptomics windows, flag suspicious middle-slice connected-component outliers, generate triad review HTML, and create conservative mid-slice rigid preview candidates. Use when Codex needs to review prev-mid-next slice consistency, avoid aggressive pairwise component transforms, mark problematic triad segments, or build non-destructive triad QC artifacts from component labels and coordinate tables.
---

# Spatial Triad Component Align

## Package paths

Run the relative commands below from this skill directory. Install this skill
beside `spateo-2d-alignment` (or keep the complete repository checkout).
The shared support tools are the fixed pre-zebrafish snapshot, located in
`../spateo-2d-alignment/pipelines/pairwise-rigid`. Verify entrypoints with its `skill.lock.yaml`.
This packaging restoration does not change alignment algorithms or defaults.


## Scope

Use this skill after per-slice component labels exist. It reviews consecutive
`prev-mid-next` windows to decide whether the middle slice or middle component
is inconsistent with both neighbors. It is conservative and does not write
final data by itself.

Triad review is useful when pairwise component edits feel too aggressive:
sometimes two adjacent slices look bad because the middle slice is damaged or
contracted, and moving both edges independently makes the 3D structure worse.

## Locked Entrypoints

Use the active release through `../spateo-2d-alignment/pipelines/pairwise-rigid/skill.lock.yaml`:

```text
../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/components/triad_component_scan.py
../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/viewers/make_triad_component_review_viewer.py
../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/components/triad_mid_slice_rigid_preview.py
```

Record `triad_component_scan`, `triad_component_review_viewer`, and
`triad_mid_slice_rigid_preview` sha256 values in provenance. Stop on hash
mismatch.

## Workflow

1. Start from the current clean slice-level or review coordinate state selected
   by `states/current_review.yaml` or `states/slice_level_baseline.yaml`.
2. Use component labels generated from the same coordinate state, proven by
   input hash.
3. Run `triad_component_scan.py` on all consecutive slice triplets.
4. Generate `review/triad_component_review.html`.
5. Inspect flagged/watch triads visually. Do not automatically correct them.
6. If the user confirms a middle-component outlier, generate a rigid preview
   that moves only that middle component.
7. Replay accepted preview recipes to full data with `spatial-alignment-compose`.

## Interpretation

- `flagged`: the middle component is a strong outlier relative to both
  neighbors.
- `watch`: possible middle outlier or one-edge issue; inspect visually.
- `normal`: no conservative evidence for a middle-component issue.

Only `mid_outlier` rows are eligible for middle-component rigid preview. Do not
treat an `edge_issue` row as permission to move the middle slice/component.

## Required Outputs

For scan/review:

```text
qc/triad_component_scan.csv
qc/triad_component_scan_summary.json
review/triad_component_review.html
review/triad_component_review_summary.json
record.json
step.yaml
decision.md
```

For a rigid preview candidate:

```text
light/triad_mid_slice_rigid_preview.csv
recipes/triad_mid_slice_rigid_preview_recipe.json
review/before_after.html
qc/triad_mid_slice_rigid_preview_metrics.json
```

## Guardrails

- Do not move prev/next slices in a triad preview.
- Do not move the whole middle slice unless the user explicitly asks for a
  slice-level operation.
- Do not inherit old aggressive component edits when the user asks to review a
  slice-level-only baseline.
- Do not include triad outputs in final h5ad without explicit acceptance and a
  replayable full-data recipe chain.
