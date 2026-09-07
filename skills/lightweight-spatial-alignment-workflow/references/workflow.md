# Lightweight Spatial Alignment Workflow v0.2.3

This workflow keeps alignment and primary manual review tied to full
coordinates. The component-balanced 300k representation is used after an edit
or for fast candidate iteration, not as the first place to decide whether a
structure is wrong.

## Step Order

1. Run full adjacent slice pairwise rigid alignment with locked code.
2. QC the pairwise result by metrics and human review.
3. Run spatial-only rigid rescue only for visibly failed edges.
4. Compose accepted slice-level transforms from source into a full baseline.
5. Detect connected components once on that full baseline and keep labels that
   the HTML viewer can use for Component color mode and in-view labels.
6. Generate a full-points review viewer from the slice-level baseline for
   human issue finding and component decisions.
7. Create `balanced300k` with `target_total_points=300000`, sampling by
   `slice x component_id` unless a project explicitly needs component rank.
8. After a candidate edit, use balanced300k viewers for fast all-slice, pair,
   triad, and component-focused structure checks. Standard viewers must include
   Component color mode when component labels are available, and should show
   component labels inside the HTML view itself.
9. Store component edits as recipes/deltas and tune them on balanced300k when
   speed matters.
10. Replay review-ready or accepted recipes to the full coordinate table.
11. Export clean user-facing coordinates from the full audit table. If the
    audit table still uses workflow row-index ids, pass an h5ad ID map so
    clean `cell_id` values are short h5ad-derived ids such as `SL14_CELL.1`.
12. Validate full replay by row count, input/output sha256, changed rows,
    changed components, and cell type preservation.

## Naming

- `review_full_points.viewer.html`: primary manual review viewer for deciding
  whether a baseline region/component is wrong.
- `review_balanced300k.viewer.html`: post-edit structure review viewer.
- `candidate_balanced300k.viewer.html`: candidate correction viewer.
- `review/review_balanced300k.viewer.html`: lightweight component map in HTML,
  with Component color mode and `Labels` overlay, used after edits or during
  fast candidate iteration.
- `qc/figures/slice_<slice>_components.png`: detector QC sidecar for checking
  whether segmentation itself looks reasonable.
- `audit_points.csv`: internal full replay table with intermediate/debug columns.
- `clean_coordinates.csv`: user-facing coordinate table with only identity and final coordinates.
- `full_replay_validation.json`: full row count, hashes, changed rows, and
  recipe-chain validation.

Do not call a 300k viewer `full_preview`.

## Project Layout

```text
spatial_alignment_projects/<sample_id>/runs/<run_id>/
  run_manifest.json
  run_manifest.yaml
  states/
    slice_level_baseline.yaml
    current_review.yaml
    current_accepted.yaml
    final.yaml
  steps/
    001_<actual_operation_slug>/
      step.yaml
      record.json
      decision.md
      provenance/
      light/
      full/
      exports/
      review/
      qc/
      recipes/
      logs/
  lineage/
    accepted_steps.jsonl
    all_operations.jsonl
  final/
```

Step names describe actual operations. QC, compose, metrics, logs, and viewers
are artifacts inside a step unless the analysis itself is the operation.

See `docs/project_layout.md` for the full v0.2.3 directory contract.

## Clean Coordinate Contract

Every full replay milestone writes an internal audit table and a clean export:

```text
full/<sample>.<step>.audit_points.csv
exports/<sample>.<step>.clean_coordinates.csv
exports/<sample>.<step>.clean_export_manifest.json
qc/<sample>.<step>.changed_cells.csv
```

The clean CSV schema is fixed:

```text
cell_id,slice_id,stage,chip_id,sl_number,celltype,x,y,z
```

In this schema, `cell_id` means `<slice>_<CellID>` from the h5ad obs metadata,
for example `SL14_CELL.1`, not a generated workflow id. Preserve the original
`CELL.n` values and gaps. Legacy ids such as `CS13_SL35_Y01636D4:0` must be
stored as `workflow_cell_id` in audit/QC sidecars and mapped during clean
export. Full h5ad `obs_names` are retained in the ID map as `h5ad_obs_name`.

Do not expose intermediate coordinate columns in the clean CSV. Store
displacements, edit labels, component ranks, and recipe ids in audit/QC
sidecars.

## State Pointers

Use `states/current_review.yaml`, `states/current_accepted.yaml`, and
`states/final.yaml` to identify the latest coordinate file. Dashboards and
handoffs must link to the clean CSV recorded in the state pointer, not to the
newest file by timestamp.

## Review Rules

- Use human judgment for suspicious edges and components.
- Start manual review from full-data coordinates. Use balanced300k only after
  an edit, for fast structure checks, or for repeated candidate iteration.
- When naming a component in a decision, include slice, rank, and component id,
  for example `SL59 rank2 id7`; component ids are not global across slices.
- Prefer conservative component edits: translation, rigid transform, or
  small-anchor transform applied to a specified component group.
- Do not propagate component transforms across long ranges unless the user
  explicitly accepts that operation.
- Do not regenerate a large full CSV for every minor candidate; replay to full
  only at review milestones.
- Do not call bulky audit files "final" or "clean".
