# Spatial Alignment Project Layout v0.2.3

BGI remote storage is the canonical project root. Local folders are mirrors for
review HTML, small manifests, and dashboards only.

## Canonical Tree

```text
spatial_alignment_projects/
  _pipeline/
    releases/v0.2.3/
  <sample_id>/
    project.yaml
    source/
      source_manifest.yaml
      source_links.txt
    runs/
      <run_id>/
        run_manifest.yaml
        states/
          slice_level_baseline.yaml
          current_review.yaml
          current_accepted.yaml
          final.yaml
        lineage/
          all_operations.jsonl
          accepted_steps.jsonl
          state_history.jsonl
          current_review.recipe_chain.json
          accepted.recipe_chain.json
        steps/
          001_<actual_operation_slug>/
        final/
          <sample>.final.clean_coordinates.csv
          <sample>.final.h5ad
          <sample>.workflow.html
          final_manifest.yaml
        logs/
```

Step numbers express order only. The slug must describe the actual operation.
QC, compose, metrics, logs, and viewers are artifacts inside an operation step
unless the analysis itself is the operation.

## Step Contract

Each step uses the same internal layout:

```text
steps/NNN_<actual_operation_slug>/
  step.yaml
  record.json
  decision.md
  provenance/
    skill_run.json
    validator_report.json
    input_manifest.json
  recipes/
  qc/
  review/
  light/
  full/
  exports/
  logs/
```

Status values are `candidate`, `review`, `accepted`, `rejected`, `not_needed`,
`superseded`, and `final`. Coordinate-changing steps must include a recipe,
provenance validator, and review artifact. Full replay milestones must include
both audit and clean coordinate artifacts.

## Coordinate Artifacts

Full replay milestones write:

```text
full/<sample>.<step>.audit_points.csv
exports/<sample>.<step>.clean_coordinates.csv
exports/<sample>.<step>.clean_export_manifest.json
qc/<sample>.<step>.changed_cells.csv
qc/<sample>.<step>.full_replay_validation.json
```

`clean_coordinates.csv` is the only user-facing coordinate table. Its schema is
fixed:

```text
cell_id,slice_id,stage,chip_id,sl_number,celltype,x,y,z
```

`cell_id` must be the short h5ad-derived id `<slice>_<CellID>`, for example
`SL14_CELL.1`. It preserves the original h5ad `CellID` values and their gaps;
it must never be renumbered. If an upstream legacy coordinate cache uses an
internal row-index id such as `slice_id:0`, build and record a sidecar map:

```text
source/<sample>.h5ad_cell_id_map.csv
```

The map schema is:

```text
workflow_cell_id,cell_id,h5ad_obs_name,slice_id,stage,chip_id,sl_number,slice_short,obs_row_index,source_h5ad,obs_cell_id,obs_CellID
```

`workflow_cell_id` may appear in `audit_points.csv`, `changed_cells.csv`, and
QC/debug sidecars, but never in `clean_coordinates.csv`. The full h5ad
`obs_names` value is retained as `h5ad_obs_name` in the map.

Intermediate columns such as `raw_*`, `manual_*`, `full_candidate_*`,
`step*`, edit labels, displacement metrics, and component ranks belong in
`audit_points.csv` or sidecar QC files, never in the clean CSV.

## State Pointers

`states/*.yaml` answer "which coordinate file is current?" Dashboard and
sharing links must read these files instead of guessing the newest CSV:

- `slice_level_baseline.yaml`: full slice-level rigid/spatial-only baseline.
- `current_review.yaml`: latest candidate under human review.
- `current_accepted.yaml`: latest accepted non-final state.
- `final.yaml`: final released state.

State files point to clean coordinates, audit coordinates, recipe chain,
reviewer HTML, validation report, row count, hashes, and validator status.

## Review Policy

Primary manual review starts from full-data coordinates. The corresponding HTML
should be named `review_full_points.viewer.html` or
`candidate_full_points.viewer.html`, not `full_preview`. Use this viewer to
decide whether a baseline region or component is truly wrong.

Component-balanced 300k viewers are for post-edit structure checks and fast
candidate iteration. A 300k viewer must not be called `full_preview`; use
`review_balanced300k.viewer.html` or `candidate_balanced300k.viewer.html`.

Connected-component labels should be visible inside the HTML viewer. The viewer
should include Component color mode and a `Labels` toggle that places labels
such as `SL59 rank2 id7` near component centroids. Full-points review can be
used for the initial decision, and balanced300k can be used after edits.

Connected-component detection steps should still keep per-slice PNGs as QC
sidecars:

```text
qc/figures/slice_<slice>_components.png
```

These figures are for checking whether the detector split or merged tissue
blocks incorrectly. They are not the main review interface. Review notes should
refer to components as `<slice> rank<rank> id<component_id>`, for example
`SL59 rank2 id7`.
