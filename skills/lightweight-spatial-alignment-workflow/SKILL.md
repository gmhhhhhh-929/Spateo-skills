---
name: lightweight-spatial-alignment-workflow
description: Run a reproducible lightweight spatial transcriptomics alignment workflow with full slice-level rigid/spatial-only alignment, connected-component detection, full-points manual review, balanced300k post-edit review, component candidate recipes, full replay validation, locked entrypoints, and workflow records. Use when Codex needs to set up, audit, or run this standardized alignment workflow for CS-stage embryo or other serial-slice spatial data.
---

# Lightweight Spatial Alignment Workflow

## Package paths

Run the relative commands below from this skill directory. Install this skill
beside `spateo-2d-alignment` (or keep the complete repository checkout).
The shared support tools are the fixed pre-zebrafish snapshot, located in
`../spateo-2d-alignment/pipelines/pairwise-rigid`. Verify entrypoints with its `skill.lock.yaml`.
This packaging restoration does not change alignment algorithms or defaults.


## Overview

Use this skill to run or organize the standardized workflow:

1. full adjacent pairwise rigid alignment;
2. pair QC and human-selected spatial-only rescue;
3. full slice-level baseline compose;
4. connected-component detection on that baseline;
5. full-points manual review for issue finding;
6. balanced300k post-edit structure review and candidate tuning;
7. full replay and validation at review milestones.

The first manual review artifact is a full-points viewer. The component-balanced
300k viewer is used after edits or during repeated candidate tuning. Full
CSV/h5ad outputs remain complete and are validated separately.

## Release Contract

Use `../spateo-2d-alignment/pipelines/pairwise-rigid` from this repository or an immutable deployed
copy of it. Resolve all scripts through `../spateo-2d-alignment/pipelines/pairwise-rigid/skill.lock.yaml`; do not call copied
entrypoints by ad hoc absolute paths unless the lock resolver validated them.

The lock supports relative `release_root: "."` inside the repository. BGI
deployments may generate a site-specific lock or record the deployed absolute
path in the run record, but entrypoint sha256 must match.

## Workflow

1. Initialize a run skeleton with `../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/workflow/init_project_run.py` and a
   project config based on `references/example.project.yaml`.
2. Verify source h5ad identity fields and, when needed, build an h5ad ID map.
   The original h5ad `obs_names`, `obs["cell_id"]`, `obs["CellID"]`, and slice
   annotations are the authority for cell identity; do not invent final ids.
3. Run pairwise rigid and spatial-only rescue through the locked pairwise
   entrypoints.
4. Compose accepted slice-level recipes on the full coordinate table.
5. Detect components once on the full baseline, preserving component id/rank
   labels for balanced300k review; keep per-slice PNGs only as detector QC
   sidecars.
6. Generate a full-points review HTML from the component-labeled full baseline
   with `../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/viewers/make_st_pointcloud_viewer.py`; use
   `--viewer-kind review_full_points`, `--component-col component_id`, and
   `--component-rank-col component_rank`.
7. Create balanced300k with `../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/sampling/create_balanced_spatial_sample.py`
   using `--target-total-points 300000`.
8. Generate post-edit or candidate review HTML from the balanced300k CSV with
   `../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/viewers/make_st_pointcloud_viewer.py`; pass source full CSV path,
   full row count, full sha256, sampling method, `--component-balanced`,
   `--component-col auto`, and `--component-rank-col auto`. The HTML viewer is
   the lightweight post-edit component map and should expose `Component` and
   `Labels`.
9. Tune component or triad candidates on balanced300k and store recipes/deltas.
10. Replay accepted or review-ready recipes to the full coordinate table, write
   `audit_points.csv`, export `clean_coordinates.csv`, and run validation.
11. If the source coordinate table uses workflow row-index ids such as
   `slice_id:0`, first build an h5ad ID map with
   `../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/core/build_h5ad_cell_id_map.py`, then export clean coordinates with
   `--cell-id-map`; user-facing `cell_id` must use the short h5ad-derived
   form `<slice>_<CellID>`, for example `SL14_CELL.1`.
12. Update `states/current_review.yaml`, `states/current_accepted.yaml`, or
   `states/final.yaml` with the clean coordinate path.
13. Record every substantive operation with `../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/workflow/write_workflow_record.py`
   and rebuild the compact dashboard with `../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/workflow/build_workflow_dashboard.py`.

## Guardrails

- Do not call a 300k viewer `full_preview`; use `review_balanced300k` or
  `candidate_balanced300k`.
- Do not use balanced300k as the first manual review surface. Start issue
  finding from `review_full_points.viewer.html`, then use balanced300k after a
  candidate edit or for repeated tuning.
- Do not create new final cell ids. Clean CSV and final h5ad identity must be
  derived from the original source h5ad annotations. Generated ids such as
  `slice_id:0` belong only in `workflow_cell_id` audit/QC sidecars.
- Do not use random viewer downsampling as the default review path. Create a
  balanced sample first, then render that exact CSV.
- Do not let component ids become ambiguous. Human-facing notes must name
  components as `<slice> rank<rank> id<component_id>`, for example
  `SL59 rank2 id7`, and should link to the review HTML that shows the
  component labels in-view.
- Do not generate review viewers without component labels when component review
  is expected. The viewer summary must record `component_col`,
  `component_rank_col`, non-empty `components`, `component_label_overlay`, and
  an appropriate `display_policy` such as `manual_full_points` or
  `balanced300k`.
- Do not redetect components for every candidate unless the full baseline
  changed materially.
- Do not mark a sampled candidate final. Final status requires full replay,
  row-count validation, hash recording, and celltype preservation checks.
- Generate full-points HTML for manual issue finding and final manual
  inspection. Do not name it `full_preview`; use `review_full_points` or
  `candidate_full_points`.
- Do not share bulky audit CSVs as final or latest coordinates. User-facing
  coordinates must be `clean_coordinates.csv` with columns
  `cell_id,slice_id,stage,chip_id,sl_number,celltype,x,y,z`.
- In clean CSVs, `cell_id` means the short h5ad-derived id
  `<slice>_<CellID>`. Internal ids such as `slice_id:zero_based_row_index`
  must be named `workflow_cell_id` and kept in audit/QC sidecars only.
- Use `states/*.yaml` to identify the latest coordinate file; do not infer it
  from modification time or a filename guess.

## References

Read `references/workflow.md` for the step order and
`references/project_layout.md` for the directory contract. Read
`references/reproducibility.md` for lock/provenance requirements.
