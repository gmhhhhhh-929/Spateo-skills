---
name: spateo-pairwise-qc
description: Quality-control and compare pairwise Spateo alignment outputs, including normal runs, spatial-only rescue runs, sigma sweeps, and ROI/drop reruns. Use when Codex needs to rank suspicious adjacent slice edges, compare candidate runs by geometry metrics, inspect cell-type discontinuities, or decide whether to hand an edge to ROI refine or compose without editing coordinates.
---

# Spateo Pairwise QC

## Package paths

Run the relative commands below from this skill directory. Install this skill
beside `spateo-2d-alignment` (or keep the complete repository checkout).
The shared support tools are the fixed pre-zebrafish snapshot, located in
`../spateo-2d-alignment/pipelines/pairwise-rigid`. Verify entrypoints with its `skill.lock.yaml`.
This packaging restoration does not change alignment algorithms or defaults.


## Overview

Use this skill after `spateo-pairwise-run` or `spateo-roi-refine` has produced
edge outputs. This is triage only: it ranks and explains candidate results but
does not submit jobs, apply recipes, or edit coordinates.

Remote/BGI remains the authoritative data-production environment. Run QC metric
generation and summaries on the remote server against remote outputs. Local
copies are for reading finished QC tables, generating visualization HTML, and
deterministic replay through `spatial-alignment-compose`; local QC must not
generate metrics from source data, create recipes, directly apply transforms, or
create adjusted h5ad files.

## Inputs

- One or more pairwise run directories containing `pairwise_edge_index.csv` or
  `edges/*/edge_transform.json`.
- Full-coordinate CSV covering the inspected slices.
- Coordinate columns. Runner defaults are `stage1_rigid_x`, `stage1_rigid_y`,
  `z_display`, `sl_number`, `celltype`, `cell_id`.

## Numeric QC

Run the shared QC script for each candidate run on the remote server:

```bash
python ../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/core/qc_pairwise_celltype_edges.py \
  --pairwise-run-dir PAIRWISE_RUN_DIR \
  --points-csv FULL_COORDINATE_CSV \
  --x-col stage1_rigid_x --y-col stage1_rigid_y --z-col z_display \
  --slice-col sl_number --celltype-col celltype --cell-id-col cell_id \
  --output-dir PAIRWISE_QC_DIR
```

Outputs:

```text
pairwise_edge_qc.csv
pairwise_celltype_edge_qc.csv
suspicious_edges.csv
suspicious_celltypes.csv
pairwise_qc_summary.json
```

Use metrics as prioritization signals, not final truth:

- final `sigma2`
- all-cell centroid distance
- bbox/shape overlap when available from downstream summaries
- PCA axis delta
- key cell-type centroid distances
- reason flags for centroid jump, axis jump, radius shift, or high sigma2

## Compare Candidates

For the same edge, compare normal, spatial-only, sigma sweep, and ROI/drop
results in one short table. Interpret them this way:

- normal good: keep normal as the candidate.
- normal bad but spatial-only good: expression likely pulled the match; send
  the spatial-only candidate to visual review and then compose if confirmed.
- spatial-only still bad: consider `spateo-roi-refine`, initialization, or a
  shape ambiguity.
- ROI/drop good: keep it as a candidate but record the drop/ROI condition in
  compose history.

## Visual QC

Generate overlays for top suspicious rows when numeric results are ambiguous:

```bash
python ../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/core/render_pairwise_visual_qc.py \
  --pairwise-run-dir PAIRWISE_RUN_DIR \
  --points-csv FULL_COORDINATE_CSV \
  --qc-csv PAIRWISE_QC_DIR/suspicious_celltypes.csv \
  --x-col stage1_rigid_x --y-col stage1_rigid_y --z-col z_display \
  --slice-col sl_number --celltype-col celltype --cell-id-col cell_id \
  --output-dir PAIRWISE_QC_DIR/visual_qc \
  --top-k 20
```

Prefer local self-contained HTML from `spatial-pointcloud-viewer` when the user
wants interactive inspection. Use a checked local replay output or copied
remote final cache when viewing corrected coordinates.

## Do Not Do

- Do not run QC scripts locally against source coordinate data.
- Do not apply transforms, compose recipes, rewrite coordinate CSVs, or create
  adjusted h5ad files directly from QC; hand confirmed candidates to
  `spatial-alignment-compose`.
- Do not treat local synchronized QC tables as authoritative if a newer remote
  run exists.

## Handoff

- If an edge needs a spatial-only rescue, hand it to `spateo-pairwise-run`.
- If an edge has a visible interfering region or small structure, hand it to
  `spateo-roi-refine`.
- If the user confirms a candidate, hand source recipes and affected slices to
  `spatial-alignment-compose`.
