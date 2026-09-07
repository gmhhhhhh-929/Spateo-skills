---
name: spateo-pairwise-run
description: Plan and run adjacent pairwise Spateo alignment jobs on a remote cluster. Use when Codex needs to create dsub payloads for normal default Spateo runs, spatial-only rescue runs, sigma sweeps, or small SL windows while keeping raw h5ad remote and syncing only pairwise outputs, audit files, provenance, logs, and coordinate CSVs.
---

# Spateo Pairwise Run

## Package paths

Run the relative commands below from this skill directory. Install this skill
beside `spateo-2d-alignment` (or keep the complete repository checkout).
The shared support tools are the fixed pre-zebrafish snapshot, located in
`../spateo-2d-alignment/pipelines/pairwise-rigid`. Verify entrypoints with its `skill.lock.yaml`.
This packaging restoration does not change alignment algorithms or defaults.


## Overview

Use this skill after `remote-workflow-intake` confirms the remote host, raw
h5ad root, workdir, conda environment, and scheduler account. Default to a
normal Spateo pairwise run; use `spatial_only` only for rescue or diagnostic
runs after QC or explicit user request.

All job planning, payload generation, Spateo execution, transform extraction,
and summaries should be produced on the remote server. Local copies are for
inspection and deterministic replay only. Local replay must use the same
remote-produced recipe and transform code, preserve the original local
coordinate cache, and pass remote/local consistency checks before being treated
as matching the remote result.

## Defaults

- Keep raw h5ad on the remote machine.
- Keep coordinate-producing outputs on the remote machine as the authoritative
  copy.
- Treat source h5ad `obs_names`, `obs["cell_id"]`, `obs["CellID"]`, and slice
  annotations as the cell identity source of truth. Pairwise runners may emit
  internal ids for audit, but these must not become final clean `cell_id`
  values.
- Submit heavy jobs through `dsub`; do not run Spateo payloads on the login
  node.
- Default run mode:
  - `EXPRESSION_MODE=normal`
  - `STAGE1_MODE=SN-S`
  - `STAGE2_MODE=none`
  - runner default sigma parameters
- Spatial-only rescue mode:
  - `EXPRESSION_MODE=spatial_only`
  - `DUMMY_REP_DIM=30`
  - dummy `X_pca` makes expression distance constant so matching is driven by
    spatial probability.
- Treat every output as a candidate until QC and visual review confirm it.

## Plan Jobs

Use the shared planner. For a normal run:

```bash
python ../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/core/make_spateo_pairwise_initial_plan.py \
  --dataset-root RAW_H5AD_ROOT \
  --runner-script SPATEO_RUNNER \
  --stage CS13 \
  --project-account PROJECT_ACCOUNT \
  --conda-env CONDA_ENV \
  --stage1-mode SN-S \
  --stage2-mode none \
  --expression-mode normal \
  --dsub-resource 'cpu=64;mem=240000;gpu=1' \
  --output-dir PAIRWISE_RUN_DIR
```

For a spatial-only rescue or sigma sweep, change only the relevant knobs:

```bash
python ../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/core/make_spateo_pairwise_initial_plan.py \
  --dataset-root RAW_H5AD_ROOT \
  --runner-script SPATEO_RUNNER \
  --stage CS13 \
  --window-start SL171 --window-end SL175 \
  --project-account PROJECT_ACCOUNT \
  --conda-env CONDA_ENV \
  --stage1-mode SN-S \
  --stage2-mode none \
  --expression-mode spatial_only \
  --dummy-rep-dim 30 \
  --sigma2-init-scale 2.0 \
  --output-dir PAIRWISE_SPATIAL_ONLY_DIR
```

Inspect `pairwise_jobs.csv`, payloads, and `submit_pairwise_spateo.sh`, then
submit on the remote machine:

```bash
bash PAIRWISE_RUN_DIR/submit_pairwise_spateo.sh
```

## Extract Outputs

After jobs finish, extract transforms and summarize:

```bash
python ../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/core/extract_pairwise_transform.py \
  --run-dir PAIRWISE_RUN_DIR \
  --stage stage1_SN-S_rigid \
  --output-edge-index PAIRWISE_RUN_DIR/pairwise_edge_index.csv

python ../spateo-2d-alignment/pipelines/pairwise-rigid/scripts/core/summarize_pairwise_initial_run.py \
  --run-dir PAIRWISE_RUN_DIR \
  --high-sigma2-threshold 0.1
```

Sync only small outputs for local review unless the user explicitly asks for a
final remote-produced cache:

- `pairwise_edge_index.csv`
- `pairwise_edge_summary.json`
- `edges/*/edge_transform.json`
- `edges/*/spateo_two_stage_aligned_points.csv`
- `edges/*/spateo_two_stage_provenance.json`
- `edges/*/spateo_pair_audit/**/spateo_pair_audit_summary.json`
- logs

Do not apply pairwise outputs directly from this skill. If a pairwise candidate
is confirmed, hand it to `spatial-alignment-compose`, which may replay the same
recipe both remotely and locally and verify consistency.

If `spateo_two_stage_aligned_points.csv` contains generated ids such as
`slice_id:0`, record them only as internal `workflow_cell_id` values. Before a
user-facing coordinate CSV or final h5ad is produced, map those rows back to
the original h5ad-derived ids with the project h5ad ID map.

## Required Checks

- Payload logs must show selected slices and run mode, especially
  `expression_mode`, `dummy_rep_dim`, `sigma2_init_scale`, and `sigma2_end`.
- Every finished edge should contain `spateo_two_stage_aligned_points.csv`,
  `spateo_two_stage_provenance.json`, and a pair audit summary.
- Pairwise output ids must either already match original h5ad-derived ids or be
  explicitly documented as internal ids that require h5ad ID map conversion.
- Hand normal and spatial-only result directories to `spateo-pairwise-qc`.
