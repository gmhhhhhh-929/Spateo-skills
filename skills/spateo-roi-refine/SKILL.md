---
name: spateo-roi-refine
description: Rerun one suspicious Spateo pair using a manually chosen ROI or drop mask, usually spatial-only SN-S rigid, to remove an interfering local structure before alignment. Use when visual QC shows a small ring, detached region, partial tissue, or other local shape that disrupts pairwise alignment and the user asks to align using only the main body or selected region.
---

# Spateo ROI Refine

## Package paths

Run the relative commands below from this skill directory. Install this skill
beside `spateo-2d-alignment` (or keep the complete repository checkout).
The shared support tools are the fixed pre-zebrafish snapshot, located in
`../spateo-2d-alignment/pipelines/pairwise-rigid`. Verify entrypoints with its `skill.lock.yaml`.
This packaging restoration does not change alignment algorithms or defaults.


## Overview

Use this skill for one edge at a time after QC or visual inspection identifies
an interfering region. The default refinement is `spatial_only + SN-S rigid`
with an explicit ROI/drop condition, for example `SL190 drop y > 3000`.

Run ROI/drop estimation and any full-coordinate candidate generation on the
remote server. Local coordinate replay, when needed, belongs to
`spatial-alignment-compose` and must use the same remote-produced recipe and
transform code with checksum or numeric consistency checks.

## Guardrails

- Keep raw h5ad on the remote machine.
- Keep ROI/drop candidate outputs and recipes on the remote machine as the
  authoritative copy.
- Treat ROI/drop output as a candidate until the user visually confirms it.
- Record the exact ROI/drop condition in provenance and in any recipe sent to
  `spatial-alignment-compose`.
- Do not apply the result to full coordinates from this skill.
- Do not create final corrected coordinate CSVs directly from this skill.
  Full-coordinate application belongs to `spatial-alignment-compose`, which may
  replay both remotely and locally from preserved original coordinate caches.
- Do not use ROI/drop to hide broad failures; if the main tissue is unstable,
  return to pairwise run/QC.

## Runner

Use a runner that supports both spatial-only mode and ROI/drop filtering. The
current local reference implementation is:

```text
../spateo-2d-alignment/pipelines/pairwise-rigid/runners/spateo_rigid_sns_alignment_spatial_only_roi.py
```

Copy or sync that runner to the remote workdir before submitting. Keep a copy
of the exact runner path in provenance.

## Required Parameters

For each ROI run, record:

- edge id, fixed slice, moving slice, and pair start order;
- `EXPRESSION_MODE=spatial_only`;
- `DUMMY_REP_DIM=30`;
- `STAGE1_MODE=SN-S`;
- `STAGE2_MODE=none`;
- sigma settings, if changed;
- ROI/drop condition, including slice, coordinate column, operator, threshold,
  and whether the condition is applied to fixed, moving, or both slices.

## Output Requirements

Each ROI candidate should produce:

```text
spateo_two_stage_aligned_points.csv
spateo_two_stage_provenance.json
spateo_pair_audit/stage1_SN-S_rigid/spateo_pair_audit_summary.json
roi_refine_recipe.json
roi_refine_summary.json
logs/
```

The recipe must be suitable for `spatial-alignment-compose` and must include:

- `operation_type: roi_drop_refine`
- source pairwise run path
- affected downstream slices
- transform steps
- expression mode and sigma settings
- exact ROI/drop condition
- user confirmation note when available

## Review

Sync recipes, provenance, QC tables, logs, and small pair CSVs for local review.
Use `spatial-alignment-compose` for full-coordinate replay and consistency
checks. Check:

- both slices are full, not sampled, unless the user explicitly requested a
  preview;
- slice labels are the intended pair;
- the interfering structure was excluded only from transform estimation, not
  silently removed from final full-coordinate review;
- metrics improve without obvious block-level distortion.

If the user confirms the ROI candidate, pass the recipe and affected slices to
`spatial-alignment-compose`.
