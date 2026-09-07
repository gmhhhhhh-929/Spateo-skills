---
name: spatial-before-after-viewer
description: Generate self-contained before/after visual comparisons for spatial transcriptomics alignment candidates. Use when Codex needs to compare original versus corrected coordinates, create HTML/summary review artifacts, inspect pairwise slice changes, verify z spacing or z exaggeration, or present non-destructive alignment previews without editing h5ad or source coordinate tables.
---

# Spatial Before/After Viewer

## Package paths

Run the commands below from this skill directory. The viewer script is bundled
locally and can be installed independently. CSV rendering requires NumPy;
the before/after viewer additionally requires Pillow. The local
`provenance/viewer.json` records the original script hash.


## Scope

Create review artifacts that show a candidate against its source coordinates.
This skill is visualization only. It must not compose recipe chains, edit h5ad
files, overwrite coordinate caches, or mark a candidate accepted.

Use it for targeted candidate review. Use `spatial-pointcloud-viewer` for
whole-embryo full-points or balanced300k pointcloud review.

## Locked Entrypoint

Use the bundled script, authenticated by `provenance/viewer.json`:

```text
scripts/make_spatial_before_after_viewer.py
```

Record `before_after_viewer` entrypoint sha256 in the review summary or
`provenance/skill_run.json`. Stop on hash mismatch.

## Workflow

1. Confirm the input table has `cell_id`, slice, before x/y/z, after x/y/z,
   and an optional color/group column such as `celltype`,
   `component_group_slot`, or `alignment_group`.
2. Keep physical z coordinates and use display z exaggeration only for review
   if the slices look visually collapsed.
3. Generate a self-contained HTML viewer plus summary JSON.
4. Link the HTML from the owning step's `review/` directory. Treat it as a
   candidate artifact, not final evidence.

## Command Shape

```bash
python scripts/make_spatial_before_after_viewer.py \
  --input candidate_points.csv \
  --cell-id-col cell_id \
  --slice-col sl_number \
  --before-x x --before-y y --before-z z \
  --after-x candidate_x --after-y candidate_y --after-z candidate_z \
  --color-col component_group_slot \
  --physical-z-step 40 \
  --z-display-multiplier 8 \
  --output-html review/before_after.html \
  --output-summary review/before_after_summary.json
```

Use actual coordinate column names from the step. Do not assume `manual_x`
exists in a clean coordinate CSV.

## Output Checks

- `total_points` matches the intended input unless explicit sampling was used.
- The summary records physical z step and display multiplier.
- Before, after, and delta modes are available.
- The review makes worsened components visible; do not crop away failures.
- The artifact is named as before/after review, not final output.

## Guardrails

- Do not overwrite original h5ad, source coordinate cache, or clean CSV.
- Do not use a before/after viewer as proof of acceptance without human review.
- Do not hide unchanged or worsened neighboring slices when the candidate could
  affect local continuity.
