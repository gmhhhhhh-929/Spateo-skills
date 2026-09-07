---
name: spatial-pointcloud-viewer
description: Generate a self-contained WebGL HTML viewer for local spatial transcriptomics 3D point clouds from full-points review tables, balanced300k post-edit samples, local coordinate caches, local replay outputs, or copied remote outputs. Use when the user wants local visualization or browser review after recipe-based coordinate updates, while preserving original coordinates and avoiding h5ad edits.
---

# Spatial Point Cloud Viewer

## Package paths

Run the commands below from this skill directory. The viewer script is bundled
locally and can be installed independently. CSV rendering requires NumPy;
the before/after viewer additionally requires Pillow. The local
`provenance/viewer.json` records the original script hash.


Generate a local, self-contained HTML viewer for 3D ST coordinates from a
full-points review table, balanced300k post-edit sample, local coordinate cache,
or local replay output generated from remote-produced recipes. The usual
workflow is to sync only
recipes, provenance, QC tables, and logs, replay the recipe locally with the
same transform code used remotely, run a remote/local consistency check, and
then build the viewer from the checked local output.

The viewer itself is for visual inspection only: it does not compose correction
chains, edit h5ad files, or load imaging atlas meshes or organ labels.

Original local coordinate caches must be treated as read-only. If a recipe is
applied locally, write the corrected coordinates to a new replay output
directory and verify it against the remote reference result before presenting it
as consistent.

Use the script in `scripts/` with explicit input paths and column names. In the
standard lightweight workflow, generate a full-points `review_full_points`
viewer first for manual issue finding. After edits, create a
component-balanced 300k CSV and render that CSV without `--max-points`.
Use `--max-points` only for an explicit random debug preview; it is
random/celltype-balanced, not component-balanced. Pass
`--slice-col` when the coordinate cache has slice labels, or `--slice-from-z`
when each slice is represented by a discrete z plane; the generated viewer then
supports cell-type and slice color modes, whole-embryo slice filtering, and a
Pair view for comparing two selected slices or stepping through adjacent slice
edges. If the CSV contains `sample_component_id/sample_component_rank` or
`component_id/component_rank`, the viewer also supports a Component color mode
with labels such as `SL59 rank2 id7`. Component mode should show both the
sidebar legend and in-view labels anchored near component centroids, with a
`Labels` toggle to hide them when they obstruct dense regions.

The viewer uses a high-contrast discrete color palette for cell types. If many
cell types are present, verify that neighboring biological categories remain
visually distinguishable in the generated summary/viewer.

Runtime requirements:

- CSV input: Python with `numpy`.
- h5ad input is supported only for explicit visualization. Do not edit local
  h5ad files.
- Remote source h5ad/raw data should stay remote; local review should use a
  local coordinate cache plus remote-produced recipes or a checked local replay
  output.

## Full Manual Review Command

Use this when deciding whether a region or component is actually wrong.

```bash
python scripts/make_st_pointcloud_viewer.py \
  --input component_labels.csv \
  --format csv \
  --x-col x --y-col y --z-col z \
  --celltype-col celltype \
  --slice-col sl_number \
  --component-col component_id \
  --component-rank-col component_rank \
  --viewer-kind review_full_points \
  --source-full-csv exports/<sample>.slice_level_baseline.clean_coordinates.csv \
  --source-full-sha256 SHA256_OF_FULL_CSV \
  --source-full-row-count FULL_ROW_COUNT \
  --output-html review_full_points.viewer.html \
  --output-summary review_full_points.viewer_summary.json
```

## Balanced300k Post-Edit Command

Use already-corrected columns if the input is a copied remote final cache or a
local replay output that passed the remote/local consistency check.

```bash
python scripts/make_st_pointcloud_viewer.py \
  --input review_balanced300k.points.csv \
  --format csv \
  --x-col x --y-col y --z-col z \
  --celltype-col celltype \
  --slice-col sl_number \
  --component-col auto \
  --component-rank-col auto \
  --viewer-kind review_balanced300k \
  --component-balanced \
  --sampling-method "slice x component_id spatial_grid_balanced" \
  --sampling-seed 13 \
  --source-full-csv exports/<sample>.<step>.clean_coordinates.csv \
  --source-full-sha256 SHA256_OF_FULL_CSV \
  --source-full-row-count FULL_ROW_COUNT \
  --output-html review_balanced300k.viewer.html \
  --output-summary review_balanced300k.viewer_summary.json
```

If the cache has no explicit slice column but each physical slice has a unique
z value, replace `--slice-col sl_number` with `--slice-from-z`.

## Optional Local h5ad Command

Use this only for visualization. Do not use this path to create or mutate
adjusted h5ad files locally.

```bash
python scripts/make_st_pointcloud_viewer.py \
  --input adjusted.h5ad \
  --format h5ad \
  --obsm-key spatial_3d \
  --celltype-col celltype \
  --slice-col sl_number \
  --output-html viewer.html \
  --output-summary viewer_summary.json
```

## Optional Debug Sampling

```bash
python scripts/make_st_pointcloud_viewer.py \
  --input corrected_points.csv --format csv \
  --x-col refined_x --y-col refined_y --z-col refined_z \
  --celltype-col celltype --slice-col sl_number \
  --max-points 300000 --seed 13 \
  --viewer-kind debug_random_preview \
  --output-html debug_random_preview.html \
  --output-summary debug_random_preview_summary.json
```

## Output Checks

After generation, verify:

- the summary `displayed_points` matches the intended full-points, balanced300k, or debug count;
- `display_policy` is `manual_full_points` for full manual review and
  `balanced300k` for post-edit/candidate review;
- `full_points` and `source_full_sha256` are populated when rendering a candidate review;
- `component_balanced` is true for standard review viewers;
- `sampled` is false when rendering an already-created balanced300k CSV without `--max-points`;
- `slices` is populated when a slice column was expected;
- `components` is populated when reviewing connected components;
- `component_col` resolves to `sample_component_id` or `component_id`, and
  `component_rank_col` resolves to `sample_component_rank` or
  `component_rank`;
- the Component mode names components as `<slice> rank<rank> id<component_id>`
  so they can be matched to the per-slice component PNGs from detection;
- component labels are visible inside the HTML viewer itself, not only in a
  separate PNG or external table;
- Pair view is available when at least two slices are present, and `Prev`/`Next`
  step through adjacent slice pairs without changing the underlying data;
- the input was either a full-points review table, balanced300k post-edit
  sample, copied remote final cache, checked local replay output, or explicit
  uncorrected local cache preview;
- all-points HTML is named `review_full_points` or `candidate_full_points`
  when it is for manual review, and `debug_all_points` only when it is a debug
  artifact;
- no source data CSV/h5ad or generated HTML is committed to the skills repo.
- standard balanced300k review viewers are never named `full_preview`.

## Do Not Do

- Do not overwrite original coordinate caches.
- Do not create or edit adjusted h5ad files locally from the viewer skill.
- Do not generate or alter transform recipes here; use remote Spateo/ROI run
  skills for recipe creation and `spatial-alignment-compose` for replay.
- Do not claim a local replay output matches the remote result unless
  `local_remote_consistency.json` or an equivalent checksum/numeric comparison
  passed.
