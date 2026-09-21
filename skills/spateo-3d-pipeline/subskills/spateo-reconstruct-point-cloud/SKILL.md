---
name: spateo-reconstruct-point-cloud
description: Validate 3D coordinates in H5AD, construct a Spateo PyVista point cloud with uniform, obs, gene, multi-gene, or external values, and deliver a previewed, manifested, round-trip-verified VTK model.
---

# Reconstruct a Spateo point cloud

Use this implemented first phase of the parent 3D pipeline when an H5AD already contains XYZ coordinates. It preserves every observation in the VTK; any preview sampling is display-only.

Read [the API contract](references/api-contract.md) before changing coloring behavior or using another Spateo revision. Use [the validation record](references/validation.md) when reporting what the result does and does not establish.

## Inspect and choose coloring

Confirm `n_obs`, unique `obs_names`, the coordinate key and unit/frame, coordinate shape/rank, requested label or genes, expression layer, and missing values. Do not repair or rescale coordinates silently.

The pinned `construct_pc` implementation copies the AnnData in memory. Review available memory before running on a very large H5AD; never downsample the scientific VTK merely to make its preview cheaper.

Choose exactly one value source:

- no source option: one uniform category;
- `--obs KEY`: categorical annotation or continuous observation value;
- one or more `--gene GENE`: one expression vector or a row-wise multi-gene sum from `X`/`--layer`;
- `--labels-file FILE`: user values joined by exact observation ID.

For stable categorical colors, supply a JSON dictionary mapping every visible label to a Matplotlib-compatible color. Hidden labels omitted from a dictionary receive a neutral internal color before their alpha is set to zero; color lists still need an entry for every sorted category. Use a scalar alpha or a dictionary covering every visible label. `--mask` is categorical only; the wrapper non-destructively sets those points' RGBA alpha to zero while preserving their original label in the VTK. This also avoids the pinned Spateo implementation's scalar-alpha mask bug. Continuous values retain the scalar and colormap; their display opacity is set later at plot time.

## Build into a new result directory

Run from this skill directory:

```bash
python scripts/build_point_cloud.py \
  --input-h5ad /data/model.h5ad \
  --output-dir /results/model_pc_v1 \
  --name model_pc \
  --coordinate-unit um --coordinate-frame aligned_xyz \
  --obs celltype \
  --key-added celltype \
  --palette-json /project/celltype_colors.json
```

For continuous gene expression:

```bash
python scripts/build_point_cloud.py \
  --input-h5ad /data/model.h5ad \
  --output-dir /results/model_gene_sum_pc_v1 \
  --gene GENE_A --gene GENE_B --layer counts \
  --key-added gene_sum --colormap viridis
```

For user-provided values, the CSV/TSV defaults to columns `obs_index` and `value`; override those names explicitly when needed. IDs must cover the H5AD exactly, irrespective of row order.

Use `--coordinate-unit` and `--coordinate-frame` to record known metadata in the manifest. Omit them rather than guessing; a missing declaration remains explicit `null` provenance.

The output directory must be new or empty. The builder writes:

- `<name>.vtk`: binary `pyvista.PolyData` saved by `st.tdr.save_model`;
- `<name>.preview.png`: XY, XZ, YZ, and isometric review views, with a visible-category legend when there are at most 20 categories;
- `<name>.manifest.json`: input/output hashes, spatial/color contracts, environment, arrays, preview sampling, and `st.tdr.read_model` round-trip results.

Downstream code reads the delivered model through Spateo, not a private wrapper:

```python
import spateo as st

pc = st.tdr.read_model("/results/model_pc_v1/model_pc.vtk")
```

Use `--skip-preview` only for a documented headless or test run. Use `--allow-planar` only when an intentionally planar dataset happens to use three coordinate columns. Never use either flag to make a failed 3D input appear valid.

## Review and handoff

Check the preview orientation, extent, slice spacing, missing or unexpectedly transparent categories, continuous value range, and whether the coordinate rank/unit matches the intended model. Confirm the manifest says `status: pass`, all observations became points, required arrays survived reload, and the VTK hash matches.

Return the three artifacts and the chosen spatial/color fields. State that a point cloud is not yet a surface, voxel, cell model, backbone, or interpolated volume. Wait for user review before moving to the next implemented 3D phase.
