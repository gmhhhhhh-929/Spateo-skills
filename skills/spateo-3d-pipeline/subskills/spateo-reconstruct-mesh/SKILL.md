---
name: spateo-reconstruct-mesh
description: Reconstruct a full-body surface and selectable annotation-specific meshes from a validated Spateo point-cloud VTK, compare native and robust marching-cubes routes, and deliver separately saved, previewed, round-trip-verified VTK meshes.
---

# Reconstruct selectable Spateo surface meshes

Use this second phase only after the point-cloud skill has produced a `pyvista.PolyData` VTK with unique `point_data["obs_index"]`. Annotation meshes also require the chosen categorical point-data array, such as `annotation`, `lineage`, or `tissue`.

Read [method selection](references/method-selection.md) before changing reconstruction parameters or choosing the native route. Use [validation](references/validation.md) when reporting evidence and limits.

For missed peripheral cells, sparse-region coverage, protrusion cleanup, tissue/body consistency or a smooth body envelope, read [coverage-aware refinement](references/coverage-and-envelopes.md). Use its reusable review runner and regional audits rather than tuning solely to global containment or a single screenshot. Own-cell coverage, tissue-in-body containment, and tissue-union occupancy are different measurements; independent tissue meshes are not automatically a space-filling anatomical partition.

## Choose the reconstruction route

Use the default `density` route for cells distributed through a tissue volume. It places the body and all requested annotations on one coordinate grid, smooths binary occupancy into a density field, chooses an isovalue against a requested source-point coverage, removes only low-support density components, extracts a marching-cubes surface, applies low-shrinkage Taubin smoothing, and repairs open edges without dropping disconnected components. It does not delete or rewrite source cells.

Use `--method spateo-marching-cube` only for a direct comparison or when the user explicitly requests Spateo's native geometry. The parameter is `mc_scale_factor`, exposed as `--mc-scale-factor`; `smooth` is a Laplacian post-processing iteration count in this mode. The wrapper uses Spateo's pinned `marching_cube_mesh` core, supplies a deterministic per-target random seed, and normalizes its output to triangle-only PolyData because the full `construct_surface` repair path can fail on real mixed connectivity.

Do not use Poisson, ball pivoting, or PyVista `reconstruct_surface` as the default for spatial-transcriptomic cell centroids: those methods assume samples lie on the object surface and the former two require meaningful oriented normals.

## Build full-body and selected meshes

Run from this skill directory. A selection can combine multiple categorical values with `|`:

```bash
python scripts/build_meshes.py \
  --input-vtk /results/model_pc/model_annotation_pc.vtk \
  --output-dir /results/model_mesh_v1 \
  --name model \
  --label-key annotation \
  --selection 'neural=CNS|CNS progenitor' \
  --selection 'gut=goblet cells|goblet cells progenitor' \
  --selection 'pharynx=pharynx|pharynx pouch'
```

The body uses every point. Each `--selection NAME=VALUE1|VALUE2` creates a mesh from only those labeled points. All outputs share voxel size and Gaussian bandwidth so overlays stay registered, but each selection has its own occupancy field, normalization, component filtering, and isovalue search. Use palette and opacity JSON objects keyed by output mesh name when stable display colors matter. Use `--no-body` to reconstruct selections alone.

The density defaults are starting values, not biological constants. Review the first result before changing `--target-long-axis` or `--voxel-size`, `--sigma`, coverage targets, component thresholds, and `--smooth`. Higher grid resolution preserves detail but costs cubic memory; larger sigma connects gaps but can bridge nearby structures. Never raise smoothing merely to hide an incorrect density field. Density-only options do not affect the native route, and native-only options do not affect density reconstruction; the manifest records the active and inactive method groups plus every value needed to reproduce either route.

For a native comparison:

```bash
python scripts/build_meshes.py \
  --input-vtk /results/model_pc/model_annotation_pc.vtk \
  --output-dir /results/model_native_mesh_v1 \
  --name model_native \
  --method spateo-marching-cube \
  --levelset 0.5 \
  --mc-scale-factor 0.8 \
  --smooth 300
```

The output directory must be new or empty. The builder stages the complete result and publishes it atomically, so a failed later tissue does not leave a partial body VTK that blocks a rerun. Every body or annotation model is saved separately as `<name>.<mesh-name>.vtk`, then reloaded with `st.tdr.read_model`. The builder also writes a four-view PNG and a JSON manifest containing exact selections, density/native parameters, component filtering, point coverage, repair topology before/after, open edges, hashes, face connectivity, label/RGBA values, and round-trip checks. It distinguishes the revision this skill was validated against from the actual runtime `mesh_methods.py` path and hash. Repair keeps valid disconnected components; only an unrepairable fragment no larger than `--repair-max-degenerate-cells` can be dropped, and every such drop changes the manifest to `pass_with_warnings`; an unexpected component split is also a warning.

## Review and handoff

For interactive inspection of these saved outputs, read the sibling [3D viewer](../spateo-render-3d-viewer/SKILL.md). Its `mesh_manifest` input can import this builder's manifest directly, with an optional traceable point cloud for grouped tissue cell overlays. Visualization does not rerun or refine reconstruction.

Inspect the isometric and orthographic views. Confirm that the body encloses the expected point cloud; internal meshes occupy plausible locations; paired or genuinely disconnected anatomy was not removed; nearby structures were not artificially bridged; z-layer gaps are not visible; and the result is not excessively shrunken or inflated.

For the default density route, prefer `status: pass`; explicitly review every message when the result is `pass_with_warnings`. Require successful Spateo round trips, triangle-only VTKs, zero open edges unless the user explicitly accepted `--skip-repair`, and final containment at or above `--minimum-final-inside-fraction` (default 0.80). Automatic isovalue search fails when its requested voxel-coverage target cannot be met. Native comparison can also finish as `pass_with_warnings` when containment is low; do not promote that result as the accepted mesh without review. Containment is still a geometric diagnostic, not proof of biological accuracy. Return the preview, manifest, and each independent VTK so the user can show or hide tissues without rebuilding the other meshes.
