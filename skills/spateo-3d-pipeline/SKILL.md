---
name: spateo-3d-pipeline
description: Build and review Spateo point clouds plus full-body or annotation-specific surface meshes from AnnData with finite XYZ coordinates; route voxel, cell, backbone, and spatial interpolation requests only when their later subskills are available.
---

# Spateo 3D pipeline

Stage 5: environment → IO → slice quality → 2D alignment → **3D reconstruction** → 4D analysis. A point cloud is the lossless model foundation for later 3D models: one input observation becomes one PyVista point and remains traceable through `obs_index`.

Use Spateo from `gmhhhhhh-929/spateo-release` commit `615644f88613bea8ceb2e2df1e2391d16de55ec1`. The library is installed separately; this skill does not vendor it.

## Route by model

| Model or operation | Status | Route |
| --- | --- | --- |
| Point cloud (`pc`) | Implemented | Read [spateo-reconstruct-point-cloud](subskills/spateo-reconstruct-point-cloud/SKILL.md) and use its validated builder. |
| Full-body or annotation surface mesh | Implemented | First create a traceable point-cloud VTK, then read [spateo-reconstruct-mesh](subskills/spateo-reconstruct-mesh/SKILL.md). |
| Voxel or reconstructed cells | Pending | Do not invent a runner or claim completion; develop and validate the next subskill with the user. |
| Backbone construction and mapping | Pending | Preserve as a later reviewed phase. |
| Spatial gene interpolation | Pending | Preserve as a later reviewed phase; do not confuse it with 4D morphogenesis GP. |

## Point-cloud gate

Require a non-empty H5AD with unique `obs_names` and numeric, finite `(n_obs, 3)` coordinates under the selected `obsm` key, normally `spatial`. Do not infer z, coordinate units, alignment, annotations, or gene semantics. Full 3D coordinate rank is the default; accept planar XYZ only when the user explicitly intends it.

Color points uniformly, from an `obs` field, from one gene, from the sum of named genes, or from an exact-ID external label/value table. Stable categorical biology colors should use an explicit full palette. Keep continuous expression as a numeric scalar plus its colormap provenance.

Every implemented single-dataset model must be saved with `st.tdr.save_model(..., "*.vtk")`, reloaded with `st.tdr.read_model`, and checked before handoff. Never overwrite an existing result. Return the VTK, manifest, and four-view preview for review; preview sampling must not alter the full VTK.

## Surface-mesh gate

Start from the validated point-cloud VTK rather than rereading coordinates through an unrelated path. Require `obs_index`; annotation selections additionally require the requested categorical point-data array. Save the body and each selected annotation as an independent `.vtk` so they can be displayed or hidden separately.

Prefer the mesh skill's density-field route for volume-filling cell centroids. It filters low-support derived components without deleting source cells, uses one coordinate grid for aligned overlays, and records coverage/topology evidence. Use the pinned Spateo marching-cubes core only as an explicit comparison or user choice; its correct parameter name is `mc_scale_factor`, and high Laplacian `smooth` values can shrink anatomy. Stop after mesh review unless the user authorizes another model phase.

Install this complete directory so the nested workflow, script, references, and tests remain together.
