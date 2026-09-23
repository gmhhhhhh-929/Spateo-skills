# Viewer input contract

Paths resolve relative to the JSON file containing them, not the shell working directory. Absolute paths also work. A dataset's meshes and points must already share one coordinate frame. Only finite nonempty PolyData `.vtk` / `.vtp` inputs are accepted. Surface models must have polygon faces; separate point clouds must have unique `point_data['obs_index']` values.

```json
{
  "title": "3D Reconstruction Review",
  "datasets": [{
    "name": "Planarian 10 dpa",
    "units": "unspecified",
    "point_cloud": "point_cloud/annotation_pc.vtk",
    "label_key": "lineage",
    "point_budget": 20000,
    "layers": [
      {"name": "Body", "role": "body", "mesh": "meshes/body.vtk", "color": "#BFC7CD"},
      {"name": "Neural", "mesh": "meshes/neural.vtk", "values": ["CNS", "CNS progenitor"], "color": "#D98CB3"},
      {"name": "Pharynx", "mesh": "meshes/pharynx.vtk", "cells": "cells/pharynx_pc.vtk", "color": "#E9C46A"}
    ]
  }]
}
```

`datasets` and each dataset's `layers` must be nonempty, with unique names within their respective scope. Each layer requires at least `mesh` or `cells` or `values` (a point-only layer is valid). `color` is a six-digit hex color. Optional `visible` and `opacity` control the initial mesh; body defaults to opacity 0.15, other meshes to 1.0. Dataset `point_budget` is a positive integer, default 20000 **per layer**, not a global bound. The HTML sidecar records actual totals. A layer's `cells` overrides shared cloud selection and cannot be combined with `values`.

For shared `point_cloud`, a body layer selects all points; a tissue requires exact `values` under `label_key`. Missing requested categories fail explicitly instead of silently rendering incomplete cells. A mesh-only tissue without `values` has no associated cells or density metric. Original selected cells, including reconstruction exclusions, are eligible for preview; cells are not clipped to the mesh. Sampling cannot be interpreted as the complete source cloud.

Alternatively replace `layers` with `"mesh_manifest": "meshes/model.manifest.json"` for the native mesh builder's list-based manifest. Its `name`, `role`, `vtk`, `selected_values` and `color` entries become layers; relative VTK paths resolve against that manifest. Still set the dataset's `point_cloud` and `label_key` explicitly if cell overlays are wanted. Refined review reports with dictionary-based `meshes` are not this schema: select their final VTK paths explicitly, not intermediate candidates. Never supply both `layers` and `mesh_manifest`.

The output sidecar includes source/config/manifest hashes, per-layer counts, mesh topology, full-source morphology, source/display selected-cell ID hashes, sampling settings, package versions, warnings, and HTML hash/size. These technical fields are not embedded as a user-facing diagnostics panel. Geometry and selected preview IDs are embedded; no H5AD expression matrix is loaded. Viewer import does not assert that geometry is contained in a body, does not recompute reconstruction coverage, and does not fill gaps.

Dependencies are taken from the active environment. Tests use the pipeline's pinned native Spateo checkout; report actual package versions at runtime. Optional static hosting does not change the single-file output format and is outside this skill's default scope.
