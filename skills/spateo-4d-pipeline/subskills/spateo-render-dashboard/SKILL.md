---
name: spateo-render-dashboard
description: Render native Spateo 4D checkpoints as a portable interactive source-target, mapping and morphometric dashboard without changing scientific results.
---

# spateo-render-dashboard

The parent pipeline generates an offline Plotly HTML from verified mapping/metric checkpoints and the manifest. It shows source/target point clouds, mapped endpoints, displacement, vector scale, biological/display groups and available scalar metrics/GLM tables. It does not render the pending 3D backbone workflow or imply GP/trajectory animations are present.

For display-only changes, update dashboard settings and create a child with the parent manifest. Offline mode requires Plotly; it fails if unavailable rather than silently loading a CDN. `dashboard.cdn=true` is an explicit network-dependent alternative. Display subsampling affects only the embedded point cloud; full scientific data stay in H5AD/CSV.

A standalone helper remains available: `python scripts/build_dashboard.py --help`. Use its explicit adata, target-adata, spatial-key, vector-key, mapped-key and manifest options. Save to a new file. Inspect stage labels, counts, finite metrics and source/target frame consistency before delivery.
