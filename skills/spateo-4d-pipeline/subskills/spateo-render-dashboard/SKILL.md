---
name: spateo-render-dashboard
description: Build an adaptive offline Spateo 4D viewer from verified registration, mapping, field, trajectory and feature artifacts, showing only available stages and actual GLM fits without recomputing science.
---

# 4 · Adaptive interactive review

Read [viewer contract](../../references/viewer-contract.md). Renderer: [viewer_payload.py](../../scripts/viewer_payload.py), template: [viewer.html](../../assets/viewer.html).

The runner builds a validated portable payload from hashed checkpoints. The viewer combines sibling trajectory and feature artifacts by explicit cell identities; it must not assume the metric H5AD also contains trajectories. Missing/skipped stages do not acquire fake data or active controls.

```bash
python ../../scripts/viewer_payload.py --payload /review/viewer_payload.json --output /review/index.html
```

For a remote run, transfer its final viewer payload and small QC artifacts, then execute this command locally using Python with Plotly. Expression matrices and source H5AD remain remote. The HTML embeds Plotly and data and works without a network connection.

The UI supports registration before/after, stage visibility, annotation filtering, shared 3D camera, point size/opacity, raw mapping versus fitted field, visible arrowheads and independent arrow-length control, mapped endpoints, trajectory model-time playback, feature coloring, GLM gene selection with fitted mean/interval, and CSV/PNG/view-setting export. Its run-evidence page records states and QC. Show total/displayed cells and trajectory counts. Full cells are the default; explicit display caps never modify scientific files.

Validate actual interactions in a browser for each supported route: alignment-only, imported aligned pair, full analysis, and display-only rerender. Check no JS errors, working offline, source/target controls, nonempty paths, and feature-to-gene agreement. Technical manifests stay in the evidence page, not the primary scene. Report any scientific QC concern visible in the final overlay; a functioning viewer does not prove correct registration.

The older `scripts/build_dashboard.py` remains a legacy standalone interface; the modular pipeline uses the renderer linked above. Do not extend the old static feature dashboard for new modular runs.
