---
name: spateo-render-3d-viewer
description: Import reconstructed Spateo surface VTKs and optional tissue point clouds into an offline interactive 3D review HTML with dataset switching, grouped Mesh/Cells controls, and source-model morphology on click. Use for inspection without rebuilding or changing geometry.
---

# Render reconstructed 3D models for review

This is the visualization phase of the 3D pipeline, not the 4D morphogenesis dashboard. Existing models are sufficient: do not require the user to rerun reconstruction or supply H5AD when VTK inputs already exist.

Read [the input contract](references/input-contract.md), then create a config naming the datasets, independent body/tissue meshes, and optional cell selections. The `mesh_manifest` shortcut accepts the list-based manifest from `spateo-reconstruct-mesh/scripts/build_meshes.py`; refined or third-party outputs use explicit `layers`. Never guess annotations from filenames.

Run from this subskill directory in the existing Spateo environment with NumPy, PyVista and Plotly installed:

```bash
python scripts/build_viewer.py --config /results/viewer.json --output-dir /results/review_v1
```

The output directory must not exist. Deliver `index.html` and `viewer_manifest.json`. HTML includes Plotly and geometry, works without a server or CDN, and can be opened directly. Keep the manifest alongside it for provenance; browser rendering does not need the source VTK files. Sharing HTML shares its embedded coordinates and sampled cell IDs. Publishing is a separate user-authorized action.

## Display contract

- All interface controls are English. Each tissue is a collapsed group containing independent Mesh and Cells toggles/opacity sliders. Cells are optional and initially hidden; mesh opacity defaults to 100% except translucent body layers. Users can rotate, zoom, reset the camera, switch datasets, inspect a mesh by clicking, or use its Inspect button.
- Surface picking selects the front-most visible layer, including translucent bodies. Hide the covering layer or use the desired tissue's Inspect button when models overlap; transparency does not imply click-through.
- Preserve the shared native coordinate frame, aspect ratio and units. Datasets are switched independently, not automatically aligned or normalized. Do not fabricate physical units.
- No smoothing, clipping, mesh decimation, coordinate rounding, or cell deletion occurs here. Point sampling is deterministic and display-only, separately budgeted per layer; report source/display counts and selected ID hashes. Mesh geometry is serialized at full precision. This conservative first version favors geometry fidelity over extremely large browser scenes.
- Morphology is computed with `st.tdr.model_morphology` on the imported source-resolution surface, not sampled cells or display approximations. For closed positive-volume meshes, cell density uses only selected source cells geometrically inside the mesh, before preview sampling. Open surfaces keep area/extents but suppress volume, volume ratios and density. Closed edges alone do not establish biological validity or absence of self-intersections.
- Technical paths, hashes, warnings and reconstruction diagnostics belong in the sidecar, not the presentation panels. Essential metric validity and unspecified units remain visible to avoid misleading interpretation. Body-envelope metrics describe the imported display shell, not anatomical growth. No tissue/body containment or coverage guarantee is newly inferred by the viewer.

## Verify before handoff

Run the builder, inspect its manifest and open the resulting HTML in a browser. Check dataset switching, camera rotation/reset, independent Mesh/Cells visibility and opacity, collapsed groups, mesh picking, and morphology display. Confirm input hashes remain unchanged. A script success is not a browser-rendering check; report any unavailable visual check honestly.

Record HTML bytes and source/display geometry counts rather than claiming a fixed browser-memory requirement. All datasets are embedded in memory, but only the selected dataset is drawn. For large inputs, first lower `point_budget`, select fewer layers/datasets, or split outputs; do not silently decimate accepted meshes. Reconstruction artifacts or containment problems route back to the mesh subskill, not to display-only fixes.
