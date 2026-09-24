---
name: spateo-align-stages
description: Register two same-species 3D timepoints using native Spateo, optionally attaching reconstruction point-cloud VTKs to paired H5AD by cell identity, and export aligned H5AD/VTK with before/after QC.
---

# 1 · Cross-timepoint registration

Read the parent [config contract](../../references/config-contract.md). Implementation: [alignment_stage.py](../../scripts/alignment_stage.py).

1. Verify unique cell/gene IDs, finite full-rank XYZ, matching coordinate units and nonnegative captured counts. Set `x_is_counts=true` only after checking X. Do not infer physical units from coordinate magnitudes.
2. Optional `inputs.pointcloud1/pointcloud2` imports previous 3D reconstruction geometry. Each VTK must have unique `obs_index` exactly matching the paired H5AD. Reorder by ID; preserve pre-import coordinates under `<spatial_key>_before_vtk`. Bare VTK lacks counts and cannot support expression-informed registration or later GLMs: locate the original H5AD.
3. Normalize counts and create a separate log layer with native `st.pp` APIs. Register common genes using `st.align.morpho_align_ref`. `SN-S` is the rigid default; `SN-N` is an explicit scientific choice that may absorb developmental deformation. Keep one declared output coordinate key for downstream analyses.
4. Export both full-cell H5ADs and `st.tdr.construct_pc` point clouds. Save/read through Spateo VTK APIs and verify coordinates and IDs after round trip. Reference sampling must not discard full-cell outputs.
5. Review the overlay, annotation continuity and geometry. Nearest-neighbor QC describes overlap, not anatomical correctness; growth may preclude perfect overlap.

```bash
python scripts/alignment_10dpa_14dpa.py --config /project/pair.json --project /project/analysis
```

The compatibility filename is stage-agnostic; it runs `--until alignment` and the adaptive viewer. No mapping, field or feature analysis is implied. Main artifacts are `stage1.h5ad`, `stage2.h5ad`, `stage1.vtk`, `stage2.vtk`, `qc.json` and the viewer payload/HTML, all indexed in the run manifest.

Mesh reconstruction is optional. Reuse the installed 3D reconstruction mesh subskill on the **aligned** VTK when needed. Do not overlay an untransformed pre-registration mesh on aligned cells.
