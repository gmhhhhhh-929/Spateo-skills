# Mesh validation record

## Automated contract tests

From the repository root, use a Python 3.10–3.12 environment containing Spateo commit `615644f88613bea8ceb2e2df1e2391d16de55ec1` and its 3D dependencies:

```bash
PYTHONPATH=/path/to/spateo-release \
  python -m pytest tests/native_skills/test_3d_mesh.py -q

python /path/to/skill-creator/scripts/quick_validate.py \
  skills/spateo-3d-pipeline/subskills/spateo-reconstruct-mesh
```

The tests cover body and multi-value annotation selection, density-grid reconstruction, component filtering, automatic coverage failure, triangle-only VTK output, disconnected-component-preserving repair, closed surfaces, preview/manifest creation, full face/array Spateo save/read round trips, deterministic native sampling, absent labels, planar inputs, case-insensitive filename collision refusal, atomic failure cleanup, and overwrite refusal.

## Real-data visual validation — 2026-09-22

The wrapper was run without modifying source data on two local H5AD-derived point clouds:

- Planarian 14 dpa: 85,403 cells, 15 z layers, 33 annotations. The final body, neural, gut, and pharynx meshes were all closed. The density mask covered 97.1–98.8% of selected source points; final mesh containment was 95.7–96.8%.
- Drosophila 16.53 h: 26,020 cells, 27 z layers, 15 tissues. The final body, CNS, muscle, and midgut outputs were closed triangle surfaces. Density-mask coverage was 97.1–99.0% and final containment was 96.9–98.0%. The run was marked `pass_with_warnings` because two recorded body fragments of at most 16 triangles were numerically degenerate and unrepairable; valid disconnected components were retained.

The deterministic direct Spateo-core Drosophila body comparison used `levelset=0.5`, `mc_scale_factor=0.8`, `dist_sample_num=100`, random seed 0, and 300 Laplacian iterations. It had 76 open edges before component-preserving repair, 0 afterward, retained 35 components, and enclosed 98.5% of the 20,000 checked points. The full public `construct_surface` path failed earlier on mixed connectivity inside MeshFix. Native output remains available for comparison; the density route is the default because it exposes explicit volume-field, coverage, component, and smoothing controls for both body and annotation meshes. This one specimen is not a universal benchmark.

These checks establish execution, topology, traceability, and a reviewable visual result. They do not establish anatomical truth, annotation validity, optimal parameters for every tissue, or publication readiness without domain review.
