# Validation and limits

## Reproduce

From the repository root, use a Python 3.10–3.12 environment containing the pinned Spateo source and its 3D dependencies:

```bash
PYTHONPATH=/path/to/spateo-release python -m pytest \
  tests/native_skills/test_3d_point_cloud.py -q
python scripts/validate_collection.py
python /path/to/skill-creator/scripts/quick_validate.py \
  skills/spateo-3d-pipeline
python /path/to/skill-creator/scripts/quick_validate.py \
  skills/spateo-3d-pipeline/subskills/spateo-reconstruct-point-cloud
```

The behavioral tests use real `st.tdr.construct_pc`, `st.tdr.save_model`, and `st.tdr.read_model` calls on synthetic AnnData. They cover categorical obs colors and masks, continuous single/multi-gene expression, exact-ID external labels, input failures, overwrite refusal, a headless four-view preview, and VTK round-trip preservation.

## Interpretation

Passing tests establish the wrapper and file contract for the pinned source revision. They do not establish biological validity, coordinate registration, physical units, meaningful annotations, suitable expression normalization, visual adequacy for every tissue, or readiness for surface/voxel/cell reconstruction. Preview sampling affects only the PNG; the VTK always retains every input observation.

The reference CS13 notebook uses a previously sampled 100,000-cell AnnData and `celltype` colors. Those are example choices, not requirements of this general point-cloud skill. Mesh reconstruction parameters in that notebook are outside this implemented phase.
