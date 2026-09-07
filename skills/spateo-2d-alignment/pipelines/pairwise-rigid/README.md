# Spatial Alignment Pipeline Release v0.2.3

Canonical locked code release for spatial alignment candidate generation and
replay. Production jobs should resolve runner entrypoints through
`skill.lock.yaml` and validate finished runs before recipe/compose.

This release intentionally freezes runner copies under `runners/`. Do not edit
files in-place. Publish a new release directory for code changes.

## Runtime Dependencies

Local review/compose scripts use standard scientific Python packages listed in
the repository root `requirements-local.txt`: `numpy`, `pandas`, `scipy`,
`Pillow`, `anndata`, `h5py`, and `PyYAML`.

The Spateo pairwise runners under `runners/` additionally require the project
Spateo environment, including `spateo`, `torch`, and `POT/ot`. Use the
cluster's validated conda environment for these runners; do not replace it with
the local review requirements file.

## Default Review Policy

Use full coordinates for pairwise rigid, spatial-only rescue, slice-level
baseline compose, full candidate replay, and validation. Full replay writes an
internal `audit_points.csv` and exports clean user-facing coordinates with
schema `cell_id,slice_id,stage,chip_id,sl_number,celltype,x,y,z`. Use
full-points HTML viewers for manual issue finding and acceptance review. Use
component-balanced 300k point samples after edits for fast structure checks and
component/triad candidate tuning.

In clean coordinates, `cell_id` is the short `<slice>_<CellID>` ids. If a legacy
source table uses generated ids such as `slice_id:0`, run
`scripts/core/build_h5ad_cell_id_map.py` and pass the resulting map to
`export_clean_coordinates.py` or `compose_alignment_recipes.py`. Generated ids
must be retained only as `workflow_cell_id` in audit/QC sidecars.

The sampling entrypoint defaults to `--target-total-points 300000` and samples
by `slice x component_id` unless overridden. Viewer summaries must record the
source full CSV, full row count, displayed row count, sampling method, seed,
coordinate columns, cell type column, component columns, and display policy.
Full-points and balanced300k viewers can color and filter by Component when
`sample_component_id/sample_component_rank` or `component_id/component_rank`
are present; labels use the form `SL59 rank2 id7`.

`states/*.yaml` are the source of truth for latest baseline, review, accepted,
and final coordinate artifacts. Dashboards must show clean coordinate links from
these state pointers.
