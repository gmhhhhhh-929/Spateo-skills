# Reproducibility Contract

Every run must bind a skill to exact code before it can produce an accepted
recipe or full replay output.

## Locked Code

- Use an immutable release directory such as `../../pipelines/pairwise-rigid`.
- Resolve entrypoints through `skill.lock.yaml`; do not hard-code runner paths
  in ad hoc payloads.
- The lock records release-relative entrypoint paths and sha256 hashes.
- Deployed BGI runs record both the release-relative path and the deployed
  absolute path.
- If any entrypoint hash differs from the lock, stop the run.

## Required Run Metadata

Each step record must include:

- skill name and version;
- pipeline release or Git tag;
- Git commit when available;
- entrypoint path and sha256;
- validator status;
- input paths and hashes for coordinate-changing steps;
- output paths and hashes for full replay outputs;
- recipe path and changed-row summary.
- clean coordinate manifest path and sha256 for any user-facing coordinate CSV.

## Viewer Metadata

Balanced300k viewers must record:

- source full CSV path and sha256;
- full row count;
- displayed row count;
- sampling unit and method;
- seed;
- coordinate columns;
- cell type column;
- whether the viewer is component-balanced.

## Validation Gates

- Spatial-only rescue must prove it used the spatial-only runner and
  `expression_mode=spatial_only`.
- A sampled component recipe must pass sampled-to-full replay validation before
  it is trusted.
- Full replay must preserve row count and cell type values.
- A sampled viewer cannot be recorded as an all-points/full viewer.
- A clean coordinate CSV must pass `validate_clean_coordinates.py`; files with
  intermediate columns such as `manual_x`, `full_candidate_x`, edit labels, or
  displacement columns cannot be linked as user-facing coordinates.
- Dashboard "latest" links must resolve through `states/*.yaml`, not file
  timestamps or ad hoc filename guesses.
