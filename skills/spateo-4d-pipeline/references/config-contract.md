# Modular config contract · spateo-4d/v4

The executable defaults and strict key validation live in `scripts/pipeline_runtime.py`; start from `assets/config.template.json`.

| Section | Contract and ownership |
| --- | --- |
| `inputs` | Required stage1/stage2 H5AD and shared coordinate_unit; optional pointcloud1/pointcloud2 VTK. Paths resolve relative to config, then persist absolute. Input contents are hashed. |
| `workflow` | entry: alignment/aligned/field; until: alignment/trajectory/features; species for documented same-species identity; frame_id required for imported coordinates/fields and must match existing frame metadata and units. |
| `labels` | Human-readable chronological stage names. No inferred elapsed biological time. |
| `alignment` | Input spatial_key, distinct aligned_key, counts_layer or explicit x_is_counts, separate normalized/log layers; SN-S/SN-N, reference sample count, native iterations/device. Export all input observations. |
| `subset` | annotation_key and one exact group. Mapping requires a group; run labels independently. |
| `mapping` | Native OT alpha/iterations, key and max_pairs guard. Source-to-target direction is fixed by input order. |
| `morphofield` | Native M, lambda_, beta, **max_iter**, tol and restart settings. Old MaxIter was ignored by the backend and is not accepted. |
| `trajectory` | enabled, key, positive t_end, interpolation_num, forward/backward/both. Model time is distinct from biological time. |
| `metrics` | enabled, selected scalar features, glm_metrics subset, glm_genes list or `"*"`, glm_min_cells, glm_top_plots, qval and optional llf thresholds. Complete test tables survive selection. |
| `gp` | Optional native spatial expression interpolation. Disabled by default; explicit genes required. |
| `dashboard` | enabled, max_points/max_target_points (0 = all), max_vectors/max_trajectories (positive explicit display caps), write_html (false for remote payload-only export). Offline renderer embeds Plotly; cdn must remain false. |
| `runtime` | cpu or native device identifier and deterministic seed. Installed library/script hashes and package versions are recorded. |

`--until` overrides workflow.until and still builds the viewer. `--stop-after` deliberately stops at an internal checkpoint. `--dry-run` validates config/imports/hashes but does not fully read matrices or prove input biology. Scientific input checks run before registration/import.

Stages: alignment → mapping → morphofield → trajectory; features (`metrics`) branch from morphofield; GP branches from mapping; viewer consumes all available branches. Display changes invalidate only viewer when implementation and prior hashes match. Source/implementation changes conservatively invalidate all stages. Completed parents are immutable.

V3 migration: set schema to v4, supply workflow and annotation group, replace morphofield.MaxIter with an explicitly chosen native max_iter (default 8, not 500), and decide full display versus caps. Previous stage filenames remain readable as data inputs, but old manifests do not silently count as equivalent new computations.
