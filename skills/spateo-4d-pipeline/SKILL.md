---
name: spateo-4d-pipeline
description: Orchestrate two-timepoint 3D Spateo analysis through four modular skills for registration, same-annotation mapping and trajectories, morphogenesis features and gene associations, and an adaptive offline viewer. Use for full analysis, alignment only, or entry from aligned H5AD or a fitted native field.
---

# Spateo cross-timepoint analysis

Use two timepoints of the same species after 3D reconstruction. This parent routes the request, defines the shared data contract and records immutable runs. Load only the subskills needed for the requested scientific stages; viewer generation accompanies each completed route.

## Four analysis subskills

| User request | Read | Scientific outputs |
| --- | --- | --- |
| Register two stages, including existing 3D point clouds | [spateo-align-stages](subskills/spateo-align-stages/SKILL.md) | Aligned full H5AD pair, ID-preserving VTK pair, before/after QC |
| Start with aligned stages; map matching cell types and infer development | [spateo-morphogenesis](subskills/spateo-morphogenesis/SKILL.md) | Subset VTKs, transport with IDs, displacement, fitted field, trajectories |
| Compute features and find associated genes | [spateo-morphogenesis-features](subskills/spateo-morphogenesis-features/SKILL.md) | Native features/Jacobian, complete GLM tables, fitted curves and figures |
| Review any available stage or change display only | [spateo-render-dashboard](subskills/spateo-render-dashboard/SKILL.md) | Self-contained interactive HTML with only available analysis tabs |

## Shared contract

Confirm temporal direction, same species, actual coordinate key and units, count semantics and biological annotation. Preserve original data and unique cell IDs. A VTK contributes geometry and labels; it cannot recreate an expression matrix. Attach it to the matching H5AD by exact `obs_index` ↔ `obs_names`, never by row order. A cleaned or sampled VTK must first be reconciled explicitly with its corresponding H5AD; do not silently drop cells.

Use native Spateo commit `615644f88613bea8ceb2e2df1e2391d16de55ec1`, audited on 2026-09-24. Dynamo is not a dependency. Read [native differences](references/protocol-migration.md) when translating notebook parameters, interpreting torsion or choosing arrow lengths. Record the actual installed implementation hash. Review annotation-specific orientation using [registration review](references/registration-review.md); whole-animal overlap can conceal a 180° anatomical reversal. Do not install or downgrade packages just to imitate historical plots.

## Run the requested route

Copy [config.template.json](assets/config.template.json) and read [config contract](references/config-contract.md). Run using the validated scientific Python environment:

```bash
python scripts/run_4d_pipeline.py --config /project/config.json --project /project/analysis --dry-run
python scripts/run_4d_pipeline.py --config /project/config.json --project /project/analysis
```

- Alignment only: `--until alignment`; includes its viewer.
- Mapping, field and trajectories: `--until trajectory`; no feature/GLM computation.
- Full analysis: `--until features` (default).
- Start at subskill 2: set `workflow.entry="aligned"`, a declared `workflow.frame_id`, and aligned H5AD inputs. No registration rerun.
- Start at subskill 3: set `workflow.entry="field"` with compatible native field H5AD and target H5AD. No mapping or field refit.
- Viewer only: use the viewer subskill directly; do not execute scientific stages.

When the user requests mapping review before downstream analysis, use `--stop-after mapping` and deliver a mapping-only review through the viewer subskill; the checkpoint command itself does not build a viewer. After approval, follow [resume after mapping review](subskills/spateo-manage-runs/SKILL.md#resume-after-mapping-review) to reuse that exact mapping in a new run. Ordinary full-analysis requests do not require this pause.

`--parent-manifest` enables hash-verified reuse; changed scientific settings invalidate dependent stages. [Run management](subskills/spateo-manage-runs/SKILL.md) and [config refinement](subskills/spateo-refine-analysis/SKILL.md) are supporting utilities, not additional scientific stages.

For local data, use the user-selected local conda environment and keep immutable run directories; record the imported library path and package versions. For remote data, verify SSH/workdir/data/environment from available configuration. Keep expression and H5AD on the server. Set `dashboard.write_html=false`, transfer only the final audited `viewer_payload.json` and small QC artifacts, then render HTML locally. Report visualization caps explicitly; the default point-cloud display includes every cell. Alignment reference sampling is recorded independently of full-cell coordinate export.

## Handoff

Deliver the four-stage route, manifest, aligned H5AD/VTK paths, mapping/trajectory/features artifacts that actually exist and a locally working viewer. Explain skipped stages and failed fits. Mapping is inferred correspondence, GLMs are association tests, and integration time is not calendar time. Uncalibrated coordinates must remain labeled as such. See [validation](references/validation.md) for test evidence and scope.
