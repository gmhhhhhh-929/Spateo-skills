---
name: spateo-morphogenesis
description: Map one observed matching annotation between aligned 3D timepoints and construct native Spateo vector fields and trajectories, preserving transport, displacement, fitted vectors and seed-cell identities. Does not perform feature GLMs.
---

# 2 · Mapping, vector fields and trajectories

Read [native migration](../../references/protocol-migration.md) and the [config contract](../../references/config-contract.md). Implementation: [morphogenesis_stage.py](../../scripts/morphogenesis_stage.py).

Use the preceding registration checkpoint, or two already aligned H5ADs with `workflow.entry="aligned"`, `workflow.frame_id` and the correct `alignment.aligned_key`. Existing point-cloud geometry may be attached by the registration input contract. An aligned-looking overlay alone does not establish a shared frame.

Set `subset.annotation_key` and a single observed `subset.group` (e.g. `anno=CNS`). Run different annotations independently; never silently allow cross-type transport or fit one smooth field through disconnected tissues. Preserve all selected cells by default. Harmonize genes expressed in both stages and normalize from their count layers after subsetting.

- `st.tdr.cell_directions` computes transport and chooses optimal mapped endpoints. Save `transport.npz` with both ID axes; verify `X_mapping - aligned = V_mapping`. It is not barycentric interpolation or cell lineage tracing.
- Rebuild subset point clouds and save both VTKs with original IDs and raw directions.
- `st.tdr.morphofield_sparsevfc` fits the native kernel field at source XYZ. Store raw directions separately from fitted velocity; save field VTK and fit QC. `morphofield.max_iter`, `beta`, `lambda_`, `M`, restart seeds and `tol` affect the fit. The old `MaxIter` parameter is rejected by the config.
- `st.tdr.morphopath` integrates all source seeds. Use an explicit model horizon and inspect extent/exit QC. Persist the field H5AD, trajectory H5AD and portable `trajectories.npz` (`cells × time × XYZ`, times, seed IDs). Do not truncate scientific paths to make a prettier picture.

```bash
python scripts/morphogenesis_10dpa_14dpa_cns.py --config /project/aligned-pair.json --project /project/analysis
```

This stage-agnostic compatibility command ends after trajectories and builds the viewer. Features/GLMs belong to the next subskill. Mapping has quadratic memory costs, including within-stage distance matrices; `max_pairs` is only a guard, not a total-RAM estimate. An actual larger requirement calls for an explicit subset/budget decision, not silent subsampling.

Optional mesh evaluation: construct surfaces with the existing 3D reconstruction mesh skill from aligned VTKs, evaluate the stored field on mesh vertices, and preserve the mesh frame and source hashes. Meshes are not required to infer the cell-based field. Arrow scale is display-only and must never enter trajectory integration.
