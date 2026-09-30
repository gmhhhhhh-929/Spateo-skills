---
name: spateo-manage-runs
description: Inspect and plan immutable Spateo 4D runs with verified input, implementation and checkpoint hashes and dependency-aware reuse.
---

# spateo-manage-runs

The parent engine owns config normalization, hashes and manifests. Use its `--dry-run --parent-manifest` to inspect reuse before a child run. `scripts/run_manager.py show --manifest PATH` prints a manifest; `plan --config PATH --project PATH --parent-manifest PATH` uses exactly the same planner as execution.

A run is completed, partial or failed. Each stage is pending, running, completed, reused, imported, skipped, failed or blocked. A stopped run stays partial. Disabled optional stages are skipped. A failed stage blocks remaining stages and saves its traceback. Never manually mark missing output as completed.

Each stage's named outputs include SHA256. Input changes at the same path, changed implementation/package versions, missing outputs or hash mismatches invalidate reuse. Reused outputs point to their immutable parent files; preserve those files with the manifest. Completed parent metadata is not edited. The DAG in [config-contract.md](../../references/config-contract.md) keeps trajectory, metrics and GP independent where scientifically appropriate.

## Resume after mapping review

Use this pause only when the user asks to review mapping before downstream analysis. Run the following from the parent skill directory with the intended full-analysis config (`workflow.until="features"`):

```bash
python scripts/run_4d_pipeline.py --config /project/config.json --project /project/analysis --run-id mapping-review --stop-after mapping
```

The manifest remains partial. Use the viewer subskill to present the available mapping, including its alpha, initialization and source/target cell identities. `mapping.initialization="uniform"` is the default; `"aligned_spatial"` supplies balanced Euclidean OT in the registered frame as native `G_init`. Retain the settings of the reviewed candidate; a particular alpha or initializer is not a universal optimum.

Once the user approves a candidate, preserve its manifest and output files unchanged. Use the same scientific config and the approved candidate's manifest as the parent; omit `--stop-after`:

```bash
python scripts/run_4d_pipeline.py --config /project/config.json --project /project/analysis --parent-manifest /project/analysis/runs/mapping-review/manifest.json --dry-run
python scripts/run_4d_pipeline.py --config /project/config.json --project /project/analysis --parent-manifest /project/analysis/runs/mapping-review/manifest.json --run-id approved-mapping-full
```

Before executing the second command, confirm that the plan reuses `alignment` and `mapping` and runs `morphofield`, `trajectory`, `metrics` and `dashboard` for the requested full route. The new run must preserve the approved mapping's source/target ID ordering and output SHA256 values, with `reused_from` pointing to its parent run. Do not silently recompute an approved mapping. If reuse fails because a config, input, implementation or output hash changed, diagnose the difference and restore the matching execution context where appropriate; never edit manifest hashes or stage status to force reuse. Keep the immutable parent available after delivering the child run's integrated HTML.
