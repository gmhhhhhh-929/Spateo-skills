---
name: spateo-manage-runs
description: Inspect and plan immutable Spateo 4D runs with verified input, implementation and checkpoint hashes and dependency-aware reuse.
---

# spateo-manage-runs

The parent engine owns config normalization, hashes and manifests. Use its `--dry-run --parent-manifest` to inspect reuse before a child run. `scripts/run_manager.py show --manifest PATH` prints a manifest; `plan --config PATH --project PATH --parent-manifest PATH` uses exactly the same planner as execution.

A run is completed, partial or failed. Each stage is pending, running, completed, reused, skipped, failed or blocked. A stopped run stays partial. Disabled optional stages are skipped. A failed stage blocks remaining stages and saves its traceback. Never manually mark missing output as completed.

Each stage's named outputs include SHA256. Input changes at the same path, changed implementation/package versions, missing outputs or hash mismatches invalidate reuse. Reused outputs point to their immutable parent files; preserve those files with the manifest. Completed parent metadata is not edited. The DAG in [config-contract.md](../../references/config-contract.md) keeps trajectory, metrics and GP independent where scientifically appropriate.
