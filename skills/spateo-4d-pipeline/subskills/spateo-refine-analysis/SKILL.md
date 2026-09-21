---
name: spateo-refine-analysis
description: Translate feedback into a new native Spateo 4D config and rerun only stages invalidated by scientific or display changes.
---

# spateo-refine-analysis

Read the parent's config.json and [config-contract.md](../../references/config-contract.md). Make a new config file; preserve unrelated settings and resolve relative inputs against the original file before relocating it. The helper `scripts/update_config.py` supports nested edits; inspect its help before use.

Classify the request by actual scientific effect: alignment→all; subset/mapping→mapping descendants; SparseVFC→field descendants; trajectory→trajectory/dashboard; metrics/GLM→metrics/dashboard; GP→GP/dashboard; dashboard→display only. Do not change scientific settings to satisfy a cosmetic request.

Run the main pipeline with `--dry-run --parent-manifest PATH`, inspect the real hash-verified plan, then execute into a new run. Report changed settings, reused checkpoints and available output comparisons. Do not describe a display resample as a new scientific fit. Old v2 configs/manifests must be migrated explicitly and cannot supply trusted v3 reuse.
