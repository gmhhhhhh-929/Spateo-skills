---
name: spateo-environment-report
description: Audit an existing Spateo interpreter and generate a health report, package snapshot, directory tree and evidence-based downstream capability matrix using isolated synthetic probes.
---

# Spateo Environment Report

Resolve `ENV_SKILL_DIR` to the parent `setup-spateo-environment` directory. Select the absolute interpreter and a new output directory outside the environment and source checkout. If the user refers to an earlier setup, use [discovery](../spateo-environment-discovery/SKILL.md) to recover its identity and provenance first.

```bash
python "$ENV_SKILL_DIR/scripts/audit_environment.py" report \
  --python "$TARGET_PYTHON" \
  --output "$AUDIT_OUTPUT" \
  --profile standard --timeout 120
```

Add `--checkout "$SPATEO_CHECKOUT"` to check an expected source checkout, or `--baseline "$PREVIOUS_REPORT_JSON"` when the request includes drift since an earlier report. Without `--checkout`, infer any editable checkout from installation metadata and check runtime imports against that source; wheel installs are checked against their installed package location. Never add a requested checkout to `sys.path` to make an import pass.

| Profile | Evidence collected |
| --- | --- |
| `metadata` | Interpreter, installed distributions, source metadata, dependency consistency and tree, without scientific imports. This cannot establish runtime readiness. |
| `standard` (default) | Dependency consistency, real scientific imports, tiny synthetic IO/preprocessing/point-cloud/native-field computations, and slice-QC/alignment/mapping/trajectory/feature API checks. |
| `full` | Standard checks plus computational mesh and GPU probes. This is not a full biological pipeline or graphical-display test. |

Use the selected interpreter for each probe from a neutral working directory, with isolated temporary/cache files and a per-probe timeout. Preserve stdout, stderr, return status and timing for success, failure and timeout. Continue independent probes after a failure so the report identifies both usable and blocked capabilities. Do not install dependencies, mutate package versions or use the user's datasets to make the audit pass.

Read [environment-audit.md](../../references/environment-audit.md) for the capability states and interpretation rules. Every downstream row must point to its actual evidence and missing prerequisites. Distinguish a successful synthetic computation from an API-only check; metadata, import success and absent GPU devices must not become claims of full scientific readiness.

Return `environment-report.md`, `environment-report.json`, `packages.json` and `environment-tree.md` with the selected executable, source identity, overall health, task-specific evidence and remaining limits. Keep `checks/` alongside them. If dependency consistency or expected source identity fails, state that in overall health even if independent probes pass. Report concrete next actions as recommendations; perform repair only when the user requests it.

The downstream handoff still requires suitable inputs: IO formats and coordinate semantics, validated slice order and representation for alignment, finite XYZ and unique cell IDs for 3D, and registered paired 3D inputs/annotations for 4D. Environment probes do not certify input data or biological accuracy.
