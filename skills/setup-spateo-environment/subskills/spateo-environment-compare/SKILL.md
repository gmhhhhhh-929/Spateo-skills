---
name: spateo-environment-compare
description: Compare two saved Spateo environment reports for interpreter, package, source and capability drift without reinstalling software or rerunning scientific analyses.
---

# Spateo Environment Compare

Use saved JSON reports produced by [spateo-environment-report](../spateo-environment-report/SKILL.md). Resolve `ENV_SKILL_DIR` to the parent `setup-spateo-environment` directory. Compare an environment before/after a change, or a selected local/server pair; different environments are a cross-environment comparison, not proof of temporal drift:

```bash
python "$ENV_SKILL_DIR/scripts/audit_environment.py" compare \
  --before "$BEFORE_REPORT_JSON" \
  --after "$AFTER_REPORT_JSON" \
  --output "$COMPARISON_OUTPUT"
```

Report changes in interpreter/platform, package additions/removals/versions, source checkout and revision, dependency health, and per-capability evidence. Preserve which report supplied each observation. Different absolute paths across machines are expected; a path change alone is not an incompatibility.

Compare like evidence with like evidence. A `metadata` report and a `full` report have different coverage: `not_tested` becoming `smoke_passed` means new evidence, not necessarily a repaired environment. Tool-hash changes may also change coverage. A missing field or probe remains unknown; unsupported report schemas fail instead of being reinterpreted. Highlight a failed check replacing the same successful check as a regression candidate, with its raw evidence; do not infer causality from a package version change alone.

Return `environment-comparison.md` and `environment-comparison.json`, with actionable changes separated from coverage differences. The comparison reads saved reports only and does not repair either environment. Its successful exit means the comparison was written, not that either environment is ready. If a fresh runtime check is required, invoke the report subskill with the intended profile and save a new report. For shared evidence conventions, read [environment-audit.md](../../references/environment-audit.md).
