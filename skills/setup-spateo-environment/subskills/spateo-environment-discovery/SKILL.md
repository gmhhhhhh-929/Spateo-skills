---
name: spateo-environment-discovery
description: Find candidate Spateo Python environments from local conda records and bounded search roots, or inspect a known SSH server, while preserving evidence of how each environment was found.
---

# Spateo Environment Discovery

Use this when an interpreter is unknown, especially when the user asks to revisit an environment prepared in an earlier session. Resolve `ENV_SKILL_DIR` to this subskill's parent `setup-spateo-environment` directory, two levels above this file's directory.

```bash
python "$ENV_SKILL_DIR/scripts/audit_environment.py" discover \
  --output "$AUDIT_OUTPUT" \
  --search-root "$KNOWN_ENV_ROOT"
```

`--search-root` is optional and repeatable. Discovery uses the current interpreter, active environment variables, local conda environment records (including `~/.conda/environments.txt`) and bounded roots supplied for this task. Explicit roots are searched to depth 4 with a shared 2,000-entry cap. Prefer known environment directories from prior work; do not recursively search an entire home directory or disk. A candidate is not a verified working environment. Retain its executable, discovery source, missing/unreachable state and any metadata errors in `environment-discovery.md` and `environment-discovery.json`.

For a server, use an existing SSH alias or host identified in the user's context and run the same discovery command there through SSH. This script does not scan hosts or guess credentials. Use noninteractive connection checks and preserve connection failures as such; an unreachable server is not evidence that its environments are absent. Transfer only the audit helper and required support files when needed, then retrieve the generated reports. Discovery and probing on the server must not install packages or read biological datasets.

Check the relevant previous setup record, manifest or session evidence before calling a candidate “created by this skill.” A matching name or the presence of Spateo alone establishes only a candidate. Record “origin unconfirmed” if creation evidence is unavailable; it can still be audited.

Select the candidate matching the user's project, source checkout and interpreter. If several remain relevant, list their differences and audit the most relevant candidate(s) within the request. Pass the exact absolute executable to [spateo-environment-tree](../spateo-environment-tree/SKILL.md) or [spateo-environment-report](../spateo-environment-report/SKILL.md); do not rely on shell activation to select the environment.

Keep host identifiers and absolute paths in the user's audit artifacts. For output and evidence conventions, read [environment-audit.md](../../references/environment-audit.md).
