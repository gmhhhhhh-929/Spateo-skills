---
name: spateo-environment-tree
description: Create a bounded, readable directory-tree document for a selected Python environment, including interpreter, package locations and editable Spateo source without scanning unrelated files.
---

# Spateo Environment Tree

Resolve `ENV_SKILL_DIR` to the parent `setup-spateo-environment` directory and select an absolute `TARGET_PYTHON` executable. Use discovery first if its location is unknown.

```bash
python "$ENV_SKILL_DIR/scripts/audit_environment.py" tree \
  --python "$TARGET_PYTHON" \
  --output "$AUDIT_OUTPUT" \
  --depth 3 --max-entries 500
```

Add `--checkout "$SPATEO_CHECKOUT"` only for an explicitly selected source tree. Omit it to document the installation selected by the interpreter. Do not infer a Spateo checkout from the Skills repository.

The tree describes the selected interpreter's prefix, `sysconfig` installation locations, package metadata and editable source when available. Keep the environment and source checkout distinguishable. Limits apply per root, with at most 60 entries listed from each directory; show symlinks as links without following them, and report inaccessible entries, depth limits and truncation rather than implying a complete inventory. The bounded tree is a navigational document, not a disk scan or a package-integrity check.

Explain the main directories in the generated document: executables, package locations, conda metadata when present, and external editable source. Return `environment-tree.md` and `environment-tree.json`; `packages.json` and `checks/inventory.json` retain supporting metadata. Markdown is a standalone document that can be opened locally. Use a separate requested document format only when the user asks for it. The tree alone does not establish import health or downstream readiness; route those questions to [spateo-environment-report](../spateo-environment-report/SKILL.md).

Read [environment-audit.md](../../references/environment-audit.md) for shared output conventions. Keep audit outputs outside the inspected environment and source checkout.
