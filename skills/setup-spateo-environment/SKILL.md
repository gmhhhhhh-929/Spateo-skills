---
name: setup-spateo-environment
description: Create, diagnose, and verify a Spateo Python environment from a separate source checkout. Use for Spateo installation, dependency or binary conflicts, and environment preparation for IO, preprocessing, alignment, and native morphogenesis.
---

# Setup Spateo Environment

Prepare a reproducible environment for the supported Spateo source profile. This skill contains setup guidance and a verifier; the Skills repository is not an installable Spateo package.

## Establish the installation target

- Resolve `SKILL_DIR` to this skill directory and `SPATEO_CHECKOUT` to a separate, complete Spateo source checkout. Never infer a package root from this skill's parent directories.
- Use the user's selected environment and source revision when supplied. Otherwise prefer a new conda environment and the fixed source revision in [references/provenance.md](references/provenance.md). Keep an existing environment intact unless an update was requested.
- Record the interpreter, platform/architecture, environment name, and source commit. This profile supports CPython 3.10–3.12; the source environment file defaults to 3.10. Another source revision may require different constraints.
- Inspect `environment.yml`, `requirements.txt`, and `setup.py` from the same checkout. `setup.py` also reads `README.md`, `dev-requirements.txt`, `3d-requirements.txt`, and `docs/requirements.txt`; do not install from an incomplete sparse checkout.

## Install or update

Use [references/installation.md](references/installation.md) for the complete-clone, conda, editable-install, and optional 3D commands. Run installation commands from `SPATEO_CHECKOUT`: its environment YAML contains `pip: - -e .`, so using another working directory installs the wrong project.

Prefer conda-forge for the compiled stack. Preserve the source constraints, particularly NumPy `<2`, AnnData `<0.12`, Shapely 2, and modern GeoPandas. Core mesh support includes `PyMCubes>=0.1.6,<0.2` (import name `mcubes`), `pymeshfix>=0.18.1,<0.19`, and `pyacvd>=0.4,<0.5`; the optional `3d` extra adds further tools.

This source profile uses native sampling, normalization, graph, and vector-field code and does not require Dynamo. Do not uninstall Dynamo merely because it is present for another workflow. Spateo metadata that still declares Dynamo signals a different installation/profile.

## Verify the actual environment

Run under the selected interpreter, from a neutral working directory rather than a directory containing another `spateo/` tree:

```bash
python -m pip check
python "$SKILL_DIR/scripts/verify_environment.py" --checkout "$SPATEO_CHECKOUT" --smoke-test --json
```

`--checkout` requires a complete source tree, matching editable-install metadata, and an imported Spateo module inside that checkout. The verifier never inserts this Skills repository or a requested checkout into `sys.path`. Omit `--checkout` to diagnose the installation already selected by the interpreter without requiring an editable source tree.

- `--metadata-only --json` diagnoses versions and installation metadata without scientific imports. A pass in this mode does not establish runtime readiness.
- `--smoke-test` computes tiny preprocessing/PCA and native vector-field examples.
- `--smoke-test-3d` additionally checks marching cubes, mesh repair, remeshing, and surface construction. It does not test every optional 3D backend or a graphical display.
- A nonzero exit, missing package, import failure, or failed `pip check` is a failed check. Report skipped smoke tests as skipped.

For development changes, run the relevant source tests in that checkout after installing its development requirements. The initial IO/preprocessing checks are `python -m pytest tests/io tests/preprocessing`; use the full suite when validating broader changes. Read [references/troubleshooting.md](references/troubleshooting.md) for ABI, cache, dependency, and notebook failures.

## Handoff and scope

Report the selected source commit, Python executable, dependency check, requested smoke-test outcomes, and remaining untested capabilities. Keep machine paths in the user's local report, not in reusable skill files.

This is the environment stage for the first three packaged skills: setup, data preprocessing, and serial alignment. Check each downstream skill's own input and runtime contracts before running it. The alignment skill's `scripts/check_runtime.py` checks its runner imports and API signatures; passing either verifier does not establish a full CPU/GPU alignment run. The two later stages of the broader pipeline are not supplied by this publication, and the native vector-field smoke test does not make them available.
