# Installation commands

Resolve `SKILL_DIR` from the installed skill's location. Choose `SPATEO_CHECKOUT` and `SPATEO_ENV` for the user's workspace and environment; quote paths because they may contain spaces. Do not set these variables to the Skills repository itself.

## Complete source checkout

When no source revision was selected, create a complete clone at a new destination and use the skill's fixed source revision:

```bash
git clone https://github.com/gmhhhhhh-929/spateo-release.git "$SPATEO_CHECKOUT"
git -C "$SPATEO_CHECKOUT" checkout --detach d6aa68addc475dd0b56f69cebe7823b1f79933a9
```

For an existing checkout, inspect its revision and local changes first; do not overwrite user changes to reach the suggested revision. This is a source-profile pin, not a lockfile for every resolved dependency.

## New conda environment

Select an unused environment name. Run from the Spateo checkout because `environment.yml` includes the relative editable install `-e .`:

```bash
cd "$SPATEO_CHECKOUT"
conda env create --name "$SPATEO_ENV" --file environment.yml
conda activate "$SPATEO_ENV"
python --version
python -m pip check
```

The supplied YAML selects Python 3.10 and a compatible conda-forge binary stack. The source supports 3.10–3.12, but choosing 3.11 or 3.12 requires an intentionally adjusted environment specification and verification on that interpreter.

## Existing environment

Activate the requested environment and confirm its interpreter. For an already compatible compiled stack, install the chosen checkout with that interpreter:

```bash
cd "$SPATEO_CHECKOUT"
python -m pip install -e .
python -m pip check
```

For an explicitly requested conda-stack update, use `conda env update --name "$SPATEO_ENV" --file environment.yml` from the checkout. That file pins Python 3.10, so account for a user's existing 3.11/3.12 environment before using it. Do not add `--prune` unless removal of unlisted packages was requested. If compiled imports conflict, a fresh environment is usually easier to diagnose than repeatedly replacing wheels.

## Optional 3D tools and development tests

```bash
cd "$SPATEO_CHECKOUT"
python -m pip install -e '.[3d]'
python -m pip check
python "$SKILL_DIR/scripts/verify_environment.py" --checkout "$SPATEO_CHECKOUT" --smoke-test --smoke-test-3d --json
```

The mesh smoke test exercises the PyMCubes route, mesh repair, remeshing, and surface construction. It does not validate Open3D, interactive rendering, every optional dependency, or a complete biological workflow. A headless machine can pass these computational checks without a working graphical display.

When source development/testing is in scope, install `dev-requirements.txt` from this same checkout and run the relevant tests. Preserve failures and interpreter/package versions in the handoff. Once the selected tests pass, additional testing should follow the scope of the change.
