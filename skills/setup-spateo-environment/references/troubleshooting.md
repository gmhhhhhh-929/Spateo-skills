# Spateo environment troubleshooting

## Binary or ABI import errors

If NumPy, SciPy, HDF5, Arrow, Shapely, GeoPandas, OpenCV, or VTK fails to import, avoid repeatedly replacing individual wheels. Create a fresh conda-forge environment from the selected Spateo checkout's `environment.yml`, following [installation.md](installation.md). Run from that checkout because its YAML installs `-e .`.

On Apple Silicon, confirm all compiled packages use the same architecture:

```bash
python -c "import platform; print(platform.machine())"
conda list 'numpy|scipy|h5py|pyarrow|shapely|geopandas'
```

## AnnData resolver conflict

This source profile supports `anndata>=0.9,<0.12`. If pip selects a newer storage API, reinstall using the selected Spateo checkout's requirements instead of forcing an untested AnnData release with `--no-deps`.

## Shapely and PyGEOS warning

A warning about GEOS versions usually means an older PyGEOS-backed GeoPandas stack is mixed with Shapely 1. Create or update to the conda-forge Shapely 2 and GeoPandas 0.14+ stack. Do not set `USE_PYGEOS=0` as a permanent substitute for correcting the environment.

## Numba cache errors

Set a writable cache directory for restricted or cluster environments:

```bash
export NUMBA_CACHE_DIR="$PWD/.cache/numba"
mkdir -p "$NUMBA_CACHE_DIR"
```

For headless test workers, a writable Matplotlib cache also avoids repeated font discovery:

```bash
export MPLCONFIGDIR="$PWD/.cache/matplotlib"
```

## `pip check` reports unrelated notebook conflicts

Treat the check as failed until resolved. For example, an old `nbformat` with a modern `nbconvert` should be repaired by installing the repository requirements, which constrain both packages to a compatible range.

## Jupyter still executes a pre-upgrade function

Python keeps imported function objects in the running kernel even after pip replaces their source files. A traceback that mentions an old argument but displays a different current source line is evidence of this stale-code state.

Before restarting, save the notebook and confirm that its file on disk is non-empty. Then restart the kernel and rerun imports and cells that create analysis objects. Do not terminate a live kernel for an unsaved or zero-byte notebook merely to refresh a package.

Confirm the refreshed implementation inside the new kernel when needed:

```python
import inspect
from spateo.tdr.models.models_individual.mesh_utils import fix_mesh

print(inspect.getsource(fix_mesh))
```

For PyMeshFix 0.18, the displayed implementation must call `meshfix.repair()` without a `verbose` keyword.

## Wrong editable checkout or missing package metadata

Use `python -m pip` with the same interpreter that runs the verifier. Inspect `python -m pip show spateo-release`, the JSON report's editable-checkout path, and the imported `spateo` path. Run the verifier from a neutral directory and inspect a user-set `PYTHONPATH` if a different source tree shadows the installation. Do not add the Skills repository to `PYTHONPATH`.

From the intended complete Spateo checkout, `python -m pip install -e .` refreshes its installation metadata. The `--checkout` option requires that editable identity; omit the option when deliberately verifying a non-editable installation. An importable source directory with missing distribution metadata is not sufficient to pass.

## Python and diagnostics before installation

`--help` needs only the Python standard library. `--metadata-only --json` requires `packaging` but avoids scientific imports and reports missing or incompatible packages without installing them. It is useful for diagnosis, but passing metadata checks does not establish runtime readiness. CPython 3.10–3.12 is the supported profile; Python 3.9 or 3.13 requires a different compatible environment.
