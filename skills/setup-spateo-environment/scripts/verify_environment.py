#!/usr/bin/env python3
"""Verify the active Spateo installation without modifying its environment.

Adapted from gmhhhhhh-929/spateo-release at
d6aa68addc475dd0b56f69cebe7823b1f79933a9 (BSD-2-Clause; see references/provenance.md).
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import platform
import subprocess
import sys
import tempfile
from importlib import import_module, metadata
from pathlib import Path
from urllib.parse import unquote, urlparse

SOURCE_COMMIT = "d6aa68addc475dd0b56f69cebe7823b1f79933a9"
CHECKOUT_FILES = (
    "environment.yml", "requirements.txt", "setup.py", "README.md",
    "dev-requirements.txt", "3d-requirements.txt", "docs/requirements.txt",
    "spateo/__init__.py",
)
REQUIRED = {
    "numpy": ">=1.23.5,<2",
    "pandas": ">=1.5.3,<2.3",
    "scipy": ">=1.10,<1.14",
    "anndata": ">=0.9,<0.12",
    "h5py": ">=3.8,<4",
    "geopandas": ">=0.14,<1.2",
    "shapely": ">=2,<3",
    "pyarrow": ">=12,<26",
    "PyMCubes": ">=0.1.6,<0.2",
    "pymeshfix": ">=0.18.1,<0.19",
    "pyacvd": ">=0.4,<0.5",
}
IMPORT_NAMES = {"PyMCubes": "mcubes"}
# Retained for the upstream mesh smoke function.
version = metadata.version


def _prepare_caches() -> None:
    # Set only this process's environment; never alter shell/global configuration.
    for variable, directory in {
        "MPLCONFIGDIR": "spateo-matplotlib-cache",
        "NUMBA_CACHE_DIR": "spateo-numba-cache",
        "XDG_CACHE_HOME": "spateo-xdg-cache",
    }.items():
        if variable not in os.environ:
            cache_path = Path(tempfile.gettempdir()) / directory
            cache_path.mkdir(parents=True, exist_ok=True)
            os.environ[variable] = str(cache_path)


def _checkout_details(checkout: Path) -> dict[str, object]:
    result = {
        "path": str(checkout),
        "missing_files": [name for name in CHECKOUT_FILES if not (checkout / name).is_file()],
    }
    # A source archive is allowed. Git identity is diagnostic when available.
    if (checkout / ".git").exists():
        try:
            head = subprocess.run(
                ["git", "-C", str(checkout), "rev-parse", "HEAD"],
                capture_output=True, text=True, timeout=10, check=True,
            ).stdout.strip()
            result["commit"] = head
            result["matches_skill_source_commit"] = head == SOURCE_COMMIT
            result["tracked_files_modified"] = bool(subprocess.run(
                ["git", "-C", str(checkout), "status", "--porcelain", "--untracked-files=no"],
                capture_output=True, text=True, timeout=10, check=True,
            ).stdout.strip())
        except (OSError, subprocess.SubprocessError) as exc:
            result["git_error"] = str(exc)
    return result


def _editable_checkout(distribution) -> Path | None:
    direct_url = distribution.read_text("direct_url.json")
    if not direct_url:
        return None
    value = json.loads(direct_url)
    parsed = urlparse(value.get("url", ""))
    if not value.get("dir_info", {}).get("editable") or parsed.scheme != "file":
        return None
    # Local editable installations only. Avoid interpreting a network URL as local.
    if parsed.netloc not in ("", "localhost"):
        return None
    return Path(unquote(parsed.path)).resolve()


def _inside(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def _run_smoke_test() -> dict[str, object]:
    import numpy as np
    from anndata import AnnData
    from scipy import sparse

    import spateo as st
    from spateo._native import sparse_vector_field

    counts = sparse.csr_matrix(np.array([[1, 0, 2], [0, 3, 1], [2, 1, 1], [1, 2, 0]]))
    adata = AnnData(counts)
    adata.obsm["spatial"] = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    result = st.pp.preprocess_spatial(
        adata,
        recipe="generic",
        min_cells=1,
        n_top_genes=2,
        n_pca_components=2,
        inplace=False,
    )
    if result is None or "counts" not in result.layers or "X_pca" not in result.obsm:
        raise RuntimeError("Spateo preprocessing smoke test did not produce the expected layers and PCA.")
    velocity = np.tile(np.asarray([0.5, -0.25]), (adata.n_obs, 1))
    vector_field = sparse_vector_field(adata.obsm["spatial"], velocity, Grid=adata.obsm["spatial"], M=4)
    if not np.isfinite(vector_field["grid_V"]).all():
        raise RuntimeError("Spateo native vector-field smoke test returned non-finite values.")
    return {
        "shape": list(result.shape),
        "pca_shape": list(result.obsm["X_pca"].shape),
        "vector_field_shape": list(vector_field["grid_V"].shape),
    }


def _run_3d_smoke_test() -> dict[str, object]:
    import itertools

    import mcubes
    import numpy as np
    import pyvista as pv

    from spateo.tdr.models.models_individual.mesh_methods import marching_cube_mesh
    from spateo.tdr.models.models_individual.mesh import construct_surface
    from spateo.tdr.models.models_individual.mesh_utils import fix_mesh, uniform_mesh

    grid = np.indices((16, 16, 16), dtype=float)
    volume = (np.square(grid - 7.5).sum(axis=0) <= 25.0).astype(float)
    vertices, triangles = mcubes.marching_cubes(volume, 0.5)
    if vertices.shape[0] == 0 or triangles.shape[0] == 0:
        raise RuntimeError("PyMCubes returned an empty surface for the 3D smoke-test volume.")

    points = np.asarray(list(itertools.product(range(5), repeat=3)), dtype=float)
    mesh = marching_cube_mesh(pv.PolyData(points), levelset=0.5, dist_sample_num=50)
    if mesh.n_points == 0 or mesh.n_cells == 0 or not np.isfinite(mesh.points).all():
        raise RuntimeError("Spateo marching_cube_mesh returned an invalid mesh.")

    open_mesh = (
        pv.Sphere(theta_resolution=30, phi_resolution=30)
        .clip(normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, 0.0), invert=False)
        .extract_surface()
        .triangulate()
        .clean()
    )

    def boundary_count(model: pv.PolyData) -> int:
        return int(
            model.extract_feature_edges(
                boundary_edges=True,
                feature_edges=False,
                manifold_edges=False,
                non_manifold_edges=False,
            ).n_cells
        )

    boundary_edges_before = boundary_count(open_mesh)
    fixed_mesh = fix_mesh(open_mesh)
    boundary_edges_after = boundary_count(fixed_mesh)
    if boundary_edges_before == 0 or boundary_edges_after != 0:
        raise RuntimeError(
            "Spateo fix_mesh did not close the 3D smoke-test surface "
            f"({boundary_edges_before} -> {boundary_edges_after} boundary edges)."
        )

    uniform = uniform_mesh(fixed_mesh.copy(), nsub=1, nclus=100)
    if uniform.n_points == 0 or uniform.n_cells == 0:
        raise RuntimeError("Spateo uniform_mesh returned an empty remeshed surface.")

    surface, inside, _ = construct_surface(
        pc=pv.PolyData(points),
        key_added="smoke_test",
        label="surface",
        cs_method="marching_cube",
        cs_args={"mc_scale_factor": 1.5, "dist_sample_num": 50},
        nsub=1,
        nclus=50,
        smooth=5,
        scale_factor=1.12,
    )
    if surface.n_points == 0 or surface.n_cells == 0 or inside.n_points == 0:
        raise RuntimeError("Spateo construct_surface returned an empty surface or point cloud.")

    return {
        "pymcubes": version("PyMCubes"),
        "pymeshfix": version("pymeshfix"),
        "pyacvd": version("pyacvd"),
        "direct_vertices": int(vertices.shape[0]),
        "direct_triangles": int(triangles.shape[0]),
        "spateo_mesh_points": int(mesh.n_points),
        "spateo_mesh_cells": int(mesh.n_cells),
        "boundary_edges_before_repair": boundary_edges_before,
        "boundary_edges_after_repair": boundary_edges_after,
        "repaired_mesh_points": int(fixed_mesh.n_points),
        "repaired_mesh_cells": int(fixed_mesh.n_cells),
        "uniform_mesh_points": int(uniform.n_points),
        "uniform_mesh_cells": int(uniform.n_cells),
        "construct_surface_points": int(surface.n_points),
        "construct_surface_cells": int(surface.n_cells),
        "construct_surface_inside_points": int(inside.n_points),
    }


def verify(args: argparse.Namespace) -> dict[str, object]:
    report: dict[str, object] = {
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "executable": sys.executable,
        "platform": platform.platform(),
        "skill_source_commit": SOURCE_COMMIT,
        "mode": "metadata-only" if args.metadata_only else "imports-and-api",
        "packages": {},
        "errors": [],
        "smoke_test": ({"status": "skipped", "reason": "prerequisite verification failed"}
                       if args.smoke_test else {"status": "not_requested"}),
        "smoke_test_3d": ({"status": "skipped", "reason": "prerequisite verification failed"}
                          if args.smoke_test_3d else {"status": "not_requested"}),
    }
    errors = report["errors"]
    packages = report["packages"]
    if platform.python_implementation() != "CPython" or not ((3, 10) <= sys.version_info[:2] <= (3, 12)):
        errors.append("This source profile supports CPython 3.10 through 3.12.")
    checkout = Path(args.checkout).expanduser().resolve() if args.checkout else None
    if checkout:
        report["checkout"] = _checkout_details(checkout)
        if report["checkout"]["missing_files"]:
            errors.append("Incomplete Spateo checkout: " + ", ".join(report["checkout"]["missing_files"]))

    # Keep --help and diagnostics usable even before scientific dependencies exist.
    try:
        from packaging.requirements import Requirement
        from packaging.specifiers import SpecifierSet
        from packaging.version import Version
    except ImportError as exc:
        errors.append(f"Cannot inspect version requirements: install packaging in this environment ({exc}).")
        report["ok"] = False
        return report

    if not args.metadata_only:
        try:
            _prepare_caches()
        except OSError as exc:
            errors.append(f"Cannot prepare temporary scientific-package caches: {exc}")
            report["ok"] = False
            return report
    for package, specifier in REQUIRED.items():
        details = {"version": "missing", "required": specifier, "compatible": False}
        try:
            details["version"] = metadata.version(package)
            details["compatible"] = Version(details["version"]) in SpecifierSet(specifier)
            if not args.metadata_only:
                import_module(IMPORT_NAMES.get(package, package.replace("-", "_")))
                details["import"] = "passed"
        except metadata.PackageNotFoundError:
            pass
        except Exception as exc:
            details["compatible"] = False
            details["error"] = f"{type(exc).__name__}: {exc}"
        packages[package] = details
        if not details["compatible"]:
            errors.append(f"{package}: found {details['version']}, expected {specifier}" +
                          (f"; {details['error']}" if "error" in details else ""))

    try:
        distribution = metadata.distribution("spateo-release")
        declared = [Requirement(item).name.lower().replace("_", "-") for item in (distribution.requires or [])]
        report["spateo_distribution"] = {
            "version": distribution.version,
            "metadata_location": str(distribution.locate_file("")),
        }
        report["dynamo_required"] = any(name in {"dynamo", "dynamo-release"} for name in declared)
        if report["dynamo_required"]:
            errors.append("Installed Spateo declares Dynamo; it does not match this native source profile.")
        editable = _editable_checkout(distribution)
        report["spateo_distribution"]["editable_checkout"] = str(editable) if editable else None
        if checkout and editable != checkout:
            errors.append("Spateo's editable-install metadata does not point to --checkout. "
                          "Install with this interpreter using python -m pip install -e . from that checkout.")
    except metadata.PackageNotFoundError:
        errors.append("spateo-release distribution metadata is missing; an importable source tree alone is not an installation.")
    except Exception as exc:
        errors.append(f"spateo-release metadata check failed: {type(exc).__name__}: {exc}")

    if not args.metadata_only:
        try:
            st = import_module("spateo")
            report["spateo"] = {"version": getattr(st, "__version__", None), "path": st.__file__}
            if checkout and (not st.__file__ or not _inside(Path(st.__file__), checkout / "spateo")):
                errors.append("Imported spateo is outside the requested checkout; inspect interpreter, PYTHONPATH and editable installs.")
            for dotted in ("io.read_atera", "io.read_visium", "pp.preprocess_spatial"):
                current = st
                for part in dotted.split("."):
                    current = getattr(current, part)
                if not callable(current):
                    raise TypeError(f"spateo.{dotted} is not callable")
        except Exception as exc:
            errors.append(f"spateo import/API check failed: {type(exc).__name__}: {exc}")

    for requested, key, action in (
        (args.smoke_test, "smoke_test", _run_smoke_test),
        (args.smoke_test_3d, "smoke_test_3d", _run_3d_smoke_test),
    ):
        if requested:
            if errors:
                report[key] = {"status": "skipped", "reason": "prerequisite verification failed"}
            else:
                try:
                    report[key] = {"status": "passed", "details": action()}
                except Exception as exc:
                    report[key] = {"status": "failed", "error": f"{type(exc).__name__}: {exc}"}
                    errors.append(f"{key} failed: {type(exc).__name__}: {exc}")
    report["ok"] = not errors
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", help="Require an editable installation from this complete Spateo checkout; does not alter sys.path.")
    parser.add_argument("--metadata-only", action="store_true", help="Check interpreter/package metadata without importing scientific packages; not a runtime readiness check.")
    parser.add_argument("--smoke-test", action="store_true", help="Run tiny preprocessing and native vector-field computations.")
    parser.add_argument("--smoke-test-3d", action="store_true", help="Run marching-cubes, mesh-repair, remeshing and surface-construction computations.")
    parser.add_argument("--json", action="store_true", help="Emit JSON on stdout; import/workflow chatter goes to stderr.")
    args = parser.parse_args(argv)
    if args.metadata_only and (args.smoke_test or args.smoke_test_3d):
        parser.error("--metadata-only cannot be combined with smoke tests")
    # Library Python prints must not corrupt machine-readable stdout.
    with contextlib.redirect_stdout(sys.stderr):
        report = verify(args)
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(f"Python: {report['python']} ({report['executable']})")
        print(f"Mode: {report['mode']}")
        for package, details in report["packages"].items():
            marker = "OK" if details["compatible"] else "FAIL"
            print(f"[{marker}] {package}: {details['version']} (required {details['required']})")
        for error in report["errors"]:
            print(f"ERROR: {error}", file=sys.stderr)
        message = "Metadata checks passed; runtime imports were not tested." if args.metadata_only else "Environment verification passed."
        print(message if report["ok"] else "Environment verification failed.")
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
