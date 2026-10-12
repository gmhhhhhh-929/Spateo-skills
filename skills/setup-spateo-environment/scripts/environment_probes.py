#!/usr/bin/env python3
"""Run one bounded environment-audit probe with the selected Python interpreter.

Only the standard library is imported at startup. Each runtime probe belongs in
its own subprocess, with a timeout enforced by the calling audit. Results go to
--result-file because native libraries may write directly to stdout. The result
file must not already exist. No checkout is injected into sys.path.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib
import importlib.util
import inspect
import json
import os
import platform
import re
import shutil
import socket
import sys
import sysconfig
import tempfile
import time
from importlib import metadata
from pathlib import Path
from urllib.parse import unquote, urlparse


def _result(status, level, details):
    return {"status": status, "level": level, "details": details}


def _safe_error(exc):
    # Errors occasionally include authenticated package-index URLs. They are
    # diagnostic text, not a reason to copy URL credentials into an audit.
    message = re.sub(r"\b(?:https?|ssh|git\+https?|git\+ssh)://[^\s]+", "<remote URL omitted>", str(exc))
    return f"{type(exc).__name__}: {message}"[:4000]


def _editable_checkout(distribution):
    """Return a local editable source path, never the raw direct_url record."""
    value = distribution.read_text("direct_url.json")
    if not value:
        return None
    try:
        data = json.loads(value)
        if not isinstance(data, dict) or not isinstance(data.get("dir_info"), dict):
            return None
        if not data["dir_info"].get("editable") or not isinstance(data.get("url"), str):
            return None
        parsed = urlparse(data["url"])
        if parsed.scheme != "file" or parsed.netloc not in ("", "localhost"):
            return None
        if parsed.query or parsed.fragment:
            return None
        path = Path(unquote(parsed.path))
        if not path.is_absolute() or "\x00" in str(path):
            return None
        return str(path.resolve())
    except (ValueError, TypeError, OSError):
        return None


def _inventory():
    packages = []
    for distribution in metadata.distributions():
        name = distribution.metadata.get("Name")
        if name:
            packages.append({"name": str(name), "version": str(distribution.version)})
    packages.sort(key=lambda item: (item["name"].casefold(), item["version"]))
    details = {
        "executable": sys.executable,
        "prefix": sys.prefix,
        "base_prefix": sys.base_prefix,
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "hostname": socket.gethostname(),
        "cpu_count": os.cpu_count(),
        "sysconfig_paths": sysconfig.get_paths(),
        "packages": packages,
        "spateo_distribution": None,
        "scope": "Metadata inventory only; no scientific package was imported or executed.",
    }
    try:
        usage = shutil.disk_usage(sys.prefix)
        details["prefix_filesystem"] = {"total_bytes": usage.total, "free_bytes": usage.free}
    except OSError as exc:
        details["prefix_filesystem_error"] = _safe_error(exc)
    try:
        distribution = metadata.distribution("spateo-release")
        details["spateo_distribution"] = {
            "version": distribution.version,
            "metadata_location": str(distribution.locate_file("")),
            "editable_checkout": _editable_checkout(distribution),
        }
    except metadata.PackageNotFoundError:
        details["spateo_distribution_missing"] = True
    return _result("passed", "metadata", details)


def _resolve(root, dotted):
    value = root
    for name in dotted.split("."):
        value = getattr(value, name)
    if not callable(value):
        raise TypeError(f"{dotted} is not callable")
    return value


def _check_apis(root, requirements):
    """Require explicit named parameters; **kwargs alone is not evidence."""
    checked = []
    for dotted, required in requirements.items():
        function = _resolve(root, dotted)
        parameters = inspect.signature(function).parameters
        missing = set(required) - set(parameters)
        if missing:
            raise RuntimeError(f"{dotted} missing explicit parameters: {sorted(missing)}")
        checked.append({"api": dotted, "required_parameters": sorted(required)})
    return checked


def _verifier_helper(name):
    # Load only the sibling helper module. Its verify() and package compatibility
    # gate are intentionally not called, so one missing optional stack does not
    # block unrelated capabilities. Never load a Spateo checkout this way.
    helper = Path(__file__).with_name("verify_environment.py")
    spec = importlib.util.spec_from_file_location("_environment_audit_verifier", helper)
    if spec is None or spec.loader is None:
        raise ImportError("Cannot load sibling verify_environment.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return getattr(module, name)()


def _io():
    import numpy as np
    import anndata as ad
    import spateo as st
    from scipy import sparse

    apis = _check_apis(st, {"io.read_atera": (), "io.read_visium": ()})
    counts = np.array([[1, 0, 2], [0, 3, 1], [2, 1, 0]], dtype=np.int32)
    data = ad.AnnData(sparse.csr_matrix(counts))
    data.obs_names = ["audit_cell_0", "audit_cell_1", "audit_cell_2"]
    data.var_names = ["audit_gene_0", "audit_gene_1", "audit_gene_2"]
    data.layers["counts"] = data.X.copy()
    data.obsm["spatial"] = np.array([[0, 0], [1, 0], [0, 1]], dtype=float)
    with tempfile.TemporaryDirectory(prefix="spateo-audit-io-") as tmp:
        path = Path(tmp) / "synthetic.h5ad"
        data.write_h5ad(path)
        loaded = ad.read_h5ad(path)
    if not np.array_equal(loaded.X.toarray(), counts):
        raise RuntimeError("H5AD round-trip changed expression values")
    if not np.array_equal(loaded.layers["counts"].toarray(), counts):
        raise RuntimeError("H5AD round-trip changed the counts layer")
    if not np.array_equal(loaded.obsm["spatial"], data.obsm["spatial"]):
        raise RuntimeError("H5AD round-trip changed spatial coordinates")
    if loaded.obs_names.tolist() != data.obs_names.tolist() or loaded.var_names.tolist() != data.var_names.tolist():
        raise RuntimeError("H5AD round-trip changed cell or gene identifiers")
    return _result("passed", "smoke", {
        "shape": list(loaded.shape), "h5ad_roundtrip": True, "checked_apis": apis,
        "scope": "Synthetic sparse H5AD values, layers, coordinates and IDs round-tripped. Atera/Visium reader APIs were checked but assay files were not read.",
    })


def _preprocess():
    details = _verifier_helper("_run_smoke_test")
    details["scope"] = "Existing verifier smoke: four cells, spatial preprocessing/PCA and a tiny native vector-field fit. No real-data preprocessing or memory-scale test."
    return _result("passed", "smoke", details)


def _slice_qc():
    # A skill-bundled runtime overlay is deliberately not activated by an audit.
    # Its existence elsewhere does not establish this installation's capability.
    module = importlib.import_module("spateo.preprocessing.slice_quality")
    apis = _check_apis(module, {
        "SliceQCConfig": ("window", "k_neighbors", "profile_genes"),
        "calculate_slice_quality": ("adata", "slice_key", "spatial_key", "config"),
        "scan_h5ad_series": ("inputs", "config", "slice_key", "spatial_key"),
        "write_slice_quality_outputs": ("result", "output_dir"),
    })
    return _result("api_only", "api", {
        "checked_apis": apis, "module_path": getattr(module, "__file__", None),
        "scope": "Installed serial-slice QC API signatures only. No assay, metric calculation, policy calibration or HTML renderer was executed; no bundled overlay was activated.",
    })


def _alignment():
    # This is the explicit keyword contract of spateo-2d-alignment/check_runtime.py.
    for name in ("numpy", "pandas", "scipy", "anndata", "h5py", "matplotlib", "torch", "ot"):
        importlib.import_module(name)
    module = importlib.import_module("spateo.alignment.morpho_alignment")
    methods = importlib.import_module("spateo.alignment.methods")
    apis = _check_apis(module, {
        "morpho_align": ("models", "rep_layer", "rep_field", "spatial_key", "key_added", "iter_key_added", "vecfld_key_added", "mode", "dissimilarity", "max_iter", "dtype", "device", "verbose"),
        "morpho_align_ref": ("models", "models_ref", "n_sampling", "sampling_method"),
    })
    apis.extend(_check_apis(methods, {
        "Morpho_pairwise.__init__": ("nn_init", "SVI_mode", "batch_size", "pre_compute_dist", "sparse_calculation_mode", "sparse_top_k", "use_chunk", "chunk_capacity", "K", "beta", "lambdaVF", "partial_robust_level", "sigma2_init_scale", "sigma2_end", "use_hvg", "init_transform", "init_layer", "init_field", "allow_flip", "nonrigid_start_iter"),
    }))
    return _result("api_only", "api", {
        "checked_apis": apis, "alignment_executed": False,
        "scope": "Alignment dependency imports and explicit keyword signatures only. No registration, GPU alignment, large-data run, or alignment accuracy was tested.",
    })


def _point_cloud():
    import numpy as np
    import anndata as ad
    import spateo as st

    data = ad.AnnData(np.ones((4, 2), dtype=np.float32))
    data.obs_names = [f"audit_cell_{i}" for i in range(4)]
    data.obs["audit_group"] = ["a", "a", "b", "b"]
    data.obsm["spatial"] = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    cloud, _ = st.tdr.construct_pc(
        adata=data, spatial_key="spatial", groupby="audit_group", key_added="audit_group",
        colormap={"a": "#336699", "b": "#cc6633"}, alphamap=1.0, mask=None,
    )
    with tempfile.TemporaryDirectory(prefix="spateo-audit-pc-") as tmp:
        path = Path(tmp) / "synthetic.vtk"
        st.tdr.save_model(cloud, str(path), binary=True)
        loaded = st.tdr.read_model(str(path))
    if cloud.n_points != data.n_obs or loaded.n_points != data.n_obs:
        raise RuntimeError("Point-cloud construction or VTK round-trip changed point count")
    if (loaded.points.shape != data.obsm["spatial"].shape
            or not np.allclose(loaded.points, data.obsm["spatial"], rtol=0, atol=1e-7)):
        raise RuntimeError("Point-cloud round-trip changed spatial coordinates")
    for name, expected in (("obs_index", data.obs_names.to_numpy()), ("audit_group", data.obs["audit_group"].to_numpy())):
        if not np.array_equal(np.asarray(loaded.point_data[name]).astype(str), np.asarray(expected).astype(str)):
            raise RuntimeError(f"Point-cloud round-trip changed point_data[{name!r}]")
    return _result("passed", "smoke", {
        "n_points": int(loaded.n_points), "vtk_roundtrip": True,
        "point_arrays": sorted(loaded.point_data.keys()),
        "scope": "Four-point native Spateo construct_pc/save_model/read_model round-trip with coordinates, cell IDs and categorical labels. No GUI, render backend or scientific surface quality was tested.",
    })


def _mesh():
    details = _verifier_helper("_run_3d_smoke_test")
    details["scope"] = "Existing headless verifier smoke: tiny marching cubes, native mesh construction, repair, remeshing and enclosing surface. No rendering or real-tissue morphology validation."
    return _result("passed", "smoke", details)


def _morphogenesis():
    import itertools
    import numpy as np
    import spateo as st
    from spateo import _native

    apis = _check_apis(st, {
        "tdr.cell_directions": ("adataA", "adataB", "spatial_key", "key_added"),
        "tdr.morphofield_sparsevfc": ("adata", "spatial_key", "V_key", "key_added", "NX", "M", "restart_num"),
        "tdr.morphopath": ("adata", "vf_key", "key_added", "t_end", "cores"),
    })
    apis.extend(_check_apis(_native, {"sparse_vector_field": ("X", "Y", "Grid", "M", "max_iter")}))
    coordinates = np.asarray(list(itertools.product((0.0, 1.0), repeat=3)))
    velocity = 0.15 * coordinates + np.asarray([0.5, -0.25, 0.125])
    fitted = _native.sparse_vector_field(coordinates, velocity, Grid=coordinates, M=8, max_iter=8)
    grid_velocity = np.asarray(fitted["grid_V"])
    if grid_velocity.shape != (8, 3) or not np.isfinite(grid_velocity).all():
        raise RuntimeError("Native 3D vector-field smoke produced invalid grid velocities")
    return _result("passed", "smoke", {
        "checked_apis": apis, "native_module_path": getattr(_native, "__file__", None),
        "vector_field_shape": list(grid_velocity.shape), "max_iter": 8,
        "scope": "Eight-point native 3D vector-field fit executed. Mapping, field wrapper and trajectory APIs were checked only; no transport, trajectory integration or complete 4D workflow was run.",
    })


def _features():
    import spateo as st

    requirements = {
        "tdr.morphofield_" + name: ("adata", "vf_key", "key_added")
        for name in ("velocity", "acceleration", "curvature", "curl", "divergence", "torsion", "jacobian")
    }
    requirements["tl.glm_degs"] = ("adata", "genes", "layer", "key_added", "fullModelFormulaStr", "reducedModelFormulaStr", "qval_threshold", "llf_threshold")
    apis = _check_apis(st, requirements)
    return _result("api_only", "api", {
        "checked_apis": apis,
        "scope": "Morphogenesis-feature and glm_degs signatures only. No derivatives, model fitting, gene associations, significance or numerical correctness was tested.",
    })


def _gpu():
    import torch

    torch.set_num_threads(1)
    x = torch.tensor([1.0, 2.0, 3.0], device="cpu")
    if float(torch.dot(x, x).item()) != 14.0:
        raise RuntimeError("Tiny PyTorch CPU tensor computation failed")
    details = {
        "torch_version": str(torch.__version__), "torch_path": getattr(torch, "__file__", None),
        "torch_cuda_version": getattr(torch.version, "cuda", None),
        "cpu_tensor_passed": True, "cuda_available": bool(torch.cuda.is_available()),
        "scope": "Tiny PyTorch CPU/CUDA dot product only. Does not establish GPU support, precision, memory capacity or performance for any Spateo pipeline.",
    }
    if not details["cuda_available"]:
        details["reason"] = "CUDA is unavailable to this interpreter/process; CPU tensor smoke passed."
        return _result("unsupported", "smoke", details)
    device = torch.device("cuda:0")
    x_cuda = x.to(device)
    answer = float(torch.dot(x_cuda, x_cuda).item())
    torch.cuda.synchronize(device)
    if answer != 14.0:
        raise RuntimeError("Tiny PyTorch CUDA tensor computation failed")
    details.update(cuda_tensor_passed=True, cuda_device_index=0, cuda_device_name=torch.cuda.get_device_name(0))
    return _result("passed", "smoke", details)


PROBES = {
    "inventory": _inventory, "io": _io, "preprocess": _preprocess,
    "slice-qc": _slice_qc, "alignment": _alignment, "point-cloud": _point_cloud,
    "mesh": _mesh, "morphogenesis": _morphogenesis, "features": _features, "gpu": _gpu,
}
LEVELS = {name: ("metadata" if name == "inventory" else "api" if name in {"slice-qc", "alignment", "features"} else "smoke") for name in PROBES}


@contextlib.contextmanager
def _temporary_runtime():
    """Keep common scientific caches and process thread counts temporary."""
    with tempfile.TemporaryDirectory(prefix="spateo-audit-runtime-") as tmp:
        settings = {
            "MPLCONFIGDIR": str(Path(tmp) / "matplotlib"),
            "NUMBA_CACHE_DIR": str(Path(tmp) / "numba"),
            "XDG_CACHE_HOME": str(Path(tmp) / "xdg"),
            "PYVISTA_USERDATA_PATH": str(Path(tmp) / "pyvista"),
            "MPLBACKEND": "Agg", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1",
        }
        old = {key: os.environ.get(key) for key in settings}
        old_bytecode = sys.dont_write_bytecode
        try:
            os.environ.update(settings)
            sys.dont_write_bytecode = True
            yield
        finally:
            sys.dont_write_bytecode = old_bytecode
            for key, value in old.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value


def run_probe(name):
    started = time.monotonic()
    try:
        if name == "inventory":
            result = PROBES[name]()
        else:
            with _temporary_runtime():
                result = PROBES[name]()
    except Exception as exc:
        result = _result("failed", LEVELS[name], {"scope": "Requested probe did not complete; no readiness claim."})
        result["error"] = _safe_error(exc)
    result["probe"] = name
    result["elapsed_seconds"] = round(time.monotonic() - started, 3)
    if name != "inventory":
        module = sys.modules.get("spateo")
        result["details"]["spateo_import_path"] = getattr(module, "__file__", None)
        result["details"]["spateo_version"] = getattr(module, "__version__", None)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", required=True, choices=tuple(PROBES))
    parser.add_argument("--result-file", required=True, type=Path, help="New JSON output file; parent directory must exist.")
    args = parser.parse_args(argv)
    # Reserve exclusively before doing work: existing reports/data cannot be
    # overwritten, and an unwritable destination fails before expensive imports.
    try:
        with args.result_file.open("x", encoding="utf-8") as output:
            with contextlib.redirect_stdout(sys.stderr):
                result = run_probe(args.probe)
            json.dump(result, output, indent=2, sort_keys=True, allow_nan=False)
            output.write("\n")
    except OSError as exc:
        parser.exit(2, _safe_error(exc) + "\n")
    return 1 if result["status"] == "failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
