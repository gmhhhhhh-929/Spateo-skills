"""Native Spateo 4D execution and immutable, content-verified checkpoints."""

from __future__ import annotations

import contextlib
import hashlib
import importlib.metadata
import json
import os
import re
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

SOURCE_COMMIT = "615644f88613bea8ceb2e2df1e2391d16de55ec1"
SCHEMA = "spateo-4d/v4"
STAGES = (
    "alignment",
    "mapping",
    "morphofield",
    "trajectory",
    "metrics",
    "gp",
    "dashboard",
)
DEPENDENCIES = {
    "alignment": (),
    "mapping": ("alignment",),
    "morphofield": ("mapping",),
    "trajectory": ("morphofield",),
    "metrics": ("morphofield",),
    "gp": ("mapping",),
    "dashboard": ("mapping", "morphofield", "trajectory", "metrics", "gp"),
}
DEFAULTS = {
    "schema_version": SCHEMA,
    "inputs": {},
    "workflow": {"entry": "alignment", "until": "features", "frame_id": None, "species": None},
    "labels": {"stage1": "source", "stage2": "target"},
    "alignment": {
        "spatial_key": "spatial_3d",
        "counts_layer": "counts",
        "x_is_counts": False,
        "log_layer": "log1p",
        "aligned_key": "aligned",
        "target_sum": 10000.0,
        "mode": "SN-S",
        "n_sampling": 2000,
        "sampling_method": "random",
        "max_iter": 200,
        "nonrigid_start_iter": 80,
    },
    "subset": {"annotation_key": "anno", "group": None},
    "mapping": {
        "key": "cells_mapping",
        "alpha": 0.001,
        "numItermax": 200,
        "numItermaxEmd": 100000,
        "target_sum": 10000.0,
        "max_pairs": 25000000,
    },
    "morphofield": {
        "key": "VecFld_morpho",
        "M": 100,
        "lambda_": 0.02,
        "restart_num": 1,
        "restart_seed": [0],
        "max_iter": 8,
        "beta": None,
        "tol": 1e-5,
    },
    "trajectory": {
        "enabled": True,
        "key": "fate_morpho",
        "t_end": 1.0,
        "interpolation_num": 50,
        "direction": "forward",
    },
    "metrics": {
        "enabled": True,
        "selected": ["speed", "acceleration", "curl", "divergence", "torsion", "curvature", "jacobian_frobenius"],
        "glm_metrics": [],
        "glm_genes": [],
        "glm_min_cells": 10,
        "glm_top_plots": 6,
        "qval_threshold": 0.05,
        "llf_threshold": None,
    },
    "gp": {
        "enabled": False,
        "genes": [],
        "training_iter": 50,
        "method": "SVGP",
        "inducing_num": 512,
    },
    "dashboard": {
        "enabled": True,
        "max_points": 0,
        "max_target_points": 0,
        "max_vectors": 600,
        "max_trajectories": 200,
        "write_html": True,
        "default_feature": "displacement",
        "cdn": False,
    },
    "runtime": {"device": "cpu", "seed": 0},
}


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name("." + path.name + "." + uuid.uuid4().hex)
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def canonical_config(path):
    path = Path(path).resolve()
    supplied = read_json(path)
    if supplied.get("schema_version") != SCHEMA:
        raise ValueError(
            f"Use schema_version={SCHEMA}; migrate old notebook configs explicitly."
        )
    if set(supplied) - set(DEFAULTS):
        raise ValueError(
            "Unknown config sections: " + str(sorted(set(supplied) - set(DEFAULTS)))
        )
    config = json.loads(json.dumps(DEFAULTS))
    for key, value in supplied.items():
        if isinstance(config[key], dict):
            if not isinstance(value, dict):
                raise ValueError(f"{key} must be an object")
            allowed = (
                set(config[key])
                if key != "inputs"
                else {"stage1", "stage2", "coordinate_unit", "pointcloud1", "pointcloud2"}
            )
            if set(value) - allowed:
                raise ValueError(
                    f"Unknown {key} settings: {sorted(set(value) - allowed)}"
                )
            config[key].update(value)
        else:
            config[key] = value
    for key in ("stage1", "stage2", "pointcloud1", "pointcloud2"):
        if key not in config["inputs"]:
            if key.startswith("pointcloud"): continue
            raise ValueError("Missing H5AD input: " + key)
        p = Path(config["inputs"][key]).expanduser()
        p = (path.parent / p).resolve() if not p.is_absolute() else p.resolve()
        if not p.is_file():
            raise FileNotFoundError(p)
        config["inputs"][key] = str(p)
    if not config["inputs"].get("coordinate_unit"):
        raise ValueError("Declare inputs.coordinate_unit for the shared XYZ frame")
    a = config["alignment"]
    if a["mode"] not in ("SN-S", "SN-N") or a["sampling_method"] not in (
        "random",
        "trn",
        "kmeans",
    ):
        raise ValueError("Unsupported alignment mode or sampling method")
    if a["nonrigid_start_iter"] < 0 or a["max_iter"] <= a["nonrigid_start_iter"] + 1:
        raise ValueError(
            "max_iter must exceed nonrigid_start_iter + 1 so reference deformation coefficients are initialized"
        )
    if a["target_sum"] is not None and a["target_sum"] <= 0:
        raise ValueError("alignment.target_sum must be positive or null")
    if a["spatial_key"] == a["aligned_key"]:
        raise ValueError("aligned_key must preserve the original spatial_key")
    if a["counts_layer"] == a["log_layer"] or a["counts_layer"] == "normalized":
        raise ValueError("Counts and derived layers must have different keys")
    m = config["mapping"]
    if (
        not 0 <= m["alpha"] <= 1
        or (m["target_sum"] is not None and m["target_sum"] <= 0)
        or m["max_pairs"] < 1
    ):
        raise ValueError("Invalid mapping alpha, target_sum or max_pairs")
    t = config["trajectory"]
    if (
        t["t_end"] <= 0
        or t["interpolation_num"] < 2
        or t["direction"] not in ("forward", "backward", "both")
    ):
        raise ValueError("Invalid trajectory duration, samples or direction")
    metrics = config["metrics"]
    supported = {"speed", "acceleration", "curl", "divergence", "torsion", "curvature", "jacobian_frobenius"} | {"jacobian_"+i+j for i in "xyz" for j in "xyz"} | {"velocity_"+i for i in "xyz"}
    if not set(metrics["selected"]) <= supported or not set(
        metrics["glm_metrics"]
    ) <= set(metrics["selected"]):
        raise ValueError("GLM features must be selected supported scalar metrics")
    if metrics["glm_metrics"] and not metrics["glm_genes"]:
        raise ValueError("Provide verified metrics.glm_genes when requesting GLMs")
    if config["gp"]["enabled"] and not config["gp"]["genes"]:
        raise ValueError("Provide verified gp.genes when enabling GP")
    for section, keys in {
        "alignment": ["n_sampling", "max_iter"],
        "morphofield": ["M", "restart_num", "max_iter"],
        "gp": ["training_iter", "inducing_num"],
        "dashboard": ["max_vectors", "max_trajectories"],
    }.items():
        if any(
            not isinstance(config[section][k], int) or config[section][k] < 1
            for k in keys
        ):
            raise ValueError(
                f"{section} sizes and iteration counts must be positive integers"
            )
    if (
        len(config["morphofield"]["restart_seed"])
        != config["morphofield"]["restart_num"]
    ):
        raise ValueError("restart_seed length must match restart_num")
    if config["dashboard"]["cdn"]:
        raise ValueError("The modular viewer is offline: dashboard.cdn must be false")
    w = config["workflow"]
    if w["entry"] not in ("alignment", "aligned", "field") or w["until"] not in ("alignment", "trajectory", "features"):
        raise ValueError("Unsupported workflow entry/until")
    if w["entry"] != "alignment" and not w["frame_id"]:
        raise ValueError("Imported aligned/field data require workflow.frame_id")
    if w["entry"] == "field" and w["until"] == "alignment":
        raise ValueError("Field entry cannot end at alignment")
    if any(config["dashboard"][k] < 0 for k in ("max_points", "max_target_points")):
        raise ValueError("Point display caps must be >=0 (0 means all)")
    return config


def fingerprints(config):
    return {
        k: {"path": config["inputs"][k], "sha256": digest(config["inputs"][k])}
        for k in ("stage1", "stage2", "pointcloud1", "pointcloud2") if k in config["inputs"]
    }


def implementation():
    import spateo as st
    import spateo.logging
    import inspect
    from spateo._native import sparse_vector_field
    if "max_iter" not in inspect.signature(sparse_vector_field).parameters:
        raise RuntimeError("Use the audited native Spateo implementation; backend signature changed")

    package = Path(st.__file__).resolve().parent
    h = hashlib.sha256()
    for p in sorted(package.rglob("*.py")):
        h.update(str(p.relative_to(package)).encode())
        h.update(p.read_bytes())
    scripts = Path(__file__).resolve().parents[1]
    sh = hashlib.sha256()
    for p in sorted(scripts.rglob("*.py")):
        sh.update(str(p.relative_to(scripts)).encode())
        sh.update(p.read_bytes())
    versions = {}
    for name in (
        "numpy",
        "scipy",
        "anndata",
        "pandas",
        "POT",
        "torch",
        "gpytorch",
        "statsmodels",
        "plotly",
    ):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    result = {
        "spateo_version": str(st.__version__),
        "spateo_python_sha256": h.hexdigest(),
        "skill_python_sha256": sh.hexdigest(),
        "viewer_asset_sha256": digest(scripts/"assets/viewer.html"),
        "python": sys.version.split()[0],
        "packages": versions,
    }
    # File hashes remain authoritative for dirty or non-Git installations.
    if (package.parent / ".git").exists():
        result["source_commit"] = subprocess.check_output(
            ["git", "-C", str(package.parent), "rev-parse", "HEAD"], text=True
        ).strip()
    return result


def descendants(names):
    invalid = set(names)
    for stage in STAGES:
        if set(DEPENDENCIES[stage]) & invalid:
            invalid.add(stage)
    return invalid


def valid_outputs(stage):
    outputs = stage.get("outputs", {})
    return bool(outputs) and all(
        Path(v["path"]).is_file() and digest(v["path"]) == v["sha256"]
        for v in outputs.values()
    )


def plan(config, parent=None, *, input_hashes=None, runtime=None):
    changed = []
    if not parent or parent.get("schema_version") != SCHEMA:
        invalid = set(STAGES)
    else:
        for section in DEFAULTS:
            if config.get(section) != parent["config"].get(section):
                changed.append(section)
        owners = {
            "inputs": "alignment",
            "workflow": "alignment",
            "labels": "alignment",
            "subset": "mapping",
            "runtime": "alignment",
            "schema_version": "alignment",
        }
        invalid = descendants(owners.get(s, s) for s in changed)
        current_runtime={k:v for k,v in (runtime or {}).items() if k != "viewer_asset_sha256"}
        previous_runtime={k:v for k,v in parent.get("implementation",{}).items() if k != "viewer_asset_sha256"}
        if input_hashes != parent.get("inputs") or current_runtime != previous_runtime:
            invalid = set(STAGES)
        elif (runtime or {}).get("viewer_asset_sha256") != parent.get("implementation",{}).get("viewer_asset_sha256"):
            invalid.add("dashboard")
        for s in STAGES:
            old = parent.get("stages", {}).get(s, {})
            enabled = config.get(s, {}).get("enabled", True)
            if old.get("status") == "skipped" and not enabled:
                continue
            if old.get("status") not in ("completed", "reused", "imported") or not valid_outputs(
                old
            ):
                invalid.update(descendants([s]))
    return {
        "changed_sections": changed,
        "recompute_stages": [s for s in STAGES if s in invalid],
        "reuse_stages": [s for s in STAGES if s not in invalid],
    }


def serializable(value):
    import numpy as np

    if isinstance(value, np.ndarray) and value.dtype.kind == 'O':
        return {"__spateo_skill_type__": "object_array", "shape": np.asarray(value.shape),
                "items": {str(i): serializable(v) for i,v in enumerate(value.ravel())}}
    if isinstance(value, dict):
        return {str(k): serializable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return {
            "__spateo_skill_type__": "tuple" if isinstance(value, tuple) else "list",
            "items": {str(i): serializable(v) for i, v in enumerate(value)},
        }
    if value is None:
        return {"__spateo_skill_type__": "none"}
    if isinstance(value, np.generic):
        return value.item()
    return value


def restore(value):
    if isinstance(value, dict):
        kind = value.get("__spateo_skill_type__")
        if kind == "object_array":
            import numpy as np
            values=[restore(value['items'][str(i)]) for i in range(len(value['items']))]
            return np.asarray(values,dtype=object).reshape(tuple(value['shape']))
        if kind == "none":
            return None
        if kind in ("list", "tuple"):
            items = [
                restore(value["items"][str(i)]) for i in range(len(value["items"]))
            ]
            return tuple(items) if kind == "tuple" else items
        return {k: restore(v) for k, v in value.items()}
    return value


def save_pair(a, b, directory):
    import anndata as ad
    import numpy as np

    outputs = {}
    for label, data in [("stage1", a), ("stage2", b)]:
        path = directory / (label + ".h5ad")
        data.uns = serializable(dict(data.uns))
        data.write_h5ad(path, compression="gzip")
        check = ad.read_h5ad(path)
        assert check.obs_names.equals(data.obs_names) and check.var_names.equals(
            data.var_names
        )
        for key in data.obsm:
            np.testing.assert_allclose(check.obsm[key], data.obsm[key], equal_nan=True)
        outputs[label + "_h5ad"] = path
    return outputs


def load_pair(outputs):
    import anndata as ad

    pair = tuple(
        ad.read_h5ad(outputs[k + "_h5ad"]["path"]) for k in ("stage1", "stage2")
    )
    for data in pair:
        data.uns = restore(dict(data.uns))
    return pair


def validate_input(data, config):
    import numpy as np
    from scipy import sparse

    a = config["alignment"]
    if (
        not data.n_obs
        or not data.n_vars
        or not data.obs_names.is_unique
        or not data.var_names.is_unique
    ):
        raise ValueError("Inputs must be nonempty with unique cell and gene IDs")
    xyz = np.asarray(data.obsm[a["spatial_key"]])
    if xyz.shape != (data.n_obs, 3) or not np.isfinite(xyz).all():
        raise ValueError("Expected finite n_cells x 3 coordinates")
    if a["counts_layer"] not in data.layers:
        if not a["x_is_counts"]:
            raise ValueError(
                "Missing counts layer; set x_is_counts only after verifying raw X"
            )
        data.layers[a["counts_layer"]] = data.X.copy()
    matrix = data.layers[a["counts_layer"]]
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.allclose(values, np.rint(values), atol=1e-6, rtol=0)
    ):
        raise ValueError(
            "Counts must be finite nonnegative integer-like captured counts"
        )
    if np.any(np.asarray(matrix.sum(axis=1)).ravel() <= 0):
        raise ValueError(
            "Zero-library cells require an explicit upstream filtering decision"
        )
    if np.linalg.matrix_rank(xyz - xyz.mean(axis=0)) < 3:
        raise ValueError("4D input must have three-dimensional spatial support")














def dashboard(config, manifest, directory):
    from viewer_payload import build_payload, render_html
    payload = build_payload(config, manifest)
    path = directory / "viewer_payload.json"
    write_json(path, payload)
    outputs = {"viewer_payload": path}
    if config["dashboard"]["write_html"]:
        html = directory / "index.html"
        html.write_text(render_html(payload), encoding='utf-8')
        outputs["dashboard_html"] = html
    return outputs


@contextlib.contextmanager
def capture_log(path):
    import logging

    previous_stderr = sys.stderr

    def handlers():
        loggers = [
            logging.getLogger(),
            *[
                v
                for v in logging.Logger.manager.loggerDict.values()
                if isinstance(v, logging.Logger)
            ],
        ]
        return {
            h
            for logger in loggers
            for h in logger.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, logging.FileHandler)
        }

    with path.open("w") as stream:
        original = {h: h.stream for h in handlers()}
        for h in original:
            h.stream = stream
        try:
            with contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
                yield
        finally:
            for h in handlers():
                if h.stream is stream:
                    h.stream = original.get(h, previous_stderr)


def execute(
    config_path,
    project,
    parent_path=None,
    run_id=None,
    dry_run=False,
    stop_after="dashboard",
    until=None,
):
    from alignment_stage import alignment
    from morphogenesis_stage import mapping, morphofield, trajectory
    from features_stage import metrics, gp
    functions = dict(alignment=alignment, mapping=mapping, morphofield=morphofield,
                     trajectory=trajectory, metrics=metrics, gp=gp)
    config = canonical_config(config_path)
    if until is not None: config["workflow"]["until"] = until
    if config["workflow"]["entry"] == "field" and config["workflow"]["until"] == "alignment":
        raise ValueError("Field entry cannot end at alignment")
    hashes, runtime = fingerprints(config), implementation()
    parent = read_json(parent_path) if parent_path else None
    change = plan(config, parent, input_hashes=hashes, runtime=runtime)
    if dry_run:
        return {
            "status": "dry-run",
            "config": config,
            "inputs": hashes,
            "implementation": runtime,
            **change,
        }
    run_id = (
        run_id
        or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ-") + uuid.uuid4().hex[:8]
    )
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", run_id):
        raise ValueError("run_id must be a safe directory basename")
    run = Path(project).expanduser().resolve() / "runs" / run_id
    run.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema_version": SCHEMA,
        "run_id": run_id,
        "run_dir": str(run),
        "config": config,
        "status": "running",
        "inputs": hashes,
        "implementation": runtime,
        "change_plan": change,
        "parent_manifest": str(Path(parent_path).resolve()) if parent_path else None,
        "parent_run_id": parent.get("run_id") if parent else None,
        "stages": {s: {"status": "pending", "outputs": {}} for s in STAGES},
    }
    write_json(run / "config.json", config)
    write_json(run / "manifest.json", manifest)
    import numpy as np

    np.random.seed(config["runtime"]["seed"])
    try:
        import torch

        torch.manual_seed(config["runtime"]["seed"])
    except ImportError:
        pass
    current = None
    try:
        for stage in STAGES:
            current = stage
            record = manifest["stages"][stage]
            until = config["workflow"]["until"]
            excluded = ((until == "alignment" and stage in ("mapping", "morphofield", "trajectory", "metrics", "gp"))
                        or (until == "trajectory" and stage in ("metrics", "gp")))
            if config["workflow"]["entry"] == "field" and stage in ("mapping", "morphofield"):
                record.update(status="imported", outputs=manifest["stages"]["alignment"]["outputs"])
            elif excluded or not config.get(stage, {}).get("enabled", True):
                record["status"] = "skipped"
            elif stage in change["reuse_stages"]:
                record.update(
                    status="reused",
                    reused_from=parent["run_id"],
                    outputs=parent["stages"][stage]["outputs"],
                )
            else:
                record["status"] = "running"
                write_json(run / "manifest.json", manifest)
                directory = run / "outputs" / stage
                directory.mkdir(parents=True)
                log = run / "logs" / (stage + ".log")
                log.parent.mkdir(exist_ok=True)
                with capture_log(log):
                    if stage == "alignment":
                        outputs = alignment(config, directory)
                    elif stage == "dashboard":
                        outputs = dashboard(config, manifest, directory)
                    else:
                        dependency = DEPENDENCIES[stage][0]
                        a, b = load_pair(manifest["stages"][dependency]["outputs"])
                        outputs = functions[stage](config, a, b, directory)
                record.update(
                    status="completed",
                    outputs={
                        k: {"path": str(v), "sha256": digest(v)}
                        for k, v in outputs.items()
                    },
                )
            write_json(run / "manifest.json", manifest)
            if stage == stop_after:
                break
        manifest["status"] = "completed" if stop_after == "dashboard" else "partial"
    except Exception as exc:
        import traceback

        error = run / "logs" / ((current or "runtime") + ".error.log")
        error.parent.mkdir(exist_ok=True)
        error.write_text(traceback.format_exc())
        manifest["status"] = "failed"
        manifest["stages"][current].update(
            status="failed", error=str(exc), error_log=str(error)
        )
        for record in manifest["stages"].values():
            if record["status"] == "pending":
                record["status"] = "blocked"
        write_json(run / "manifest.json", manifest)
        raise
    write_json(run / "manifest.json", manifest)
    return {
        "status": manifest["status"],
        "run_id": run_id,
        "manifest": str(run / "manifest.json"),
        "reused_stages": [
            s for s, v in manifest["stages"].items() if v["status"] == "reused"
        ],
        "viewer_payload": manifest["stages"]["dashboard"]["outputs"].get("viewer_payload", {}).get("path"),
        "dashboard": manifest["stages"]["dashboard"]["outputs"]
        .get("dashboard_html", {})
        .get("path"),
    }
