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
SCHEMA = "spateo-4d/v3"
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
        "MaxIter": 500,
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
        "selected": ["acceleration", "curl", "divergence", "torsion", "curvature"],
        "glm_metrics": [],
        "glm_genes": [],
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
        "max_points": 6000,
        "max_target_points": 4000,
        "max_vectors": 1200,
        "default_feature": "displacement",
        "cdn": False,
    },
    "runtime": {"device": "cpu", "seed": 0},
}


def read_json(path):
    return json.loads(Path(path).read_text())


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
                else {"stage1", "stage2", "coordinate_unit"}
            )
            if set(value) - allowed:
                raise ValueError(
                    f"Unknown {key} settings: {sorted(set(value) - allowed)}"
                )
            config[key].update(value)
        else:
            config[key] = value
    for key in ("stage1", "stage2"):
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
    supported = {"acceleration", "curl", "divergence", "torsion", "curvature"}
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
        "morphofield": ["M", "restart_num", "MaxIter"],
        "gp": ["training_iter", "inducing_num"],
        "dashboard": ["max_points", "max_target_points", "max_vectors"],
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
    return config


def fingerprints(config):
    return {
        k: {"path": config["inputs"][k], "sha256": digest(config["inputs"][k])}
        for k in ("stage1", "stage2")
    }


def implementation():
    import spateo as st
    import spateo.logging

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
            "labels": "alignment",
            "subset": "mapping",
            "runtime": "alignment",
            "schema_version": "alignment",
        }
        invalid = descendants(owners.get(s, s) for s in changed)
        if input_hashes != parent.get("inputs") or runtime != parent.get(
            "implementation"
        ):
            invalid = set(STAGES)
        for s in STAGES:
            old = parent.get("stages", {}).get(s, {})
            enabled = config.get(s, {}).get("enabled", True)
            if old.get("status") == "skipped" and not enabled:
                continue
            if old.get("status") not in ("completed", "reused") or not valid_outputs(
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


def alignment(config, directory):
    import anndata as ad
    import numpy as np
    import spateo as st

    a, b = (ad.read_h5ad(config["inputs"][k]) for k in ("stage1", "stage2"))
    p = config["alignment"]
    for data in (a, b):
        validate_input(data, config)
        st.pp.normalize_total(
            data,
            layer=p["counts_layer"],
            out_layer="normalized",
            target_sum=p["target_sum"],
            inplace=True,
        )
        st.pp.log1p_layer(
            data,
            layer="normalized",
            out_layer=p["log_layer"],
            set_X=False,
            inplace=True,
        )
    common = [g for g in a.var_names if g in b.var_names]
    if not common:
        raise ValueError("No shared genes for cross-stage alignment")
    aligned, _, _, _ = st.align.morpho_align_ref(
        models=[a, b],
        rep_layer=p["log_layer"],
        rep_field="layer",
        spatial_key=p["spatial_key"],
        key_added=p["aligned_key"],
        genes=common,
        mode=p["mode"],
        n_sampling=min(p["n_sampling"], a.n_obs, b.n_obs),
        sampling_method=p["sampling_method"],
        max_iter=p["max_iter"],
        nonrigid_start_iter=p["nonrigid_start_iter"],
        device=config["runtime"]["device"],
        verbose=False,
    )
    qc = {}
    for name, original, data in zip(("stage1", "stage2"), (a, b), aligned):
        np.testing.assert_array_equal(
            data.obsm[p["spatial_key"]], original.obsm[p["spatial_key"]]
        )
        coords = np.asarray(data.obsm[p["aligned_key"]])
        if not np.isfinite(coords).all():
            raise ValueError("Alignment produced nonfinite coordinates")
        qc[name] = {
            "n_cells": data.n_obs,
            "min": coords.min(0).tolist(),
            "max": coords.max(0).tolist(),
        }
    write_json(directory / "qc.json", qc)
    return {**save_pair(*aligned, directory), "qc": directory / "qc.json"}


def mapping(config, a, b, directory):
    import numpy as np
    import pandas as pd
    import spateo as st

    s, p, al = config["subset"], config["mapping"], config["alignment"]
    if s["group"] is not None:
        pair = []
        for data in (a, b):
            labels = data.obs[s["annotation_key"]]
            if labels.isna().any():
                raise ValueError("Missing subset annotation")
            part = data[labels.astype(str) == str(s["group"])].copy()
            if not part.n_obs:
                raise ValueError("Requested subset has no cells in one stage")
            pair.append(part)
        a, b = pair
    common = [g for g in a.var_names if g in b.var_names]
    totals_a = dict(
        zip(a.var_names, np.asarray(a.layers[al["counts_layer"]].sum(0)).ravel())
    )
    totals_b = dict(
        zip(b.var_names, np.asarray(b.layers[al["counts_layer"]].sum(0)).ravel())
    )
    common = [g for g in common if totals_a[g] > 0 and totals_b[g] > 0]
    if not common:
        raise ValueError("No expressed shared genes in selected subsets")
    a, b = a[:, common].copy(), b[:, common].copy()
    if a.n_obs * b.n_obs > p["max_pairs"]:
        raise ValueError(
            "Mapping exceeds mapping.max_pairs; choose an explicit subset or increase the reviewed memory budget"
        )
    for data in (a, b):
        if np.any(np.asarray(data.layers[al["counts_layer"]].sum(1)).ravel() <= 0):
            raise ValueError("Gene harmonization left zero-library cells")
        st.pp.normalize_total(
            data,
            layer=al["counts_layer"],
            out_layer="normalized",
            target_sum=p["target_sum"],
            inplace=True,
        )
        st.pp.log1p_layer(
            data,
            layer="normalized",
            out_layer=al["log_layer"],
            set_X=True,
            inplace=True,
        )
        if s["annotation_key"] not in data.obs:
            data.obs[s["annotation_key"]] = "all cells (display group)"
    _, pi = st.tdr.cell_directions(
        adataA=a,
        adataB=b,
        layer=al["log_layer"],
        spatial_key=al["aligned_key"],
        key_added=p["key"],
        alpha=p["alpha"],
        numItermax=p["numItermax"],
        numItermaxEmd=p["numItermaxEmd"],
        device=config["runtime"]["device"],
        inplace=True,
    )
    x, v = a.obsm["X_" + p["key"]], a.obsm["V_" + p["key"]]
    if not np.isfinite(v).all() or not np.isfinite(pi).all():
        raise ValueError("Mapping produced nonfinite results")
    np.testing.assert_allclose(x - a.obsm[al["aligned_key"]], v)
    np.savez_compressed(
        directory / "transport.npz",
        pi=pi,
        source_ids=a.obs_names.to_numpy(dtype=str),
        target_ids=b.obs_names.to_numpy(dtype=str),
    )
    table = pd.DataFrame(
        {
            "cell_id": a.obs_names,
            "group": a.obs[s["annotation_key"]].astype(str).to_numpy(),
            "displacement": np.linalg.norm(v, axis=1),
        }
    )
    table.to_csv(directory / "mapping.csv", index=False)
    summary = table.groupby("group")["displacement"].agg(
        n="size",
        mean="mean",
        median="median",
        p05=lambda x: x.quantile(0.05),
        p95=lambda x: x.quantile(0.95),
    )
    summary.to_csv(directory / "mapping_summary.csv")
    return {
        **save_pair(a, b, directory),
        "transport": directory / "transport.npz",
        "mapping_table": directory / "mapping.csv",
        "mapping_summary": directory / "mapping_summary.csv",
    }


def morphofield(config, a, b, directory):
    import numpy as np
    import spateo as st

    p = config["morphofield"].copy()
    key = p.pop("key")
    p["M"] = min(p["M"], a.n_obs)
    xyz = a.obsm[config["alignment"]["aligned_key"]]
    st.tdr.morphofield_sparsevfc(
        a,
        spatial_key=config["alignment"]["aligned_key"],
        V_key="V_" + config["mapping"]["key"],
        key_added=key,
        NX=xyz.copy(),
        **p,
    )
    if not np.isfinite(a.uns[key]["V"]).all():
        raise ValueError("Vector-field fit produced nonfinite values")
    return save_pair(a, b, directory)


def trajectory(config, a, b, directory):
    import numpy as np
    import spateo as st

    p = config["trajectory"]
    st.tdr.morphopath(
        a,
        vf_key=config["morphofield"]["key"],
        key_added=p["key"],
        t_end=p["t_end"],
        interpolation_num=p["interpolation_num"],
        direction=p["direction"],
        cores=1,
    )
    for values in a.uns[p["key"]]["prediction"].values():
        if not np.isfinite(values).all():
            raise ValueError("Trajectory contains nonfinite coordinates")
    return save_pair(a, b, directory)


def metrics(config, a, b, directory):
    import numpy as np
    import spateo as st

    p = config["metrics"]
    extra = {}
    st.tdr.morphofield_velocity(
        a, vf_key=config["morphofield"]["key"], key_added="velocity"
    )
    for metric in p["selected"]:
        getattr(st.tdr, "morphofield_" + metric)(
            a, vf_key=config["morphofield"]["key"], key_added=metric
        )
        if not np.isfinite(a.obs[metric].to_numpy()).all():
            raise ValueError("Nonfinite metric: " + metric)
        if metric in p["glm_metrics"]:
            missing = set(p["glm_genes"]) - set(a.var_names)
            if missing:
                raise ValueError("Unknown GLM genes: " + str(sorted(missing)))
            key = "glm_degs_" + metric
            if a.obs[metric].nunique() < 4:
                raise ValueError(
                    "Spline GLM needs at least four distinct feature values"
                )
            st.tl.glm_degs(
                a,
                layer="normalized",
                genes=p["glm_genes"],
                key_added=key,
                fullModelFormulaStr=f"~cr({metric}, df=3)",
                qval_threshold=p["qval_threshold"],
                llf_threshold=p["llf_threshold"],
            )
            table = a.uns[key]["glm_result"]
            table.to_csv(directory / (key + ".csv"))
            extra[key] = directory / (key + ".csv")
    a.obs[p["selected"]].to_csv(directory / "metrics.csv", index_label="cell_id")
    return {
        **save_pair(a, b, directory),
        **extra,
        "metrics_table": directory / "metrics.csv",
    }


def gp(config, a, b, directory):
    import numpy as np
    import spateo as st

    p = config["gp"]
    missing = set(p["genes"]) - set(a.var_names)
    if missing:
        raise ValueError("Unknown GP genes: " + str(sorted(missing)))
    if set(p["genes"]) & set(a.obs.columns):
        raise ValueError("GP gene IDs conflict with observation field names")
    result = st.tdr.gp_interpolation(
        a,
        target_points=a.obsm[config["alignment"]["aligned_key"]].copy(),
        keys=p["genes"],
        spatial_key=config["alignment"]["aligned_key"],
        layer="normalized",
        training_iter=p["training_iter"],
        method=p["method"],
        inducing_num=min(p["inducing_num"], a.n_obs),
        device=config["runtime"]["device"],
        verbose=False,
    )
    if not np.isfinite(np.asarray(result.X)).all():
        raise ValueError("GP produced nonfinite expression")
    result.write_h5ad(directory / "gene_interpolation.h5ad")
    return {"gene_interpolation": directory / "gene_interpolation.h5ad"}


def dashboard(config, manifest, directory):
    import importlib.util
    from argparse import Namespace

    script = (
        Path(__file__).resolve().parents[1]
        / "subskills/spateo-render-dashboard/scripts/build_dashboard.py"
    )
    spec = importlib.util.spec_from_file_location("spateo_dashboard", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Metrics and trajectories are sibling checkpoints; display the metric checkpoint.
    source = (
        manifest["stages"]["metrics"]
        if manifest["stages"]["metrics"]["status"] != "skipped"
        else manifest["stages"]["morphofield"]
    )
    pair = source["outputs"]
    p = config["dashboard"]
    args = Namespace(
        adata=Path(pair["stage1_h5ad"]["path"]),
        target_adata=Path(pair["stage2_h5ad"]["path"]),
        manifest=None,
        spatial_key=config["alignment"]["aligned_key"],
        target_spatial_key=None,
        groupby=config["subset"]["annotation_key"],
        max_points=p["max_points"],
        max_target_points=p["max_target_points"],
        seed=config["runtime"]["seed"],
        features=["displacement", "velocity", *config["metrics"]["selected"]],
        vector_key="V_" + config["mapping"]["key"],
        mapped_key="X_" + config["mapping"]["key"],
        mapping_summary=Path(
            manifest["stages"]["mapping"]["outputs"]["mapping_summary"]["path"]
        ),
        metric_summary=None,
        glm_dir=None,
        max_vectors=p["max_vectors"],
        default_feature=p["default_feature"],
    )
    metric_outputs = manifest["stages"]["metrics"].get("outputs", {})
    if "metrics_table" in metric_outputs:
        args.glm_dir = Path(metric_outputs["metrics_table"]["path"]).parent
    payload = module.build_payload(args)
    payload.update(
        runId=manifest["run_id"],
        runStatus="completed",
        stages=[
            {
                "name": k,
                "status": ("completed" if k == "dashboard" else v["status"]),
                "reusedFrom": v.get("reused_from"),
            }
            for k, v in manifest["stages"].items()
        ],
    )
    output = directory / "index.html"
    output.write_text(module.render_html(payload, p["cdn"], None))
    return {"dashboard_html": output}


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
):
    config = canonical_config(config_path)
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
            if not config.get(stage, {}).get("enabled", True):
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
                        outputs = globals()[stage](config, a, b, directory)
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
        "dashboard": manifest["stages"]["dashboard"]["outputs"]
        .get("dashboard_html", {})
        .get("path"),
    }
