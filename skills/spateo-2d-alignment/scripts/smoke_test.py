#!/usr/bin/env python3
"""Temporary synthetic CPU checks for packaging and input contracts, without Spateo/GPU."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import shutil
import sys
import tempfile
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    sys.modules[name] = value
    spec.loader.exec_module(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    skill = Path(__file__).resolve().parents[1]
    pair = skill / "pipelines/pairwise-rigid"
    continuity = skill / "pipelines/continuity-guided"
    results = []
    with tempfile.TemporaryDirectory(prefix="spateo-2d-smoke-") as temporary:
        root = Path(temporary)
        inputs = root / "slices"
        inputs.mkdir()
        rng = np.random.default_rng(41)
        points = rng.normal(size=(80, 2))
        labels = np.asarray(["a"] * 40 + ["b"] * 40)
        onehot = np.eye(2, dtype=np.float32)[np.repeat([0, 1], 40)]
        for index in range(2):
            ids = [f"cell_{index}_{cell}" for cell in range(80)]
            model = ad.AnnData(X=sparse.csr_matrix(onehot), obs=pd.DataFrame({"anno": labels, "physical_z": index * 20.0}, index=ids))
            model.obsm["spatial"] = (points + np.asarray([3 * index, -index])).astype(np.float32)
            model.obsm["X_pca"] = onehot
            model.write_h5ad(inputs / f"CS01_SL{index + 1:03d}_YTEST.Spatial.h5ad")
        validator = module("packaging_validator", pair / "validate_inputs.py")
        summary, _ = validator.validate_inputs(inputs, "annotation-onehot")
        assert summary["n_cells"] == 160 and summary["global_ids_unique"] and summary["physical_z_order"] == "strictly_increasing"
        results.append("shared annotation one-hot and globally unique input IDs")
        try:
            validator.validate_inputs(inputs, "expression-pca")
        except ValueError:
            results.append("one-hot correctly rejected as expression PCA")
        else:
            raise AssertionError("one-hot was accepted as expression PCA")
        nested = inputs / "nested_archive"
        nested.mkdir()
        nested_file = nested / "CS01_SL003_YTEST.Spatial.h5ad"
        shutil.copyfile(next(inputs.glob("*.h5ad")), nested_file)
        try:
            validator.validate_inputs(inputs, "annotation-onehot")
        except ValueError as exc:
            assert "flat directory" in str(exc)
        else:
            raise AssertionError("Nested unreviewed H5AD was accepted")
        nested_file.unlink()
        nested.rmdir()
        blank = root / "blank-labels"
        blank.mkdir()
        for path in inputs.glob("*.h5ad"):
            obj = ad.read_h5ad(path)
            obj.obs["anno"] = " "
            obj.obsm["X_pca"] = np.tile([1, 0], (obj.n_obs, 1)).astype(np.float32)
            obj.write_h5ad(blank / path.name)
        try:
            validator.validate_inputs(blank, "annotation-onehot")
        except ValueError as exc:
            assert "nonblank" in str(exc)
        else:
            raise AssertionError("Blank annotation labels were accepted")
        results.append("nested unreviewed H5AD and blank annotation values are rejected")
        os.environ.update(DATASET_ROOT=str(inputs), OUTDIR=str(root / "not-created"), EXPRESSION_MODE="spatial_only")
        normal = module("pairwise_normal", pair / "runners/spateo_rigid_sns_alignment.py")
        spatial = module("pairwise_spatial", pair / "runners/spateo_rigid_sns_alignment_spatial_only.py")
        order = normal.collect_global_order()
        loaded = normal.load_slice(order.iloc[0], rng)
        constant = spatial.load_slice(order.iloc[0], rng)
        assert np.array_equal(loaded.obsm["X_pca"], onehot)
        assert constant.obsm["X_pca"].shape == (80, 30) and np.all(constant.obsm["X_pca"] == 1)
        assert loaded.obs_names[0] == "CS01_SL001_YTEST:0"
        results.append("unchanged pairwise loaders preserve features and spatial-only ones30 semantics")
        sys.path.insert(0, str(continuity))
        core = module("_core", continuity / "_core.py")
        scoring = module("_continuity_scoring", continuity / "_continuity_scoring.py")
        engine = module("engine", continuity / "engine.py")
        args = SimpleNamespace(slice_dir=inputs, annotation_key="anno", spatial_key="spatial", representation_key="X_pca", sample_cap=60, seed=42)
        slices = core.load_slices(args)
        models = engine.load_models(slices, args)
        assert len(models) == 2 and sum(model.n_obs for model in models) == 160
        expected_ascii_hash = core.sha256_arrays(slices[0].cell_ids.astype("U").astype("S"),
                            slices[0].annotations.astype("U").astype("S"), slices[0].xy.astype("<f8"))
        assert slices[0].content_sha256 == expected_ascii_hash
        unicode_dir = root / "unicode"
        unicode_dir.mkdir()
        for index, path in enumerate(sorted(inputs.glob("*.h5ad"))):
            obj = ad.read_h5ad(path)
            obj.obs_names = [f"细胞_{index}_{cell}" for cell in range(obj.n_obs)]
            obj.obs["anno"] = ["类型甲"] * 40 + ["类型乙"] * 40
            obj.write_h5ad(unicode_dir / path.name)
        validator.validate_inputs(unicode_dir, "annotation-onehot", filename_style="continuity")
        unicode_args = SimpleNamespace(**{**vars(args), "slice_dir": unicode_dir})
        unicode_slices = core.load_slices(unicode_args)
        unicode_models = engine.load_models(unicode_slices, unicode_args)
        assert unicode_slices[0].cell_ids[0] == "细胞_0_0" and unicode_models[0].obs["anno"].iloc[0] == "类型甲"
        results.append("UTF-8 IDs/labels survive continuity hashing; ASCII content hashes remain unchanged")
        transform = core.rigid_matrix(33, np.asarray([4, -7]))
        target = core.apply_matrix(points, transform)
        fitted = core.fit_proper_rigid(points, target)
        np.testing.assert_allclose(core.apply_matrix(points, fitted), target, atol=1e-10)
        assert np.linalg.det(fitted[:2, :2]) > 0
        baseline = scoring.score_alignment_v2(points, labels, points, labels)
        displaced = scoring.score_alignment_v2(points + 20, labels, points, labels)
        assert baseline["total"] < displaced["total"]
        results.append("continuity dependency closure loads, proper rigid fit recovers a synthetic transform, score responds to displacement")
        for script in [pair / "run.py", continuity / "run.py", skill / "scripts/check_runtime.py"]:
            subprocess.run([sys.executable, str(script), "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        dry = subprocess.run([sys.executable, str(pair / "run.py"), "--slice-dir", str(inputs),
                              "--output-dir", str(root / "dry-run"), "--representation", "annotation-onehot", "--dry-run"],
                             text=True, stdout=subprocess.PIPE, check=True)
        assert json.loads(dry.stdout)["algorithm_executed"] is False and not (root / "dry-run").exists()
        results.append("both public CLIs and runtime-check help work; dry-run has no output-directory side effect")
        cli = module("pairwise_public_cli", pair / "run.py")
        os.environ["MAX_ITER"] = "999"
        os.environ["STAGE1_MODE"] = "SN-N"
        os.environ["INITIAL_COORDINATE_X_COL"] = "hidden_coordinate"
        cleaned, known_keys, removed = cli.clean_environment(pair / "runners/spateo_rigid_sns_alignment.py", {"STAGE1_MODE": "SN-S"})
        assert "MAX_ITER" not in cleaned and "INITIAL_COORDINATE_X_COL" not in cleaned and cleaned["STAGE1_MODE"] == "SN-S"
        assert {"MAX_ITER", "STAGE1_MODE", "INITIAL_COORDINATE_X_COL"} <= set(removed)
        results.append("ambient recognized runner variables are removed before explicit CLI settings")
    print(json.dumps({"status": "passed", "checks": results, "spateo_or_gpu_alignment_executed": False,
                      "accuracy_benchmark_repeated": False}, indent=2))


if __name__ == "__main__":
    main()
