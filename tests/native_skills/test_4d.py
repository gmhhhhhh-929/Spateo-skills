import copy
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import pipeline_runtime as rt


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    root = tmp_path_factory.mktemp("native-4d")
    rng = np.random.default_rng(731)
    n = 40
    g = 12
    xyz = rng.normal(size=(n, 3)) * 10
    counts = rng.poisson(8, size=(n, g)) + 1
    for i in (1, 2):
        data = ad.AnnData(
            sparse.csr_matrix(
                counts if i == 1 else counts + rng.poisson(1, size=(n, g))
            ),
            obs=pd.DataFrame(
                {"anno": ["CNS"] * n}, index=[f"c{i}_{j}" for j in range(n)]
            ),
            var=pd.DataFrame(index=[f"g{j}" for j in range(g)]),
        )
        data.layers["counts"] = data.X.copy()
        data.obsm["spatial_3d"] = (
            xyz.copy()
            if i == 1
            else xyz * 1.1
            + np.column_stack([np.sin(xyz[:, 0]), np.cos(xyz[:, 1]), np.sin(xyz[:, 2])])
            + [1, 2, 3]
        )
        data.write_h5ad(root / f"stage{i}.h5ad")
    config = copy.deepcopy(rt.DEFAULTS)
    config["inputs"] = {
        "stage1": str(root / "stage1.h5ad"),
        "stage2": str(root / "stage2.h5ad"),
        "coordinate_unit": "um",
    }
    config["alignment"].update(n_sampling=30, max_iter=5, nonrigid_start_iter=1)
    config["subset"]["group"] = "CNS"
    config["morphofield"].update(M=15, MaxIter=20)
    config["trajectory"].update(t_end=0.1, interpolation_num=6)
    config["metrics"].update(
        glm_metrics=["acceleration"], glm_genes=["g0", "g1"], qval_threshold=1.0
    )
    config["gp"].update(enabled=True, genes=["g0"], training_iter=2, inducing_num=10)
    path = root / "config.json"
    rt.write_json(path, config)
    return root, path, config


@pytest.fixture(scope="module")
def completed(inputs):
    root, path, config = inputs
    result = rt.execute(path, root / "project", run_id="base")
    return result, rt.read_json(result["manifest"])


def test_complete_native_pipeline_and_persisted_outputs(inputs, completed):
    root, path, config = inputs
    result, m = completed
    assert result["status"] == "completed"
    assert all(v["status"] == "completed" for v in m["stages"].values())
    assert all(rt.valid_outputs(v) for v in m["stages"].values())
    a, b = rt.load_pair(m["stages"]["alignment"]["outputs"])
    original = ad.read_h5ad(root / "stage1.h5ad")
    np.testing.assert_array_equal(a.obsm["spatial_3d"], original.obsm["spatial_3d"])
    np.testing.assert_array_equal(
        a.layers["counts"].toarray(), original.layers["counts"].toarray()
    )
    mapped, _ = rt.load_pair(m["stages"]["mapping"]["outputs"])
    np.testing.assert_allclose(
        mapped.obsm["X_cells_mapping"] - mapped.obsm["aligned"],
        mapped.obsm["V_cells_mapping"],
    )
    transport = np.load(m["stages"]["mapping"]["outputs"]["transport"]["path"])
    assert transport["pi"].shape == (40, 40)
    assert transport["source_ids"].tolist() == mapped.obs_names.tolist()
    trajectories, _ = rt.load_pair(m["stages"]["trajectory"]["outputs"])
    assert len(trajectories.uns["fate_morpho"]["prediction"]) == 40
    geom, _ = rt.load_pair(m["stages"]["metrics"]["outputs"])
    assert np.isfinite(geom.obs[config["metrics"]["selected"]].to_numpy()).all()
    assert "glm_degs_acceleration" in m["stages"]["metrics"]["outputs"]
    glm = pd.read_csv(
        m["stages"]["metrics"]["outputs"]["glm_degs_acceleration"]["path"]
    )
    assert len(glm) == 2 and glm["status"].eq("ok").all()
    summary = pd.read_csv(m["stages"]["mapping"]["outputs"]["mapping_summary"]["path"])
    assert summary["n"].sum() == 40 and summary["group"].tolist() == ["CNS"]
    gp = ad.read_h5ad(m["stages"]["gp"]["outputs"]["gene_interpolation"]["path"])
    assert gp.shape == (40, 1) and gp.var_names.tolist() == ["g0"]
    html = Path(result["dashboard"]).read_text()
    assert (
        "cdn.plot.ly" not in html.split("</head>")[0]
        or '<script src="https://cdn.plot.ly' not in html
    )
    assert '"runStatus":"completed"' in html.replace(" ", "")
    assert '"status":"running"' not in html.replace(" ", "")


def test_display_only_child_reuses_verified_science(inputs, completed):
    root, path, config = inputs
    result, parent = completed
    frozen = Path(result["manifest"]).read_bytes()
    new = copy.deepcopy(config)
    new["dashboard"]["max_points"] = 20
    child = root / "display.json"
    rt.write_json(child, new)
    dry = rt.execute(child, root / "project", result["manifest"], dry_run=True)
    assert dry["recompute_stages"] == ["dashboard"]
    follow = rt.execute(child, root / "project", result["manifest"], run_id="display")
    assert follow["reused_stages"] == list(rt.STAGES[:-1])
    assert Path(result["manifest"]).read_bytes() == frozen


@pytest.mark.parametrize(
    "section,expected",
    [
        ("trajectory", ["trajectory", "dashboard"]),
        ("metrics", ["metrics", "dashboard"]),
        ("gp", ["gp", "dashboard"]),
        ("morphofield", ["morphofield", "trajectory", "metrics", "dashboard"]),
        ("mapping", list(rt.STAGES[1:])),
    ],
)
def test_dag_invalidation(inputs, completed, section, expected):
    root, path, config = inputs
    _, parent = completed
    config = copy.deepcopy(config)
    keys = {
        "trajectory": ("t_end", 0.2),
        "metrics": ("qval_threshold", 0.5),
        "gp": ("training_iter", 3),
        "morphofield": ("M", 12),
        "mapping": ("alpha", 0.1),
    }
    key, value = keys[section]
    config[section][key] = value
    plan = rt.plan(
        config, parent, input_hashes=parent["inputs"], runtime=parent["implementation"]
    )
    assert plan["recompute_stages"] == expected


def test_input_and_output_content_changes_invalidate(inputs, completed, tmp_path):
    _, _, config = inputs
    _, parent = completed
    modified = copy.deepcopy(parent["inputs"])
    modified["stage1"]["sha256"] = "changed"
    assert rt.plan(
        config, parent, input_hashes=modified, runtime=parent["implementation"]
    )["recompute_stages"] == list(rt.STAGES)
    p = copy.deepcopy(parent)
    f = tmp_path / "tampered"
    f.write_text("changed")
    p["stages"]["trajectory"]["outputs"]["stage1_h5ad"]["path"] = str(f)
    plan = rt.plan(
        config, p, input_hashes=parent["inputs"], runtime=parent["implementation"]
    )
    assert plan["recompute_stages"] == ["trajectory", "dashboard"]


def test_failure_blocks_downstream_and_retains_parent(inputs, completed):
    root, _, config = inputs
    result, parent = completed
    config = copy.deepcopy(config)
    config["subset"]["group"] = "absent"
    path = root / "bad.json"
    rt.write_json(path, config)
    with pytest.raises(ValueError, match="no cells"):
        rt.execute(path, root / "project", result["manifest"], run_id="bad")
    m = rt.read_json(root / "project/runs/bad/manifest.json")
    assert m["status"] == "failed" and m["stages"]["mapping"]["status"] == "failed"
    assert m["stages"]["alignment"]["status"] == "reused"
    assert all(m["stages"][s]["status"] == "blocked" for s in rt.STAGES[2:])
    assert Path(m["stages"]["mapping"]["error_log"]).is_file()


def test_partial_run_resume_and_disabled_stages(inputs):
    root, _, config = inputs
    config = copy.deepcopy(config)
    config["gp"]["enabled"] = False
    config["metrics"]["enabled"] = False
    config["trajectory"]["enabled"] = False
    path = root / "partial.json"
    rt.write_json(path, config)
    a = rt.execute(
        path, root / "partial-project", run_id="aligned", stop_after="alignment"
    )
    assert a["status"] == "partial"
    b = rt.execute(path, root / "partial-project", a["manifest"], run_id="continued")
    assert b["reused_stages"] == ["alignment"]
    m = rt.read_json(b["manifest"])
    assert b["status"] == "completed" and all(
        m["stages"][s]["status"] == "skipped" for s in ("gp", "metrics", "trajectory")
    )


def test_config_rejects_unknown_old_and_uninitialized_alignment(inputs, tmp_path):
    _, _, config = inputs
    for change in (
        {"old_unused_setting": True},
        {"schema_version": "spateo-4d-run/v2"},
    ):
        bad = copy.deepcopy(config)
        bad.update(change)
        p = tmp_path / "bad.json"
        rt.write_json(p, bad)
        with pytest.raises(ValueError):
            rt.canonical_config(p)
    bad = copy.deepcopy(config)
    bad["alignment"]["nonrigid_start_iter"] = 80
    p = tmp_path / "iter.json"
    rt.write_json(p, bad)
    with pytest.raises(ValueError, match="coefficients"):
        rt.canonical_config(p)


def test_original_counts_require_explicit_assertion(inputs):
    root, _, config = inputs
    a = ad.read_h5ad(root / "stage1.h5ad")
    del a.layers["counts"]
    with pytest.raises(ValueError, match="Missing counts"):
        rt.validate_input(a, config)
    config = copy.deepcopy(config)
    config["alignment"]["x_is_counts"] = True
    rt.validate_input(a, config)
    a.layers["counts"] = a.layers["counts"].astype(float)
    a.layers["counts"].data[0] = 0.5
    with pytest.raises(ValueError, match="integer"):
        rt.validate_input(a, config)
