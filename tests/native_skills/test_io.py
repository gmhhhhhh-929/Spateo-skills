import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import spateo_io as cli
from validate_anndata import validate


def visium(root, bad=False):
    root.mkdir(parents=True)
    X = sparse.csc_matrix([[1, 3], [2, 4]])
    with h5py.File(root / "filtered_feature_bc_matrix.h5", "w") as f:
        g = f.create_group("matrix")
        for key, data in [
            ("data", X.data),
            ("indices", X.indices),
            ("indptr", X.indptr),
            ("shape", X.shape),
            ("barcodes", np.array(["c1", "c2"], dtype="S")),
        ]:
            g.create_dataset(key, data=data)
        feats = g.create_group("features")
        for key, data in [
            ("id", ["g1", "g2"]),
            ("name", ["G1", "G2"]),
            ("feature_type", ["Gene Expression"] * 2),
        ]:
            feats.create_dataset(key, data=np.array(data, dtype="S"))
    (root / "spatial").mkdir()
    pd.DataFrame(
        {
            "barcode": ["c2", "missing" if bad else "c1"],
            "in_tissue": [1, 1],
            "array_row": [1, 0],
            "array_col": [1, 0],
            "pxl_row_in_fullres": [11, 10],
            "pxl_col_in_fullres": [21, 20],
        }
    ).to_csv(root / "spatial/tissue_positions.csv", index=False)
    return root


def test_discovery_is_deferred_not_fake_anndata(tmp_path, capsys):
    p = visium(tmp_path / "input")
    assert cli.main(["discover", str(p)]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "pending"
    assert all(e["status"] == "deferred" for e in report["datasets"].values())
    assert not list(p.glob("*.h5ad"))


def test_real_read_export_and_roundtrip(tmp_path, capsys):
    p = visium(tmp_path / "input")
    out = tmp_path / "out"
    assert (
        cli.main(
            [
                "read",
                str(p),
                "--output-dir",
                str(out),
                "--matrix-semantics",
                "counts",
                "--coordinate-unit",
                "fullres_pixel",
            ]
        )
        == 0
    )
    report = json.loads(capsys.readouterr().out)
    saved = next(iter(report["exports"].values()))
    a = ad.read_h5ad(saved["output"])
    np.testing.assert_array_equal(a.X.toarray(), [[1, 2], [3, 4]])
    np.testing.assert_array_equal(a.obsm["spatial"], [[20, 10], [21, 11]])
    assert a.obs_names.tolist() == ["c1", "c2"]
    assert json.loads(Path(saved["manifest"]).read_text())["roundtrip"] == "pass"
    with pytest.raises(FileExistsError):
        cli.main(["read", str(p), "--output-dir", str(out)])


def test_partial_scope_preserves_good_output_and_bad_diagnostics(tmp_path, capsys):
    visium(tmp_path / "inputs/good")
    visium(tmp_path / "inputs/bad", bad=True)
    out = tmp_path / "export"
    assert cli.main(["read", str(tmp_path / "inputs"), "--output-dir", str(out)]) == 2
    r = json.loads(capsys.readouterr().out)
    assert r["status"] == "partial" and len(r["exports"]) == 1
    assert sorted(e["status"] for e in r["datasets"].values()) == ["failed", "ready"]
    assert r["export_status"] == "incomplete"
    assert (out / "read_report.json").is_file()


def test_h5ad_preserves_layers_raw_xyz_and_images(tmp_path, capsys):
    a = ad.AnnData(
        sparse.csr_matrix([[1, 2], [3, 4]]),
        obs=pd.DataFrame({"label": ["A", "B"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    a.layers["counts"] = a.X.copy()
    a.raw = a.copy()
    a.obsm["spatial"] = np.arange(6).reshape(2, 3)
    a.uns["spatial"] = {
        "library": {
            "images": {"hires": np.arange(16).reshape(4, 4)},
            "scalefactors": {"x": 0.5},
        }
    }
    a.uns["__type"] = "original_type"
    src = tmp_path / "source.h5ad"
    dst = tmp_path / "new.h5ad"
    a.write_h5ad(src)
    original = src.read_bytes()
    assert cli.main(["convert", str(src), str(dst), "--reader", "read_h5ad"]) == 0
    b = ad.read_h5ad(dst)
    assert b.uns["__type"] == "original_type"
    cli.assert_equal(a.uns, b.uns)
    cli.assert_equal(a.layers, b.layers)
    cli.assert_equal(a.raw.X, b.raw.X)
    assert src.read_bytes() == original
    with pytest.raises(FileExistsError):
        cli.main(["convert", str(src), str(dst), "--reader", "read_h5ad"])


def test_bad_counts_coordinates_and_agg_rejected():
    a = ad.AnnData(np.array([[1.0, 2.0], [3.0, 4.0]]))
    a.obsm["spatial"] = np.array([[0.0, 1.0], [2.0, 3.0]])
    assert validate(a, matrix_semantics="counts")["status"] == "pass"
    a.X[0, 0] = 0.5
    assert validate(a, matrix_semantics="counts")["status"] == "fail"
    a.X[0, 0] = 1
    a.obsm["spatial"][0, 0] = np.nan
    assert validate(a)["status"] == "fail"
    a.obsm["spatial"][0, 0] = 0
    a.uns["__type"] = "AGG"
    assert validate(a)["status"] == "fail"


def test_unknown_dataset_selection_is_not_success(tmp_path, capsys):
    p = visium(tmp_path / "input")
    assert (
        cli.main(
            [
                "read",
                str(p),
                "--output-dir",
                str(tmp_path / "out"),
                "--dataset-key",
                "nonexistent",
            ]
        )
        == 2
    )
    report = json.loads(capsys.readouterr().out)
    assert report["unknown_selected_keys"] == ["nonexistent"] and not report["exports"]


def test_explicit_stable_gene_ids_preserve_symbols_and_counts(tmp_path, capsys):
    p = visium(tmp_path / "input")
    with h5py.File(p / "filtered_feature_bc_matrix.h5", "r+") as f:
        del f["matrix/features/name"]
        f["matrix/features"].create_dataset(
            "name", data=np.array(["same", "same"], dtype="S")
        )
    rejected = tmp_path / "rejected"
    assert cli.main(["read", str(p), "--output-dir", str(rejected)]) == 2
    capsys.readouterr()
    out = tmp_path / "stable"
    assert (
        cli.main(
            [
                "read",
                str(p),
                "--output-dir",
                str(out),
                "--feature-id-column",
                "gene_ids",
            ]
        )
        == 0
    )
    report = json.loads(capsys.readouterr().out)
    record = next(iter(report["exports"].values()))
    a = ad.read_h5ad(record["output"])
    assert a.var_names.tolist() == ["g1", "g2"]
    assert a.var.source_var_name.tolist() == ["same", "same"]
    np.testing.assert_array_equal(a.X.toarray(), [[1, 2], [3, 4]])
    assert (
        json.loads(Path(record["manifest"]).read_text())["feature_id_column"]
        == "gene_ids"
    )
