import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyvista as pv
import pytest
import spateo as st

import build_point_cloud as point_cloud


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_input(path: Path, *, spatial=None) -> ad.AnnData:
    n_obs = 12
    matrix = np.arange(n_obs * 4, dtype=float).reshape(n_obs, 4)
    obs = pd.DataFrame(
        {
            "celltype": ["A", "B", "C"] * 4,
            "score": np.linspace(0.0, 1.0, n_obs),
        },
        index=[f"cell-{i:02d}" for i in range(n_obs)],
    )
    adata = ad.AnnData(
        matrix, obs=obs, var=pd.DataFrame(index=["g1", "g2", "g3", "g4"])
    )
    adata.layers["counts"] = matrix + 1.0
    if spatial is None:
        spatial = np.column_stack(
            [
                np.linspace(0.0, 11.0, n_obs),
                np.array([0, 2, 1, 4, 3, 7, 8, 5, 10, 6, 9, 12], dtype=float),
                np.array([0, 1, 4, 2, 8, 3, 7, 11, 5, 13, 6, 10], dtype=float),
            ]
        )
    adata.obsm["spatial"] = spatial
    adata.write_h5ad(path)
    return adata


def run_builder(*arguments: str) -> None:
    assert point_cloud.main(list(arguments)) == 0


def test_categorical_mask_and_vtk_roundtrip(tmp_path):
    input_path = tmp_path / "input.h5ad"
    write_input(input_path)
    before_hash = sha256(input_path)
    palette = tmp_path / "palette.json"
    palette.write_text(json.dumps({"A": "#e41a1c", "B": "#377eb8", "C": "#4daf4a"}))
    output = tmp_path / "categorical"

    run_builder(
        "--input-h5ad",
        str(input_path),
        "--output-dir",
        str(output),
        "--name",
        "cells",
        "--obs",
        "celltype",
        "--key-added",
        "celltype",
        "--palette-json",
        str(palette),
        "--mask",
        "C",
        "--skip-preview",
    )

    assert sha256(input_path) == before_hash
    model = st.tdr.read_model(str(output / "cells.vtk"))
    assert isinstance(model, pv.PolyData)
    assert model.n_points == 12
    assert set(model.point_data) == {"celltype", "celltype_rgba", "obs_index"}
    assert np.asarray(model.point_data["obs_index"]).astype(str).tolist() == [
        f"cell-{i:02d}" for i in range(12)
    ]
    labels = np.asarray(model.point_data["celltype"]).astype(str)
    rgba = np.asarray(model.point_data["celltype_rgba"])
    assert set(labels) == {"A", "B", "C"}
    assert np.all(rgba[labels == "C", 3] == 0)
    manifest = json.loads((output / "cells.manifest.json").read_text())
    assert manifest["status"] == "pass"
    assert manifest["model"]["n_points"] == manifest["input"]["n_obs"] == 12
    assert manifest["roundtrip"]["point_data_equal"] is True
    assert manifest["output"]["preview"] is None


def test_multi_gene_sum_is_continuous_and_manifested(tmp_path):
    input_path = tmp_path / "input.h5ad"
    adata = write_input(input_path)
    output = tmp_path / "genes"

    run_builder(
        "--input-h5ad",
        str(input_path),
        "--output-dir",
        str(output),
        "--gene",
        "g1",
        "--gene",
        "g2",
        "--layer",
        "counts",
        "--key-added",
        "gene_sum",
        "--colormap",
        "viridis",
        "--skip-preview",
    )

    model = st.tdr.read_model(str(output / "point_cloud.vtk"))
    expected = np.asarray(adata.layers["counts"][:, :2].sum(axis=1)).reshape(-1)
    np.testing.assert_allclose(model.point_data["gene_sum"], expected)
    assert "gene_sum_rgba" not in model.point_data
    manifest = json.loads((output / "point_cloud.manifest.json").read_text())
    assert manifest["coloring"]["mode"] == "gene_sum"
    assert manifest["coloring"]["plot_cmap_returned_by_spateo"] == "viridis"
    assert manifest["coloring"]["summary"]["kind"] == "continuous"


def test_external_labels_are_joined_by_id_and_preview_renders(tmp_path):
    input_path = tmp_path / "input.h5ad"
    adata = write_input(input_path)
    labels = pd.DataFrame(
        {
            "cell_id": list(reversed(adata.obs_names.tolist())),
            "region": (["left", "right"] * 6)[::-1],
        }
    )
    label_path = tmp_path / "labels.tsv"
    labels.to_csv(label_path, sep="\t", index=False)
    output = tmp_path / "external"

    run_builder(
        "--input-h5ad",
        str(input_path),
        "--output-dir",
        str(output),
        "--labels-file",
        str(label_path),
        "--labels-id-column",
        "cell_id",
        "--labels-value-column",
        "region",
        "--key-added",
        "region",
        "--preview-max-points",
        "7",
    )

    preview = output / "point_cloud.preview.png"
    assert preview.is_file() and preview.stat().st_size > 0
    model = st.tdr.read_model(str(output / "point_cloud.vtk"))
    assert (
        np.asarray(model.point_data["region"]).astype(str).tolist()
        == ["left", "right"] * 6
    )
    manifest = json.loads((output / "point_cloud.manifest.json").read_text())
    assert manifest["coloring"]["mode"] == "external"
    assert manifest["output"]["preview"]["n_points"] == 7
    assert manifest["output"]["preview"]["sampling"] == "uniform_without_replacement"
    assert manifest["output"]["preview"]["legend"] == {
        "included": True,
        "categories": ["left", "right"],
    }


@pytest.mark.parametrize(
    "spatial, message",
    [
        (np.ones((12, 2)), "must have shape"),
        (
            np.column_stack([np.arange(12), np.arange(12) ** 2, np.full(12, np.nan)]),
            "non-finite",
        ),
        (
            np.column_stack([np.arange(12), np.arange(12) ** 2, np.zeros(12)]),
            "coordinate rank",
        ),
    ],
)
def test_invalid_spatial_coordinates_fail(tmp_path, spatial, message):
    input_path = tmp_path / "invalid.h5ad"
    write_input(input_path, spatial=spatial)
    with pytest.raises(ValueError, match=message):
        run_builder(
            "--input-h5ad",
            str(input_path),
            "--output-dir",
            str(tmp_path / "output"),
            "--skip-preview",
        )


def test_planar_override_external_id_mismatch_and_overwrite_refusal(tmp_path):
    input_path = tmp_path / "planar.h5ad"
    planar = np.column_stack([np.arange(12), np.arange(12) ** 2, np.zeros(12)])
    write_input(input_path, spatial=planar)
    good_output = tmp_path / "planar_output"
    run_builder(
        "--input-h5ad",
        str(input_path),
        "--output-dir",
        str(good_output),
        "--allow-planar",
        "--ascii",
        "--skip-preview",
    )
    planar_manifest = json.loads(
        (good_output / "point_cloud.manifest.json").read_text()
    )
    assert planar_manifest["output"]["vtk_binary"] is False
    with pytest.raises(FileExistsError, match="non-empty"):
        run_builder(
            "--input-h5ad",
            str(input_path),
            "--output-dir",
            str(good_output),
            "--allow-planar",
            "--skip-preview",
        )

    bad_labels = tmp_path / "bad.csv"
    pd.DataFrame({"obs_index": ["cell-00"], "value": ["x"]}).to_csv(
        bad_labels, index=False
    )
    with pytest.raises(ValueError, match="cover exactly"):
        run_builder(
            "--input-h5ad",
            str(input_path),
            "--output-dir",
            str(tmp_path / "bad_external"),
            "--labels-file",
            str(bad_labels),
            "--allow-planar",
            "--skip-preview",
        )
