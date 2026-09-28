"""Display source contracts must reproduce detector coordinates/capture, not guess."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "referee_contract_renderer",
    ROOT / "skills/spatial-slice-quality-qc/scripts/slice_quality_visualization.py",
)
viewer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(viewer)


@pytest.fixture
def cs8_source(tmp_path):
    a = AnnData(sparse.csr_matrix(np.full((8, 3), 0.5)), obs=pd.DataFrame({
        "sample_order": ["001"] * 4 + ["002"] * 4,
        "spatial_X": np.arange(8, dtype=float),
        "spatial_Y": [0, 1, 0, 1, 0, 1, 0, 1],
        "total_counts": [999] * 8,
        "n_genes_by_counts": [99] * 8,
    }, index=[f"c{i}" for i in range(8)]))
    # An unrelated 3D reconstruction must never override manifest obs coordinates.
    a.obsm["spatial3D"] = np.full((8, 3), 900.0)
    a.layers["counts"] = sparse.csr_matrix(np.tile([1, 0, 4], (8, 1)))
    path = tmp_path / "cs8.h5ad"
    a.write_h5ad(path)
    record = dict(path=str(path), slice_source="obs:sample_order",
                  coordinate_source="obs:spatial_X,spatial_Y", count_layer="counts",
                  count_like=True, expression_capture_available=True,
                  capture_loss_evidence_available=True)
    return path, record


def test_cs8_obs_xy_contract_and_selected_layer(cs8_source):
    path, record = cs8_source
    resolved = viewer._manifest_source_contract({"sources": [record]}, None, None, None)
    assert resolved == (path, "sample_order", "obs:spatial_X,spatial_Y")
    points, coordinates, metadata = viewer._prepare_source_evidence_points(
        path, ["001", "002"], 1.0, 100,
        slice_key=resolved[1], spatial_key=resolved[2], source_record=record,
    )
    np.testing.assert_array_equal(coordinates["001"][:, 0], np.arange(4))
    assert points["002"]["captured_counts"] == [5.0] * 4
    assert points["002"]["n_genes"] == [2.0] * 4
    assert metadata["total_counts_key"] == "counts:row_sum"
    assert metadata["absolute_expression_available"] is True


@pytest.mark.parametrize("source", ["obsm:spatial[:,0:2]", "obsm:spatial[:,:2]", "obsm:spatial"])
def test_obsm_coordinate_contract_formats(source):
    record = dict(path="source.h5ad", slice_source="obs:slice", coordinate_source=source)
    assert viewer._manifest_source_contract({"sources": [record]}, None, None, None)[2] == "spatial"


def test_annotation_or_normalized_unmeasured_capture_stays_unavailable(cs8_source):
    path, record = cs8_source
    record.update(count_layer="X", count_like=False, expression_capture_available=False,
                  capture_loss_evidence_available=False)
    points, _, metadata = viewer._prepare_source_evidence_points(
        path, ["001", "002"], 1.0, 100, slice_key="sample_order",
        spatial_key="obs:spatial_X,spatial_Y", source_record=record,
    )
    assert points["001"]["captured_counts"] == []
    assert points["001"]["n_genes"] == []
    assert metadata["absolute_expression_available"] is False


def test_explicit_bad_slice_key_never_guesses_z(cs8_source):
    path, record = cs8_source
    with pytest.raises(KeyError, match="Explicit slice key"):
        viewer._prepare_source_evidence_points(
            path, ["001"], 1.0, 100, slice_key="misspelled_slice",
            spatial_key="spatial3D", source_record=record,
        )


def test_sample_fallback_keeps_absolute_counts_missingness_and_total():
    samples = {"001": dict(x=[0, 1, 2], y=[0, 1, 0], counts=[12, None, 3], n_total=600)}
    points, _ = viewer._prepare_kde_points(samples, ["001"], 1.0, 100)
    assert points["001"]["captured_counts"] == [12.0, None, 3.0]
    assert points["001"]["log_captured_counts"][1] is None
    assert points["001"]["available"] == 600
    assert points["001"]["displayed"] == 3
    assert points["001"]["evidence_geometry_scope"] == "display sample only"


def test_source_missing_requested_slice_is_not_silently_empty(cs8_source):
    path, record = cs8_source
    with pytest.raises(ValueError, match="does not resolve requested slice"):
        viewer._prepare_source_evidence_points(
            path, ["absent"], 1.0, 100, slice_key="sample_order",
            spatial_key="obs:spatial_X,spatial_Y", source_record=record,
        )
