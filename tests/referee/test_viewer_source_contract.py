"""Display source contracts must reproduce detector coordinates/capture, not guess."""
import importlib.util
import json
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


@pytest.mark.parametrize("route", [None, np.nan, pd.NA, "", " "])
@pytest.mark.parametrize("language,label", [("en", "No confirmed route"), ("zh", "无确认缺失路径")])
def test_missing_loss_route_has_explicit_display_label(route, language, label):
    markup = viewer._loss_diagnostic_html({
        "loss_evidence_version": "bilateral-v3", "tissue_loss_route": route,
        "tissue_loss_reason": np.nan,
    }, language)
    assert f"<p>{label}</p>" in markup
    assert "nan" not in markup
    assert "None" not in markup


def test_presentation_refresh_preserves_all_other_payload_fields_and_html(tmp_path):
    report = tmp_path / "index.html"
    payload = {"language": "en", "order": ["s1"], "points": {"s1": {"x": [1.25]}},
               "regions": {"s1": {"note": "</script> must remain escaped"}}, "records": [{
                   "slice_id": "s1", "final_call": "keep", "quality_anomaly_score": 0.123,
                   "loss_evidence_version": "bilateral-v3", "tissue_loss_route": None,
                   "loss_diagnostic_html": "nan · stale markup",
               }]}
    start = '<html>unchanged<script id="payload" type="application/json">'
    end = '</script><p>unchanged footer</p></html>'
    report.write_text(start + json.dumps(payload).replace("</", "<\\/") + end)
    result = viewer.refresh_loss_diagnostic_html(report)
    content = report.read_text()
    assert content.startswith(start) and content.endswith(end)
    refreshed = json.loads(content[len(start):-len(end)])
    refreshed["records"][0].pop("loss_diagnostic_html")
    payload["records"][0].pop("loss_diagnostic_html")
    assert refreshed == payload
    assert result["updated_records"] == 1
    assert result["scientific_payload_unchanged"] is True
    assert viewer.refresh_loss_diagnostic_html(report)["updated_records"] == 0


def test_presentation_refresh_preserves_omitted_window_explanations_from_audit(tmp_path):
    details = json.dumps([{"window": 3, "geometry_loss_candidate": True, "capture_loss_candidate": False}])
    full_record = {"slice_id": "s1", "loss_evidence_version": "bilateral-v3",
                   "tissue_loss_route": None, "adaptive_window_details": details}
    expected = viewer._loss_diagnostic_html(full_record, "en")
    record = {key: value for key, value in full_record.items() if key != "adaptive_window_details"}
    record["loss_diagnostic_html"] = expected.replace("No confirmed route", "nan")
    report = tmp_path / "index.html"
    report.write_text('<script id="payload" type="application/json">' +
                      json.dumps({"language": "en", "records": [record]}) + '</script>')
    with pytest.raises(ValueError, match="corresponding binary audit"):
        viewer.refresh_loss_diagnostic_html(report)
    pd.DataFrame([{"slice_id": "s1", "adaptive_window_details": details}]).to_csv(
        tmp_path / "slice_quality_binary_audit.csv", index=False)
    result = viewer.refresh_loss_diagnostic_html(report)
    payload = json.loads(report.read_text().split('application/json">')[1].split('</script>')[0])
    assert payload["records"][0]["loss_diagnostic_html"] == expected
    assert "adaptive_window_details" not in payload["records"][0]
    assert result["window_details_audit_sha256"]
