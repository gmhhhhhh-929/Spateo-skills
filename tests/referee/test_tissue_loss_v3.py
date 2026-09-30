"""Independent adversarial regression tests for the experimental loss resolver.

The synthetic flat series are identifiability/unit-test fixtures, not biological
validation or an estimate of sensitivity. No injected dose is supplied to the
detector. Run with anndata/numpy/scipy/pandas/pytest installed; importing all of
Spateo is deliberately unnecessary.
"""
from __future__ import annotations

import copy
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from numpy.testing import assert_allclose, assert_array_equal
from pandas.testing import assert_frame_equal


MODULE_PATH = Path(__file__).resolve().parents[2] / "skills/spatial-slice-quality-qc/runtime/spateo/preprocessing/slice_quality.py"


@pytest.fixture(scope="session")
def qc():
    spec = importlib.util.spec_from_file_location("_referee_v3_independent_tests", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _policy(qc, **changes):
    values = dict(
        keep_max_score=0.129,
        exclude_min_score=0.700,
        calibration_id="independent-adversarial-fixture-not-certification",
        unresolved_action="keep",
        enable_tissue_loss_resolver=True,
        loss_min_confirming_windows=2,
    )
    values.update(changes)
    return qc.HighConfidencePolicy.from_mapping(values)


def _publish(qc, metrics, **changes):
    evidence = qc.add_multiscale_exclusion_evidence(metrics, config=qc.SliceQCConfig())
    return qc.apply_high_confidence_policy(evidence, _policy(qc, **changes))


def _raw_series(n_slices=7):
    side = np.arange(31, dtype=float)
    xx, yy = np.meshgrid(side, side)
    xy = np.column_stack((xx.ravel(), yy.ravel()))
    # Repeated measured integer matrices ensure a known flat, intact reference.
    counts = np.random.default_rng(491).poisson(1.5, (len(xy), 40)).astype(np.int32)
    counts[:, 0] += 1
    matrix = np.tile(counts, (n_slices, 1))
    labels = np.repeat([f"s{i}" for i in range(n_slices)], len(xy))
    obs = pd.DataFrame({"section": labels}, index=[f"cell{i}" for i in range(len(matrix))])
    adata = AnnData(matrix, obs=obs)
    adata.obsm["spatial"] = np.tile(xy, (n_slices, 1))
    return adata


def _inject(adata, kind, fraction=0.9, target="s3"):
    labels = np.asarray(adata.obs["section"].astype(str))
    indices = np.flatnonzero(labels == target)
    xy = np.asarray(adata.obsm["spatial"])[indices]
    rng = np.random.default_rng(1042)
    if kind == "capture":
        out = adata.copy()
        out.X[indices] = rng.binomial(out.X[indices], 1.0 - fraction)
        return out
    if kind == "edge":
        # Cutting from one edge retains a dense, compact residual strip.
        rank = np.lexsort((xy[:, 1], xy[:, 0]))
    elif kind == "tear":
        # Largest distance from the central line preserves both outer strips.
        rank = np.argsort(-np.abs(xy[:, 0] - np.mean(xy[:, 0])), kind="stable")
    elif kind == "random":
        rank = rng.permutation(len(xy))
    elif kind == "hole":
        rank = np.argsort(-np.linalg.norm(xy - xy.mean(axis=0), axis=1), kind="stable")
    else:
        raise ValueError(kind)
    keep = np.ones(adata.n_obs, dtype=bool)
    keep[indices] = False
    keep[indices[rank[:max(1, int(round(len(indices) * (1-fraction))) )]]] = True
    return adata[keep].copy()


def _calculate(qc, adata, *, enabled=True):
    return qc.calculate_slice_quality(
        adata, slice_key="section", spatial_key="spatial", layer="X",
        config=qc.SliceQCConfig(tissue_loss_enabled=enabled),
    )


def _row(frame, section="s3"):
    return frame.set_index("slice_id").loc[section]


def _metric_frame(n_slices=7):
    return pd.DataFrame({
        "slice_id": [f"s{i}" for i in range(n_slices)],
        "n_locations": [1000.0]*n_slices,
        "hull_area": [100.0]*n_slices,
        "cell_density": [10.0]*n_slices,
        "median_total_counts": [100.0]*n_slices,
        "median_n_genes": [50.0]*n_slices,
        "expression_capture_available": [True]*n_slices,
        "coordinate_valid_fraction": [1.0]*n_slices,
        "hole_fraction": [0.0]*n_slices,
        "fragmentation": [0.0]*n_slices,
        "largest_component_fraction": [1.0]*n_slices,
        "knn_tail_ratio": [0.2]*n_slices,
        "expression_profile_anomaly": [0.0]*n_slices,
        "celltype_composition_anomaly": [0.0]*n_slices,
    })


def _score(qc, frame, *, enabled=True, window=3):
    return qc._score_metrics(frame, None, {}, qc.SliceQCConfig(tissue_loss_enabled=enabled), window)


@pytest.mark.parametrize("kind", ["edge", "tear", "random", "hole"])
@pytest.mark.parametrize("fraction", [0.7, 0.9, 0.95])
def test_raw_geometry_losses_excluded_on_flat_identifiable_series(qc, kind, fraction):
    adata = _inject(_raw_series(), kind, fraction)
    scored = _calculate(qc, adata)
    target = _row(scored)
    assert target["loss_evidence_version"] == "bilateral-v3"
    assert bool(target["geometry_loss_candidate"])
    assert not bool(target["capture_loss_candidate"])
    published = _publish(qc, scored)
    assert _row(published)["final_call"] == "exclude"
    assert published.loc[published.slice_id.ne("s3"), "final_call"].eq("keep").all()


def test_raw_capture_loss_has_independent_route(qc):
    scored = _calculate(qc, _inject(_raw_series(), "capture", 0.9))
    assert bool(_row(scored)["capture_loss_candidate"])
    assert not bool(_row(scored)["geometry_loss_candidate"])
    assert _row(_publish(qc, scored))["final_call"] == "exclude"


def test_unchanged_raw_series_has_no_new_exclusions(qc):
    scored = _calculate(qc, _raw_series())
    assert not scored["geometry_loss_candidate"].any()
    assert not scored["capture_loss_candidate"].any()
    assert _publish(qc, scored)["final_call"].eq("keep").all()


def test_uniform_raw_subsampling_is_not_focal_loss(qc):
    adata = _raw_series()
    # Same sparse sampling in every section; no reference to the pristine input.
    keep_per_slice = np.sort(np.random.default_rng(100).choice(31*31, 96, replace=False))
    keep = np.concatenate([keep_per_slice+i*31*31 for i in range(7)])
    scored = _calculate(qc, adata[keep].copy())
    assert not scored["geometry_loss_candidate"].any()
    assert _publish(qc, scored)["final_call"].eq("keep").all()


def test_raw_monotone_taper_not_mistaken_for_focal_deletion(qc):
    adata = _raw_series()
    xy = np.asarray(adata.obsm["spatial"])
    labels = np.asarray(adata.obs.section.astype(str))
    keep = np.concatenate([np.flatnonzero((labels == f"s{i}") & (xy[:, 0] <= width))
                           for i, width in enumerate([30, 26, 22, 18, 14, 10, 6])])
    scored = _calculate(qc, adata[keep].copy())
    assert not scored["geometry_loss_candidate"].any()
    assert _publish(qc, scored)["final_call"].eq("keep").all()


@pytest.mark.parametrize("target", ["s0", "s6"])
def test_deleted_endpoint_remains_evidence_limited(qc, target):
    scored = _calculate(qc, _inject(_raw_series(), "edge", 0.95, target=target))
    assert not bool(_row(scored, target)["geometry_loss_candidate"])
    assert _row(_publish(qc, scored), target)["final_call"] == "keep"


@pytest.mark.parametrize("n_slices", [1, 2])
def test_singleton_and_pair_do_not_fabricate_bilateral_context(qc, n_slices):
    scored = _calculate(qc, _inject(_raw_series(n_slices), "edge", 0.95, target="s0"))
    assert not scored["geometry_loss_candidate"].any()
    assert not scored["capture_loss_candidate"].any()


def test_annotation_one_hot_does_not_fabricate_capture_evidence(qc):
    adata = _inject(_raw_series(), "capture", 0.95)
    adata.uns["blind_input_contract"] = {"expression_representation": "annotation one-hot"}
    scored = _calculate(qc, adata)
    assert not scored["capture_loss_candidate"].any()
    assert not scored["expression_capture_available"].any()


def test_normalized_noncount_matrix_is_not_measured_capture(qc):
    adata = _raw_series()
    adata.X = np.log1p(adata.X.astype(float))
    adata.X[np.asarray(adata.obs.section.astype(str)) == "s3"] *= 0.05
    scored = _calculate(qc, adata)
    assert not scored["capture_loss_candidate"].any()


def test_geometry_loss_does_not_require_expression_measurement(qc):
    adata = _inject(_raw_series(), "edge", 0.9)
    adata.uns["blind_input_contract"] = {"expression_representation": "annotation one-hot"}
    scored = _calculate(qc, adata)
    assert not scored["capture_loss_candidate"].any()
    assert bool(_row(scored)["geometry_loss_candidate"])
    assert _row(_publish(qc, scored))["final_call"] == "exclude"


def test_normalized_matrix_with_both_raw_summaries_can_supply_capture(qc):
    adata = _inject(_raw_series(), "capture", 0.9)
    adata.obs["total_counts"] = adata.X.sum(axis=1)
    adata.obs["n_genes_by_counts"] = np.count_nonzero(adata.X, axis=1)
    adata.X = np.log1p(adata.X.astype(float))
    scored = _calculate(qc, adata)
    assert bool(_row(scored)["capture_loss_candidate"])
    assert _row(_publish(qc, scored))["final_call"] == "exclude"


@pytest.mark.parametrize("bad_area", [0.0, -1.0, np.nan])
def test_invalid_geometry_does_not_become_correlated_corroboration(qc, bad_area):
    frame = _metric_frame()
    frame.loc[3, ["n_locations", "hull_area", "cell_density"]] = [0.0, bad_area, np.nan]
    scored = _score(qc, frame)
    assert not bool(_row(scored)["geometry_loss_candidate"])


@pytest.mark.parametrize("orientation", ["vertical", "diagonal"])
def test_collinear_survivor_does_not_invent_finite_area(qc, orientation):
    adata = _raw_series()
    labels = np.asarray(adata.obs.section.astype(str))
    xy = adata.obsm["spatial"]
    line = xy[:, 0] == (0 if orientation == "vertical" else xy[:, 1])
    keep = (labels != "s3") | line
    scored = _calculate(qc, adata[keep].copy())
    assert not bool(_row(scored)["geometry_loss_candidate"])


def test_legacy_disabled_configuration_matches_locked_v2_golden(qc):
    # Golden recorded from dfe1fbc... v2 on this exact sufficient-statistic fixture.
    frame = _metric_frame()
    frame.loc[3, ["n_locations", "hull_area"]] = [200.0, 20.0]
    scored = _score(qc, frame, enabled=False)
    assert_allclose(scored.quality_anomaly_score, [0, 0, 0, 0.4566, 0, 0, 0], atol=1e-12)
    assert_allclose(scored.density_domain_score, [0, 0, 0, 0.65, 0, 0, 0], atol=1e-12)
    assert_array_equal(scored.partial_structure_protection, [False, False, False, True, False, False, False])
    assert scored.recommendation.tolist() == ["keep", "keep", "keep", "review", "keep", "keep", "keep"]


def test_loss_publisher_does_not_early_return_at_old_low_score(qc):
    frame = _metric_frame()
    frame.loc[3, ["n_locations", "hull_area"]] = [200.0, 20.0]
    scored = _score(qc, frame)
    scored.loc[3, "quality_anomaly_score"] = 0.0
    scored.loc[3, "recommendation"] = "keep"
    scored.loc[3, "partial_structure_protection"] = True
    assert _row(_publish(qc, scored))["final_call"] == "exclude"


@pytest.mark.parametrize("missing", ["loss_evidence_version", "loss_evidence_config", "geometry_loss_candidate", "capture_loss_candidate"])
def test_missing_loss_contract_fails_closed_with_explicit_error(qc, missing):
    scored = _score(qc, _metric_frame())
    evidence = qc.add_multiscale_exclusion_evidence(scored, config=qc.SliceQCConfig())
    assert missing in evidence
    with pytest.raises((ValueError, KeyError), match="(?i)(loss|evidence|geometry|capture|missing|rescor)"):
        qc.apply_high_confidence_policy(evidence.drop(columns=[missing]), _policy(qc))


def test_wrong_version_cannot_be_published_as_v3(qc):
    scored = _score(qc, _metric_frame())
    evidence = qc.add_multiscale_exclusion_evidence(scored, config=qc.SliceQCConfig())
    evidence["loss_evidence_version"] = "unverified-old-cache"
    with pytest.raises((ValueError, KeyError), match="(?i)(loss|evidence|version|rescor)"):
        qc.apply_high_confidence_policy(evidence, _policy(qc))


def test_metric_input_not_mutated(qc):
    frame = _metric_frame()
    frame.loc[3, ["n_locations", "hull_area"]] = [200.0, 20.0]
    original = frame.copy(deep=True)
    scored = _score(qc, frame)
    _publish(qc, scored)
    assert_frame_equal(frame, original)


def test_raw_input_not_mutated(qc):
    adata = _inject(_raw_series(), "tear", 0.9)
    matrix = adata.X.copy()
    xy = adata.obsm["spatial"].copy()
    obs = adata.obs.copy(deep=True)
    uns = copy.deepcopy(adata.uns)
    _calculate(qc, adata)
    assert_array_equal(adata.X, matrix)
    assert_array_equal(adata.obsm["spatial"], xy)
    assert_frame_equal(adata.obs, obs)
    assert adata.uns == uns


def test_metric_physical_scale_invariance(qc):
    frame = _metric_frame()
    frame.loc[3, ["n_locations", "hull_area"]] = [200.0, 20.0]
    base = _publish(qc, _score(qc, frame))
    for scale in [0.01, 0.5, 10.0, 1000.0]:
        scaled = frame.copy()
        scaled["hull_area"] *= scale**2
        scaled["cell_density"] /= scale**2
        result = _publish(qc, _score(qc, scaled))
        assert_array_equal(result.geometry_loss_candidate, base.geometry_loss_candidate)
        assert_array_equal(result.final_call, base.final_call)


def test_raw_rigid_and_scale_invariance(qc):
    adata = _inject(_raw_series(), "edge", 0.9)
    base = _publish(qc, _calculate(qc, adata))
    transformed = adata.copy()
    angle = 0.719
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    transformed.obsm["spatial"] = (transformed.obsm["spatial"] @ rotation.T)*100.0 + [391.0, -81.0]
    result = _publish(qc, _calculate(qc, transformed))
    assert_array_equal(result.geometry_loss_candidate, base.geometry_loss_candidate)
    assert_array_equal(result.final_call, base.final_call)


def test_widening_tear_loss_route_cannot_disappear_on_flat_fixture(qc):
    outcomes = []
    for fraction in [0.7, 0.9, 0.95]:
        scored = _calculate(qc, _inject(_raw_series(), "tear", fraction))
        outcomes.append(_row(_publish(qc, scored))["final_call"])
    assert outcomes == ["exclude", "exclude", "exclude"]


def test_multiscale_details_preserve_typed_evidence(qc):
    frame = _metric_frame()
    frame.loc[3, ["n_locations", "hull_area"]] = [200.0, 20.0]
    evidence = qc.add_multiscale_exclusion_evidence(_score(qc, frame), config=qc.SliceQCConfig())
    details = json.loads(_row(evidence)["adaptive_window_details"])
    assert [item["window"] for item in details] == [3, 5, 7]
    assert all(item["geometry_loss_candidate"] for item in details)
    assert all(item["loss_evidence_version"] == "bilateral-v3" for item in details)


def test_short_series_confirms_only_available_window(qc):
    frame = _metric_frame(3)
    frame.loc[1, ["n_locations", "hull_area"]] = [200.0, 20.0]
    scored = _score(qc, frame)
    published = _publish(qc, scored)
    assert _row(published, "s1")["final_call"] == "exclude"
    details = json.loads(_row(published, "s1")["adaptive_window_details"])
    assert [item["window"] for item in details] == [3]


def test_missing_new_fields_in_window_cache_raises(qc):
    evidence = qc.add_multiscale_exclusion_evidence(_score(qc, _metric_frame()), config=qc.SliceQCConfig())
    details = json.loads(evidence.loc[3, "adaptive_window_details"])
    for item in details:
        item.pop("geometry_loss_candidate", None)
    evidence.loc[3, "adaptive_window_details"] = json.dumps(details)
    with pytest.raises((ValueError, KeyError), match="(?i)(loss|evidence|geometry|window|rescor)"):
        qc.apply_high_confidence_policy(evidence, _policy(qc))


def test_new_resolver_cannot_claim_certified_scope(qc, tmp_path):
    evidence = qc.add_multiscale_exclusion_evidence(_score(qc, _metric_frame()), config=qc.SliceQCConfig())
    with pytest.raises(ValueError, match="(?i)(experimental|certif)"):
        qc.write_high_confidence_outputs(evidence, _policy(qc), tmp_path / "wrong_scope")
    output = qc.write_high_confidence_outputs(evidence, _policy(qc), tmp_path / "experimental",
                                               application_scope="experimental_policy")
    audit = pd.read_csv(output["audit"])
    summary = json.loads(Path(output["summary"]).read_text())
    assert audit.certified_call.isna().all()
    assert not summary["new_input_independently_validated"]
    assert summary["certified"] == 0


def test_custom_threshold_cache_cannot_be_published_under_frozen_defaults(qc):
    frame = _metric_frame()
    frame.loc[3, ["n_locations", "hull_area"]] = [600.0, 60.0]
    custom = qc.SliceQCConfig(tissue_loss_min_fraction=0.2, tissue_loss_min_geometry_fraction=0.2)
    scored = qc._score_metrics(frame, None, {}, custom, 3)
    assert bool(_row(scored)["geometry_loss_candidate"])
    evidence = qc.add_multiscale_exclusion_evidence(scored, config=custom)
    with pytest.raises(ValueError, match="(?i)(frozen|configuration|threshold|rescor)"):
        qc.apply_high_confidence_policy(evidence, _policy(qc))


def test_stale_window_thresholds_cannot_mix_with_new_primary_evidence(qc):
    evidence = qc.add_multiscale_exclusion_evidence(_score(qc, _metric_frame()), config=qc.SliceQCConfig())
    details = json.loads(evidence.loc[3, "adaptive_window_details"])
    stale = json.loads(details[1]["loss_evidence_config"])
    stale["tissue_loss_min_fraction"] = 0.2
    details[1]["loss_evidence_config"] = json.dumps(stale)
    evidence.loc[3, "adaptive_window_details"] = json.dumps(details)
    with pytest.raises(ValueError, match="(?i)(frozen|configuration|threshold|rescor)"):
        qc.apply_high_confidence_policy(evidence, _policy(qc))


@pytest.mark.parametrize("minimum", [0, -1, 1.5, True])
def test_confirmation_count_must_be_a_positive_integer(qc, minimum):
    evidence = qc.add_multiscale_exclusion_evidence(_score(qc, _metric_frame()), config=qc.SliceQCConfig())
    with pytest.raises(ValueError, match="(?i)(integer|window|confirm)"):
        qc.apply_high_confidence_policy(evidence, _policy(qc, loss_min_confirming_windows=minimum))
