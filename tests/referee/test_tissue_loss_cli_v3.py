"""Independent CLI publication tests; no installed Spateo package required.

These
tests call the real publication command function and write temporary artifacts.
They do not claim a full shell/bootstrap/viewer integration test.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import types

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[2]
SKILL = ROOT / "skills/spatial-slice-quality-qc"
RUNTIME = SKILL / "runtime/spateo/preprocessing/slice_quality.py"
if not RUNTIME.is_file():
    raise FileNotFoundError("Packaged QC runtime is required for CLI tests")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def cli_context(tmp_path, monkeypatch):
    for name in ["spateo", "spateo.preprocessing"]:
        package = types.ModuleType(name)
        package.__path__ = []
        monkeypatch.setitem(sys.modules, name, package)
    qc = _load("_v3_cli_independent_runtime", RUNTIME)
    monkeypatch.setitem(sys.modules, "spateo.preprocessing.slice_quality", qc)
    cli = _load("_v3_cli_independent_command", SKILL / "scripts/run_slice_quality_qc.py")
    source = tmp_path / "old_source"
    source.mkdir()
    frame = pd.DataFrame({
        "slice_id": [f"s{i}" for i in range(7)],
        "n_locations": [1000.0]*7, "hull_area": [100.0]*7, "cell_density": [10.0]*7,
        "median_total_counts": [100.0]*7, "median_n_genes": [50.0]*7,
        "expression_capture_available": [True]*7, "coordinate_valid_fraction": [1.0]*7,
        "hole_fraction": [0.0]*7, "fragmentation": [0.0]*7,
        "largest_component_fraction": [1.0]*7, "knn_tail_ratio": [0.2]*7,
        "expression_profile_anomaly": [0.0]*7, "celltype_composition_anomaly": [0.0]*7,
    })
    frame.loc[3, ["n_locations", "hull_area"]] = [200.0, 20.0]
    legacy = qc.SliceQCConfig(tissue_loss_enabled=False)
    old = qc._score_metrics(frame, None, {}, legacy, 3)
    old = qc.add_multiscale_exclusion_evidence(old, config=legacy)
    old.to_csv(source / "slice_quality_metrics.csv", index=False)
    (source / "slice_quality_manifest.json").write_text(json.dumps({"config": asdict(legacy), "selected_window": 3}))
    (source / "slice_quality_display_payload.json").write_text('{"fixture":"immutable display cache"}')
    config_keys = ["tissue_loss_min_fraction", "tissue_loss_min_geometry_fraction",
                   "tissue_loss_min_reference_points", "capture_loss_min_counts_fraction", "capture_loss_min_genes_fraction"]
    default = asdict(qc.SliceQCConfig())
    policy = {
        "keep_max_score": 0.129, "exclude_min_score": 0.7,
        "calibration_id": "cli-fixture-not-certification", "unresolved_action": "keep",
        "enable_tissue_loss_resolver": True, "loss_min_confirming_windows": 2,
        "benchmark_summary": {"base_policy_validation_passed": False,
                              "raw_matrix_stress_validation_passed": False,
                              "adaptive_extension_validation_passed": False,
                              "loss_config": {key: default[key] for key in config_keys}},
    }
    policy_path = tmp_path / "policy.json"
    policy_path.write_text(json.dumps({"policy": policy}))
    args = argparse.Namespace(input_dir=str(source), output_dir=str(tmp_path / "new_output"),
                              policy=str(policy_path), application_scope="experimental_policy",
                              complete_binary=True, allow_unvalidated_policy=True)
    return cli, source, policy_path, args


@pytest.mark.parametrize("same", [None, "exact", "normalized"])
def test_new_loss_publication_never_overwrites_input(cli_context, same):
    cli, source, _, args = cli_context
    before = {p.name: p.read_bytes() for p in source.iterdir()}
    args.output_dir = None if same is None else str(source if same == "exact" else source / ".." / source.name)
    with pytest.raises(SystemExit, match="(?i)(distinct|preserve|output)"):
        cli.run_publish(args)
    assert before == {p.name: p.read_bytes() for p in source.iterdir()}


def test_new_loss_publication_requires_experimental_scope(cli_context):
    cli, _, _, args = cli_context
    args.application_scope = "certified"
    with pytest.raises(SystemExit, match="(?i)(experimental|certif)"):
        cli.run_publish(args)
    assert not Path(args.output_dir).exists()


def test_unvalidated_policy_requires_explicit_override(cli_context):
    cli, _, _, args = cli_context
    args.allow_unvalidated_policy = False
    with pytest.raises(SystemExit, match="(?i)(validation|evaluation|validat)"):
        cli.run_publish(args)
    assert not Path(args.output_dir).exists()


def test_new_loss_publication_rescores_caches_and_preserves_source(cli_context):
    cli, source, policy_path, args = cli_context
    before = {p.name: p.read_bytes() for p in source.iterdir()}
    cli.run_publish(args)
    output = Path(args.output_dir)
    assert before == {p.name: p.read_bytes() for p in source.iterdir()}
    assert (output / "slice_quality_display_payload.json").read_bytes() == before["slice_quality_display_payload.json"]
    new = pd.read_csv(output / "slice_quality_metrics.csv")
    audit = pd.read_csv(output / "slice_quality_binary_audit.csv")
    manifest = json.loads((output / "slice_quality_manifest.json").read_text())
    summary = json.loads((output / "binary_policy_application.json").read_text())
    assert new.loss_evidence_version.eq("bilateral-v3").all()
    assert bool(new.set_index("slice_id").loc["s3", "geometry_loss_candidate"])
    assert audit.set_index("slice_id").loc["s3", "final_call"] == "exclude"
    assert not new.capture_loss_candidate.any()  # Legacy cache had no measured-count provenance.
    assert audit.certified_call.isna().all()
    assert manifest["loss_rescore"]["raw_matrix_reextracted"] is False
    assert manifest["loss_rescore"]["parent_metrics_sha256"] == hashlib.sha256(before["slice_quality_metrics.csv"]).hexdigest()
    assert summary["policy_sha256"] == hashlib.sha256(policy_path.read_bytes()).hexdigest()
    assert summary["application_scope"] == "experimental_policy"
    assert not summary["new_input_independently_validated"]
    assert manifest["config"]["tissue_loss_enabled"] is True


def test_incomplete_geometry_cache_requires_scan(cli_context):
    cli, source, _, args = cli_context
    path = source / "slice_quality_metrics.csv"
    frame = pd.read_csv(path).drop(columns=["hull_area"])
    frame.to_csv(path, index=False)
    with pytest.raises(SystemExit, match="(?i)(geometry|scan)"):
        cli.run_publish(args)
    assert not Path(args.output_dir).exists()
