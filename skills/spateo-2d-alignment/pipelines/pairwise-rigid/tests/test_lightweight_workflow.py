#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import shutil
import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


RELEASE = Path(__file__).resolve().parents[1]
LOCK = RELEASE / "skill.lock.yaml"


def run_json(cmd: list[str]) -> dict:
    proc = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE)
    return json.loads(proc.stdout)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_fixture(tmp: Path) -> tuple[Path, Path]:
    coord_rows: list[dict[str, object]] = []
    comp_rows: list[dict[str, object]] = []
    spec = [
        (1, "c1_small", 1, 12),
        (1, "c1_large", 2, 180),
        (2, "c2_small", 1, 8),
        (2, "c2_large", 2, 160),
    ]
    for sl, comp_id, rank, count in spec:
        for i in range(count):
            cid = f"SL{sl}_{comp_id}_{i:04d}"
            coord_rows.append(
                {
                    "cell_id": cid,
                    "sl_number": sl,
                    "celltype": "type_a" if i % 2 == 0 else "type_b",
                    "x": sl * 1000 + rank * 100 + (i % 30),
                    "y": rank * 1000 + (i // 30),
                    "z": sl * 40,
                    "manual_x": sl * 1000 + rank * 100 + (i % 30),
                    "manual_y": rank * 1000 + (i // 30),
                    "manual_z": sl * 40,
                }
            )
            comp_rows.append(
                {
                    "cell_id": cid,
                    "sl_number": sl,
                    "component_rank": rank,
                    "component_id": comp_id,
                    "component_points": count,
                }
            )
    coords = tmp / "coords.csv"
    comps = tmp / "components.csv"
    write_csv(coords, coord_rows)
    write_csv(comps, comp_rows)
    return coords, comps


class LightweightWorkflowTests(unittest.TestCase):
    def test_lock_resolves_relative_release_root(self) -> None:
        report = run_json(
            [
                sys.executable,
                str(RELEASE / "spalign/lock_resolver.py"),
                "--lock",
                str(LOCK),
                "--mode",
                "balanced_spatial_sample",
                "--json",
            ]
        )
        self.assertEqual(report["pipeline_release"], "v0.2.3")
        self.assertTrue(report["entrypoint_path"].endswith("scripts/sampling/create_balanced_spatial_sample.py"))

    def test_hash_mismatch_fails(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            bad_lock = Path(td) / "skill.lock.yaml"
            data = json.loads(LOCK.read_text(encoding="utf-8"))
            data["release_root"] = str(RELEASE)
            data["allowed_entrypoints"]["balanced_spatial_sample"]["sha256"] = "0" * 64
            bad_lock.write_text(json.dumps(data), encoding="utf-8")
            proc = subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "spalign/lock_resolver.py"),
                    "--lock",
                    str(bad_lock),
                    "--mode",
                    "balanced_spatial_sample",
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("hash mismatch", proc.stderr.lower() + proc.stdout.lower())

    def test_spatial_only_validator_pass_and_fail(self) -> None:
        resolved = run_json(
            [
                sys.executable,
                str(RELEASE / "spalign/lock_resolver.py"),
                "--lock",
                str(LOCK),
                "--mode",
                "spatial_only_rigid",
                "--json",
            ]
        )
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            edge = root / "edge"
            edge.mkdir()
            payload = root / "payload.sh"
            payload.write_text(
                "\n".join(
                    [
                        f"export RUNNER_SCRIPT={shlex.quote(resolved['entrypoint_path'])}",
                        "export EXPRESSION_MODE=spatial_only",
                        "export SIGMA2_INIT_SCALE=2.0",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            provenance = {
                "runner_script": resolved["entrypoint_path"],
                "expression_mode": "spatial_only",
                "sigma2_init_scale": 2.0,
            }
            edge_transform = {
                "runner_script": resolved["entrypoint_path"],
                "runner_parameters": {"expression_mode": "spatial_only", "sigma2_init_scale": 2.0},
            }
            (edge / "spateo_two_stage_provenance.json").write_text(json.dumps(provenance), encoding="utf-8")
            (edge / "edge_transform.json").write_text(json.dumps(edge_transform), encoding="utf-8")
            report = run_json(
                [
                    sys.executable,
                    str(RELEASE / "spalign/validate_run.py"),
                    "--lock",
                    str(LOCK),
                    "--mode",
                    "spatial_only_rigid",
                    "--payload",
                    str(payload),
                    "--edge-dir",
                    str(edge),
                ]
            )
            self.assertEqual(report["status"], "pass")

            payload.write_text(payload.read_text(encoding="utf-8").replace("spatial_only", "normal"), encoding="utf-8")
            report = run_json(
                [
                    sys.executable,
                    str(RELEASE / "spalign/validate_run.py"),
                    "--lock",
                    str(LOCK),
                    "--mode",
                    "spatial_only_rigid",
                    "--payload",
                    str(payload),
                    "--edge-dir",
                    str(edge),
                ]
            )
            self.assertEqual(report["status"], "fail")

    def test_balanced300k_target_sampling(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            coords, comps = make_fixture(tmp)
            out_csv = tmp / "balanced.csv"
            out_summary = tmp / "summary.csv"
            out_manifest = tmp / "manifest.json"
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/sampling/create_balanced_spatial_sample.py"),
                    "--coordinates",
                    str(coords),
                    "--components",
                    str(comps),
                    "--output-csv",
                    str(out_csv),
                    "--output-summary-csv",
                    str(out_summary),
                    "--output-manifest",
                    str(out_manifest),
                    "--target-total-points",
                    "120",
                    "--keep-all-below",
                    "20",
                    "--min-per-component",
                    "10",
                    "--component-key",
                    "component_id",
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            manifest = json.loads(out_manifest.read_text(encoding="utf-8"))
            self.assertEqual(manifest["sampled_points"], 120)
            self.assertEqual(manifest["target_mode"], "target_total_points")
            with out_csv.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 120)
            self.assertEqual({row["sl_number"] for row in rows}, {"1", "2"})
            small_counts = {
                key: sum(1 for row in rows if row["sample_component_id"] == key)
                for key in ["c1_small", "c2_small"]
            }
            self.assertEqual(small_counts, {"c1_small": 12, "c2_small": 8})
            self.assertTrue({"celltype", "sample_component_id", "sample_component_key"} <= set(rows[0]))

    @unittest.skipUnless(shutil.which(sys.executable), "Python executable is required")
    def test_viewer_metadata_when_numpy_available(self) -> None:
        try:
            import numpy  # noqa: F401
        except Exception:
            self.skipTest("viewer test requires numpy")
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            coords, comps = make_fixture(tmp)
            sample = tmp / "balanced.csv"
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/sampling/create_balanced_spatial_sample.py"),
                    "--coordinates",
                    str(coords),
                    "--components",
                    str(comps),
                    "--output-csv",
                    str(sample),
                    "--output-summary-csv",
                    str(tmp / "summary.csv"),
                    "--output-manifest",
                    str(tmp / "manifest.json"),
                    "--target-total-points",
                    "120",
                    "--keep-all-below",
                    "20",
                    "--min-per-component",
                    "10",
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            summary = tmp / "viewer_summary.json"
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/viewers/make_st_pointcloud_viewer.py"),
                    "--input",
                    str(sample),
                    "--format",
                    "csv",
                    "--x-col",
                    "manual_x",
                    "--y-col",
                    "manual_y",
                    "--z-col",
                    "manual_z",
                    "--celltype-col",
                    "celltype",
                    "--slice-col",
                    "sl_number",
                    "--viewer-kind",
                    "review_balanced300k",
                    "--component-balanced",
                    "--sampling-method",
                    "slice x component_id spatial_grid_balanced",
                    "--sampling-seed",
                    "13",
                    "--source-full-csv",
                    str(coords),
                    "--source-full-row-count",
                    "360",
                    "--output-html",
                    str(tmp / "viewer.html"),
                    "--output-summary",
                    str(summary),
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            meta = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(meta["displayed_points"], 120)
            self.assertEqual(meta["full_points"], 360)
            self.assertTrue(meta["component_balanced"])
            self.assertEqual(meta["viewer_kind"], "review_balanced300k")
            self.assertEqual(meta["celltype_col"], "celltype")
            self.assertEqual(meta["component_col"], "sample_component_id")
            self.assertEqual(meta["component_rank_col"], "sample_component_rank")
            self.assertTrue(meta["component_label_overlay"])
            self.assertTrue(any(group["name"].startswith("SL1 rank1 id") for group in meta["components"]))
            self.assertTrue(all("centroid" in group for group in meta["components"]))
            html = (tmp / "viewer.html").read_text(encoding="utf-8")
            self.assertIn("Component", html)
            self.assertIn("Labels", html)

    @unittest.skipUnless(shutil.which(sys.executable), "Python executable is required")
    def test_full_points_review_viewer_policy(self) -> None:
        try:
            import numpy  # noqa: F401
        except Exception:
            self.skipTest("viewer test requires numpy")
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            coords, comps = make_fixture(tmp)
            with coords.open(encoding="utf-8", newline="") as handle:
                coord_rows = {row["cell_id"]: row for row in csv.DictReader(handle)}
            with comps.open(encoding="utf-8", newline="") as handle:
                comp_rows = list(csv.DictReader(handle))
            full_review = tmp / "full_component_labeled.csv"
            merged_rows = []
            for row in comp_rows:
                merged = coord_rows[row["cell_id"]].copy()
                merged.update(row)
                merged_rows.append(merged)
            write_csv(full_review, merged_rows)
            summary = tmp / "full_viewer_summary.json"
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/viewers/make_st_pointcloud_viewer.py"),
                    "--input",
                    str(full_review),
                    "--format",
                    "csv",
                    "--x-col",
                    "manual_x",
                    "--y-col",
                    "manual_y",
                    "--z-col",
                    "manual_z",
                    "--celltype-col",
                    "celltype",
                    "--slice-col",
                    "sl_number",
                    "--component-col",
                    "component_id",
                    "--component-rank-col",
                    "component_rank",
                    "--viewer-kind",
                    "review_full_points",
                    "--source-full-csv",
                    str(coords),
                    "--source-full-row-count",
                    "360",
                    "--output-html",
                    str(tmp / "full_viewer.html"),
                    "--output-summary",
                    str(summary),
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            meta = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(meta["display_policy"], "manual_full_points")
            self.assertEqual(meta["displayed_points"], 360)
            self.assertFalse(meta["sampled"])
            self.assertEqual(meta["component_col"], "component_id")

    def test_clean_export_and_validation_rejects_pollution(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            audit = tmp / "audit.csv"
            rows = [
                {
                    "cell_id": "a",
                    "slice_id": "s1",
                    "stage": "CS",
                    "chip_id": "chip",
                    "sl_number": "1",
                    "celltype": "neural",
                    "manual_x": "0",
                    "manual_y": "0",
                    "manual_z": "40",
                    "final_x": "1",
                    "final_y": "2",
                    "final_z": "40",
                    "full_candidate_edit_label": "demo",
                }
            ]
            write_csv(audit, rows)
            clean = tmp / "clean.csv"
            manifest = tmp / "manifest.json"
            changed = tmp / "changed.csv"
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/core/export_clean_coordinates.py"),
                    "--input",
                    str(audit),
                    "--output",
                    str(clean),
                    "--manifest",
                    str(manifest),
                    "--changed-cells-output",
                    str(changed),
                    "--x-col",
                    "final_x",
                    "--y-col",
                    "final_y",
                    "--z-col",
                    "final_z",
                    "--baseline-x-col",
                    "manual_x",
                    "--baseline-y-col",
                    "manual_y",
                    "--baseline-z-col",
                    "manual_z",
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            with clean.open(encoding="utf-8") as handle:
                clean_rows = list(csv.DictReader(handle))
            self.assertEqual(list(clean_rows[0].keys()), ["cell_id", "slice_id", "stage", "chip_id", "sl_number", "celltype", "x", "y", "z"])
            self.assertEqual(clean_rows[0]["x"], "1")
            clean_manifest = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertEqual(clean_manifest["row_count"], 1)
            self.assertEqual(clean_manifest["changed_row_count"], 1)
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/core/validate_clean_coordinates.py"),
                    "--input",
                    str(clean),
                    "--expected-row-count",
                    "1",
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            polluted = tmp / "polluted.csv"
            write_csv(polluted, [{"cell_id": "a", "slice_id": "s1", "stage": "CS", "chip_id": "chip", "sl_number": "1", "celltype": "neural", "x": "1", "y": "2", "z": "40", "manual_x": "0"}])
            proc = subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/core/validate_clean_coordinates.py"),
                    "--input",
                    str(polluted),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertNotEqual(proc.returncode, 0)

    def test_clean_export_can_remap_workflow_ids_to_short_slice_cell_ids(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            audit = tmp / "audit.csv"
            write_csv(
                audit,
                [
                    {
                        "cell_id": "CS13_SL35_Y01636D4:0",
                        "slice_id": "CS13_SL35_Y01636D4",
                        "stage": "CS13",
                        "chip_id": "Y01636D4",
                        "sl_number": "35",
                        "celltype": "neural",
                        "before_x": "0",
                        "before_y": "0",
                        "before_z": "280",
                        "candidate_x": "1",
                        "candidate_y": "2",
                        "candidate_z": "280",
                        "candidate_edit_label": "demo_edit",
                    }
                ],
            )
            id_map = tmp / "cell_id_map.csv"
            write_csv(
                id_map,
                [
                    {
                        "workflow_cell_id": "CS13_SL35_Y01636D4:0",
                        "cell_id": "SL35_CELL.1",
                    }
                ],
            )
            clean = tmp / "clean.csv"
            manifest = tmp / "manifest.json"
            changed = tmp / "changed.csv"
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/core/export_clean_coordinates.py"),
                    "--input",
                    str(audit),
                    "--output",
                    str(clean),
                    "--manifest",
                    str(manifest),
                    "--changed-cells-output",
                    str(changed),
                    "--x-col",
                    "candidate_x",
                    "--y-col",
                    "candidate_y",
                    "--z-col",
                    "candidate_z",
                    "--baseline-x-col",
                    "before_x",
                    "--baseline-y-col",
                    "before_y",
                    "--baseline-z-col",
                    "before_z",
                    "--cell-id-map",
                    str(id_map),
                    "--clean-cell-id-source",
                    "short_slice_cellid",
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            with clean.open(encoding="utf-8") as handle:
                clean_rows = list(csv.DictReader(handle))
            self.assertEqual(clean_rows[0]["cell_id"], "SL35_CELL.1")
            self.assertNotIn("workflow_cell_id", clean_rows[0])
            with changed.open(encoding="utf-8") as handle:
                changed_rows = list(csv.DictReader(handle))
            self.assertEqual(changed_rows[0]["cell_id"], "SL35_CELL.1")
            self.assertEqual(changed_rows[0]["workflow_cell_id"], "CS13_SL35_Y01636D4:0")
            clean_manifest = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertEqual(clean_manifest["mapped_cell_id_count"], 1)
            self.assertEqual(clean_manifest["identity_mapping"]["clean_cell_id_source"], "short_slice_cellid")
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/core/validate_clean_coordinates.py"),
                    "--input",
                    str(clean),
                    "--expected-row-count",
                    "1",
                    "--forbid-row-index-cell-id-pattern",
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            old_id_clean = tmp / "old_id_clean.csv"
            write_csv(
                old_id_clean,
                [
                    {
                        "cell_id": "CS13_SL35_Y01636D4:0",
                        "slice_id": "CS13_SL35_Y01636D4",
                        "stage": "CS13",
                        "chip_id": "Y01636D4",
                        "sl_number": "35",
                        "celltype": "neural",
                        "x": "1",
                        "y": "2",
                        "z": "280",
                    }
                ],
            )
            proc = subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/core/validate_clean_coordinates.py"),
                    "--input",
                    str(old_id_clean),
                    "--forbid-row-index-cell-id-pattern",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertNotEqual(proc.returncode, 0)

    def test_v2_layout_init_and_dashboard_state(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            config = tmp / "project.yaml"
            config.write_text(
                "\n".join(
                    [
                        "sample_id: TEST",
                        "source_h5ad: /remote/source.h5ad",
                        "columns:",
                        "  cell_id: cell_id",
                        "  slice: sl_number",
                        "  celltype: celltype",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            run_root = tmp / "TEST" / "runs" / "run1"
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/workflow/init_project_run.py"),
                    "--project-config",
                    str(config),
                    "--run-root",
                    str(run_root),
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            for rel in ["states/current_review.yaml", "lineage/state_history.jsonl", "run_manifest.yaml"]:
                self.assertTrue((run_root / rel).exists(), rel)
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/workflow/validate_project_layout.py"),
                    "--run-root",
                    str(run_root),
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            clean = run_root / "steps/001_demo/exports/TEST.step001.clean_coordinates.csv"
            write_csv(clean, [{"cell_id": "a", "slice_id": "s1", "stage": "CS", "chip_id": "chip", "sl_number": "1", "celltype": "neural", "x": "1", "y": "2", "z": "40"}])
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/workflow/update_state_pointer.py"),
                    "--run-root",
                    str(run_root),
                    "--state",
                    "current_review",
                    "--step",
                    "steps/001_demo",
                    "--status",
                    "review",
                    "--summary",
                    "demo state",
                    "--clean-coordinates",
                    str(clean),
                    "--row-count",
                    "1",
                    "--validator-status",
                    "pass",
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            html = run_root / "final/TEST.workflow.html"
            subprocess.run(
                [
                    sys.executable,
                    str(RELEASE / "scripts/workflow/build_workflow_dashboard.py"),
                    "--project-root",
                    str(run_root),
                    "--output-html",
                    str(html),
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            text = html.read_text(encoding="utf-8")
            self.assertIn("Current States", text)
            self.assertIn("TEST.step001.clean_coordinates.csv", text)


if __name__ == "__main__":
    unittest.main()
