#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


REQUIRED_RUN_DIRS = ["states", "lineage", "steps", "final", "logs"]
REQUIRED_LINEAGE_FILES = ["all_operations.jsonl", "accepted_steps.jsonl", "state_history.jsonl"]
REQUIRED_STATE_FILES = ["slice_level_baseline.yaml", "current_review.yaml", "current_accepted.yaml", "final.yaml"]
STEP_SUBDIRS = ["provenance", "recipes", "qc", "review", "light", "full", "exports", "logs"]


def add_error(report: dict[str, Any], message: str) -> None:
    report["status"] = "fail"
    report.setdefault("errors", []).append(message)


def load_jsonish(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def validate_run(root: Path, strict_steps: bool) -> dict[str, Any]:
    report: dict[str, Any] = {
        "version": "spatial_project_layout_validation_v1",
        "run_root": str(root),
        "status": "pass",
        "errors": [],
        "warnings": [],
        "steps": [],
    }
    if not root.exists():
        add_error(report, f"missing run root: {root}")
        return report
    for dirname in REQUIRED_RUN_DIRS:
        if not (root / dirname).is_dir():
            add_error(report, f"missing run directory: {dirname}")
    for fname in REQUIRED_LINEAGE_FILES:
        if not (root / "lineage" / fname).exists():
            add_error(report, f"missing lineage file: lineage/{fname}")
    for fname in REQUIRED_STATE_FILES:
        path = root / "states" / fname
        if not path.exists():
            add_error(report, f"missing state pointer: states/{fname}")
        else:
            state = load_jsonish(path)
            artifacts = state.get("artifacts") or {}
            clean = (artifacts.get("clean_coordinates") or {}).get("path")
            if clean and not (root / clean).exists() and not Path(clean).is_absolute():
                add_error(report, f"state {fname} clean_coordinates path does not exist: {clean}")
    for step in sorted((root / "steps").glob("[0-9][0-9][0-9]_*")):
        if not step.is_dir():
            continue
        item = {"step": step.name, "status": "pass", "warnings": [], "errors": []}
        for filename in ["step.yaml", "record.json", "decision.md"]:
            if not (step / filename).exists():
                target = item["errors"] if strict_steps and filename in {"step.yaml", "record.json"} else item["warnings"]
                target.append(f"missing {filename}")
        for dirname in STEP_SUBDIRS:
            if not (step / dirname).is_dir():
                target = item["errors"] if strict_steps else item["warnings"]
                target.append(f"missing {dirname}/")
        if item["errors"]:
            report["status"] = "fail"
        report["steps"].append(item)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate v0.2 spatial alignment project/run layout.")
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--strict-steps", action="store_true")
    parser.add_argument("--output-report", type=Path)
    args = parser.parse_args()
    report = validate_run(args.run_root.expanduser().resolve(), args.strict_steps)
    text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output_report:
        args.output_report.parent.mkdir(parents=True, exist_ok=True)
        args.output_report.write_text(text, encoding="utf-8")
    print(text, end="")
    if report["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
