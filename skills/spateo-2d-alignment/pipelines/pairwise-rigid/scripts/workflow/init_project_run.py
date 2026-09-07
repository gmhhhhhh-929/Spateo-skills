#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import json
import shutil
from pathlib import Path
from typing import Any


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def load_config(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    try:
        obj = json.loads(text)
    except json.JSONDecodeError:
        try:
            import yaml  # type: ignore
        except Exception as exc:
            obj = load_minimal_yaml(text)
        else:
            obj = yaml.safe_load(text)
    if not isinstance(obj, dict):
        raise SystemExit(f"Project config must be an object: {path}")
    return obj


def parse_scalar(value: str) -> Any:
    value = value.strip()
    if value in {"true", "True"}:
        return True
    if value in {"false", "False"}:
        return False
    if value in {"null", "None", "~"}:
        return None
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def load_minimal_yaml(text: str) -> dict[str, Any]:
    root: dict[str, Any] = {}
    stack: list[tuple[int, dict[str, Any]]] = [(-1, root)]
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip(" "))
        line = raw.strip()
        if ":" not in line or line.startswith("-"):
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]
        if value == "":
            child: dict[str, Any] = {}
            parent[key] = child
            stack.append((indent, child))
        else:
            parent[key] = parse_scalar(value)
    return root


def main() -> None:
    parser = argparse.ArgumentParser(description="Initialize a canonical lightweight spatial alignment run directory.")
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--pipeline-release", default="v0.2.3")
    parser.add_argument("--code-commit", default="")
    parser.add_argument("--allow-existing", action="store_true")
    args = parser.parse_args()

    if args.run_root.exists() and any(args.run_root.iterdir()) and not args.allow_existing:
        raise SystemExit(f"Refusing to initialize non-empty run root without --allow-existing: {args.run_root}")

    config = load_config(args.project_config)
    for name in ["states", "steps", "lineage", "final", "logs"]:
        (args.run_root / name).mkdir(parents=True, exist_ok=True)

    manifest = {
        "version": "lightweight_spatial_alignment_run_v2",
        "created_at": now_iso(),
        "sample_id": config.get("sample_id", ""),
        "run_id": args.run_root.name,
        "run_root": str(args.run_root),
        "project_config": str(args.project_config),
        "pipeline_release": args.pipeline_release,
        "code_commit": args.code_commit,
        "viewer_policy": {
            "manual_review_default": "review_full_points",
            "post_edit_default": "balanced300k",
            "all_points_html": "manual_review_and_final_debug",
        },
        "coordinate_columns": config.get("coordinate_columns", {}),
        "columns": {
            "cell_id": (config.get("columns") or {}).get("cell_id", config.get("cell_id_column", "cell_id")),
            "slice": (config.get("columns") or {}).get("slice", config.get("slice_column", "sl_number")),
            "celltype": (config.get("columns") or {}).get("celltype", config.get("celltype_column", "celltype")),
        },
        "clean_coordinate_schema": ["cell_id", "slice_id", "stage", "chip_id", "sl_number", "celltype", "x", "y", "z"],
        "artifact_policy": {
            "clean_coordinates": "user_facing_and_downstream",
            "audit_points": "reproducibility_and_debug_only",
            "all_points_html": "manual_review_and_final_debug",
        },
    }
    (args.run_root / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (args.run_root / "run_manifest.yaml").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    shutil.copy2(args.project_config, args.run_root / "project_config.snapshot.yaml")
    for rel in ["lineage/accepted_steps.jsonl", "lineage/all_operations.jsonl", "lineage/state_history.jsonl"]:
        path = args.run_root / rel
        if not path.exists():
            path.write_text("", encoding="utf-8")
    for state_name in ["slice_level_baseline", "current_review", "current_accepted", "final"]:
        state_path = args.run_root / "states" / f"{state_name}.yaml"
        if not state_path.exists():
            state = {
                "version": "spatial_alignment_state_pointer_v1",
                "updated_at": now_iso(),
                "state": state_name,
                "status": "empty",
                "step": "",
                "summary": "",
                "artifacts": {},
            }
            state_path.write_text(json.dumps(state, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
