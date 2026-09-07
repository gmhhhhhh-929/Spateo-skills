#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

from pairwise_core import (
    json_dump,
    load_json,
    natural_key,
    now_iso,
    parse_sl_number,
    read_point_table,
    require_columns,
    to_jsonable,
    write_csv_rows,
)


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(block_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def file_manifest(path: Path, role: str, hash_file: bool = True) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    item: dict[str, Any] = {
        "role": role,
        "path": str(resolved),
        "exists": resolved.exists(),
    }
    if resolved.exists():
        stat = resolved.stat()
        item.update(
            {
                "size_bytes": stat.st_size,
                "mtime": stat.st_mtime,
                "mtime_iso": now_iso_from_timestamp(stat.st_mtime),
            }
        )
        if hash_file:
            item["sha256"] = sha256_file(resolved)
    return item


def now_iso_from_timestamp(timestamp: float) -> str:
    import datetime as dt

    return dt.datetime.fromtimestamp(timestamp, dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def hash_jsonable(obj: Any) -> str:
    payload = json.dumps(to_jsonable(obj), ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def coord_hash(coords: np.ndarray) -> str:
    arr = np.asarray(coords, dtype=np.float64)
    return hashlib.sha256(arr.tobytes(order="C")).hexdigest()


def load_operations(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[Path]]:
    if args.recipe_chain:
        chain_path = Path(args.recipe_chain).expanduser().resolve()
        chain = load_json(chain_path)
        ops = chain.get("operations")
        if not isinstance(ops, list):
            raise SystemExit("recipe_chain.json must contain an operations list.")
        return [dict(op) for op in ops], [chain_path]

    operation_paths = [Path(p).expanduser().resolve() for p in args.operation]
    operations = []
    for i, path in enumerate(operation_paths, start=1):
        recipe = load_json(path)
        operations.append(operation_from_recipe(recipe, path, i))
    return operations, operation_paths


def operation_from_recipe(recipe: dict[str, Any], path: Path, index: int) -> dict[str, Any]:
    steps = recipe.get("steps")
    if not isinstance(steps, list) or not steps:
        raise SystemExit(f"Recipe has no steps: {path}")
    op_type = str(recipe.get("operation_type") or infer_operation_type(recipe, path))
    return {
        "operation_id": f"op_{index:03d}",
        "operation_type": op_type,
        "source_recipe_path": str(path),
        "source_recipe_sha256": sha256_file(path),
        "edge_id": recipe.get("edge_id") or infer_edge_id(recipe, path),
        "fixed_slice": recipe.get("fixed_slice") or recipe.get("fixed_sl_number") or "",
        "moving_slice": recipe.get("moving_slice") or recipe.get("moving_sl_number") or "",
        "affected_slices": infer_affected_slices(recipe),
        "coordinate_input_columns": recipe.get("coordinate_source") or {},
        "coordinate_output_columns": {"x": "manual_x", "y": "manual_y", "z": "manual_z"},
        "transform": {"steps": copy.deepcopy(steps), "transform_hash": hash_jsonable(steps)},
        "expression_mode": recipe.get("expression_mode") or infer_expression_mode(recipe),
        "sigma": {
            "sigma2_init_scale": recipe.get("sigma2_init_scale", ""),
            "sigma2_end": recipe.get("sigma2_end", ""),
            "sigma2": recipe.get("sigma2", ""),
        },
        "roi_or_drop": recipe.get("roi_or_drop") or infer_roi_or_drop(recipe),
        "timestamp": recipe.get("created_at") or now_iso(),
        "user_confirmation_note": recipe.get("user_confirmation_note") or recipe.get("selected_candidate") or "",
        "source_run_directory": str(Path(str(recipe.get("selected_pairwise_run") or path)).expanduser().resolve().parent)
        if recipe.get("selected_pairwise_run")
        else str(path.parent),
        "raw_recipe": recipe,
    }


def infer_operation_type(recipe: dict[str, Any], path: Path) -> str:
    text = json.dumps(recipe, ensure_ascii=False).lower() + " " + str(path).lower()
    if "drop" in text or "roi" in text:
        return "roi_drop_refine"
    if "spatial_only" in text:
        return "spatial_only_rescue"
    if "edge_transform" in text or "pairwise" in text:
        return "pairwise_transform"
    return "manual_recipe"


def infer_edge_id(recipe: dict[str, Any], path: Path) -> str:
    for key in ["edge_id", "selected_candidate", "selected_pairwise_run", "transform_derivation"]:
        value = str(recipe.get(key, ""))
        match = __import__("re").search(r"SL\s*0*(\d+)[_-]+SL\s*0*(\d+)", value, flags=__import__("re").IGNORECASE)
        if match:
            return f"SL{int(match.group(1))}__SL{int(match.group(2))}"
    match = __import__("re").search(r"SL\s*0*(\d+)[_-]+SL\s*0*(\d+)", str(path), flags=__import__("re").IGNORECASE)
    if match:
        return f"SL{int(match.group(1))}__SL{int(match.group(2))}"
    return ""


def infer_affected_slices(recipe: dict[str, Any]) -> list[str]:
    labels: set[str] = set()
    for step in recipe.get("steps") or []:
        scope = step.get("scope") or {}
        for label in scope.get("slice_labels") or []:
            labels.add(str(label))
    return sorted(labels, key=natural_key)


def infer_expression_mode(recipe: dict[str, Any]) -> str:
    text = json.dumps(recipe, ensure_ascii=False).lower()
    if "spatial_only" in text:
        return "spatial_only"
    return str(recipe.get("expression_mode") or "normal")


def infer_roi_or_drop(recipe: dict[str, Any]) -> dict[str, Any]:
    text = json.dumps(recipe, ensure_ascii=False).lower()
    out: dict[str, Any] = {}
    if "drop" in text:
        out["mode"] = "drop"
    if "roi" in text:
        out["uses_roi"] = True
    derivation = recipe.get("transform_derivation")
    if derivation:
        out["description"] = derivation
    return out


def mask_for_scope(table, scope: dict[str, Any]) -> np.ndarray:
    mask = np.zeros(table.n_obs, dtype=bool)
    labels = {str(x) for x in scope.get("slice_labels") or []}
    if labels:
        numeric = {parse_sl_number(x) for x in labels}
        numeric.discard(None)
        for i, value in enumerate(table.slices):
            sl = parse_sl_number(value)
            if str(value) in labels or (sl is not None and sl in numeric):
                mask[i] = True
    celltype = scope.get("celltype_filter")
    if celltype:
        ct_mask = table.celltypes == str(celltype)
        mask = mask & ct_mask if mask.any() else ct_mask
    return mask


def apply_transform_step(coords: np.ndarray, step: dict[str, Any], mask: np.ndarray) -> np.ndarray:
    out = coords.copy()
    if not mask.any():
        return out
    op = str(step.get("operation", ""))
    params = step.get("parameters") or {}
    center = step.get("center") or {}
    affected = out[mask].copy()
    if op == "translate":
        affected[:, 0] += float(params.get("dx", 0))
        affected[:, 1] += float(params.get("dy", 0))
        affected[:, 2] += float(params.get("dz", 0))
    elif op == "rotate_z":
        cx = float(center.get("x", 0))
        cy = float(center.get("y", 0))
        theta = math.radians(float(params.get("degrees", 0)))
        rot = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]], dtype=float)
        affected[:, :2] = (affected[:, :2] - np.array([cx, cy])) @ rot.T + np.array([cx, cy])
    elif op == "scale":
        c = np.array([float(center.get("x", 0)), float(center.get("y", 0)), float(center.get("z", 0))], dtype=float)
        affected = (affected - c) * float(params.get("scale", 1.0)) + c
    elif op in {"affine_2d", "matrix_2d"}:
        matrix = np.asarray(params.get("matrix"), dtype=float)
        if matrix.shape != (2, 2):
            raise SystemExit("affine_2d/matrix_2d requires parameters.matrix with shape 2x2.")
        t = np.asarray(params.get("translation", [0, 0]), dtype=float)
        affected[:, :2] = affected[:, :2] @ matrix.T + t
    else:
        raise SystemExit(f"Unsupported recipe operation: {op}")
    out[mask] = affected
    return out


def apply_operation(table, coords: np.ndarray, operation: dict[str, Any]) -> tuple[np.ndarray, dict[str, Any]]:
    before = coords
    after = coords.copy()
    step_summaries = []
    for step in ((operation.get("transform") or {}).get("steps") or []):
        mask = mask_for_scope(table, step.get("scope") or {})
        before_step = after
        after = apply_transform_step(after, step, mask)
        changed = np.linalg.norm(after - before_step, axis=1) > 1e-9
        step_summaries.append(
            {
                "step_id": step.get("id", ""),
                "operation": step.get("operation", ""),
                "scope": step.get("scope") or {},
                "points_in_scope": int(mask.sum()),
                "points_changed": int(changed.sum()),
                "reason": step.get("reason", ""),
            }
        )
    changed = np.linalg.norm(after - before, axis=1) > 1e-9
    disp = np.linalg.norm(after - before, axis=1)
    return after, {
        "points_affected": int(changed.sum()),
        "mean_displacement": float(disp[changed].mean()) if changed.any() else 0.0,
        "max_displacement": float(disp[changed].max()) if changed.any() else 0.0,
        "step_summaries": step_summaries,
    }


def normalize_operations(operations: list[dict[str, Any]], args: argparse.Namespace) -> list[dict[str, Any]]:
    out = []
    for i, op in enumerate(operations, start=1):
        item = copy.deepcopy(op)
        item.setdefault("operation_id", f"op_{i:03d}")
        item.setdefault("operation_type", "manual_recipe")
        item.setdefault("coordinate_output_columns", {"x": args.out_x_col, "y": args.out_y_col, "z": args.out_z_col})
        item["step_order"] = i
        out.append(item)
    return out


def write_outputs(args: argparse.Namespace, table, operations: list[dict[str, Any]], source_paths: list[Path]) -> None:
    outdir = Path(args.output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    coords = table.coords.copy()
    initial_hash = coord_hash(coords)
    history_rows = []
    per_step = []
    for i, op in enumerate(operations, start=1):
        coords, summary = apply_operation(table, coords, op)
        affected_slices = op.get("affected_slices") or []
        if not affected_slices:
            affected_slices = sorted(
                {
                    str(label)
                    for step in ((op.get("transform") or {}).get("steps") or [])
                    for label in ((step.get("scope") or {}).get("slice_labels") or [])
                },
                key=natural_key,
            )
        history_rows.append(
            {
                "step_order": i,
                "operation_id": op.get("operation_id", ""),
                "edge_id": op.get("edge_id", ""),
                "operation_type": op.get("operation_type", ""),
                "affected_slices": ";".join(map(str, affected_slices)),
                "points_affected": summary["points_affected"],
                "coordinate_columns_written": f"{args.out_x_col};{args.out_y_col};{args.out_z_col}",
                "source_recipe_path": op.get("source_recipe_path", ""),
                "source_run_directory": op.get("source_run_directory", ""),
                "expression_mode": op.get("expression_mode", ""),
                "sigma2_init_scale": (op.get("sigma") or {}).get("sigma2_init_scale", ""),
                "sigma2_end": (op.get("sigma") or {}).get("sigma2_end", ""),
                "sigma2": (op.get("sigma") or {}).get("sigma2", ""),
                "mean_displacement": summary["mean_displacement"],
                "max_displacement": summary["max_displacement"],
                "user_confirmation_note": op.get("user_confirmation_note", ""),
            }
        )
        per_step.append({"operation_id": op.get("operation_id", ""), **summary})

    fields = list(table.rows[0].keys()) if table.rows else []
    require_columns(fields, [args.x_col, args.y_col, args.z_col, args.slice_col])
    for col in [args.out_x_col, args.out_y_col, args.out_z_col]:
        if col not in fields:
            fields.append(col)
    rows = []
    for i, row in enumerate(table.rows):
        item = dict(row)
        item[args.out_x_col] = f"{coords[i, 0]:.8g}"
        item[args.out_y_col] = f"{coords[i, 1]:.8g}"
        item[args.out_z_col] = f"{coords[i, 2]:.8g}"
        rows.append(item)
    audit_path = outdir / "audit_points.csv"
    legacy_path = outdir / "corrected_points.csv"
    write_csv_rows(audit_path, rows, fields)
    if legacy_path.exists():
        legacy_path.unlink()
    try:
        os.link(audit_path, legacy_path)
    except OSError:
        import shutil

        shutil.copy2(audit_path, legacy_path)
    write_csv_rows(outdir / "correction_history.csv", history_rows)

    final_hash = coord_hash(coords)
    chain = {
        "version": "spatial_alignment_recipe_chain_v1",
        "created_at": now_iso(),
        "source_input": str(Path(args.input).expanduser().resolve()),
        "source_input_sha256": sha256_file(Path(args.input).expanduser().resolve()),
        "source_coordinate_columns": {"x": args.x_col, "y": args.y_col, "z": args.z_col},
        "final_coordinate_columns": {"x": args.out_x_col, "y": args.out_y_col, "z": args.out_z_col},
        "slice_column": args.slice_col,
        "celltype_column": args.celltype_col or "",
        "cell_id_column": args.cell_id_col or "",
        "operations": operations,
    }
    chain_hash = hash_jsonable(chain)
    chain["recipe_chain_sha256"] = chain_hash
    json_dump(outdir / "recipe_chain.json", chain)

    manifest = {
        "version": "spatial_alignment_input_manifest_v1",
        "created_at": now_iso(),
        "files": [file_manifest(Path(args.input), "source_coordinates")]
        + [file_manifest(path, "operation_recipe") for path in source_paths],
    }
    json_dump(outdir / "input_file_manifest.json", manifest)

    summary = {
        "version": "spatial_alignment_compose_summary_v1",
        "created_at": now_iso(),
        "total_points": table.n_obs,
        "final_slices": table.slice_labels,
        "n_operations": len(operations),
        "points_modified_per_step": per_step,
        "final_coordinate_columns": {"x": args.out_x_col, "y": args.out_y_col, "z": args.out_z_col},
        "source_input_hash": sha256_file(Path(args.input).expanduser().resolve()),
        "initial_coordinate_hash": initial_hash,
        "final_coordinate_hash": final_hash,
        "recipe_chain_hash": chain_hash,
        "replay_checked": False,
    }
    if args.replay_check:
        replay_coords = replay(table.coords.copy(), table, operations)
        max_abs_delta = float(np.max(np.abs(replay_coords - coords))) if coords.size else 0.0
        summary["replay_checked"] = True
        summary["replay_max_abs_delta"] = max_abs_delta
        summary["replay_coordinate_hash"] = coord_hash(replay_coords)
        if max_abs_delta > args.replay_atol:
            json_dump(outdir / "apply_summary.json", summary)
            raise SystemExit(f"Replay check failed: max_abs_delta={max_abs_delta} > {args.replay_atol}")
    clean_output = Path(args.clean_output).expanduser().resolve() if args.clean_output else outdir / "clean_coordinates.csv"
    clean_manifest = Path(args.clean_manifest).expanduser().resolve() if args.clean_manifest else outdir / "clean_export_manifest.json"
    changed_cells = Path(args.changed_cells_output).expanduser().resolve() if args.changed_cells_output else outdir / "changed_cells.csv"
    exporter = Path(__file__).with_name("export_clean_coordinates.py")
    export_cmd = [
        sys.executable,
        str(exporter),
        "--input",
        str(audit_path),
        "--output",
        str(clean_output),
        "--manifest",
        str(clean_manifest),
        "--changed-cells-output",
        str(changed_cells),
        "--cell-id-col",
        args.cell_id_col,
        "--sl-number-col",
        args.slice_col,
        "--celltype-col",
        args.celltype_col or "celltype",
        "--x-col",
        args.out_x_col,
        "--y-col",
        args.out_y_col,
        "--z-col",
        args.out_z_col,
        "--baseline-x-col",
        args.x_col,
        "--baseline-y-col",
        args.y_col,
        "--baseline-z-col",
        args.z_col,
        "--code-release",
        args.code_release,
        "--entrypoint-sha256",
        args.entrypoint_sha256,
    ]
    if args.cell_id_map:
        export_cmd.extend(
            [
                "--cell-id-map",
                args.cell_id_map,
                "--cell-id-map-key-col",
                args.cell_id_map_key_col,
                "--cell-id-map-value-col",
                args.cell_id_map_value_col,
                "--clean-cell-id-source",
                args.clean_cell_id_source or "short_slice_cellid",
            ]
        )
    subprocess.run(export_cmd, check=True, stdout=subprocess.PIPE, text=True)
    summary["audit_points_csv"] = str(audit_path)
    summary["legacy_corrected_points_csv"] = str(legacy_path)
    summary["clean_coordinates_csv"] = str(clean_output)
    summary["clean_export_manifest"] = str(clean_manifest)
    summary["changed_cells_csv"] = str(changed_cells)
    json_dump(outdir / "apply_summary.json", summary)


def replay(initial_coords: np.ndarray, table, operations: list[dict[str, Any]]) -> np.ndarray:
    coords = initial_coords.copy()
    for op in operations:
        coords, _ = apply_operation(table, coords, op)
    return coords


def main() -> None:
    parser = argparse.ArgumentParser(description="Compose confirmed spatial alignment recipes with replayable edit history.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--operation", action="append", default=[], help="Recipe JSON path. Repeat in execution order.")
    parser.add_argument("--recipe-chain", help="Replay an existing recipe_chain.json instead of --operation.")
    parser.add_argument("--x-col", required=True)
    parser.add_argument("--y-col", required=True)
    parser.add_argument("--z-col", required=True)
    parser.add_argument("--slice-col", required=True)
    parser.add_argument("--celltype-col")
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--cell-id-map", default="", help="Optional CSV mapping source workflow cell ids to clean output cell ids.")
    parser.add_argument("--cell-id-map-key-col", default="workflow_cell_id")
    parser.add_argument("--cell-id-map-value-col", default="cell_id")
    parser.add_argument("--clean-cell-id-source", default="")
    parser.add_argument("--out-x-col", default="manual_x")
    parser.add_argument("--out-y-col", default="manual_y")
    parser.add_argument("--out-z-col", default="manual_z")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--clean-output", help="Clean user-facing CSV path. Defaults to OUTPUT_DIR/clean_coordinates.csv.")
    parser.add_argument("--clean-manifest", help="Clean export manifest path. Defaults to OUTPUT_DIR/clean_export_manifest.json.")
    parser.add_argument("--changed-cells-output", help="Changed-cell audit CSV path. Defaults to OUTPUT_DIR/changed_cells.csv.")
    parser.add_argument("--code-release", default="v0.2.3")
    parser.add_argument("--entrypoint-sha256", default="")
    parser.add_argument("--replay-check", action="store_true")
    parser.add_argument("--replay-atol", type=float, default=1e-8)
    args = parser.parse_args()

    if bool(args.recipe_chain) == bool(args.operation):
        raise SystemExit("Use exactly one of --recipe-chain or one or more --operation.")

    table = read_point_table(
        Path(args.input).expanduser().resolve(),
        args.x_col,
        args.y_col,
        args.z_col,
        args.slice_col,
        args.celltype_col,
        args.cell_id_col,
    )
    operations, source_paths = load_operations(args)
    operations = normalize_operations(operations, args)
    write_outputs(args, table, operations, source_paths)
    print(Path(args.output_dir).expanduser().resolve() / "audit_points.csv")


if __name__ == "__main__":
    main()
