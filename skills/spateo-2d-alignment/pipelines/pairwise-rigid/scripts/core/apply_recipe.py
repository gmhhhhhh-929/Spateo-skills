#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from pairwise_core import json_dump, load_json, now_iso, read_point_table, require_columns, write_csv_rows


def mask_for_scope(table, scope: dict[str, object]) -> np.ndarray:
    mask = np.zeros(table.n_obs, dtype=bool)
    labels = {str(x) for x in scope.get("slice_labels") or []}
    if labels:
        mask |= np.array([str(x) in labels for x in table.slices], dtype=bool)
        numeric = {int(x) for x in labels if str(x).isdigit()}
        if numeric:
            mask |= np.array([int(str(x)) in numeric if str(x).isdigit() else False for x in table.slices], dtype=bool)
    celltype = scope.get("celltype_filter")
    if celltype:
        ct_mask = table.celltypes == str(celltype)
        mask = mask & ct_mask if mask.any() else ct_mask
    return mask


def apply_step(coords: np.ndarray, step: dict[str, object], mask: np.ndarray) -> np.ndarray:
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
    else:
        raise SystemExit(f"Unsupported recipe operation: {op}")
    out[mask] = affected
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply a confirmed spatial alignment recipe to a coordinate CSV.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--recipe", required=True)
    parser.add_argument("--x-col", required=True)
    parser.add_argument("--y-col", required=True)
    parser.add_argument("--z-col", required=True)
    parser.add_argument("--slice-col", required=True)
    parser.add_argument("--celltype-col")
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    table = read_point_table(
        Path(args.input).expanduser().resolve(),
        args.x_col,
        args.y_col,
        args.z_col,
        args.slice_col,
        args.celltype_col,
        args.cell_id_col,
    )
    recipe = load_json(Path(args.recipe).expanduser().resolve())
    coords = table.coords.copy()
    for step in recipe.get("steps", []):
        mask = mask_for_scope(table, step.get("scope") or {})
        coords = apply_step(coords, step, mask)
    fields = list(table.rows[0].keys()) if table.rows else []
    require_columns(fields, [args.x_col, args.y_col, args.z_col, args.slice_col])
    for col in ["manual_x", "manual_y", "manual_z"]:
        if col not in fields:
            fields.append(col)
    rows = []
    for i, row in enumerate(table.rows):
        item = dict(row)
        item["manual_x"] = f"{coords[i, 0]:.8g}"
        item["manual_y"] = f"{coords[i, 1]:.8g}"
        item["manual_z"] = f"{coords[i, 2]:.8g}"
        rows.append(item)
    outdir = Path(args.output_dir).expanduser().resolve()
    write_csv_rows(outdir / "corrected_points.csv", rows, fields)
    disp = np.linalg.norm(coords - table.coords, axis=1)
    qc_rows = []
    for label in table.slice_labels:
        mask = table.slices == label
        qc_rows.append({"slice": label, "n": int(mask.sum()), "mean_displacement": float(disp[mask].mean()), "max_displacement": float(disp[mask].max())})
    write_csv_rows(outdir / "qc_by_slice.csv", qc_rows)
    json_dump(outdir / "apply_summary.json", {"version": "spatial_alignment_apply_v1", "created_at": now_iso(), "recipe": str(Path(args.recipe).expanduser().resolve())})
    print(outdir / "corrected_points.csv")


if __name__ == "__main__":
    main()
