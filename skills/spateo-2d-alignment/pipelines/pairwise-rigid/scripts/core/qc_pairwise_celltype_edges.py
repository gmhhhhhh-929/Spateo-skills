#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from pairwise_core import (
    edge_masks,
    json_dump,
    load_json_or_none,
    maybe_float,
    now_iso,
    pair_metrics,
    read_csv_rows,
    read_point_table,
    write_csv_rows,
)


def flag_reasons(
    row: dict[str, object],
    min_cells: int,
    normalized_jump_threshold: float,
    axis_delta_threshold: float,
    radius_ratio_min: float,
    radius_ratio_max: float,
) -> list[str]:
    reasons = []
    n_fixed = int(row.get("n_fixed") or 0)
    n_moving = int(row.get("n_moving") or 0)
    if n_fixed < min_cells or n_moving < min_cells:
        return ["low_cell_count"]
    jump = maybe_float(row.get("normalized_centroid_jump"))
    angle = maybe_float(row.get("axis_angle_delta_deg"))
    ratio = maybe_float(row.get("radius_ratio"))
    sigma2 = maybe_float(row.get("sigma2"))
    if jump is not None and jump >= normalized_jump_threshold:
        reasons.append("centroid_jump")
    if angle is not None and angle >= axis_delta_threshold:
        reasons.append("axis_angle_jump")
    if ratio is not None and (ratio <= radius_ratio_min or ratio >= radius_ratio_max):
        reasons.append("radius_ratio_shift")
    if sigma2 is not None and bool(row.get("high_sigma2")):
        reasons.append("high_sigma2")
    return reasons


def load_edges(run_dir: Path) -> list[dict[str, object]]:
    index_path = run_dir / "pairwise_edge_index.csv"
    rows = read_csv_rows(index_path) if index_path.exists() else []
    edges = []
    if rows:
        for row in rows:
            edge_json = Path(row.get("edge_dir", "")) / "edge_transform.json" if row.get("edge_dir") else None
            obj = load_json_or_none(edge_json) if edge_json else None
            edge = dict(row)
            if isinstance(obj, dict):
                edge.update(obj)
            edges.append(edge)
        return edges
    edges_root = run_dir / "edges"
    for edge_dir in sorted(edges_root.iterdir()) if edges_root.exists() else []:
        if not edge_dir.is_dir():
            continue
        obj = load_json_or_none(edge_dir / "edge_transform.json")
        if isinstance(obj, dict):
            edges.append(obj)
    return edges


def main() -> None:
    parser = argparse.ArgumentParser(description="QC adjacent pairwise Spateo edges by all-cell and cell-type geometry.")
    parser.add_argument("--pairwise-run-dir", required=True)
    parser.add_argument("--points-csv", required=True)
    parser.add_argument("--x-col", default="stage1_rigid_x")
    parser.add_argument("--y-col", default="stage1_rigid_y")
    parser.add_argument("--z-col", default="z_display")
    parser.add_argument("--slice-col", default="sl_number")
    parser.add_argument("--celltype-col", default="celltype")
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--min-cells-per-celltype-per-slice", type=int, default=30)
    parser.add_argument("--high-sigma2-threshold", type=float, default=0.1)
    parser.add_argument("--normalized-centroid-jump-threshold", type=float, default=1.5)
    parser.add_argument("--axis-angle-delta-threshold", type=float, default=35.0)
    parser.add_argument("--radius-ratio-min", type=float, default=0.5)
    parser.add_argument("--radius-ratio-max", type=float, default=2.0)
    parser.add_argument("--top-k-suspicious-celltypes", type=int, default=50)
    args = parser.parse_args()

    run_dir = Path(args.pairwise_run_dir).expanduser().resolve()
    outdir = Path(args.output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    table = read_point_table(
        Path(args.points_csv).expanduser().resolve(),
        args.x_col,
        args.y_col,
        args.z_col,
        args.slice_col,
        args.celltype_col,
        args.cell_id_col,
    )
    edges = load_edges(run_dir)
    if not edges:
        raise SystemExit(f"No pairwise edges found in {run_dir}")

    edge_rows = []
    celltype_rows = []
    for edge in edges:
        fixed_mask, moving_mask = edge_masks(table, edge)
        sigma2 = maybe_float(edge.get("sigma2"))
        all_metrics = pair_metrics(table.coords[fixed_mask], table.coords[moving_mask])
        high_sigma2 = bool(sigma2 is not None and sigma2 > args.high_sigma2_threshold)
        edge_row = {
            "edge_id": edge.get("edge_id", ""),
            "fixed_slice_id": edge.get("fixed_slice_id", ""),
            "moving_slice_id": edge.get("moving_slice_id", ""),
            "fixed_sl_number": edge.get("fixed_sl_number", ""),
            "moving_sl_number": edge.get("moving_sl_number", ""),
            "status": edge.get("status", ""),
            "sigma2": edge.get("sigma2", ""),
            "gamma": edge.get("gamma", ""),
            "high_sigma2": high_sigma2,
            "n_fixed": int(fixed_mask.sum()),
            "n_moving": int(moving_mask.sum()),
            **all_metrics,
        }
        edge_reasons = flag_reasons(
            edge_row,
            1,
            args.normalized_centroid_jump_threshold,
            args.axis_angle_delta_threshold,
            args.radius_ratio_min,
            args.radius_ratio_max,
        )
        edge_row["reason_flags"] = ";".join([r for r in edge_reasons if r != "low_cell_count"])
        edge_row["priority_score"] = priority_score(edge_row)
        edge_rows.append(edge_row)

        for celltype in table.celltype_labels:
            ct = table.celltypes == celltype
            fm = fixed_mask & ct
            mm = moving_mask & ct
            metrics = pair_metrics(table.coords[fm], table.coords[mm])
            row = {
                "edge_id": edge.get("edge_id", ""),
                "fixed_slice_id": edge.get("fixed_slice_id", ""),
                "moving_slice_id": edge.get("moving_slice_id", ""),
                "fixed_sl_number": edge.get("fixed_sl_number", ""),
                "moving_sl_number": edge.get("moving_sl_number", ""),
                "celltype": celltype,
                "n_fixed": int(fm.sum()),
                "n_moving": int(mm.sum()),
                "sigma2": edge.get("sigma2", ""),
                "gamma": edge.get("gamma", ""),
                "high_sigma2": high_sigma2,
                **metrics,
            }
            reasons = flag_reasons(
                row,
                args.min_cells_per_celltype_per_slice,
                args.normalized_centroid_jump_threshold,
                args.axis_angle_delta_threshold,
                args.radius_ratio_min,
                args.radius_ratio_max,
            )
            row["reason_flags"] = ";".join(reasons)
            row["priority_score"] = priority_score(row)
            celltype_rows.append(row)

    suspicious_edges = [
        row
        for row in sorted(edge_rows, key=lambda r: float(r.get("priority_score") or 0), reverse=True)
        if row.get("reason_flags")
    ]
    suspicious_celltypes = [
        row
        for row in sorted(celltype_rows, key=lambda r: float(r.get("priority_score") or 0), reverse=True)
        if row.get("reason_flags") and "low_cell_count" not in str(row.get("reason_flags"))
    ][: args.top_k_suspicious_celltypes]

    write_csv_rows(outdir / "pairwise_edge_qc.csv", edge_rows)
    write_csv_rows(outdir / "pairwise_celltype_edge_qc.csv", celltype_rows)
    write_csv_rows(outdir / "suspicious_edges.csv", suspicious_edges)
    write_csv_rows(outdir / "suspicious_celltypes.csv", suspicious_celltypes)
    summary = {
        "version": "spateo_pairwise_qc_v1",
        "created_at": now_iso(),
        "pairwise_run_dir": str(run_dir),
        "points_csv": str(Path(args.points_csv).expanduser().resolve()),
        "n_edges": len(edge_rows),
        "n_celltype_edge_rows": len(celltype_rows),
        "n_suspicious_edges": len(suspicious_edges),
        "n_suspicious_celltype_rows": len(suspicious_celltypes),
        "thresholds_are_triage_not_final_decisions": True,
        "thresholds": {
            "min_cells_per_celltype_per_slice": args.min_cells_per_celltype_per_slice,
            "high_sigma2_threshold": args.high_sigma2_threshold,
            "normalized_centroid_jump_threshold": args.normalized_centroid_jump_threshold,
            "axis_angle_delta_threshold": args.axis_angle_delta_threshold,
            "radius_ratio_min": args.radius_ratio_min,
            "radius_ratio_max": args.radius_ratio_max,
        },
        "top_suspicious_edges": [row.get("edge_id") for row in suspicious_edges[:20]],
    }
    json_dump(outdir / "pairwise_qc_summary.json", summary)
    print(outdir / "pairwise_qc_summary.json")


def priority_score(row: dict[str, object]) -> float:
    score = 0.0
    jump = maybe_float(row.get("normalized_centroid_jump"))
    angle = maybe_float(row.get("axis_angle_delta_deg"))
    ratio = maybe_float(row.get("radius_ratio"))
    sigma2 = maybe_float(row.get("sigma2"))
    n_fixed = maybe_float(row.get("n_fixed")) or 0.0
    n_moving = maybe_float(row.get("n_moving")) or 0.0
    if jump is not None and np.isfinite(jump):
        score += min(jump, 10.0) * 2.0
    if angle is not None and np.isfinite(angle):
        score += min(angle / 15.0, 6.0)
    if ratio is not None and np.isfinite(ratio) and ratio > 0:
        score += abs(np.log2(ratio))
    if sigma2 is not None and np.isfinite(sigma2):
        score += min(sigma2 * 20.0, 10.0)
    score += min(np.log10(max(min(n_fixed, n_moving), 1.0)), 4.0) * 0.2
    if "low_cell_count" in str(row.get("reason_flags", "")):
        score *= 0.1
    return float(score)


if __name__ == "__main__":
    main()
