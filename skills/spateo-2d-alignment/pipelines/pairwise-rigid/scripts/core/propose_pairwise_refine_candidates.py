#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from pairwise_core import (
    edge_masks,
    estimate_cloud_axis_rigid_2d,
    json_dump,
    load_json_or_none,
    now_iso,
    pair_metrics,
    read_point_table,
    rigid_transform_about_center,
    sample_points,
    sanitize_token,
    write_csv_rows,
)


def load_edge(run_dir: Path, edge_id: str) -> dict[str, object]:
    path = run_dir / "edges" / edge_id / "edge_transform.json"
    obj = load_json_or_none(path)
    if isinstance(obj, dict):
        return obj
    provenance = load_json_or_none(run_dir / "edges" / edge_id / "provenance.json")
    if isinstance(provenance, dict):
        return provenance
    raise SystemExit(f"Cannot read edge metadata for {edge_id}")


def parse_grid(text: str | None, default: list[float]) -> list[float]:
    if not text:
        return default
    out = []
    for token in text.replace(",", " ").split():
        out.append(float(token))
    return out


def apply_candidate(table, moving_mask: np.ndarray, degrees: float, dx: float, dy: float, center: np.ndarray) -> np.ndarray:
    coords = table.coords.copy()
    coords[moving_mask] = rigid_transform_about_center(coords[moving_mask], degrees, dx, dy, center)
    return coords


def write_corrected_pair_csv(table, coords: np.ndarray, pair_mask: np.ndarray, out: Path) -> None:
    fields = list(table.rows[0].keys()) if table.rows else []
    for col in ["manual_x", "manual_y", "manual_z"]:
        if col not in fields:
            fields.append(col)
    rows = []
    for i, row in enumerate(table.rows):
        if not pair_mask[i]:
            continue
        item = dict(row)
        item["manual_x"] = f"{coords[i, 0]:.8g}"
        item["manual_y"] = f"{coords[i, 1]:.8g}"
        item["manual_z"] = f"{coords[i, 2]:.8g}"
        rows.append(item)
    write_csv_rows(out, rows, fields)


def candidate_recipe(edge: dict[str, object], degrees: float, dx: float, dy: float, center: np.ndarray, reason: str) -> dict[str, object]:
    moving_sl = edge.get("moving_sl_number") or edge.get("moving_slice_id")
    scope = {"slice_column": "sl_number", "slice_labels": [str(moving_sl)]}
    return {
        "version": "spateo_pairwise_refine_candidate_v1",
        "created_at": now_iso(),
        "edge_id": edge.get("edge_id", ""),
        "fixed_slice_id": edge.get("fixed_slice_id", ""),
        "moving_slice_id": edge.get("moving_slice_id", ""),
        "steps": [
            {
                "id": "step_001",
                "type": "pairwise_moving_slice_rigid_transform",
                "operation": "rotate_z",
                "scope": scope,
                "parameters": {"degrees": degrees},
                "center": {"mode": "manual", "x": float(center[0]), "y": float(center[1]), "z": 0.0},
                "reason": reason,
            },
            {
                "id": "step_002",
                "type": "pairwise_moving_slice_rigid_transform",
                "operation": "translate",
                "scope": scope,
                "parameters": {"dx": dx, "dy": dy, "dz": 0.0},
                "center": {"mode": "manual", "x": float(center[0]), "y": float(center[1]), "z": 0.0},
                "reason": reason,
            }
        ],
    }


def score_candidate(table, edge, fixed_mask, moving_mask, target_mask, corrected: np.ndarray, min_guardrail_cells: int) -> dict[str, object]:
    before_target = pair_metrics(table.coords[fixed_mask & target_mask], table.coords[moving_mask & target_mask])
    after_target = pair_metrics(corrected[fixed_mask & target_mask], corrected[moving_mask & target_mask])
    before_all = pair_metrics(table.coords[fixed_mask], table.coords[moving_mask])
    after_all = pair_metrics(corrected[fixed_mask], corrected[moving_mask])

    target_before = float(before_target.get("normalized_centroid_jump") or math.nan)
    target_after = float(after_target.get("normalized_centroid_jump") or math.nan)
    target_improvement = target_before - target_after if np.isfinite(target_before) and np.isfinite(target_after) else 0.0
    all_before = float(before_all.get("normalized_centroid_jump") or math.nan)
    all_after = float(after_all.get("normalized_centroid_jump") or math.nan)
    all_degradation = all_after - all_before if np.isfinite(all_before) and np.isfinite(all_after) else 0.0

    worst_other_degradation = 0.0
    affected_other = []
    for ct in table.celltype_labels:
        ct_mask = table.celltypes == ct
        if ct_mask[target_mask].any():
            continue
        if int((fixed_mask & ct_mask).sum()) < min_guardrail_cells or int((moving_mask & ct_mask).sum()) < min_guardrail_cells:
            continue
        before = pair_metrics(table.coords[fixed_mask & ct_mask], table.coords[moving_mask & ct_mask])
        after = pair_metrics(corrected[fixed_mask & ct_mask], corrected[moving_mask & ct_mask])
        b = float(before.get("normalized_centroid_jump") or math.nan)
        a = float(after.get("normalized_centroid_jump") or math.nan)
        degradation = a - b if np.isfinite(a) and np.isfinite(b) else 0.0
        worst_other_degradation = max(worst_other_degradation, degradation)
        if degradation > 0.5:
            affected_other.append({"celltype": ct, "degradation": degradation})

    rank_score = -target_improvement + max(0.0, all_degradation) * 0.5 + max(0.0, worst_other_degradation) * 1.5
    if target_improvement <= 0:
        recommendation = "reject_no_target_improvement"
    elif all_degradation > 1.0 or worst_other_degradation > 1.0:
        recommendation = "review_guardrail_degradation"
    elif target_improvement > 0.75:
        recommendation = "promising_candidate"
    else:
        recommendation = "weak_candidate"
    return {
        "target_before_normalized_jump": target_before,
        "target_after_normalized_jump": target_after,
        "target_improvement": target_improvement,
        "all_before_normalized_jump": all_before,
        "all_after_normalized_jump": all_after,
        "all_cell_degradation": all_degradation,
        "other_celltype_worst_degradation": worst_other_degradation,
        "affected_other_celltypes": affected_other[:20],
        "rank_score": rank_score,
        "recommendation": recommendation,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate rigid refine candidates for one suspect pairwise edge.")
    parser.add_argument("--pairwise-run-dir", required=True)
    parser.add_argument("--points-csv", required=True)
    parser.add_argument("--edge-id", required=True)
    parser.add_argument("--target-celltype", required=True)
    parser.add_argument("--x-col", default="stage1_rigid_x")
    parser.add_argument("--y-col", default="stage1_rigid_y")
    parser.add_argument("--z-col", default="z_display")
    parser.add_argument("--slice-col", default="sl_number")
    parser.add_argument("--celltype-col", default="celltype")
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--rotation-offset-grid", default="-10 -5 0 5 10")
    parser.add_argument("--dx-offset-grid", default="-1500 0 1500")
    parser.add_argument("--dy-offset-grid", default="-1500 0 1500")
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--max-target-points", type=int, default=5000)
    parser.add_argument("--min-target-cells-per-slice", type=int, default=20)
    parser.add_argument("--min-guardrail-cells-per-slice", type=int, default=30)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--rng-seed", type=int, default=20260616)
    args = parser.parse_args()

    run_dir = Path(args.pairwise_run_dir).expanduser().resolve()
    outdir = Path(args.output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    edge = load_edge(run_dir, args.edge_id)
    table = read_point_table(
        Path(args.points_csv).expanduser().resolve(),
        args.x_col,
        args.y_col,
        args.z_col,
        args.slice_col,
        args.celltype_col,
        args.cell_id_col,
    )
    fixed_mask, moving_mask = edge_masks(table, edge)
    target_mask = table.celltypes == args.target_celltype
    if int((fixed_mask & target_mask).sum()) < args.min_target_cells_per_slice:
        raise SystemExit("Target cell type has too few cells in fixed slice.")
    if int((moving_mask & target_mask).sum()) < args.min_target_cells_per_slice:
        raise SystemExit("Target cell type has too few cells in moving slice.")

    rng = np.random.default_rng(args.rng_seed)
    fixed_target = sample_points(table.coords[fixed_mask & target_mask], args.max_target_points, rng)
    moving_target = sample_points(table.coords[moving_mask & target_mask], args.max_target_points, rng)
    base_deg, base_dx, base_dy = estimate_cloud_axis_rigid_2d(moving_target, fixed_target)
    center = moving_target.mean(axis=0)
    rot_offsets = parse_grid(args.rotation_offset_grid, [0.0])
    dx_offsets = parse_grid(args.dx_offset_grid, [0.0])
    dy_offsets = parse_grid(args.dy_offset_grid, [0.0])
    pair_mask = fixed_mask | moving_mask

    candidates = []
    candidate_payloads = []
    seen = set()
    for drot in rot_offsets:
        for ddx in dx_offsets:
            for ddy in dy_offsets:
                deg = base_deg + drot
                dx = base_dx + ddx
                dy = base_dy + ddy
                key = (round(deg, 6), round(dx, 3), round(dy, 3))
                if key in seen:
                    continue
                seen.add(key)
                corrected = apply_candidate(table, moving_mask, deg, dx, dy, center)
                score = score_candidate(table, edge, fixed_mask, moving_mask, target_mask, corrected, args.min_guardrail_cells_per_slice)
                candidate_payloads.append((score["rank_score"], deg, dx, dy, corrected, score))

    candidate_payloads.sort(key=lambda item: float(item[0]))
    for i, (_, deg, dx, dy, corrected, score) in enumerate(candidate_payloads[: args.top_k], start=1):
        candidate_id = f"candidate_{i:03d}"
        cand_dir = outdir / "candidates" / candidate_id
        cand_dir.mkdir(parents=True, exist_ok=True)
        recipe = candidate_recipe(edge, deg, dx, dy, center, f"pairwise refine guided by {args.target_celltype}")
        json_dump(cand_dir / "transform_recipe.json", recipe)
        json_dump(cand_dir / "candidate_score.json", score)
        write_corrected_pair_csv(table, corrected, pair_mask, cand_dir / "corrected_pair_points.csv")
        row = {
            "candidate_id": candidate_id,
            "edge_id": args.edge_id,
            "fixed_slice_id": edge.get("fixed_slice_id", ""),
            "moving_slice_id": edge.get("moving_slice_id", ""),
            "fixed_sl_number": edge.get("fixed_sl_number", ""),
            "moving_sl_number": edge.get("moving_sl_number", ""),
            "pair_start_order": edge.get("pair_start_order", ""),
            "target_celltype": args.target_celltype,
            "rotation_deg": deg,
            "dx": dx,
            "dy": dy,
            "center_x": float(center[0]),
            "center_y": float(center[1]),
            **{k: v for k, v in score.items() if k != "affected_other_celltypes"},
            "candidate_dir": str(cand_dir),
            "transform_recipe": str(cand_dir / "transform_recipe.json"),
            "corrected_pair_points": str(cand_dir / "corrected_pair_points.csv"),
        }
        candidates.append(row)

    write_csv_rows(outdir / "candidate_ranking.csv", candidates)
    json_dump(
        outdir / "inputs.json",
        {
            "version": "spateo_pairwise_refine_candidates_v1",
            "created_at": now_iso(),
            "pairwise_run_dir": str(run_dir),
            "points_csv": str(Path(args.points_csv).expanduser().resolve()),
            "edge_id": args.edge_id,
            "target_celltype": args.target_celltype,
            "base_estimate": {"rotation_deg": base_deg, "dx": base_dx, "dy": base_dy, "center": center.tolist()},
            "notes": [
                "Candidates are not final coordinates.",
                "The base estimate uses target-celltype cloud centroid and principal-axis continuity; it does not assume one-to-one cell correspondence.",
                "Each candidate applies the rigid transform to the whole moving slice for guardrail scoring.",
            ],
        },
    )
    if candidates:
        json_dump(outdir / "recommended_candidate.json", candidates[0])
    print(outdir / "candidate_ranking.csv")


if __name__ == "__main__":
    main()
