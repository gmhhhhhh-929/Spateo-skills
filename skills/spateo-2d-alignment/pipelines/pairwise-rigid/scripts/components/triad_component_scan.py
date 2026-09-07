#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


SLICE_SORT_KEY = lambda x: int(float(str(x).replace("SL", "")))


def read_component_summary(path: Path, min_points: int) -> dict[str, list[dict[str, Any]]]:
    by_slice: dict[str, list[dict[str, Any]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "slice",
            "component_id",
            "component_rank",
            "points",
            "fraction",
            "centroid_x",
            "centroid_y",
            "min_x",
            "max_x",
            "min_y",
            "max_y",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Component summary missing columns: {sorted(missing)}")
        for row in reader:
            points = int(float(row["points"]))
            if int(float(row["component_rank"])) <= 0 or points < min_points:
                continue
            sx = str(int(float(row["slice"])))
            width = float(row["max_x"]) - float(row["min_x"])
            height = float(row["max_y"]) - float(row["min_y"])
            item = {
                "slice": sx,
                "component_id": str(row["component_id"]),
                "component_rank": int(float(row["component_rank"])),
                "points": points,
                "fraction": float(row["fraction"]),
                "centroid_x": float(row["centroid_x"]),
                "centroid_y": float(row["centroid_y"]),
                "min_x": float(row["min_x"]),
                "max_x": float(row["max_x"]),
                "min_y": float(row["min_y"]),
                "max_y": float(row["max_y"]),
                "width": width,
                "height": height,
                "area": max(width * height, 1.0),
                "diag": math.hypot(width, height),
                "aspect": width / max(height, 1e-6),
                "top_celltypes": row.get("top_celltypes", ""),
            }
            by_slice[sx].append(item)
    for items in by_slice.values():
        items.sort(key=lambda v: (v["component_rank"], -v["points"]))
    return by_slice


def dist(a: dict[str, Any], b: dict[str, Any]) -> float:
    return math.hypot(float(a["centroid_x"]) - float(b["centroid_x"]), float(a["centroid_y"]) - float(b["centroid_y"]))


def signed_delta(a: dict[str, Any], b: dict[str, Any]) -> tuple[float, float]:
    return float(b["centroid_x"]) - float(a["centroid_x"]), float(b["centroid_y"]) - float(a["centroid_y"])


def scale_for(a: dict[str, Any], b: dict[str, Any], c: dict[str, Any]) -> float:
    return max(float(a["diag"]), float(b["diag"]), float(c["diag"]), 100.0)


def size_similarity(a: dict[str, Any], b: dict[str, Any]) -> float:
    ratio = min(float(a["points"]), float(b["points"])) / max(float(a["points"]), float(b["points"]), 1.0)
    area_ratio = min(float(a["area"]), float(b["area"])) / max(float(a["area"]), float(b["area"]), 1.0)
    return 0.65 * ratio + 0.35 * area_ratio


def shape_delta(a: dict[str, Any], b: dict[str, Any]) -> float:
    aspect = abs(math.log(max(float(a["aspect"]), 1e-6) / max(float(b["aspect"]), 1e-6)))
    diag = abs(math.log(max(float(a["diag"]), 1e-6) / max(float(b["diag"]), 1e-6)))
    return aspect + diag


def choose_rank(slice_items: list[dict[str, Any]], rank: int) -> dict[str, Any] | None:
    for item in slice_items:
        if int(item["component_rank"]) == rank:
            return item
    return None


def classify(row: dict[str, Any], args: argparse.Namespace) -> tuple[str, str]:
    score = float(row["triad_outlier_score"])
    mid_outlier = row["issue_type"] == "mid_outlier"
    edge_issue = row["issue_type"] == "edge_issue"
    if mid_outlier and score >= args.flag_threshold and int(row["mid_points"]) >= args.min_flag_points:
        return "flagged", "mid component deviates from the prev/next consensus"
    if mid_outlier and score >= args.watch_threshold:
        return "watch", "possible mid component deviation"
    if edge_issue and score >= args.watch_threshold:
        return "watch", "one adjacent edge changes more than the triad consensus"
    return "normal", "triad geometry is within conservative thresholds"


def scan(by_slice: dict[str, list[dict[str, Any]]], args: argparse.Namespace) -> list[dict[str, Any]]:
    slices = sorted(by_slice, key=SLICE_SORT_KEY)
    rows: list[dict[str, Any]] = []
    for i in range(1, len(slices) - 1):
        prev_sl, mid_sl, next_sl = slices[i - 1], slices[i], slices[i + 1]
        max_rank = max(
            [item["component_rank"] for item in by_slice[prev_sl][: args.max_rank]]
            + [item["component_rank"] for item in by_slice[mid_sl][: args.max_rank]]
            + [item["component_rank"] for item in by_slice[next_sl][: args.max_rank]]
        )
        for rank in range(1, min(max_rank, args.max_rank) + 1):
            prev = choose_rank(by_slice[prev_sl], rank)
            mid = choose_rank(by_slice[mid_sl], rank)
            nxt = choose_rank(by_slice[next_sl], rank)
            if prev is None or mid is None or nxt is None:
                continue
            scale = scale_for(prev, mid, nxt)
            prev_next = dist(prev, nxt)
            mid_expected_x = (float(prev["centroid_x"]) + float(nxt["centroid_x"])) / 2.0
            mid_expected_y = (float(prev["centroid_y"]) + float(nxt["centroid_y"])) / 2.0
            mid_consensus_dist = math.hypot(float(mid["centroid_x"]) - mid_expected_x, float(mid["centroid_y"]) - mid_expected_y)
            prev_mid = dist(prev, mid)
            mid_next = dist(mid, nxt)
            imbalance = abs(prev_mid - mid_next)
            jump_sum = prev_mid + mid_next
            continuity_gap = max(0.0, jump_sum - max(prev_next, 1.0))
            prev_mid_dx, prev_mid_dy = signed_delta(prev, mid)
            mid_next_dx, mid_next_dy = signed_delta(mid, nxt)
            reversal = 0.0
            norm_a = math.hypot(prev_mid_dx, prev_mid_dy)
            norm_b = math.hypot(mid_next_dx, mid_next_dy)
            if norm_a > 1e-6 and norm_b > 1e-6:
                reversal = max(0.0, -((prev_mid_dx * mid_next_dx + prev_mid_dy * mid_next_dy) / (norm_a * norm_b)))
            size_stability = min(size_similarity(prev, nxt), size_similarity(prev, mid), size_similarity(mid, nxt))
            shape_change = max(shape_delta(prev, mid), shape_delta(mid, nxt))
            points_weight = min(1.0, math.log10(max(float(mid["points"]), 10.0)) / 4.0)
            normalized_mid_offset = mid_consensus_dist / scale
            normalized_gap = continuity_gap / scale
            normalized_prev_next = prev_next / scale
            issue_type = "normal"
            if normalized_mid_offset >= args.mid_offset_threshold and normalized_prev_next <= args.neighbor_stability_threshold:
                issue_type = "mid_outlier"
            elif max(prev_mid, mid_next) / scale >= args.edge_issue_threshold and abs(prev_mid - mid_next) / scale >= args.edge_imbalance_threshold:
                issue_type = "edge_issue"
            score = (
                70.0 * normalized_mid_offset
                + 20.0 * normalized_gap
                + 16.0 * reversal
                + 8.0 * max(0.0, 1.0 - size_stability)
                + 4.0 * min(shape_change, 2.5)
            ) * (0.45 + 0.55 * points_weight)
            row = {
                "triad_id": f"SL{prev_sl}_SL{mid_sl}_SL{next_sl}",
                "prev_slice": prev_sl,
                "mid_slice": mid_sl,
                "next_slice": next_sl,
                "component_rank": rank,
                "prev_component_id": prev["component_id"],
                "mid_component_id": mid["component_id"],
                "next_component_id": nxt["component_id"],
                "prev_points": prev["points"],
                "mid_points": mid["points"],
                "next_points": nxt["points"],
                "prev_mid_distance": prev_mid,
                "mid_next_distance": mid_next,
                "prev_next_distance": prev_next,
                "mid_consensus_distance": mid_consensus_dist,
                "mid_expected_x": mid_expected_x,
                "mid_expected_y": mid_expected_y,
                "mid_centroid_x": mid["centroid_x"],
                "mid_centroid_y": mid["centroid_y"],
                "scale": scale,
                "normalized_mid_offset": normalized_mid_offset,
                "normalized_gap": normalized_gap,
                "edge_imbalance": imbalance / scale,
                "direction_reversal": reversal,
                "size_stability": size_stability,
                "shape_change": shape_change,
                "issue_type": issue_type,
                "triad_outlier_score": score,
            }
            status, reason = classify(row, args)
            row["status"] = status
            row["reason"] = reason
            rows.append(row)
    rows.sort(key=lambda r: (-float(r["triad_outlier_score"]), r["triad_id"], int(r["component_rank"])))
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--component-summary", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--min-component-points", type=int, default=50)
    parser.add_argument("--min-flag-points", type=int, default=100)
    parser.add_argument("--max-rank", type=int, default=8)
    parser.add_argument("--mid-offset-threshold", type=float, default=0.13)
    parser.add_argument("--neighbor-stability-threshold", type=float, default=0.55)
    parser.add_argument("--edge-issue-threshold", type=float, default=0.42)
    parser.add_argument("--edge-imbalance-threshold", type=float, default=0.15)
    parser.add_argument("--flag-threshold", type=float, default=28.0)
    parser.add_argument("--watch-threshold", type=float, default=16.0)
    args = parser.parse_args()

    by_slice = read_component_summary(args.component_summary, args.min_component_points)
    rows = scan(by_slice, args)
    write_csv(args.output_csv, rows)
    status_counts: dict[str, int] = defaultdict(int)
    issue_counts: dict[str, int] = defaultdict(int)
    flagged_triads: dict[str, int] = defaultdict(int)
    for row in rows:
        status_counts[str(row["status"])] += 1
        issue_counts[str(row["issue_type"])] += 1
        if row["status"] in {"flagged", "watch"}:
            flagged_triads[str(row["triad_id"])] += 1
    summary = {
        "component_summary": str(args.component_summary),
        "triad_windows": len({row["triad_id"] for row in rows}),
        "component_rows": len(rows),
        "status_counts": dict(status_counts),
        "issue_counts": dict(issue_counts),
        "flagged_or_watch_triads": len(flagged_triads),
        "top_rows": rows[:30],
        "parameters": {
            "min_component_points": args.min_component_points,
            "min_flag_points": args.min_flag_points,
            "max_rank": args.max_rank,
            "mid_offset_threshold": args.mid_offset_threshold,
            "neighbor_stability_threshold": args.neighbor_stability_threshold,
            "edge_issue_threshold": args.edge_issue_threshold,
            "edge_imbalance_threshold": args.edge_imbalance_threshold,
            "flag_threshold": args.flag_threshold,
            "watch_threshold": args.watch_threshold,
        },
    }
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
