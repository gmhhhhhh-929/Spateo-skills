#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_int(text: str, seed: int) -> int:
    data = f"{seed}|{text}".encode("utf-8")
    return int(hashlib.sha1(data).hexdigest()[:16], 16)


def norm_int_text(value: str) -> str:
    return str(int(float(value)))


def component_group_key(label: dict[str, str], component_key: str) -> tuple[str, str]:
    return (label["sl_number"], label[component_key])


def first_pass(
    coordinates: Path,
    components: Path,
    args: argparse.Namespace,
) -> tuple[dict[str, dict[str, str]], dict[tuple[str, str], dict[str, Any]], int]:
    comp_by_cell: dict[str, dict[str, str]] = {}
    stats: dict[tuple[str, str], dict[str, Any]] = {}
    with components.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            args.component_cell_id_col,
            args.component_slice_col,
            args.component_rank_col,
            args.component_id_col,
            args.component_points_col,
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Component labels missing columns: {sorted(missing)}")
        for row in reader:
            cell_id = row[args.component_cell_id_col]
            sl = norm_int_text(row[args.component_slice_col])
            rank = norm_int_text(row[args.component_rank_col])
            comp_by_cell[cell_id] = {
                "sl_number": sl,
                "component_rank": rank,
                "component_id": str(row[args.component_id_col]),
                "component_points": norm_int_text(row[args.component_points_col]),
            }

    coordinate_rows = 0
    with coordinates.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {args.cell_id_col, args.slice_col, args.x_col, args.y_col}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Coordinates missing columns: {sorted(missing)}")
        for row in reader:
            coordinate_rows += 1
            label = comp_by_cell.get(row[args.cell_id_col])
            if label is None:
                continue
            key = component_group_key(label, args.component_key)
            item = stats.setdefault(
                key,
                {
                    "slice": label["sl_number"],
                    "component_rank": label["component_rank"],
                    "component_id": label["component_id"],
                    "component_key": args.component_key,
                    "component_key_value": key[1],
                    "component_points": int(label["component_points"]),
                    "count": 0,
                    "min_x": float("inf"),
                    "max_x": float("-inf"),
                    "min_y": float("inf"),
                    "max_y": float("-inf"),
                },
            )
            x = float(row[args.x_col])
            y = float(row[args.y_col])
            item["count"] += 1
            item["min_x"] = min(item["min_x"], x)
            item["max_x"] = max(item["max_x"], x)
            item["min_y"] = min(item["min_y"], y)
            item["max_y"] = max(item["max_y"], y)
    return comp_by_cell, stats, coordinate_rows


def fraction_target(count: int, args: argparse.Namespace) -> int:
    if count <= args.keep_all_below:
        return count
    target = max(args.min_per_component, int(math.ceil(count * args.fraction)))
    if args.max_per_component > 0:
        target = min(target, args.max_per_component)
    return min(count, target)


def compute_targets(stats: dict[tuple[str, str], dict[str, Any]], args: argparse.Namespace) -> tuple[dict[tuple[str, str], int], dict[str, Any]]:
    source_points = sum(int(item["count"]) for item in stats.values())
    if args.target_total_points <= 0 or source_points <= args.target_total_points:
        targets = {key: (int(item["count"]) if source_points <= args.target_total_points and args.target_total_points > 0 else fraction_target(int(item["count"]), args)) for key, item in stats.items()}
        return targets, {
            "target_mode": "all_points" if source_points <= args.target_total_points and args.target_total_points > 0 else "fraction",
            "target_total_points": int(args.target_total_points),
            "effective_target_total": int(sum(targets.values())),
            "target_exceeded_due_to_minimums": False,
        }

    targets: dict[tuple[str, str], int] = {}
    caps: dict[tuple[str, str], int] = {}
    for key, item in stats.items():
        count = int(item["count"])
        if count <= args.keep_all_below:
            base = count
        else:
            base = min(count, args.min_per_component)
        cap = count if args.max_per_component <= 0 else min(count, args.max_per_component)
        cap = max(cap, base)
        targets[key] = base
        caps[key] = cap

    requested = int(args.target_total_points)
    base_total = sum(targets.values())
    if base_total >= requested:
        return targets, {
            "target_mode": "target_total_points",
            "target_total_points": requested,
            "effective_target_total": int(base_total),
            "target_exceeded_due_to_minimums": True,
        }

    remaining = requested - base_total
    while remaining > 0:
        open_keys = [key for key in stats if targets[key] < caps[key]]
        if not open_keys:
            break
        weights = {key: math.sqrt(max(int(stats[key]["count"]) - targets[key], 1)) for key in open_keys}
        weight_sum = sum(weights.values())
        raw = {key: remaining * weights[key] / weight_sum for key in open_keys}
        increments = {key: min(caps[key] - targets[key], int(math.floor(raw[key]))) for key in open_keys}
        assigned = sum(increments.values())
        leftover = remaining - assigned
        if leftover > 0:
            for key in sorted(open_keys, key=lambda k: (raw[k] - math.floor(raw[k]), int(stats[k]["count"])), reverse=True):
                if leftover <= 0:
                    break
                room = caps[key] - targets[key] - increments[key]
                if room <= 0:
                    continue
                increments[key] += 1
                leftover -= 1
        progress = sum(increments.values())
        if progress <= 0:
            break
        for key, inc in increments.items():
            targets[key] += inc
        remaining -= progress

    return targets, {
        "target_mode": "target_total_points",
        "target_total_points": requested,
        "effective_target_total": int(sum(targets.values())),
        "target_exceeded_due_to_minimums": False,
    }


def select_grid_balanced(
    bins: dict[tuple[int, int], list[tuple[int, dict[str, str], str]]],
    target: int,
) -> list[tuple[int, dict[str, str], str]]:
    for rows in bins.values():
        rows.sort(key=lambda v: v[0])
    bin_keys = sorted(bins)
    positions = {key: 0 for key in bin_keys}
    selected: list[tuple[int, dict[str, str], str]] = []
    while len(selected) < target:
        progressed = False
        for key in bin_keys:
            pos = positions[key]
            rows = bins[key]
            if pos >= len(rows):
                continue
            selected.append(rows[pos])
            positions[key] = pos + 1
            progressed = True
            if len(selected) >= target:
                break
        if not progressed:
            break
    return selected


def second_pass(
    coordinates: Path,
    comp_by_cell: dict[str, dict[str, str]],
    stats: dict[tuple[str, str], dict[str, Any]],
    out_csv: Path,
    summary_csv: Path,
    args: argparse.Namespace,
) -> dict[str, Any]:
    targets, target_meta = compute_targets(stats, args)
    selected_by_key: dict[tuple[str, str], list[tuple[int, dict[str, str], str]]] = defaultdict(list)
    bins_by_key: dict[tuple[str, str], dict[tuple[int, int], list[tuple[int, dict[str, str], str]]]] = defaultdict(lambda: defaultdict(list))
    fields: list[str] | None = None

    with coordinates.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit("Coordinates missing header")
        fields = list(reader.fieldnames)
        for col in [
            "sample_component_id",
            "sample_component_rank",
            "sample_component_key",
            "sample_component_key_value",
            "sample_strategy",
            "sample_key",
        ]:
            if col not in fields:
                fields.append(col)
        for row in reader:
            label = comp_by_cell.get(row[args.cell_id_col])
            if label is None:
                continue
            key = component_group_key(label, args.component_key)
            item = stats[key]
            target = targets[key]
            row["sample_component_id"] = label["component_id"]
            row["sample_component_rank"] = label["component_rank"]
            row["sample_component_key"] = args.component_key
            row["sample_component_key_value"] = key[1]
            row["sample_key"] = f"SL{label['sl_number']}_{args.component_key}_{key[1]}"
            score = stable_int(row[args.cell_id_col], args.seed)
            if int(item["count"]) <= target:
                row["sample_strategy"] = "keep_all_component"
                selected_by_key[key].append((score, row, row["sample_strategy"]))
                continue
            x = float(row[args.x_col])
            y = float(row[args.y_col])
            span_x = max(float(item["max_x"]) - float(item["min_x"]), 1e-6)
            span_y = max(float(item["max_y"]) - float(item["min_y"]), 1e-6)
            grid = max(1, int(math.sqrt(max(target, 1))))
            bx = min(grid - 1, max(0, int((x - float(item["min_x"])) / span_x * grid)))
            by = min(grid - 1, max(0, int((y - float(item["min_y"])) / span_y * grid)))
            row["sample_strategy"] = "spatial_grid_balanced"
            bins_by_key[key][(bx, by)].append((score, row, row["sample_strategy"]))

    for key, bins in bins_by_key.items():
        selected_by_key[key].extend(select_grid_balanced(bins, targets[key]))

    final_rows: list[dict[str, str]] = []
    for key, rows in selected_by_key.items():
        target = targets[key]
        rows.sort(key=lambda v: v[0])
        final_rows.extend([row for _score, row, _strategy in rows[:target]])
    final_rows.sort(
        key=lambda r: (
            int(float(r[args.slice_col])),
            int(float(r.get("sample_component_rank", 999999))),
            stable_int(r[args.cell_id_col], args.seed),
        )
    )

    if fields is None:
        raise SystemExit("Coordinates missing header")
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(final_rows)

    selected_counts: dict[tuple[str, str], int] = defaultdict(int)
    for row in final_rows:
        label = comp_by_cell[row[args.cell_id_col]]
        selected_counts[component_group_key(label, args.component_key)] += 1

    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        fields_summary = [
            "slice",
            "component_rank",
            "component_id",
            "component_key",
            "component_key_value",
            "full_count",
            "target_count",
            "sampled_count",
            "sample_fraction",
            "strategy",
            "min_x",
            "max_x",
            "min_y",
            "max_y",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields_summary)
        writer.writeheader()
        for key, item in sorted(stats.items(), key=lambda kv: (int(kv[1]["slice"]), int(kv[1]["component_rank"]), str(kv[1]["component_id"]))):
            sampled = selected_counts.get(key, 0)
            writer.writerow(
                {
                    "slice": item["slice"],
                    "component_rank": item["component_rank"],
                    "component_id": item["component_id"],
                    "component_key": item["component_key"],
                    "component_key_value": item["component_key_value"],
                    "full_count": item["count"],
                    "target_count": targets[key],
                    "sampled_count": sampled,
                    "sample_fraction": sampled / max(int(item["count"]), 1),
                    "strategy": "keep_all_component" if int(item["count"]) <= targets[key] else "spatial_grid_balanced",
                    "min_x": item["min_x"],
                    "max_x": item["max_x"],
                    "min_y": item["min_y"],
                    "max_y": item["max_y"],
                }
            )

    by_slice: dict[str, int] = defaultdict(int)
    for row in final_rows:
        by_slice[norm_int_text(row[args.slice_col])] += 1

    source_points = sum(int(item["count"]) for item in stats.values())
    return {
        "source_points": source_points,
        "sampled_points": len(final_rows),
        "components": len(stats),
        "slices": len({key[0] for key in stats}),
        "by_slice": dict(sorted(by_slice.items(), key=lambda kv: int(kv[0]))),
        "sampling_unit": f"{args.slice_col} x {args.component_key}",
        "component_key": args.component_key,
        "parameters": {
            "target_total_points": args.target_total_points,
            "fraction": args.fraction,
            "min_per_component": args.min_per_component,
            "max_per_component": args.max_per_component,
            "keep_all_below": args.keep_all_below,
            "seed": args.seed,
        },
        **target_meta,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a component-balanced spatial review sample.")
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--components", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--output-summary-csv", required=True, type=Path)
    parser.add_argument("--output-manifest", required=True, type=Path)
    parser.add_argument("--x-col", default="x")
    parser.add_argument("--y-col", default="y")
    parser.add_argument("--z-col", default="z")
    parser.add_argument("--slice-col", default="sl_number")
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--component-cell-id-col", default="cell_id")
    parser.add_argument("--component-slice-col", default="sl_number")
    parser.add_argument("--component-rank-col", default="component_rank")
    parser.add_argument("--component-id-col", default="component_id")
    parser.add_argument("--component-points-col", default="component_points")
    parser.add_argument("--component-key", choices=["component_id", "component_rank"], default="component_id")
    parser.add_argument("--target-total-points", type=int, default=300000, help="0 disables target mode and uses --fraction.")
    parser.add_argument("--fraction", type=float, default=0.08)
    parser.add_argument("--min-per-component", type=int, default=80)
    parser.add_argument("--max-per-component", type=int, default=0, help="0 means no hard per-component cap.")
    parser.add_argument("--keep-all-below", type=int, default=300)
    parser.add_argument("--seed", type=int, default=13)
    args = parser.parse_args()

    comp_by_cell, stats, coordinate_rows = first_pass(args.coordinates, args.components, args)
    manifest = second_pass(args.coordinates, comp_by_cell, stats, args.output_csv, args.output_summary_csv, args)
    manifest.update(
        {
            "coordinates": str(args.coordinates),
            "components": str(args.components),
            "coordinate_rows": coordinate_rows,
            "output_csv": str(args.output_csv),
            "output_summary_csv": str(args.output_summary_csv),
            "coordinate_columns": {"x": args.x_col, "y": args.y_col, "z": args.z_col},
            "coordinate_schema": "clean_coordinates_v1" if [args.x_col, args.y_col, args.z_col] == ["x", "y", "z"] else "custom",
            "id_columns": {"coordinates_cell_id": args.cell_id_col, "components_cell_id": args.component_cell_id_col},
            "input_hashes": {
                "coordinates": sha256(args.coordinates),
                "components": sha256(args.components),
            },
            "output_hashes": {
                "sample_csv": sha256(args.output_csv),
                "summary_csv": sha256(args.output_summary_csv),
            },
        }
    )
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.output_manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
