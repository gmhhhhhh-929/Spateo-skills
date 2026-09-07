#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_scan_row(path: Path, triad_id: str, component_rank: int) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("triad_id") == triad_id and int(float(row.get("component_rank", "0"))) == component_rank:
                return row
    raise SystemExit(f"No scan row for {triad_id} rank{component_rank}")


def read_component_labels(path: Path, slices: set[str]) -> dict[str, dict[str, str]]:
    labels: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sl = str(int(float(row["sl_number"])))
            if sl in slices:
                labels[row["cell_id"]] = {
                    "sl_number": sl,
                    "component_rank": str(int(float(row["component_rank"]))),
                    "component_id": str(row["component_id"]),
                }
    return labels


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--components", required=True, type=Path)
    parser.add_argument("--scan", required=True, type=Path)
    parser.add_argument("--triad-id", required=True)
    parser.add_argument("--component-rank", required=True, type=int)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--output-recipe", required=True, type=Path)
    parser.add_argument("--x-col", default="manual_x")
    parser.add_argument("--y-col", default="manual_y")
    parser.add_argument("--z-col", default="manual_z")
    parser.add_argument("--slice-col", default="sl_number")
    parser.add_argument("--celltype-col", default="celltype")
    parser.add_argument("--out-x-col", default="triad_preview_x")
    parser.add_argument("--out-y-col", default="triad_preview_y")
    parser.add_argument("--out-z-col", default="triad_preview_z")
    args = parser.parse_args()

    scan = read_scan_row(args.scan, args.triad_id, args.component_rank)
    mid_slice = str(int(float(scan["mid_slice"])))
    prev_slice = str(int(float(scan["prev_slice"])))
    next_slice = str(int(float(scan["next_slice"])))
    dx = float(scan["mid_expected_x"]) - float(scan["mid_centroid_x"])
    dy = float(scan["mid_expected_y"]) - float(scan["mid_centroid_y"])
    labels = read_component_labels(args.components, {prev_slice, mid_slice, next_slice})

    fields_written = False
    moved = 0
    total = 0
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.coordinates.open(newline="", encoding="utf-8") as src, args.output_csv.open("w", newline="", encoding="utf-8") as dst:
        reader = csv.DictReader(src)
        if reader.fieldnames is None:
            raise SystemExit("Coordinates have no header")
        fields = list(reader.fieldnames)
        for col in [args.out_x_col, args.out_y_col, args.out_z_col, "triad_preview_action"]:
            if col not in fields:
                fields.append(col)
        writer = csv.DictWriter(dst, fieldnames=fields)
        writer.writeheader()
        for row in reader:
            total += 1
            sl = str(int(float(row[args.slice_col])))
            label = labels.get(row["cell_id"], {})
            is_target = sl == mid_slice and label.get("component_rank") == str(args.component_rank)
            row[args.out_x_col] = str(float(row[args.x_col]) + (dx if is_target else 0.0))
            row[args.out_y_col] = str(float(row[args.y_col]) + (dy if is_target else 0.0))
            row[args.out_z_col] = row[args.z_col]
            row["triad_preview_action"] = "move_mid_component_to_triad_consensus" if is_target else "unchanged"
            moved += int(is_target)
            writer.writerow(row)
    recipe = {
        "recipe_type": "triad_mid_slice_rigid_preview",
        "triad_id": args.triad_id,
        "component_rank": args.component_rank,
        "prev_slice": prev_slice,
        "mid_slice": mid_slice,
        "next_slice": next_slice,
        "transform": {"operation": "translate", "dx": dx, "dy": dy, "dz": 0.0},
        "scope": {"slice": mid_slice, "component_rank": args.component_rank},
        "points_moved": moved,
        "total_points": total,
        "guardrails": ["preview only", "moves only the flagged mid component", "does not modify final h5ad"],
    }
    args.output_recipe.parent.mkdir(parents=True, exist_ok=True)
    args.output_recipe.write_text(json.dumps(recipe, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(recipe, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
