#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def load_sample_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            ids.add(row["cell_id"])
    return ids


def load_coords(path: Path, ids: set[str], x_col: str, y_col: str, z_col: str) -> dict[str, tuple[float, float, float]]:
    out: dict[str, tuple[float, float, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"cell_id", x_col, y_col, z_col}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{path} missing columns: {sorted(missing)}")
        for row in reader:
            cid = row["cell_id"]
            if cid in ids:
                out[cid] = (float(row[x_col]), float(row[y_col]), float(row[z_col]))
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-input", required=True, type=Path)
    parser.add_argument("--sample-replay", required=True, type=Path)
    parser.add_argument("--full-replay", required=True, type=Path)
    parser.add_argument("--x-col", default="test_x")
    parser.add_argument("--y-col", default="test_y")
    parser.add_argument("--z-col", default="test_z")
    parser.add_argument("--atol", type=float, default=1e-6)
    parser.add_argument("--output-json", required=True, type=Path)
    args = parser.parse_args()

    ids = load_sample_ids(args.sample_input)
    sample = load_coords(args.sample_replay, ids, args.x_col, args.y_col, args.z_col)
    full = load_coords(args.full_replay, ids, args.x_col, args.y_col, args.z_col)
    missing_sample = sorted(ids - set(sample))
    missing_full = sorted(ids - set(full))
    max_abs_delta = 0.0
    max_l2_delta = 0.0
    mismatches = 0
    for cid in ids & set(sample) & set(full):
        ds = [abs(sample[cid][i] - full[cid][i]) for i in range(3)]
        l2 = math.sqrt(sum(d * d for d in ds))
        max_abs_delta = max(max_abs_delta, max(ds))
        max_l2_delta = max(max_l2_delta, l2)
        mismatches += int(max(ds) > args.atol)
    report = {
        "status": "pass" if not missing_sample and not missing_full and mismatches == 0 else "fail",
        "sample_ids": len(ids),
        "sample_replay_rows_matched": len(sample),
        "full_replay_rows_matched": len(full),
        "missing_sample_count": len(missing_sample),
        "missing_full_count": len(missing_full),
        "mismatch_count": mismatches,
        "max_abs_delta": max_abs_delta,
        "max_l2_delta": max_l2_delta,
        "atol": args.atol,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if report["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
