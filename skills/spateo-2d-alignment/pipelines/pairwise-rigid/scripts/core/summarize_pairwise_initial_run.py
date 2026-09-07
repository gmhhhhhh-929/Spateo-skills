#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from pairwise_core import edge_transform_row, find_edge_dirs, json_dump, load_json_or_none, maybe_float, now_iso, write_csv_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize pairwise initial Spateo run edge status and native metrics.")
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--edges-root")
    parser.add_argument("--output-dir")
    parser.add_argument("--high-sigma2-threshold", type=float, default=0.1)
    args = parser.parse_args()

    run_dir = Path(args.run_dir).expanduser().resolve()
    outdir = Path(args.output_dir).expanduser().resolve() if args.output_dir else run_dir
    edge_dirs = find_edge_dirs(run_dir, Path(args.edges_root) if args.edges_root else None)
    rows = []
    for edge_dir in edge_dirs:
        obj = load_json_or_none(edge_dir / "edge_transform.json")
        if not isinstance(obj, dict):
            provenance = load_json_or_none(edge_dir / "provenance.json")
            obj = {
                "edge_id": edge_dir.name,
                "status": "missing_edge_transform",
                "warnings": ["missing_edge_transform"],
                "fixed_slice_id": (provenance or {}).get("fixed_slice_id", "") if isinstance(provenance, dict) else "",
                "moving_slice_id": (provenance or {}).get("moving_slice_id", "") if isinstance(provenance, dict) else "",
                "fixed_sl_number": (provenance or {}).get("fixed_sl_number", "") if isinstance(provenance, dict) else "",
                "moving_sl_number": (provenance or {}).get("moving_sl_number", "") if isinstance(provenance, dict) else "",
                "pair_start_order": (provenance or {}).get("pair_start_order", "") if isinstance(provenance, dict) else "",
                "paths": {"edge_dir": str(edge_dir)},
            }
        row = edge_transform_row(obj)
        sigma2 = maybe_float(row.get("sigma2"))
        row["high_sigma2"] = bool(sigma2 is not None and sigma2 > args.high_sigma2_threshold)
        rows.append(row)
    outdir.mkdir(parents=True, exist_ok=True)
    write_csv_rows(outdir / "pairwise_edge_index.csv", rows)
    status_counts: dict[str, int] = {}
    for row in rows:
        key = str(row.get("status", ""))
        status_counts[key] = status_counts.get(key, 0) + 1
    sigma_values = [maybe_float(row.get("sigma2")) for row in rows]
    sigma_values = [x for x in sigma_values if x is not None]
    summary = {
        "version": "spateo_pairwise_edge_summary_v2_clean",
        "created_at": now_iso(),
        "run_dir": str(run_dir),
        "n_edges": len(rows),
        "status_counts": status_counts,
        "high_sigma2_threshold": args.high_sigma2_threshold,
        "n_high_sigma2": sum(1 for row in rows if row.get("high_sigma2")),
        "high_sigma2_edges": [row["edge_id"] for row in rows if row.get("high_sigma2")],
        "missing_edges": [row["edge_id"] for row in rows if str(row.get("status", "")).startswith("missing")],
        "mean_sigma2": float(np.mean(sigma_values)) if sigma_values else None,
        "median_sigma2": float(np.median(sigma_values)) if sigma_values else None,
        "max_sigma2": float(np.max(sigma_values)) if sigma_values else None,
    }
    json_dump(outdir / "pairwise_edge_summary.json", summary)
    print(outdir / "pairwise_edge_index.csv")


if __name__ == "__main__":
    main()
