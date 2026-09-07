#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from pairwise_core import edge_transform_row, find_edge_dirs, json_dump, now_iso, read_csv_rows, write_csv_rows, write_edge_transform


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract and summarize pairwise Spateo refine reruns.")
    parser.add_argument("--refine-run-dir", required=True)
    parser.add_argument("--stage", default="stage1_SN-S_rigid")
    parser.add_argument("--jobs-csv")
    parser.add_argument("--output-dir")
    args = parser.parse_args()

    run_dir = Path(args.refine_run_dir).expanduser().resolve()
    outdir = Path(args.output_dir).expanduser().resolve() if args.output_dir else run_dir
    reruns_root = run_dir / "reruns"
    edge_dirs = find_edge_dirs(edges_root=reruns_root)
    rows = []
    for edge_dir in edge_dirs:
        item = write_edge_transform(edge_dir, args.stage)
        row = edge_transform_row(item)
        row["candidate_id"] = edge_dir.name
        rows.append(row)
    jobs = read_csv_rows(Path(args.jobs_csv).expanduser().resolve()) if args.jobs_csv else []
    job_by_candidate = {row.get("candidate_id"): row for row in jobs}
    for row in rows:
        job = job_by_candidate.get(row.get("candidate_id"), {})
        for key in ["edge_id", "target_celltype", "rank_score", "initial_pretransform_recipe"]:
            if job.get(key):
                row[key] = job[key]
    outdir.mkdir(parents=True, exist_ok=True)
    write_csv_rows(outdir / "pairwise_refine_edge_index.csv", rows)
    summary = {
        "version": "spateo_pairwise_refine_summary_v1",
        "created_at": now_iso(),
        "refine_run_dir": str(run_dir),
        "n_reruns": len(rows),
        "ok_candidates": [row.get("candidate_id") for row in rows if row.get("status") == "ok"],
        "failed_or_missing_candidates": [row.get("candidate_id") for row in rows if row.get("status") != "ok"],
    }
    json_dump(outdir / "refine_summary.json", summary)
    print(outdir / "refine_summary.json")


if __name__ == "__main__":
    main()
