#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from pairwise_core import edge_transform_row, find_edge_dirs, write_csv_rows, write_edge_transform


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract per-edge Spateo transform audit into edge_transform.json files.")
    parser.add_argument("--run-dir")
    parser.add_argument("--edges-root")
    parser.add_argument("--edge-dir")
    parser.add_argument("--stage", default="stage1_SN-S_rigid")
    parser.add_argument("--output-edge-index")
    args = parser.parse_args()

    edge_dirs = find_edge_dirs(
        Path(args.run_dir) if args.run_dir else None,
        Path(args.edges_root) if args.edges_root else None,
        Path(args.edge_dir) if args.edge_dir else None,
    )
    if not edge_dirs:
        raise SystemExit("No edge directories found.")
    transforms = [write_edge_transform(edge_dir, args.stage) for edge_dir in edge_dirs]
    if args.output_edge_index:
        write_csv_rows(Path(args.output_edge_index).expanduser().resolve(), [edge_transform_row(item) for item in transforms])
    print(f"extracted {len(transforms)} edge transform(s)")


if __name__ == "__main__":
    main()
