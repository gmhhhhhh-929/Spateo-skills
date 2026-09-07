#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
from pathlib import Path
from typing import Any


FIELDNAMES = [
    "workflow_cell_id",
    "cell_id",
    "h5ad_obs_name",
    "slice_id",
    "stage",
    "chip_id",
    "sl_number",
    "slice_short",
    "obs_row_index",
    "source_h5ad",
    "obs_cell_id",
    "obs_CellID",
]


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def natural_key(text: str) -> list[Any]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(block_size), b""):
            h.update(chunk)
    return h.hexdigest()


def file_meta(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size_bytes": stat.st_size,
        "mtime": stat.st_mtime,
        "sha256": sha256_file(resolved),
    }


def parse_slice(path: Path) -> dict[str, Any] | None:
    match = re.search(r"(?P<stage>CS\d+)_SL(?P<sl>\d+)_(?P<chip>Y\w+)\.Spatial", path.name)
    if not match:
        return None
    return {
        "slice_id": path.name.split(".Spatial")[0],
        "stage": match.group("stage"),
        "chip_id": match.group("chip"),
        "sl_number": int(match.group("sl")),
    }


def build_map(args: argparse.Namespace) -> dict[str, Any]:
    try:
        import anndata as ad
    except Exception as exc:  # pragma: no cover - depends on runtime env
        raise SystemExit("build_h5ad_cell_id_map.py requires anndata in the active environment") from exc

    root = args.h5ad_root.expanduser().resolve()
    output = args.output_csv.expanduser().resolve()
    manifest_path = args.output_manifest.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    h5ad_paths = []
    for path in sorted(root.rglob("*.h5ad"), key=lambda p: natural_key(str(p))):
        meta = parse_slice(path)
        if meta is None:
            continue
        if args.sample_id and meta["stage"] != args.sample_id:
            continue
        h5ad_paths.append(path)
    if not h5ad_paths:
        raise SystemExit(f"No matching h5ad files found under {root}")

    row_count = 0
    seen_workflow: set[str] = set()
    seen_cell_id: set[str] = set()
    slice_counts: dict[str, int] = {}
    h5ad_meta: list[dict[str, Any]] = []

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        for path in h5ad_paths:
            meta = parse_slice(path)
            if meta is None:
                continue
            adata = ad.read_h5ad(path, backed="r")
            try:
                obs_names = [str(value) for value in adata.obs_names.tolist()]
                obs_CellID = (
                    [str(value) for value in adata.obs["CellID"].tolist()]
                    if "CellID" in adata.obs
                    else [""] * len(obs_names)
                )
                obs_cell_id = (
                    [str(value) for value in adata.obs["cell_id"].tolist()]
                    if "cell_id" in adata.obs
                    else [value.replace("CELL.", "") if value.startswith("CELL.") else value for value in obs_CellID]
                )
                obs_slice = (
                    [str(value) for value in adata.obs["slice"].tolist()]
                    if "slice" in adata.obs
                    else [f"SL{meta['sl_number']}"] * len(obs_names)
                )
                for i, obs_name in enumerate(obs_names):
                    workflow_cell_id = f"{meta['slice_id']}:{i}"
                    short_cell_id = f"{obs_slice[i]}_{obs_CellID[i]}"
                    if workflow_cell_id in seen_workflow:
                        raise SystemExit(f"Duplicate workflow_cell_id: {workflow_cell_id}")
                    if short_cell_id in seen_cell_id:
                        raise SystemExit(f"Duplicate clean short cell_id: {short_cell_id}")
                    seen_workflow.add(workflow_cell_id)
                    seen_cell_id.add(short_cell_id)
                    writer.writerow(
                        {
                            "workflow_cell_id": workflow_cell_id,
                            "cell_id": short_cell_id,
                            "h5ad_obs_name": obs_name,
                            "slice_id": meta["slice_id"],
                            "stage": meta["stage"],
                            "chip_id": meta["chip_id"],
                            "sl_number": meta["sl_number"],
                            "slice_short": obs_slice[i],
                            "obs_row_index": i,
                            "source_h5ad": str(path),
                            "obs_cell_id": obs_cell_id[i],
                            "obs_CellID": obs_CellID[i],
                        }
                    )
                row_count += len(obs_names)
                slice_counts[str(meta["sl_number"])] = len(obs_names)
                h5ad_meta.append({**meta, "path": str(path), "n_obs": int(adata.n_obs)})
            finally:
                adata.file.close()

    manifest = {
        "version": "spatial_h5ad_cell_id_map_manifest_v1",
        "created_at": now_iso(),
        "status": "pass",
        "h5ad_root": str(root),
        "sample_id": args.sample_id,
        "workflow_cell_id_rule": "slice_id + ':' + zero_based_h5ad_obs_row_index",
        "clean_cell_id_source": "short_slice_cellid",
        "clean_cell_id_rule": "h5ad obs['slice'] + '_' + h5ad obs['CellID']",
        "h5ad_obs_name_column": "h5ad_obs_name",
        "output": file_meta(output),
        "row_count": row_count,
        "unique_workflow_cell_ids": len(seen_workflow),
        "unique_clean_cell_ids": len(seen_cell_id),
        "slice_count": len(slice_counts),
        "slice_counts": slice_counts,
        "h5ad_files": h5ad_meta,
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a workflow-cell-id to original h5ad obs_name map.")
    parser.add_argument("--h5ad-root", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--output-manifest", required=True, type=Path)
    parser.add_argument("--sample-id", default="", help="Optional stage/sample id such as CS13.")
    args = parser.parse_args()
    print(json.dumps(build_map(args), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
