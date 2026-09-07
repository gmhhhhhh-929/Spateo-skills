#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
from typing import Any


CLEAN_COLUMNS = ["cell_id", "slice_id", "stage", "chip_id", "sl_number", "celltype", "x", "y", "z"]


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


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


def load_cell_id_map(path: Path, key_col: str, value_col: str) -> tuple[dict[str, str], dict[str, Any]]:
    resolved = path.expanduser().resolve()
    mapping: dict[str, str] = {}
    duplicate_keys: list[str] = []
    duplicate_values: list[str] = []
    seen_values: set[str] = set()
    with resolved.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []
        missing = [col for col in [key_col, value_col] if col not in header]
        if missing:
            raise SystemExit(f"Cell ID map missing required columns {missing}: {resolved}")
        for i, row in enumerate(reader, start=2):
            key = str(row.get(key_col, "")).strip()
            value = str(row.get(value_col, "")).strip()
            if not key or not value:
                raise SystemExit(f"Cell ID map has empty key/value at line {i}: {resolved}")
            if key in mapping:
                duplicate_keys.append(key)
            if value in seen_values:
                duplicate_values.append(value)
            mapping[key] = value
            seen_values.add(value)
            if duplicate_keys or duplicate_values:
                break
    if duplicate_keys:
        raise SystemExit(f"Cell ID map duplicate key: {duplicate_keys[0]}")
    if duplicate_values:
        raise SystemExit(f"Cell ID map duplicate output cell_id: {duplicate_values[0]}")
    meta = file_meta(resolved)
    meta.update({"key_col": key_col, "value_col": value_col, "row_count": len(mapping)})
    return mapping, meta


def fnum(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def changed_by_coords(row: dict[str, str], args: argparse.Namespace) -> tuple[bool, float]:
    if not (args.baseline_x_col and args.baseline_y_col and args.baseline_z_col):
        return False, 0.0
    needed = [args.baseline_x_col, args.baseline_y_col, args.baseline_z_col, args.x_col, args.y_col, args.z_col]
    if any(col not in row for col in needed):
        return False, 0.0
    dx = fnum(row[args.x_col]) - fnum(row[args.baseline_x_col])
    dy = fnum(row[args.y_col]) - fnum(row[args.baseline_y_col])
    dz = fnum(row[args.z_col]) - fnum(row[args.baseline_z_col])
    disp = math.sqrt(dx * dx + dy * dy + dz * dz)
    return disp > args.changed_atol, disp


def changed_by_audit(row: dict[str, str]) -> bool:
    for key, value in row.items():
        lk = key.lower()
        text = str(value).strip()
        if "edit_label" in lk:
            if text and text.lower() not in {"baseline_unchanged", "unchanged", "none", "na", "null", "0"}:
                return True
    return False


def selected_audit_columns(fieldnames: list[str]) -> list[str]:
    keep = []
    for col in fieldnames:
        lk = col.lower()
        if any(token in lk for token in ["edit", "displacement", "component", "rank", "recipe", "operation"]):
            keep.append(col)
    return keep


def export_clean(args: argparse.Namespace) -> dict[str, Any]:
    source = args.input.expanduser().resolve()
    output = args.output.expanduser().resolve()
    manifest_path = args.manifest.expanduser().resolve()
    changed_path = args.changed_cells_output.expanduser().resolve() if args.changed_cells_output else None
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if changed_path:
        changed_path.parent.mkdir(parents=True, exist_ok=True)
    cell_id_map: dict[str, str] | None = None
    cell_id_map_meta: dict[str, Any] | None = None
    if args.cell_id_map:
        cell_id_map, cell_id_map_meta = load_cell_id_map(
            args.cell_id_map,
            args.cell_id_map_key_col,
            args.cell_id_map_value_col,
        )

    required = [args.cell_id_col, args.sl_number_col, args.celltype_col, args.x_col, args.y_col, args.z_col]
    optional_map = {
        "slice_id": args.slice_id_col,
        "stage": args.stage_col,
        "chip_id": args.chip_id_col,
    }

    row_count = 0
    changed_count = 0
    mapped_cell_id_count = 0
    clean_cell_ids: set[str] = set()
    audit_keep: list[str] = []
    source_fields: list[str] = []
    with source.open(newline="", encoding="utf-8") as inp:
        reader = csv.DictReader(inp)
        if reader.fieldnames is None:
            raise SystemExit(f"Input CSV has no header: {source}")
        source_fields = list(reader.fieldnames)
        missing = [col for col in required if col not in source_fields]
        if missing:
            raise SystemExit(f"Input CSV missing required columns: {missing}")
        audit_keep = selected_audit_columns(source_fields)
        changed_fields = [
            "cell_id",
            *(["workflow_cell_id"] if cell_id_map is not None else []),
            "sl_number",
            "changed_by_coordinate_delta",
            "displacement",
            "baseline_x",
            "baseline_y",
            "baseline_z",
            "x",
            "y",
            "z",
            *audit_keep,
        ]
        with output.open("w", newline="", encoding="utf-8") as out_handle:
            writer = csv.DictWriter(out_handle, fieldnames=CLEAN_COLUMNS)
            writer.writeheader()
            changed_handle = changed_path.open("w", newline="", encoding="utf-8") if changed_path else None
            try:
                changed_writer = csv.DictWriter(changed_handle, fieldnames=changed_fields) if changed_handle else None
                if changed_writer:
                    changed_writer.writeheader()
                for row in reader:
                    row_count += 1
                    source_cell_id = row[args.cell_id_col]
                    clean_cell_id = source_cell_id
                    if cell_id_map is not None:
                        clean_cell_id = cell_id_map.get(source_cell_id, "")
                        if not clean_cell_id:
                            raise SystemExit(
                                f"Cell ID map has no entry for {args.cell_id_col}={source_cell_id!r} "
                                f"at input row {row_count}"
                            )
                        mapped_cell_id_count += 1
                    if clean_cell_id in clean_cell_ids:
                        raise SystemExit(f"Duplicate output clean cell_id: {clean_cell_id!r}")
                    clean_cell_ids.add(clean_cell_id)
                    clean = {
                        "cell_id": clean_cell_id,
                        "slice_id": row.get(optional_map["slice_id"], ""),
                        "stage": row.get(optional_map["stage"], ""),
                        "chip_id": row.get(optional_map["chip_id"], ""),
                        "sl_number": row[args.sl_number_col],
                        "celltype": row[args.celltype_col],
                        "x": row[args.x_col],
                        "y": row[args.y_col],
                        "z": row[args.z_col],
                    }
                    writer.writerow(clean)
                    coord_changed, disp = changed_by_coords(row, args)
                    audit_changed = changed_by_audit(row)
                    if changed_writer and (coord_changed or audit_changed):
                        changed_count += 1
                        changed_row = {
                                "cell_id": clean_cell_id,
                                "sl_number": row[args.sl_number_col],
                                "changed_by_coordinate_delta": str(bool(coord_changed)).lower(),
                                "displacement": f"{disp:.12g}",
                                "baseline_x": row.get(args.baseline_x_col, "") if args.baseline_x_col else "",
                                "baseline_y": row.get(args.baseline_y_col, "") if args.baseline_y_col else "",
                                "baseline_z": row.get(args.baseline_z_col, "") if args.baseline_z_col else "",
                                "x": row[args.x_col],
                                "y": row[args.y_col],
                                "z": row[args.z_col],
                                **{col: row.get(col, "") for col in audit_keep},
                        }
                        if cell_id_map is not None:
                            changed_row["workflow_cell_id"] = source_cell_id
                        changed_writer.writerow(changed_row)
            finally:
                if changed_handle:
                    changed_handle.close()

    source_meta = file_meta(source)
    output_meta = file_meta(output)
    changed_meta = file_meta(changed_path) if changed_path and changed_path.exists() else None
    manifest = {
        "version": "spatial_clean_coordinates_manifest_v1",
        "created_at": now_iso(),
        "status": "pass",
        "source": source_meta,
        "output": output_meta,
        "changed_cells": changed_meta,
        "row_count": row_count,
        "changed_row_count": changed_count,
        "mapped_cell_id_count": mapped_cell_id_count,
        "clean_columns": CLEAN_COLUMNS,
        "source_columns": source_fields,
        "coordinate_mapping": {"x": args.x_col, "y": args.y_col, "z": args.z_col},
        "identity_mapping": {
            "source_cell_id": args.cell_id_col,
            "clean_cell_id_source": args.clean_cell_id_source
            or ("short_slice_cellid" if cell_id_map is not None else args.cell_id_col),
            "cell_id_map": cell_id_map_meta,
            "slice_id": args.slice_id_col,
            "stage": args.stage_col,
            "chip_id": args.chip_id_col,
            "sl_number": args.sl_number_col,
            "celltype": args.celltype_col,
        },
        "baseline_coordinate_mapping": {
            "x": args.baseline_x_col,
            "y": args.baseline_y_col,
            "z": args.baseline_z_col,
        },
        "code_release": args.code_release,
        "entrypoint_sha256": args.entrypoint_sha256,
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Export a clean user-facing coordinate CSV from a full audit coordinate table.")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--changed-cells-output", type=Path)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--cell-id-map", type=Path, help="Optional CSV mapping source workflow cell ids to clean output cell ids.")
    parser.add_argument("--cell-id-map-key-col", default="workflow_cell_id")
    parser.add_argument("--cell-id-map-value-col", default="cell_id")
    parser.add_argument("--clean-cell-id-source", default="")
    parser.add_argument("--slice-id-col", default="slice_id")
    parser.add_argument("--stage-col", default="stage")
    parser.add_argument("--chip-id-col", default="chip_id")
    parser.add_argument("--sl-number-col", default="sl_number")
    parser.add_argument("--celltype-col", default="celltype")
    parser.add_argument("--x-col", required=True)
    parser.add_argument("--y-col", required=True)
    parser.add_argument("--z-col", required=True)
    parser.add_argument("--baseline-x-col", default="")
    parser.add_argument("--baseline-y-col", default="")
    parser.add_argument("--baseline-z-col", default="")
    parser.add_argument("--changed-atol", type=float, default=1e-6)
    parser.add_argument("--code-release", default="")
    parser.add_argument("--entrypoint-sha256", default="")
    args = parser.parse_args()
    manifest = export_clean(args)
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
