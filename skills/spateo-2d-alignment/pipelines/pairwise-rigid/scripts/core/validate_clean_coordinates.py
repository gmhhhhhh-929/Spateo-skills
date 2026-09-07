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


CLEAN_COLUMNS = ["cell_id", "slice_id", "stage", "chip_id", "sl_number", "celltype", "x", "y", "z"]
FORBIDDEN_PREFIXES = ("raw_", "initial_", "stage1_", "stage2_", "manual_", "full_candidate_", "step")
FORBIDDEN_TOKENS = ("edit_label", "displacement", "component_rank", "component_id", "recipe", "operation")


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fail(report: dict[str, Any], reason: str) -> None:
    report["status"] = "fail"
    report.setdefault("errors", []).append(reason)


def validate(args: argparse.Namespace) -> dict[str, Any]:
    path = args.input.expanduser().resolve()
    report: dict[str, Any] = {
        "version": "spatial_clean_coordinates_validation_v1",
        "created_at": now_iso(),
        "input": str(path),
        "status": "pass",
        "errors": [],
        "warnings": [],
    }
    if not path.exists():
        fail(report, f"missing input: {path}")
        return report

    seen: set[str] = set()
    row_count = 0
    celltype_values: set[str] = set()
    row_index_id_re = re.compile(args.row_index_cell_id_regex) if args.forbid_row_index_cell_id_pattern else None
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []
        report["columns"] = header
        if header != CLEAN_COLUMNS:
            fail(report, f"clean CSV header must be exactly {CLEAN_COLUMNS}; got {header}")
        for col in header:
            lower = col.lower()
            if lower.startswith(FORBIDDEN_PREFIXES) or any(token in lower for token in FORBIDDEN_TOKENS):
                fail(report, f"forbidden intermediate/audit column in clean CSV: {col}")
        for row in reader:
            row_count += 1
            cid = row.get("cell_id", "")
            if not cid:
                fail(report, f"empty cell_id at row {row_count}")
                if len(report["errors"]) > 20:
                    break
            elif row_index_id_re and row_index_id_re.match(cid):
                fail(report, f"cell_id looks like workflow row-index id, not original h5ad id: {cid}")
                if len(report["errors"]) > 20:
                    break
            elif not args.skip_unique_cell_id:
                if cid in seen:
                    fail(report, f"duplicate cell_id: {cid}")
                    if len(report["errors"]) > 20:
                        break
                seen.add(cid)
            for coord in ("x", "y", "z"):
                try:
                    float(row.get(coord, ""))
                except Exception:
                    fail(report, f"non-numeric {coord} at row {row_count}: {row.get(coord, '')!r}")
                    if len(report["errors"]) > 20:
                        break
            if row.get("celltype", ""):
                celltype_values.add(str(row["celltype"]))
    report["row_count"] = row_count
    report["sha256"] = sha256_file(path) if path.exists() else ""
    report["unique_cell_ids"] = len(seen) if not args.skip_unique_cell_id else None
    report["celltype_value_count"] = len(celltype_values)
    if args.expected_row_count >= 0 and row_count != args.expected_row_count:
        fail(report, f"row_count mismatch: expected {args.expected_row_count}, got {row_count}")
    suspicious = {"baseline", "step", "review", "candidate", "full_candidate"}
    if celltype_values and celltype_values <= suspicious:
        fail(report, f"celltype values look like workflow labels, not biology: {sorted(celltype_values)}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate clean spatial coordinate CSV schema and basic invariants.")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--expected-row-count", type=int, default=-1)
    parser.add_argument("--skip-unique-cell-id", action="store_true")
    parser.add_argument("--forbid-row-index-cell-id-pattern", action="store_true")
    parser.add_argument("--row-index-cell-id-regex", default=r"^.+:\d+$")
    parser.add_argument("--output-report", type=Path)
    args = parser.parse_args()
    report = validate(args)
    text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output_report:
        args.output_report.parent.mkdir(parents=True, exist_ok=True)
        args.output_report.write_text(text, encoding="utf-8")
    print(text, end="")
    if report["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
