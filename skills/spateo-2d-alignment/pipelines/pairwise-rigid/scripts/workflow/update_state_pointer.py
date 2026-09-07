#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
from typing import Any


VALID_STATES = {"slice_level_baseline", "current_review", "current_accepted", "final"}


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rel_or_abs(path: str, root: Path) -> str:
    if not path:
        return ""
    p = Path(path).expanduser()
    try:
        return os.path.relpath(p.resolve(), root.resolve())
    except Exception:
        return str(path)


def path_info(path: str, root: Path) -> dict[str, Any]:
    if not path:
        return {}
    p = (root / path).resolve() if not Path(path).is_absolute() else Path(path).expanduser().resolve()
    info: dict[str, Any] = {"path": rel_or_abs(str(p), root), "exists": p.exists()}
    if p.exists() and p.is_file():
        stat = p.stat()
        info.update({"size_bytes": stat.st_size, "sha256": sha256_file(p)})
    return info


def main() -> None:
    parser = argparse.ArgumentParser(description="Update a run states/*.yaml pointer and append lineage/state_history.jsonl.")
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--state", required=True, choices=sorted(VALID_STATES))
    parser.add_argument("--step", required=True)
    parser.add_argument("--status", required=True)
    parser.add_argument("--summary", default="")
    parser.add_argument("--clean-coordinates", default="")
    parser.add_argument("--audit-points", default="")
    parser.add_argument("--viewer", default="")
    parser.add_argument("--recipe-chain", default="")
    parser.add_argument("--validation", default="")
    parser.add_argument("--manifest", default="")
    parser.add_argument("--row-count", type=int, default=-1)
    parser.add_argument("--code-release", default="")
    parser.add_argument("--validator-status", default="")
    args = parser.parse_args()

    root = args.run_root.expanduser().resolve()
    (root / "states").mkdir(parents=True, exist_ok=True)
    (root / "lineage").mkdir(parents=True, exist_ok=True)
    record = {
        "version": "spatial_alignment_state_pointer_v1",
        "updated_at": now_iso(),
        "state": args.state,
        "step": args.step,
        "status": args.status,
        "summary": args.summary,
        "row_count": args.row_count if args.row_count >= 0 else None,
        "code_release": args.code_release,
        "validator_status": args.validator_status,
        "artifacts": {
            "clean_coordinates": path_info(args.clean_coordinates, root),
            "audit_points": path_info(args.audit_points, root),
            "viewer": path_info(args.viewer, root),
            "recipe_chain": path_info(args.recipe_chain, root),
            "validation": path_info(args.validation, root),
            "manifest": path_info(args.manifest, root),
        },
    }
    state_path = root / "states" / f"{args.state}.yaml"
    state_path.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    history = root / "lineage" / "state_history.jsonl"
    with history.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
    print(json.dumps(record, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
