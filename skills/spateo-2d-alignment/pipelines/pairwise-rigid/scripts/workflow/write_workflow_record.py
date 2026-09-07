#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
from typing import Any


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def parse_kv(values: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise SystemExit(f"Expected key=value: {value}")
        key, item = value.split("=", 1)
        out[key] = item
    return out


def append_jsonl(path: Path, obj: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(obj, ensure_ascii=False, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Write a compact spatial workflow record.")
    parser.add_argument("--record-dir", required=True, type=Path)
    parser.add_argument("--record-id", required=True)
    parser.add_argument("--operation-type", required=True)
    parser.add_argument("--status", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--skill-name", default="")
    parser.add_argument("--skill-version", default="")
    parser.add_argument("--code-release", default="")
    parser.add_argument("--entrypoint-sha256", default="")
    parser.add_argument("--validator-status", default="")
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--output", action="append", default=[], help="key=relative/path")
    parser.add_argument("--artifact", action="append", default=[])
    parser.add_argument("--metric", action="append", default=[], help="key=value")
    parser.add_argument("--guardrail", action="append", default=[])
    parser.add_argument("--index", type=Path)
    args = parser.parse_args()

    record = {
        "record_id": args.record_id,
        "created_at": now_iso(),
        "operation_type": args.operation_type,
        "status": args.status,
        "summary": args.summary,
        "skill_name": args.skill_name,
        "skill_version": args.skill_version,
        "code_release": args.code_release,
        "entrypoint_sha256": args.entrypoint_sha256,
        "validator_status": args.validator_status,
        "inputs": args.input,
        "outputs": parse_kv(args.output),
        "metrics": parse_kv(args.metric),
        "artifacts": args.artifact,
        "guardrails": args.guardrail,
    }
    args.record_dir.mkdir(parents=True, exist_ok=True)
    record_path = args.record_dir / "record.json"
    record_path.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if args.index:
        append_jsonl(args.index, {"record_path": str(record_path), **record})
    print(str(record_path))


if __name__ == "__main__":
    main()
