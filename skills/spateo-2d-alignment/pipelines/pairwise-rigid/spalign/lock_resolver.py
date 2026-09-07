#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any


def load_yaml(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    try:
        obj = json.loads(text)
        if isinstance(obj, dict):
            return obj
    except json.JSONDecodeError:
        pass
    try:
        import yaml  # type: ignore
    except Exception as exc:
        raise SystemExit("PyYAML is required for non-JSON YAML lock files. Use JSON-compatible skill.lock.yaml.") from exc
    obj = yaml.safe_load(text)
    if not isinstance(obj, dict):
        raise SystemExit(f"Invalid YAML object: {path}")
    return obj


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_release_root(lock_path: Path, release_root_value: str) -> Path:
    release_root = Path(str(release_root_value)).expanduser()
    if not release_root.is_absolute():
        release_root = lock_path.parent / release_root
    return release_root.resolve()


def resolve_entrypoint(lock_path: Path, mode: str) -> dict[str, Any]:
    lock = load_yaml(lock_path)
    entrypoints = lock.get("allowed_entrypoints") or {}
    if mode not in entrypoints:
        raise SystemExit(f"Mode {mode!r} is not declared in {lock_path}")
    entry = dict(entrypoints[mode])
    release_root = resolve_release_root(lock_path, str(lock["release_root"]))
    path = (release_root / str(entry["relative_path"])).resolve()
    if not path.exists():
        raise SystemExit(f"Locked entrypoint does not exist: {path}")
    actual_sha = sha256_file(path)
    expected_sha = str(entry["sha256"])
    if actual_sha != expected_sha:
        raise SystemExit(
            f"Locked entrypoint hash mismatch for {mode}: expected {expected_sha}, got {actual_sha} at {path}"
        )
    return {
        "skill_name": lock.get("skill_name"),
        "skill_version": lock.get("skill_version"),
        "pipeline_release": lock.get("pipeline_release"),
        "release_root": str(release_root),
        "mode": mode,
        "entrypoint_path": str(path),
        "entrypoint_sha256": actual_sha,
        "conda_env": lock.get("conda_env", ""),
        "required_payload_env": entry.get("required_payload_env") or {},
        "required_provenance": entry.get("required_provenance") or {},
        "forbidden_null_fields": entry.get("forbidden_null_fields") or [],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Resolve a locked skill entrypoint.")
    parser.add_argument("--lock", required=True)
    parser.add_argument("--mode", required=True)
    parser.add_argument("--json", action="store_true", help="Print full JSON resolution.")
    parser.add_argument("--field", choices=["entrypoint_path", "entrypoint_sha256", "conda_env"])
    args = parser.parse_args()

    resolved = resolve_entrypoint(Path(args.lock).expanduser().resolve(), args.mode)
    if args.field:
        print(resolved[args.field])
    elif args.json:
        print(json.dumps(resolved, indent=2, ensure_ascii=False))
    else:
        print(resolved["entrypoint_path"])


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(1)
