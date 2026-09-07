#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shlex
import sys
from pathlib import Path
from typing import Any

from lock_resolver import load_yaml, resolve_entrypoint, sha256_file


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        obj = json.load(handle)
    if not isinstance(obj, dict):
        raise SystemExit(f"Invalid JSON object: {path}")
    return obj


def get_nested(obj: dict[str, Any], dotted: str) -> Any:
    cur: Any = obj
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


def parse_payload_env(path: Path) -> dict[str, str]:
    env: dict[str, str] = {}
    if not path.exists():
        return env
    export_re = re.compile(r"^\s*export\s+([A-Za-z_][A-Za-z0-9_]*)=(.*)$")
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = export_re.match(line)
        if not match:
            continue
        key, raw = match.groups()
        raw = raw.strip()
        try:
            parts = shlex.split(raw, posix=True)
            env[key] = parts[0] if parts else ""
        except ValueError:
            env[key] = raw.strip("\"'")
    return env


def runner_matches_locked_code(actual: str, expected_path: str, expected_sha: str) -> tuple[bool, str]:
    if actual == expected_path:
        return True, "exact locked path"
    actual_path = Path(actual).expanduser()
    if not actual_path.exists():
        return False, f"{actual} does not exist"
    actual_sha = sha256_file(actual_path)
    if actual_sha == expected_sha:
        return True, f"equivalent locked code hash at {actual}"
    return False, f"hash mismatch at {actual}: {actual_sha} vs {expected_sha}"


def values_match(actual: Any, expected: Any) -> bool:
    if isinstance(expected, float) or isinstance(actual, float):
        try:
            return math.isclose(float(actual), float(expected), rel_tol=0, abs_tol=1e-9)
        except Exception:
            return False
    return str(actual) == str(expected)


def validate_run(
    *,
    lock_path: Path,
    mode: str,
    payload_path: Path | None,
    edge_dir: Path | None,
    provenance_path: Path | None,
    edge_transform_path: Path | None,
) -> dict[str, Any]:
    resolved = resolve_entrypoint(lock_path, mode)
    checks: list[dict[str, Any]] = []

    def add(name: str, passed: bool, detail: str = "", **extra: Any) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail, **extra})

    payload_env: dict[str, str] = {}
    if payload_path:
        payload_env = parse_payload_env(payload_path)
        add("payload_exists", payload_path.exists(), str(payload_path))
        runner = payload_env.get("RUNNER_SCRIPT")
        if runner:
            runner_ok, runner_detail = runner_matches_locked_code(
                runner, resolved["entrypoint_path"], resolved["entrypoint_sha256"]
            )
        else:
            runner_ok, runner_detail = False, "RUNNER_SCRIPT missing"
        add("payload_runner_matches_locked_code", runner_ok, runner_detail)
        for key, expected in resolved["required_payload_env"].items():
            add(f"payload_env_{key}", payload_env.get(key) == str(expected), f"{payload_env.get(key)} vs {expected}")
    else:
        add("payload_provided", False, "No payload path was provided.")

    if edge_dir:
        provenance_path = provenance_path or edge_dir / "spateo_two_stage_provenance.json"
        edge_transform_path = edge_transform_path or edge_dir / "edge_transform.json"

    provenance: dict[str, Any] = {}
    if provenance_path and provenance_path.exists():
        provenance = read_json(provenance_path)
        add("provenance_exists", True, str(provenance_path))
    else:
        add("provenance_exists", False, str(provenance_path or ""))

    if provenance:
        runner_script = str(provenance.get("runner_script") or "")
        expected_base = resolved["required_provenance"].get("runner_script_basename")
        if expected_base:
            if runner_script:
                add("provenance_runner_basename", Path(runner_script).name == str(expected_base), f"{Path(runner_script).name} vs {expected_base}")
            else:
                add(
                    "provenance_runner_basename",
                    True,
                    "not present in provenance; payload and edge_transform runner checks are authoritative",
                )
        for dotted, expected in resolved["required_provenance"].items():
            if dotted == "runner_script_basename":
                continue
            actual = get_nested(provenance, dotted)
            add(f"provenance_{dotted}", values_match(actual, expected), f"{actual} vs {expected}")
        for dotted in resolved["forbidden_null_fields"]:
            actual = get_nested(provenance, dotted)
            add(f"provenance_non_null_{dotted}", actual is not None and actual != "", f"{actual}")

    if edge_transform_path and edge_transform_path.exists():
        edge = read_json(edge_transform_path)
        add("edge_transform_exists", True, str(edge_transform_path))
        edge_runner = str(edge.get("runner_script") or "")
        expected_base = resolved["required_provenance"].get("runner_script_basename")
        if expected_base:
            add("edge_transform_runner_basename", Path(edge_runner).name == str(expected_base), f"{Path(edge_runner).name} vs {expected_base}")
        params = edge.get("runner_parameters") or {}
        if mode == "spatial_only_rigid":
            add("edge_transform_expression_mode", params.get("expression_mode") == "spatial_only", f"{params.get('expression_mode')}")
            add("edge_transform_sigma2_init_scale", str(params.get("sigma2_init_scale")) == "2.0", f"{params.get('sigma2_init_scale')}")
        if mode == "normal_rigid":
            add("edge_transform_expression_mode", params.get("expression_mode") == "normal", f"{params.get('expression_mode')}")
    else:
        add("edge_transform_exists", False, str(edge_transform_path or ""))

    actual_entrypoint_path = Path(resolved["entrypoint_path"])
    add("locked_entrypoint_hash", sha256_file(actual_entrypoint_path) == resolved["entrypoint_sha256"], resolved["entrypoint_sha256"])

    status = "pass" if all(c["passed"] for c in checks) else "fail"
    return {
        "validator": "spalign.validate_run",
        "validator_version": "0.1.0",
        "status": status,
        "skill_name": resolved["skill_name"],
        "skill_version": resolved["skill_version"],
        "pipeline_release": resolved["pipeline_release"],
        "mode": mode,
        "entrypoint_path": resolved["entrypoint_path"],
        "entrypoint_sha256": resolved["entrypoint_sha256"],
        "payload_path": str(payload_path) if payload_path else "",
        "provenance_path": str(provenance_path) if provenance_path else "",
        "edge_transform_path": str(edge_transform_path) if edge_transform_path else "",
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate a locked spatial alignment run.")
    parser.add_argument("--lock", required=True)
    parser.add_argument("--mode", required=True)
    parser.add_argument("--payload")
    parser.add_argument("--edge-dir")
    parser.add_argument("--provenance")
    parser.add_argument("--edge-transform")
    parser.add_argument("--output")
    parser.add_argument("--fail-on-error", action="store_true")
    args = parser.parse_args()

    report = validate_run(
        lock_path=Path(args.lock).expanduser().resolve(),
        mode=args.mode,
        payload_path=Path(args.payload).expanduser().resolve() if args.payload else None,
        edge_dir=Path(args.edge_dir).expanduser().resolve() if args.edge_dir else None,
        provenance_path=Path(args.provenance).expanduser().resolve() if args.provenance else None,
        edge_transform_path=Path(args.edge_transform).expanduser().resolve() if args.edge_transform else None,
    )
    payload = json.dumps(report, indent=2, ensure_ascii=False)
    if args.output:
        Path(args.output).expanduser().resolve().write_text(payload + "\n", encoding="utf-8")
    print(payload)
    if args.fail_on_error and report["status"] != "pass":
        raise SystemExit(2)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(1)
