#!/usr/bin/env python3
"""Create a child Spateo config using dotted JSON assignments."""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path
from typing import Any


def read_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.expanduser().read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def load_base(path: Path) -> dict[str, Any]:
    value = read_object(path)
    config = value.get("config", value)
    if not isinstance(config, dict):
        raise ValueError(f"Manifest config is not an object: {path}")
    config = copy.deepcopy(config)
    for key in ("stage1", "stage2"):
        raw = config.get("inputs", {}).get(key)
        if raw:
            source = Path(raw).expanduser()
            config["inputs"][key] = str(
                (path.resolve().parent / source).resolve()
                if not source.is_absolute()
                else source.resolve()
            )
    return config


def parse_assignment(raw: str) -> tuple[list[str], Any]:
    key, separator, encoded = raw.partition("=")
    parts = [part for part in key.split(".") if part]
    if not separator or not parts:
        raise ValueError(f"Expected dotted.path=JSON, got {raw!r}")
    try:
        value = json.loads(encoded)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Value for {key!r} must be JSON; strings need quotes"
        ) from exc
    return parts, value


def set_path(config: dict[str, Any], parts: list[str], value: Any) -> None:
    cursor = config
    for part in parts[:-1]:
        existing = cursor.get(part)
        if existing is None:
            existing = {}
            cursor[part] = existing
        if not isinstance(existing, dict):
            raise ValueError(
                f"Cannot descend through non-object config path: {'.'.join(parts)}"
            )
        cursor = existing
    cursor[parts[-1]] = value


def delete_path(config: dict[str, Any], dotted: str) -> None:
    parts = [part for part in dotted.split(".") if part]
    if not parts:
        raise ValueError("Delete path cannot be empty")
    cursor = config
    for part in parts[:-1]:
        child = cursor.get(part)
        if not isinstance(child, dict):
            return
        cursor = child
    cursor.pop(parts[-1], None)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--set",
        dest="assignments",
        action="append",
        default=[],
        metavar="PATH=JSON",
    )
    parser.add_argument("--delete", action="append", default=[], metavar="PATH")
    return parser


def main(argv: list[str] | None = None) -> int:
    try:
        args = build_parser().parse_args(argv)
        config = load_base(args.base)
        for raw in args.assignments:
            parts, value = parse_assignment(raw)
            set_path(config, parts, value)
        for path in args.delete:
            delete_path(config, path)
        output = args.output.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite config: {output}")
        output.write_text(
            json.dumps(config, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(json.dumps({"output": str(output), "config": config}, ensure_ascii=False))
        return 0
    except Exception as exc:  # noqa: BLE001 - CLI boundary
        print(
            json.dumps(
                {"error": type(exc).__name__, "message": str(exc)}, ensure_ascii=False
            ),
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
