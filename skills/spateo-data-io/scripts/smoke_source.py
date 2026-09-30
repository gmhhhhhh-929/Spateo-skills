#!/usr/bin/env python3
"""Run pinned native-source IO/runtime tests plus skill IO behavior tests."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

COMMIT = json.loads(
    (Path(__file__).resolve().parents[1] / "references/source_manifest.json").read_text()
)["commit"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-root", type=Path, required=True)
    args = p.parse_args()
    source = args.source_root.resolve()
    actual = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    if actual != COMMIT:
        raise ValueError("Source revision differs from verified skill: " + actual)
    subprocess.check_call([sys.executable, str(Path(__file__).with_name("refresh_source_manifest.py")),
                           "--source-root", str(source), "--check"])
    env = dict(
        os.environ,
        PYTHONPATH=str(source),
        MPLBACKEND="Agg",
        PYTHONDONTWRITEBYTECODE="1",
    )
    tests = [
        source / "tests/io",
        source / "tests/test_native_runtime.py",
        Path(__file__).with_name("test_spateo_io.py"),
    ]
    packaged = Path(__file__).resolve().parents[3] / "tests/native_skills/test_io.py"
    if packaged.is_file():
        tests.append(packaged)
    return subprocess.call(
        [sys.executable, "-m", "pytest", *[str(t) for t in tests], "-q"],
        env=env,
        cwd=source,
    )


if __name__ == "__main__":
    raise SystemExit(main())
