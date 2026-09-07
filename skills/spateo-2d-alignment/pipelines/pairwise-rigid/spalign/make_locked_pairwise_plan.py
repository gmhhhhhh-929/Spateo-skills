#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from lock_resolver import resolve_entrypoint


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a pairwise plan using a locked skill entrypoint.")
    parser.add_argument("--lock", required=True)
    parser.add_argument("--mode", required=True, choices=["normal_rigid", "spatial_only_rigid"])
    parser.add_argument("--planner", required=True, help="Path to make_spateo_pairwise_initial_plan.py")
    parser.add_argument("planner_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    resolved = resolve_entrypoint(Path(args.lock).expanduser().resolve(), args.mode)
    planner_args = list(args.planner_args)
    if planner_args and planner_args[0] == "--":
        planner_args = planner_args[1:]
    forbidden = {"--runner-script", "--expression-mode", "--dummy-rep-dim", "--sigma2-init-scale"}
    present_forbidden = [arg for arg in planner_args if arg in forbidden or any(arg.startswith(f"{flag}=") for flag in forbidden)]
    if present_forbidden:
        raise SystemExit(
            "Locked pairwise plan forbids overriding runner/expression parameters: "
            + ", ".join(present_forbidden)
        )
    cmd = [
        sys.executable,
        str(Path(args.planner).expanduser().resolve()),
        "--runner-script",
        resolved["entrypoint_path"],
        *planner_args,
    ]
    if args.mode == "spatial_only_rigid":
        cmd.extend(["--expression-mode", "spatial_only", "--dummy-rep-dim", "30", "--sigma2-init-scale", "2.0"])
    elif args.mode == "normal_rigid":
        cmd.extend(["--expression-mode", "normal"])
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
