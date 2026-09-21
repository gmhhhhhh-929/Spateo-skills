#!/usr/bin/env python3
"""Inspect or plan native 4D runs; stage mutation belongs to the runner."""
import argparse, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
from pipeline_runtime import read_json, execute

p = argparse.ArgumentParser(description=__doc__)
s = p.add_subparsers(dest="command", required=True)
a = s.add_parser("show")
a.add_argument("--manifest", required=True)
b = s.add_parser("plan")
b.add_argument("--config", required=True)
b.add_argument("--project", required=True)
b.add_argument("--parent-manifest")
if __name__ == "__main__":
    args = p.parse_args()
    value = (
        read_json(args.manifest)
        if args.command == "show"
        else execute(args.config, args.project, args.parent_manifest, dry_run=True)
    )
    print(json.dumps(value, indent=2))
