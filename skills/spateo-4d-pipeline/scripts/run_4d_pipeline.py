#!/usr/bin/env python3
"""Execute the native, versioned Spateo 4D pipeline from a v3 JSON config."""
import argparse
import json
import sys
from pipeline_runtime import STAGES, execute


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--parent-manifest")
    parser.add_argument("--run-id")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--stop-after", choices=STAGES, default="dashboard")
    args = parser.parse_args(argv)
    try:
        result = execute(
            args.config,
            args.project,
            args.parent_manifest,
            args.run_id,
            args.dry_run,
            args.stop_after,
        )
        print(json.dumps(result, indent=2))
        return 0
    except Exception as exc:
        print(json.dumps({"status": "error", "error": str(exc)}), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
