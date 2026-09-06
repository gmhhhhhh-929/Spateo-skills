#!/usr/bin/env python3
"""Verify the packaged source snapshot hashes without importing scientific libraries."""
import hashlib
import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "provenance/source_migration.json").read_text())
    mismatches = []
    for entry in manifest["files"]:
        path = root / entry["path"]
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        if actual != entry["migrated_sha256"]:
            mismatches.append({"path": entry["path"], "expected": entry["migrated_sha256"], "actual": actual})
    print(json.dumps({"status": "failed" if mismatches else "passed", "source_files_checked": len(manifest["files"]),
                      "mismatches": mismatches, "scope": "Migrated source code files; no data or GPU execution"}, indent=2))
    return 1 if mismatches else 0


if __name__ == "__main__":
    raise SystemExit(main())
