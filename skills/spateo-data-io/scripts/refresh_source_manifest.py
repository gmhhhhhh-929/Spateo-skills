#!/usr/bin/env python3
"""Regenerate source signatures and IO hashes after reviewing a Spateo revision."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess


def collect(source):
    records = []
    for path in sorted((source / "spateo/io").rglob("*.py")):
        payload = path.read_bytes()
        tree = ast.parse(payload, filename=str(path))
        signatures = []
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and not node.name.startswith("_"):
                signatures.append("def " + node.name + "(" + ast.unparse(node.args) + ")")
            elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
                for method in node.body:
                    if isinstance(method, ast.FunctionDef) and not method.name.startswith("_"):
                        signatures.append("def " + node.name + "." + method.name + "(" + ast.unparse(method.args) + ")")
        records.append({"path": path.relative_to(source).as_posix(),
                        "sha256": hashlib.sha256(payload).hexdigest(), "signatures": signatures})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--check", action="store_true", help="Verify saved commit and every IO file; write nothing")
    args = parser.parse_args()
    source = args.source_root.resolve()
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    manifest = {"repository": "https://github.com/gmhhhhhh-929/spateo-release",
                "commit": revision, "files": collect(source)}
    target = Path(__file__).resolve().parents[1] / "references"
    if args.check:
        saved = json.loads((target / "source_manifest.json").read_text())
        if manifest != saved:
            raise SystemExit("Source revision or IO file hashes differ; review API changes before refreshing")
        print(json.dumps({"status": "pass", "commit": revision, "files_checked": len(manifest["files"])}))
        return 0
    # A clean committed snapshot prevents a dirty API being labelled as a Git commit.
    dirty = subprocess.check_output(["git", "-C", str(source), "status", "--porcelain", "--", "spateo/io"], text=True)
    if dirty.strip():
        raise SystemExit("Commit reviewed IO changes before regenerating the source manifest")
    (target / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    lines = ["# Current IO source API", "", f"Source commit: `{revision}`. Automatic APIs return SpatialReadResult.",
             "Hashes cover every Python file under `spateo/io`; signatures below are generated from source AST.", ""]
    for record in manifest["files"]:
        if record["signatures"]:
            lines += ["## " + record["path"], "", "```python", *record["signatures"], "```", ""]
    (target / "source-api.md").write_text("\n".join(lines))
    print(json.dumps({"status": "updated", "commit": revision, "files": len(manifest["files"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
