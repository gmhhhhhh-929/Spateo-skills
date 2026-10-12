#!/usr/bin/env python3
"""Discover, document and compare existing Spateo environments without installing packages.

The controller and inventory path use only the Python standard library. Scientific
checks run separately in the selected interpreter with bounded execution time.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
SCHEMA = "spateo-environment-audit/v1"
CAPABILITIES = {
    "io": ("H5AD IO", "H5AD round-trip only; platform readers require matching source files."),
    "preprocess": ("Spatial preprocessing", "Tiny preprocessing/PCA and native vector field; dataset scale is untested."),
    "slice-qc": ("Serial slice quality", "Installed API/signature check only; real ordered slices and keys must be validated."),
    "alignment": ("2D / 3D alignment", "API/signature check only; no full registration, accuracy or resource benchmark."),
    "point-cloud": ("3D point cloud", "Tiny point-cloud round-trip; real 3D coordinates and cell identity remain input requirements."),
    "mesh": ("3D surface mesh", "Tiny native mesh computation; no display, full reconstruction or optional Open3D guarantee."),
    "morphogenesis": ("Native morphogenesis", "Native field/API check only; transport, trajectories and full 4D pipeline remain untested."),
    "features": ("Morphogenesis features / GLM", "API check only; no real feature model or gene association fit."),
    "gpu": ("CUDA GPU tensor runtime", "Tiny CUDA tensor check only; MPS, Spateo GPU pipelines and throughput are untested."),
}
STANDARD = ("io", "preprocess", "slice-qc", "alignment", "point-cloud", "morphogenesis", "features")


def save_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def new_output(path):
    path = Path(path).expanduser().resolve()
    if path.exists() and any(path.iterdir()):
        raise ValueError(f"Output must be a new or empty directory: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def stamp():
    return datetime.now(timezone.utc).isoformat()


def executable(value):
    # Preserve the selected venv's bin/python symlink: resolving it can select base Python.
    found = shutil.which(os.path.expanduser(value))
    return os.path.abspath(found or os.path.expanduser(value))


def run_command(command, timeout, cwd, env, log):
    start = time.monotonic()
    timed_out = False
    try:
        # File-backed logs prevent package/native chatter from filling controller memory.
        with Path(log).open("wb") as handle:
            process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT,
                                       cwd=cwd, env=env, start_new_session=(os.name == "posix"))
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                if os.name == "posix":
                    os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
                code = process.wait()
        return {"returncode": code, "timed_out": timed_out,
                "elapsed_seconds": round(time.monotonic() - start, 3), "log": str(Path(log).name)}
    except OSError as exc:
        Path(log).write_text(f"{type(exc).__name__}: {exc}\n", encoding="utf-8")
        return {"returncode": None, "timed_out": False, "error": str(exc),
                "elapsed_seconds": round(time.monotonic() - start, 3), "log": str(Path(log).name)}


def child_environment(scratch):
    env = dict(os.environ)
    for name in ("MPLCONFIGDIR", "NUMBA_CACHE_DIR", "XDG_CACHE_HOME"):
        path = Path(scratch) / name.lower()
        path.mkdir()
        env[name] = str(path)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"
    env["MPLBACKEND"] = "Agg"
    # Bound CPU oversubscription without changing user's persistent shell settings.
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[name] = "1"
    return env


def probe(name, python, output, scratch, env, timeout):
    result_path = output / "checks" / f"{name}.json"
    log_path = output / "checks" / f"{name}.log"
    execution = run_command([python, "-B", str(HERE / "environment_probes.py"),
                             "--probe", name, "--result-file", str(result_path)],
                            timeout, scratch, env, log_path)
    try:
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if not isinstance(result, dict) or "status" not in result:
            raise ValueError("probe result has no status")
    except (OSError, ValueError) as exc:
        result = {"status": "failed", "error": f"No valid probe result: {exc}"}
    if execution["timed_out"]:
        result.update(status="failed", error=f"Timed out after {timeout} seconds")
    elif execution["returncode"] != 0:
        result.update(status="failed", error=result.get("error") or "Probe process failed")
    result["execution"] = execution
    save_json(result_path, result)
    return result


def bounded_tree(root, depth=3, max_entries=500, per_directory=60, excludes=()):
    """Never recurse through symlinks; record listing failures and truncation."""
    root = Path(root)
    lines = [str(root)]
    state = {"entries": 0, "truncated": False, "errors": [], "depth_limited": False}
    excluded = {str(Path(p).resolve()) for p in excludes}

    def visit(folder, level, prefix):
        if level >= depth:
            state["depth_limited"] = True
            lines.append(prefix + "... [depth limit]")
            return
        try:
            with os.scandir(folder) as iterator:
                entries = list(itertools.islice(iterator, per_directory + 1))
        except OSError as exc:
            state["errors"].append(f"{folder}: {exc}")
            lines.append(prefix + f"[unreadable: {type(exc).__name__}]")
            return
        overflow = len(entries) > per_directory
        for entry in sorted(entries[:per_directory], key=lambda e: e.name.lower()):
            if state["entries"] >= max_entries:
                state["truncated"] = True
                lines.append(prefix + "... [total entry limit]")
                return
            state["entries"] += 1
            try:
                if entry.is_symlink():
                    lines.append(prefix + entry.name + " -> " + os.readlink(entry.path) + " [not followed]")
                elif entry.is_dir(follow_symlinks=False):
                    lines.append(prefix + entry.name + "/")
                    if entry.name in {".git", "__pycache__", ".cache"} or str(Path(entry.path).resolve()) in excluded:
                        lines.append(prefix + "  ... [excluded]")
                    else:
                        visit(entry.path, level + 1, prefix + "  ")
                else:
                    lines.append(prefix + entry.name)
            except OSError as exc:
                state["errors"].append(f"{entry.path}: {exc}")
                lines.append(prefix + entry.name + " [unreadable]")
        if overflow:
            state["truncated"] = True
            lines.append(prefix + f"... [directory limit: first {per_directory} entries, sorted for display]")

    if root.is_symlink():
        lines[0] += " -> " + os.readlink(root) + " [not followed]"
    elif root.is_dir():
        visit(root, 0, "  ")
    else:
        state["errors"].append(f"Not an accessible directory: {root}")
    return {"root": str(root), "text": "\n".join(lines), "limits": {
        "depth": depth, "max_entries": max_entries, "per_directory": per_directory}, **state}


def inventory_details(result):
    return result.get("details", {}) if result.get("status") == "passed" else {}


def write_trees(output, inventory, checkout, depth, max_entries):
    roots = []
    prefix = inventory.get("prefix")
    if prefix:
        roots.append(("Environment prefix", prefix))
    # Support both conventional sysconfig field names while preserving raw inventory.
    paths = inventory.get("sysconfig_paths", inventory.get("paths", {}))
    for name in ("purelib", "platlib", "scripts", "include"):
        if paths.get(name):
            roots.append((name, paths[name]))
    distribution = inventory.get("spateo_distribution") or {}
    source = checkout or distribution.get("editable_checkout")
    if source:
        roots.append(("Spateo source (separate from environment)", source))
    seen, trees = set(), []
    for label, root in roots:
        if root in seen:
            continue
        seen.add(root)
        trees.append({"label": label, **bounded_tree(root, depth, max_entries, excludes=(output,))})
    save_json(output / "environment-tree.json", trees)
    lines = ["# Environment directory map", "", "Generated: " + stamp(), "",
             "Bounded directory listing; file contents are not read and symlinks are not followed.",
             "Limits apply per root. Depth, entry caps and unreadable directories are recorded explicitly.", ""]
    for tree in trees:
        lines += ["## " + tree["label"], "", "```text", tree["text"], "```", "",
                  f"Entries: {tree['entries']}; truncated: {tree['truncated']}; depth limited: {tree['depth_limited']}.", ""]
    if not trees:
        lines += ["No directory roots could be obtained. See checks/inventory.json."]
    (output / "environment-tree.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return trees


def get_source_identity(inventory, checkout, output, scratch, env):
    distribution = inventory.get("spateo_distribution") or {}
    path = checkout or distribution.get("editable_checkout")
    result = {"requested_checkout": checkout, "editable_checkout": distribution.get("editable_checkout")}
    if path:
        source = Path(path).expanduser().resolve()
        result["path"] = str(source)
        if (source / ".git").exists():
            for key, arguments in (("commit", ["rev-parse", "HEAD"]),
                                   ("tracked_changes", ["status", "--porcelain", "--untracked-files=no"])):
                log = output / "checks" / f"source-{key}.log"
                run = run_command(["git", "--no-optional-locks", "-C", str(source)] + arguments, 10, scratch, env, log)
                if run["returncode"] == 0:
                    result[key] = log.read_text(encoding="utf-8", errors="replace").strip()
                else:
                    result[key + "_error"] = run
    return result


def capability_status(result):
    if result is None:
        return "not_tested"
    if result.get("status") == "passed" and result.get("level") == "smoke":
        return "smoke_passed"
    if result.get("status") in ("passed", "api_only") and result.get("level") == "api":
        return "api_only"
    if result.get("status") == "passed" and result.get("level") == "metadata":
        return "not_tested"
    return "blocked"


def escape(value):
    return str(value).replace("|", "\\|").replace("\n", " ")


def render_report(report):
    identity = report["identity"]
    lines = ["# Spateo environment report", "", f"Generated (UTC): {report['generated_at']}", "",
             f"Overall: **{report['overall']}** · Profile: `{report['profile']}`", "",
             f"Selected interpreter: `{report['selected_python']}`", "",
             f"Actual interpreter: `{identity.get('executable', 'unavailable')}`", "",
             f"Python: `{identity.get('python', 'unavailable')}` · Prefix: `{identity.get('prefix', 'unavailable')}`", "",
             f"Source: `{report['source'].get('path', 'not resolved')}`", "",
             f"Source commit: `{report['source'].get('commit', 'not recorded')}`", "",
             f"Spateo distribution: `{(identity.get('spateo_distribution') or {}).get('version', 'unavailable')}` · Host: `{identity.get('hostname', 'unavailable')}` · CPU threads: `{identity.get('cpu_count', 'unavailable')}`", "",
             "## Health checks", "", "| Check | Status | Evidence |", "| --- | --- | --- |"]
    for name, item in report["health"].items():
        evidence = escape(item.get("reason", ""))
        if item.get("log"):
            evidence = f"[Raw log](checks/{item['log']})"
        elif name == "inventory":
            evidence = "[Raw inventory](checks/inventory.json)"
        lines.append(f"| {escape(name)} | {item['status']} | {evidence} |")
    lines += ["", "## Downstream capability evidence", "",
              "`smoke_passed` = bounded synthetic calculation passed; `api_only` = interface check only; "
              "`blocked` = missing/failed prerequisite or probe; `not_tested` = no execution evidence.", "",
              "| Task | Evidence level | Boundary / remaining input |", "| --- | --- | --- |"]
    for name, item in report["capabilities"].items():
        detail = item.get("error") or item["boundary"]
        label = f"[{CAPABILITIES[name][0]}]({item['evidence']})" if item.get('evidence') else CAPABILITIES[name][0]
        lines.append(f"| {label} | {item['status']} | {escape(detail)} |")
    lines += ["", "A passing task probe does not override a failed dependency or source-identity check. "
              "API presence does not establish a successful full pipeline. GPU/display, dataset inputs, "
              "memory, biological accuracy and pipeline-scale performance need their own validation.", "",
              "## Next actions", ""]
    lines += ["- " + action for action in report["next_actions"]]
    lines += ["", "## Files", "", "- [Directory map](environment-tree.md)",
              "- [Machine-readable report](environment-report.json)", "- [Package snapshot](packages.json)",
              "- Raw results and logs: `checks/`", "",
              "No package install/uninstall/update is performed. Temporary caches and synthetic fixtures are isolated. "
              "Package versions are a snapshot, not a cross-platform lockfile. Machine paths are local diagnostics; "
              "review before sharing publicly.", ""]
    return "\n".join(lines)


def tool_hashes():
    return {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in ("audit_environment.py", "environment_probes.py", "verify_environment.py")}


def audit(args):
    implementation_at_start = tool_hashes()
    output = new_output(args.output)
    (output / "checks").mkdir()
    python = executable(args.python)
    with tempfile.TemporaryDirectory(prefix="spateo-audit-") as scratch:
        env = child_environment(scratch)
        inventory_result = probe("inventory", python, output, scratch, env, args.timeout)
        inventory = inventory_details(inventory_result)
        save_json(output / "packages.json", inventory.get("packages", []))
        trees = write_trees(output, inventory, args.checkout, args.depth, args.max_entries)
        if args.command == "tree":
            print(output / "environment-tree.md")
            return 0 if inventory and trees and not any(t["errors"] for t in trees) else 1
        source = get_source_identity(inventory, args.checkout, output, scratch, env)
        health = {"inventory": {"status": "passed" if inventory else "failed",
                                "reason": inventory_result.get("error", "checks/inventory.json")}}
        pip_execution = run_command([python, "-B", "-m", "pip", "check"], args.timeout, scratch, env,
                                    output / "checks" / "pip-check.log")
        health["pip_check"] = {"status": "passed" if pip_execution["returncode"] == 0 else "failed", **pip_execution}
        verify_command = [python, "-B", str(HERE / "verify_environment.py"), "--json"]
        if args.profile == "metadata":
            verify_command.append("--metadata-only")
        expected_checkout = args.checkout or (inventory.get("spateo_distribution") or {}).get("editable_checkout")
        if expected_checkout:
            verify_command += ["--checkout", str(Path(expected_checkout).expanduser().resolve())]
        verification = run_command(verify_command, args.timeout, scratch, env, output / "checks" / "verification.log")
        # The verifier sends import chatter to stderr; raw mixed log is retained and its exit code is authoritative.
        health["source_profile_verification"] = {"status": "passed" if verification["returncode"] == 0 else "failed", **verification}
        checks = {}
        requested = () if args.profile == "metadata" else STANDARD
        if args.profile == "full":
            requested += ("mesh", "gpu")
        for name in requested:
            print(f"Checking {name}...", file=sys.stderr, flush=True)
            checks[name] = probe(name, python, output, scratch, env, args.timeout)
        capabilities = {}
        for name, (_, boundary) in CAPABILITIES.items():
            result = checks.get(name)
            capabilities[name] = {"status": capability_status(result), "boundary": boundary}
            if result:
                capabilities[name]["evidence"] = f"checks/{name}.json"
                if result.get("error"):
                    capabilities[name]["error"] = result["error"]
                elif result.get("status") == "unsupported":
                    capabilities[name]["error"] = result.get("details", {}).get("reason", "Capability unsupported")
        source["imported_paths"] = sorted({r.get("details", {}).get("spateo_import_path") for r in checks.values()
                                           if r.get("details", {}).get("spateo_import_path")})
        if args.profile != "metadata" and source["imported_paths"]:
            metadata_root = (inventory.get("spateo_distribution") or {}).get("metadata_location")
            expected_root = Path(expected_checkout or metadata_root).resolve() / "spateo" if (expected_checkout or metadata_root) else None
            mismatches = []
            for imported in source["imported_paths"]:
                try:
                    Path(imported).resolve().relative_to(expected_root)
                except (ValueError, TypeError):
                    mismatches.append(imported)
            health["import_identity"] = {"status": "failed" if mismatches else "passed",
                                         "reason": "Imported modules outside installed source: " + ", ".join(mismatches) if mismatches else "All measured Spateo imports match installed source metadata."}
        health["audit_implementation"] = {
            "status": "passed" if implementation_at_start == tool_hashes() else "failed",
            "reason": "Audit scripts must remain unchanged during the run."}
        failed_health = [k for k, v in health.items() if v["status"] != "passed"]
        blocked = [k for k, v in capabilities.items() if v["status"] == "blocked"]
        overall = "blocked" if failed_health else ("metadata_only" if args.profile == "metadata" else
                                                   "partial" if blocked else "checked_with_limits")
        actions = []
        if failed_health:
            actions.append("Inspect failed health checks before starting a new pipeline: " + ", ".join(failed_health) + ".")
        unsupported = [name for name in blocked if checks[name].get("status") == "unsupported"]
        failed_probes = [name for name in blocked if name not in unsupported]
        if failed_probes:
            actions.append("Inspect the raw error and resolve failed capability probes: " + ", ".join(failed_probes) + ".")
        if unsupported:
            actions.append("Unavailable optional capability: " + ", ".join(unsupported) + ". Use a supported CPU path or select an appropriate GPU host; CPU evidence is reported separately.")
        if args.profile == "metadata":
            actions.append("Run the standard profile to collect runtime evidence; metadata alone cannot establish readiness.")
        actions.append("Validate the selected downstream skill's source revision, data keys and compute requirements on real inputs.")
        if any(v["status"] == "not_tested" for v in capabilities.values()):
            actions.append("Use the full profile only when mesh and GPU evidence are needed; display rendering remains untested.")
        report = {"schema_version": SCHEMA, "generated_at": stamp(), "profile": args.profile,
                  "selected_python": python, "identity": inventory, "source": source,
                  "health": health, "capabilities": capabilities, "overall": overall,
                  "next_actions": actions, "execution_context": {
                      "neutral_working_directory": True, "pythonpath_present": bool(os.environ.get("PYTHONPATH")),
                      "cache_policy": "temporary", "timeout_seconds_per_check": args.timeout},
                  "tool_sha256": implementation_at_start}
        save_json(output / "environment-report.json", report)
        (output / "environment-report.md").write_text(render_report(report), encoding="utf-8")
        if args.baseline:
            before = json.loads(Path(args.baseline).read_text(encoding="utf-8"))
            write_comparison(before, report, output)
        print(output / "environment-report.md")
        return 1 if overall in ("blocked", "partial") else 0


def package_map(report):
    packages = report.get("identity", {}).get("packages", [])
    return {p["name"].lower().replace("_", "-"): p["version"] for p in packages}


def compare_reports(before, after):
    for report in (before, after):
        if report.get("schema_version") != SCHEMA:
            raise ValueError(f"Expected report schema {SCHEMA}")
    complete_inventory = all(r.get("health", {}).get("inventory", {}).get("status") == "passed"
                             and isinstance(r.get("identity", {}).get("packages"), list)
                             and r.get("identity", {}).get("prefix")
                             and r.get("identity", {}).get("executable") for r in (before, after))
    a, b = (package_map(before), package_map(after)) if complete_inventory else ({}, {})
    changes = {"schema_version": "spateo-environment-comparison/v1", "generated_at": stamp(),
               "before_generated_at": before.get("generated_at"), "after_generated_at": after.get("generated_at"),
               "packages_added": {k: b[k] for k in b.keys() - a.keys()},
               "packages_removed": {k: a[k] for k in a.keys() - b.keys()},
               "packages_changed": {k: {"before": a[k], "after": b[k]} for k in a.keys() & b.keys() if a[k] != b[k]},
               "identity_changes": {}, "source_changes": {}, "health_changes": {}, "capability_changes": {},
               "profile_changed": before.get("profile") != after.get("profile"),
               "comparison_status": "complete" if complete_inventory else "incomplete_inventory",
               "tool_changed": before.get("tool_sha256") != after.get("tool_sha256")}
    for section, keys in (("identity", ("hostname", "executable", "prefix", "python", "platform")),
                          ("source", ("path", "commit", "tracked_changes", "editable_checkout"))):
        if section == "identity" and not complete_inventory:
            continue
        for key in keys:
            x, y = before.get(section, {}).get(key), after.get(section, {}).get(key)
            if x != y:
                changes[section + "_changes"][key] = {"before": x, "after": y}
    for section, field in (("health", "health_changes"), ("capabilities", "capability_changes")):
        for key in before.get(section, {}).keys() | after.get(section, {}).keys():
            x = before.get(section, {}).get(key, {}).get("status")
            y = after.get(section, {}).get(key, {}).get("status")
            if x != y:
                changes[field][key] = {"before": x, "after": y}
    changes["same_environment_identity"] = not changes["identity_changes"] if complete_inventory else None
    if not complete_inventory:
        for key in ("packages_added", "packages_removed", "packages_changed"):
            changes[key] = None
    return changes


def write_comparison(before, after, output):
    changes = compare_reports(before, after)
    save_json(output / "environment-comparison.json", changes)
    lines = ["# Environment comparison", "", "Version and evidence comparison; no packages or environments were modified.", "",
             "Different interpreters/hosts are a cross-environment comparison, not proven temporal drift. "
             "Different profiles or tool revisions can change evidence without a dependency regression.", "",
             f"Comparison coverage: {changes['comparison_status']} (null means unknown). Same environment identity: {changes['same_environment_identity']}; profile changed: {changes['profile_changed']}; tool changed: {changes['tool_changed']}.", ""]
    for section in ("identity_changes", "source_changes", "packages_added", "packages_removed", "packages_changed", "health_changes", "capability_changes"):
        lines += ["## " + section.replace("_", " "), "", "```json", json.dumps(changes[section], ensure_ascii=False, indent=2, sort_keys=True), "```", ""]
    (output / "environment-comparison.md").write_text("\n".join(lines), encoding="utf-8")
    return changes


def discover(args):
    output = new_output(args.output)
    candidates, errors = {}, []

    def add(prefix, evidence, python_override=None):
        prefix = Path(prefix).expanduser().absolute()
        py = Path(python_override) if python_override else prefix / ("python.exe" if os.name == "nt" else "bin/python")
        if not python_override and not py.exists() and (prefix / "bin/python3").is_file():
            py = prefix / "bin/python3"
        if not py.is_file():
            # stale registry entries remain useful diagnostics
            present = False
        else:
            present = True
        key = str(prefix)
        record = candidates.setdefault(key, {"prefix": key, "python": str(py), "interpreter_exists": present,
                                             "evidence": [], "provenance": "setup-skill origin not established", "readiness": "not_tested"})
        if evidence not in record["evidence"]:
            record["evidence"].append(evidence)

    add(sys.prefix, "controller interpreter", sys.executable)
    for key in ("CONDA_PREFIX", "VIRTUAL_ENV"):
        if os.environ.get(key):
            add(os.environ[key], key)
    registry = Path.home() / ".conda/environments.txt"
    if registry.is_file():
        try:
            for line in registry.read_text(encoding="utf-8").splitlines():
                if line.strip():
                    add(line.strip(), "conda environment registry")
        except OSError as exc:
            errors.append(str(exc))
    scanned = 0
    # Traverse only explicitly selected roots, bounded by depth and a global entry cap.
    def scan(root, depth):
        nonlocal scanned
        if scanned >= 2000:
            return
        root = Path(root).expanduser()
        if root.is_symlink():
            return
        if any((root / name).is_file() for name in ("bin/python", "bin/python3", "python.exe")):
            add(root, "explicit bounded search root")
            # Conda base can own envs; avoid recursively scanning its packages.
            if (root / "envs").is_dir() and depth > 0:
                scan(root / "envs", depth - 1)
            return
        if depth <= 0:
            return
        try:
            with os.scandir(root) as iterator:
                for entry in iterator:
                    if scanned >= 2000:
                        break
                    scanned += 1
                    if entry.is_dir(follow_symlinks=False) and entry.name not in {".git", "site-packages", "pkgs", "node_modules"}:
                        scan(Path(entry.path), depth - 1)
        except OSError as exc:
            errors.append(f"{root}: {exc}")
    for root in args.search_root:
        scan(root, 4)
    result = {"schema_version": "spateo-environment-discovery/v1", "generated_at": stamp(),
              "scope": "local registry, active interpreter and explicit roots only; no SSH or filesystem-wide scan",
              "candidates": sorted(candidates.values(), key=lambda p: p["prefix"]),
              "errors": errors, "search_entries": scanned, "search_entry_limit": 2000,
              "search_depth_limit": 4, "search_truncated": scanned >= 2000}
    save_json(output / "environment-discovery.json", result)
    lines = ["# Environment candidates", "", result["scope"], "",
             "Candidates are not proof of Spateo installation or prior skill use. Verify selected interpreters with `report`.", "",
             "| Python | Exists | Discovery evidence |", "| --- | --- | --- |"]
    for item in result["candidates"]:
        lines.append(f"| `{escape(item['python'])}` | {item['interpreter_exists']} | {escape(', '.join(item['evidence']))} |")
    lines += ["", f"Search truncated: {result['search_truncated']}; errors: {len(errors)}.", ""]
    (output / "environment-discovery.md").write_text("\n".join(lines), encoding="utf-8")
    print(output / "environment-discovery.md")
    return 0


def positive(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    find = commands.add_parser("discover", help="List candidates without importing scientific packages or using SSH")
    find.add_argument("--search-root", action="append", default=[])
    find.add_argument("--output", required=True)
    for command in ("tree", "report"):
        sub = commands.add_parser(command)
        sub.add_argument("--python", required=True, help="Selected environment's interpreter (venv symlinks preserved)")
        sub.add_argument("--output", required=True)
        sub.add_argument("--checkout")
        sub.add_argument("--timeout", type=positive, default=120, help="Timeout in seconds for each isolated process")
        sub.add_argument("--depth", type=positive, default=3)
        sub.add_argument("--max-entries", type=positive, default=500)
        if command == "report":
            sub.add_argument("--profile", choices=("metadata", "standard", "full"), default="standard")
            sub.add_argument("--baseline", help="Earlier environment-report.json for comparison")
    diff = commands.add_parser("compare")
    diff.add_argument("--before", required=True)
    diff.add_argument("--after", required=True)
    diff.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "discover":
            return discover(args)
        if args.command in ("tree", "report"):
            return audit(args)
        before = json.loads(Path(args.before).read_text(encoding="utf-8"))
        after = json.loads(Path(args.after).read_text(encoding="utf-8"))
        comparison_output = new_output(args.output)
        write_comparison(before, after, comparison_output)
        print(comparison_output / "environment-comparison.md")
        return 0
    except (OSError, ValueError) as exc:
        print(f"Audit error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
