#!/usr/bin/env python3
"""Controller invariants with temporary fixtures; no scientific dependencies."""
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("audit", Path(__file__).with_name("audit_environment.py"))
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class AuditTests(unittest.TestCase):
    def test_tree_does_not_follow_symlink_cycle_or_external_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "env"
            root.mkdir()
            (root / "cycle").symlink_to(root, target_is_directory=True)
            (root / "external").symlink_to(Path(tmp), target_is_directory=True)
            result = audit.bounded_tree(root)
            self.assertEqual(result["entries"], 2)
            self.assertEqual(result["text"].count("[not followed]"), 2)

    def test_tree_has_global_entry_budget_and_explicit_truncation(self):
        with tempfile.TemporaryDirectory() as tmp:
            for i in range(30):
                (Path(tmp) / str(i)).mkdir()
                (Path(tmp) / str(i) / "data").write_text("do not read")
            result = audit.bounded_tree(tmp, depth=5, max_entries=7)
            self.assertEqual(result["entries"], 7)
            self.assertTrue(result["truncated"])
            self.assertNotIn("do not read", result["text"])

    def test_existing_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            file = Path(tmp) / "keep"
            file.write_text("user data")
            with self.assertRaises(ValueError):
                audit.new_output(tmp)
            self.assertEqual(file.read_text(), "user data")

    def test_venv_python_symlink_identity_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            py = Path(tmp) / "python"
            py.symlink_to(sys.executable)
            self.assertEqual(audit.executable(str(py)), str(py))

    def test_timeout_retains_log_and_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "log"
            result = audit.run_command([sys.executable, "-c", "import time; print('started', flush=True); time.sleep(30)"],
                                       0.1, tmp, os.environ.copy(), log)
            self.assertTrue(result["timed_out"])
            self.assertNotEqual(result["returncode"], 0)
            self.assertIn("started", log.read_text())

    def test_missing_interpreter_is_structured_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = audit.run_command([str(Path(tmp) / "missing-python")], 1, tmp, os.environ.copy(), Path(tmp) / "log")
            self.assertIsNone(result["returncode"])
            self.assertIn("error", result)

    def test_api_and_metadata_never_become_smoke_passes(self):
        self.assertEqual(audit.capability_status(None), "not_tested")
        self.assertEqual(audit.capability_status({"status": "api_only", "level": "api"}), "api_only")
        self.assertEqual(audit.capability_status({"status": "passed", "level": "metadata"}), "not_tested")
        self.assertEqual(audit.capability_status({"status": "failed", "level": "smoke"}), "blocked")
        self.assertEqual(audit.capability_status({"status": "passed", "level": "smoke"}), "smoke_passed")

    def test_comparison_reports_changes_without_mutation(self):
        before = {"schema_version": audit.SCHEMA, "profile": "full", "identity": {"prefix": "/a", "executable": "/a/bin/python", "packages": [
            {"name": "NumPy", "version": "1"}, {"name": "removed", "version": "1"}]},
            "health": {"inventory": {"status": "passed"}, "pip": {"status": "passed"}}, "capabilities": {"io": {"status": "smoke_passed"}}}
        after = {"schema_version": audit.SCHEMA, "profile": "metadata", "identity": {"prefix": "/b", "executable": "/b/bin/python", "packages": [
            {"name": "numpy", "version": "2"}, {"name": "added", "version": "3"}]},
            "health": {"inventory": {"status": "passed"}, "pip": {"status": "failed"}}, "capabilities": {"io": {"status": "not_tested"}}}
        original = json.dumps(before, sort_keys=True)
        result = audit.compare_reports(before, after)
        self.assertEqual(result["packages_changed"]["numpy"], {"before": "1", "after": "2"})
        self.assertEqual(result["packages_added"], {"added": "3"})
        self.assertEqual(result["packages_removed"], {"removed": "1"})
        self.assertFalse(result["same_environment_identity"])
        self.assertTrue(result["profile_changed"])
        self.assertEqual(result["health_changes"]["pip"]["after"], "failed")
        self.assertEqual(json.dumps(before, sort_keys=True), original)

    def test_failed_inventory_is_unknown_not_mass_package_removal(self):
        before = {"schema_version": audit.SCHEMA, "identity": {"prefix": "/a", "executable": "/a/bin/python", "packages": [{"name": "numpy", "version": "1"}]}, "health": {"inventory": {"status": "passed"}}}
        failed = {"schema_version": audit.SCHEMA, "identity": {}, "health": {"inventory": {"status": "failed"}}}
        result = audit.compare_reports(before, failed)
        self.assertIsNone(result["packages_removed"])
        self.assertIsNone(result["same_environment_identity"])
        self.assertEqual(result["comparison_status"], "incomplete_inventory")
        self.assertIsNone(audit.compare_reports(failed, failed)["same_environment_identity"])

    def test_inferred_editable_checkout_and_shadowed_runtime_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "audit"
            source = Path(tmp) / "installed-source"
            source.mkdir()
            inventory = {"executable": sys.executable, "prefix": tmp, "packages": [],
                         "spateo_distribution": {"editable_checkout": str(source)}}
            calls = []
            def command(argv, timeout, cwd, env, log):
                calls.append(argv)
                Path(log).write_text("fixture")
                return {"returncode": 0, "log": Path(log).name}
            def check(name, *args):
                if name == "inventory":
                    return {"status": "passed", "level": "metadata", "details": inventory}
                return {"status": "api_only", "level": "api", "details": {"spateo_import_path": str(Path(tmp) / "shadow/spateo/__init__.py")}}
            args = SimpleNamespace(output=str(output), python=sys.executable, checkout=None,
                                   timeout=1, depth=1, max_entries=10, command="report", profile="standard", baseline=None)
            with patch.object(audit, "probe", side_effect=check), patch.object(audit, "run_command", side_effect=command):
                self.assertEqual(audit.audit(args), 1)
            verifier = next(argv for argv in calls if "--json" in argv)
            self.assertEqual(verifier[verifier.index("--checkout") + 1], str(source.resolve()))
            report = json.loads((output / "environment-report.json").read_text())
            self.assertEqual(report["health"]["import_identity"]["status"], "failed")
            self.assertEqual(report["overall"], "blocked")

    def test_unknown_schema_is_not_silently_compared(self):
        with self.assertRaises(ValueError):
            audit.compare_reports({}, {})


if __name__ == "__main__":
    unittest.main(verbosity=2)
