#!/usr/bin/env python3
"""Probe isolation, evidence and failure contracts; standard library only.

These fixture tests do not replace the separate real-interpreter smoke probes.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


SCRIPT = Path(__file__).with_name("environment_probes.py")
before_path = list(sys.path)
spec = importlib.util.spec_from_file_location("environment_probes", SCRIPT)
probes = importlib.util.module_from_spec(spec)
spec.loader.exec_module(probes)
PATH_UNCHANGED = before_path == sys.path


class DistributionFixture:
    def __init__(self, url, editable=True, name="spateo-release", version="1.1.2"):
        self.direct_url = json.dumps({"url": url, "dir_info": {"editable": editable}})
        self.metadata = {"Name": name}
        self.version = version

    def read_text(self, name):
        return self.direct_url if name == "direct_url.json" else None

    def locate_file(self, name):
        return Path("/fixture/site-packages") / name


class EnvironmentProbesTests(unittest.TestCase):
    def test_help_and_inventory_work_without_site_packages(self):
        self.assertTrue(PATH_UNCHANGED)
        with tempfile.TemporaryDirectory(prefix="audit tests ") as tmp:
            help_run = subprocess.run([sys.executable, "-S", str(SCRIPT), "--help"], capture_output=True, text=True)
            self.assertEqual(help_run.returncode, 0, help_run.stderr)
            result_file = Path(tmp) / "metadata.json"
            run = subprocess.run([sys.executable, "-S", str(SCRIPT), "--probe", "inventory", "--result-file", str(result_file)], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            result = json.loads(result_file.read_text())
            self.assertEqual((result["status"], result["level"]), ("passed", "metadata"))
            self.assertEqual(result["details"]["packages"], [])
            self.assertIsNone(result["details"]["spateo_distribution"])
            self.assertEqual(result["details"]["executable"], sys.executable)

    def test_inventory_avoids_scientific_imports_and_sorts_packages(self):
        fixtures = [DistributionFixture("file:///tmp/source%20tree"), DistributionFixture("https://unused.invalid", name="AnnData", version="0.10")]
        with patch.object(probes.metadata, "distributions", return_value=fixtures), \
             patch.object(probes.metadata, "distribution", return_value=fixtures[0]), \
             patch.object(probes.importlib, "import_module", side_effect=AssertionError("unexpected scientific import")):
            result = probes.run_probe("inventory")
        self.assertEqual(result["status"], "passed")
        self.assertEqual([item["name"] for item in result["details"]["packages"]], ["AnnData", "spateo-release"])
        self.assertEqual(result["details"]["spateo_distribution"]["editable_checkout"], str(Path("/tmp/source tree").resolve()))

    def test_editable_urls_never_export_credentials_or_remote_locations(self):
        for url in (
            "https://user:supersecret@example.invalid/repo",
            "git+ssh://user:supersecret@example.invalid/repo",
            "file://user:supersecret@example.invalid/repo",
            "file://server.example.invalid/repo",
            "file:///tmp/source?token=supersecret",
            "file:///tmp/source#supersecret",
            "file:relative/source",
            "file:///tmp/bad%00path",
        ):
            with self.subTest(url=url):
                fixture = DistributionFixture(url)
                self.assertIsNone(probes._editable_checkout(fixture))
                with patch.object(probes.metadata, "distributions", return_value=[fixture]), \
                     patch.object(probes.metadata, "distribution", return_value=fixture):
                    serialized = json.dumps(probes.run_probe("inventory"))
                self.assertNotIn("supersecret", serialized)
                self.assertNotIn("example.invalid", serialized)

    def test_malformed_direct_url_is_ignored_without_echo(self):
        fixture = DistributionFixture("file:///tmp/source")
        for value in ("bad secret JSON", "[]", "null", '{"dir_info": true}', '{"dir_info":{"editable":true},"url":12}'):
            fixture.direct_url = value
            self.assertIsNone(probes._editable_checkout(fixture))
        self.assertIsNone(probes._editable_checkout(DistributionFixture("file:///tmp/source", editable=False)))

    def test_missing_dependency_yields_json_failure_in_bare_interpreter(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "failure.json"
            process = subprocess.run([sys.executable, "-S", str(SCRIPT), "--probe", "io", "--result-file", str(path)], capture_output=True, text=True)
            self.assertEqual(process.returncode, 1, process.stderr)
            result = json.loads(path.read_text())
            self.assertEqual(result["status"], "failed")
            self.assertIn("ModuleNotFoundError", result["error"])
            self.assertIsNone(result["details"]["spateo_import_path"])

    def test_explicit_signature_required_not_only_kwargs(self):
        root = SimpleNamespace(api=lambda **kwargs: None)
        with self.assertRaisesRegex(RuntimeError, "missing explicit parameters.*spatial_key"):
            probes._check_apis(root, {"api": ["spatial_key"]})
        root.api = lambda spatial_key=None, **kwargs: None
        self.assertEqual(probes._check_apis(root, {"api": ["spatial_key"]})[0]["required_parameters"], ["spatial_key"])
        root.api = None
        with self.assertRaises(TypeError):
            probes._check_apis(root, {"api": []})

    def test_api_only_is_preserved_and_error_does_not_claim_readiness(self):
        with patch.dict(probes.PROBES, {"alignment": lambda: probes._result("api_only", "api", {"scope": "Only signatures"})}):
            result = probes.run_probe("alignment")
        self.assertEqual((result["status"], result["level"]), ("api_only", "api"))
        def broken():
            raise RuntimeError("fetch https://u:password@example.invalid/test failed")
        with patch.dict(probes.PROBES, {"alignment": broken}):
            result = probes.run_probe("alignment")
        self.assertEqual(result["status"], "failed")
        self.assertNotIn("password", json.dumps(result))
        self.assertIn("no readiness claim", result["details"]["scope"])

    def test_runtime_environment_and_bytecode_flag_are_restored_after_failure(self):
        original_bytecode = sys.dont_write_bytecode
        observed = {}
        def fail():
            observed["cache"] = os.environ["MPLCONFIGDIR"]
            self.assertTrue(sys.dont_write_bytecode)
            self.assertEqual(os.environ["OMP_NUM_THREADS"], "1")
            raise RuntimeError("fixture")
        with patch.dict(os.environ, {"MPLCONFIGDIR": "/do/not/write/here", "OMP_NUM_THREADS": "7"}), \
             patch.dict(probes.PROBES, {"mesh": fail}):
            self.assertEqual(probes.run_probe("mesh")["status"], "failed")
            self.assertEqual(os.environ["MPLCONFIGDIR"], "/do/not/write/here")
            self.assertEqual(os.environ["OMP_NUM_THREADS"], "7")
        self.assertFalse(Path(observed["cache"]).parent.exists())
        self.assertEqual(sys.dont_write_bytecode, original_bytecode)

    def test_result_file_cannot_overwrite_existing_user_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "existing.json"
            destination.write_text("preserve me")
            with patch.object(probes, "run_probe") as action, contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as outcome:
                    probes.main(["--probe", "inventory", "--result-file", str(destination)])
            self.assertEqual(outcome.exception.code, 2)
            action.assert_not_called()
            self.assertEqual(destination.read_text(), "preserve me")

    def test_json_file_is_separate_from_native_stdout(self):
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "result.json"
            source = """
import importlib.util, os, sys
spec = importlib.util.spec_from_file_location('probes', sys.argv[1])
probes = importlib.util.module_from_spec(spec)
spec.loader.exec_module(probes)
def chatter():
    print('python import chatter')
    os.write(1, b'native stdout chatter\\n')
    return probes._result('passed', 'smoke', {'scope': 'fixture'})
probes.PROBES['io'] = chatter
raise SystemExit(probes.main(['--probe', 'io', '--result-file', sys.argv[2]]))
"""
            process = subprocess.run([sys.executable, "-S", "-c", source, str(SCRIPT), str(destination)], capture_output=True, text=True)
            self.assertEqual(process.returncode, 0, process.stderr)
            self.assertIn("native stdout chatter", process.stdout)
            self.assertIn("python import chatter", process.stderr)
            self.assertEqual(json.loads(destination.read_text())["status"], "passed")


if __name__ == "__main__":
    unittest.main()
