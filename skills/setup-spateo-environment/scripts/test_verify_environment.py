#!/usr/bin/env python3
"""Small verifier contract tests; no Spateo installation or scientific imports.

Requires packaging, as does the verifier's metadata path. Scientific-package and
distribution fixtures test diagnostic logic, not the actual Spateo environment.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

# Load packaging before simulating a supported interpreter. In particular,
# Python 3.9 must not enter packaging's real Python 3.10 typing branch.
import packaging.requirements
import packaging.specifiers
import packaging.version


SCRIPT = Path(__file__).with_name("verify_environment.py")
original_path = list(sys.path)
spec = importlib.util.spec_from_file_location("spateo_environment_verifier", SCRIPT)
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)
PATH_UNCHANGED = original_path == sys.path

VERSIONS = {
    "numpy": "1.26.4", "pandas": "2.2.3", "scipy": "1.13.1",
    "anndata": "0.10.9", "h5py": "3.11.0", "geopandas": "0.14.4",
    "shapely": "2.0.6", "pyarrow": "18.0.0", "PyMCubes": "0.1.6",
    "pymeshfix": "0.18.1", "pyacvd": "0.4.0",
}


class DistributionFixture:
    version = "1.1.2"
    requires = ["numpy>=1.23.5,<2"]

    def __init__(self, root, editable=True):
        self.root = root
        self.editable = editable

    def read_text(self, name):
        if name == "direct_url.json":
            return json.dumps({"url": self.root.as_uri(), "dir_info": {"editable": self.editable}})
        return None

    def locate_file(self, name):
        return self.root / name


class EnvironmentVerifierTests(unittest.TestCase):
    def setUp(self):
        self.stack = contextlib.ExitStack()
        self.addCleanup(self.stack.close)
        self.temporary = tempfile.TemporaryDirectory(prefix="spateo verifier test ")
        self.addCleanup(self.temporary.cleanup)
        self.checkout = Path(self.temporary.name).resolve() / "source with spaces"
        for name in verifier.CHECKOUT_FILES:
            destination = self.checkout / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text("# fixture\n")
        for target, value in (
            ("version", VERSIONS.__getitem__),
            ("distribution", lambda _: DistributionFixture(self.checkout)),
        ):
            self.stack.enter_context(patch.object(verifier.metadata, target, side_effect=value))
        self.stack.enter_context(patch.object(verifier.sys, "version_info", (3, 10, 9)))
        self.stack.enter_context(patch.object(verifier.platform, "python_implementation", return_value="CPython"))
        self.stack.enter_context(patch.object(verifier, "_prepare_caches"))

    def args(self, **overrides):
        values = dict(checkout=str(self.checkout), metadata_only=True, smoke_test=False, smoke_test_3d=False)
        values.update(overrides)
        return argparse.Namespace(**values)

    def fake_import(self, name):
        if name == "spateo":
            return SimpleNamespace(
                __file__=str(self.checkout / "spateo/__init__.py"), __version__="1.1.2",
                io=SimpleNamespace(read_atera=lambda: None, read_visium=lambda: None),
                pp=SimpleNamespace(preprocess_spatial=lambda: None),
            )
        return SimpleNamespace()

    def test_portable_import_has_no_path_override_and_help_needs_only_stdlib(self):
        self.assertTrue(PATH_UNCHANGED)
        process = subprocess.run([sys.executable, "-S", str(SCRIPT), "--help"], capture_output=True, text=True)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertIn("--checkout", process.stdout)

    def test_metadata_success_with_space_in_checkout_avoids_imports(self):
        with patch.object(verifier, "import_module", side_effect=AssertionError("unexpected import")):
            result = verifier.verify(self.args())
        self.assertTrue(result["ok"], result["errors"])
        self.assertEqual(result["mode"], "metadata-only")
        self.assertEqual(result["spateo_distribution"]["editable_checkout"], str(self.checkout))
        verifier._prepare_caches.assert_not_called()

    def test_wrong_editable_identity_is_a_failure(self):
        other = self.checkout.parent / "other source"
        with patch.object(verifier.metadata, "distribution", return_value=DistributionFixture(other)):
            result = verifier.verify(self.args())
        self.assertFalse(result["ok"])
        self.assertTrue(any("metadata does not point" in error for error in result["errors"]))

    def test_missing_source_installation_file_is_a_failure(self):
        (self.checkout / "docs/requirements.txt").unlink()
        result = verifier.verify(self.args())
        self.assertFalse(result["ok"])
        self.assertEqual(result["checkout"]["missing_files"], ["docs/requirements.txt"])

    def test_importable_source_without_distribution_is_not_a_pass(self):
        with patch.object(verifier.metadata, "distribution", side_effect=verifier.metadata.PackageNotFoundError), \
             patch.object(verifier, "import_module", side_effect=self.fake_import):
            result = verifier.verify(self.args(metadata_only=False))
        self.assertFalse(result["ok"])
        self.assertTrue(any("distribution metadata is missing" in error for error in result["errors"]))

    def test_shadowed_module_fails_and_smoke_is_skipped(self):
        module = self.fake_import("spateo")
        module.__file__ = str(self.checkout.parent / "wrong/spateo/__init__.py")
        with patch.object(verifier, "import_module", return_value=module), \
             patch.object(verifier, "_run_smoke_test") as smoke:
            result = verifier.verify(self.args(metadata_only=False, smoke_test=True))
        self.assertFalse(result["ok"])
        self.assertEqual(result["smoke_test"]["status"], "skipped")
        smoke.assert_not_called()

    def test_json_output_separates_python_library_chatter(self):
        def noisy_verify(_):
            print("library diagnostic")
            return {"ok": True, "errors": []}

        stdout, stderr = io.StringIO(), io.StringIO()
        with patch.object(verifier, "verify", side_effect=noisy_verify), \
             contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            code = verifier.main(["--json"])
        self.assertEqual(code, 0)
        self.assertTrue(json.loads(stdout.getvalue())["ok"])
        self.assertIn("library diagnostic", stderr.getvalue())

    def test_unwritable_cache_is_a_structured_failure_and_skips_smoke(self):
        with patch.object(verifier, "_prepare_caches", side_effect=PermissionError("fixture cache is read-only")):
            result = verifier.verify(self.args(metadata_only=False, smoke_test=True))
        self.assertFalse(result["ok"])
        self.assertEqual(result["smoke_test"]["status"], "skipped")
        self.assertTrue(any("temporary scientific-package caches" in error for error in result["errors"]))

    def test_bare_environment_has_structured_failure(self):
        process = subprocess.run(
            [sys.executable, "-S", str(SCRIPT), "--metadata-only", "--json"], capture_output=True, text=True,
        )
        self.assertEqual(process.returncode, 1, process.stderr)
        result = json.loads(process.stdout)
        self.assertFalse(result["ok"])
        self.assertTrue(any("packaging" in error for error in result["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
