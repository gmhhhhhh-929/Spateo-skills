"""Run real Spateo APIs while forbidding any Dynamo import."""

import importlib.abc
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "skills/spateo-4d-pipeline/scripts"))
sys.path.insert(0, str(ROOT / "skills/spateo-data-io/scripts"))
sys.path.insert(
    0,
    str(
        ROOT
        / "skills/spateo-3d-pipeline/subskills/spateo-reconstruct-point-cloud/scripts"
    ),
)
sys.path.insert(
    0,
    str(ROOT / "skills/spateo-3d-pipeline/subskills/spateo-reconstruct-mesh/scripts"),
)


class NoDynamo(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] == "dynamo":
            raise AssertionError("Native skill attempted a Dynamo import: " + fullname)


@pytest.fixture(autouse=True, scope="session")
def native_only():
    assert not any(k.split(".")[0] == "dynamo" for k in sys.modules)
    guard = NoDynamo()
    sys.meta_path.insert(0, guard)
    yield
    sys.meta_path.remove(guard)
    assert not any(k.split(".")[0] == "dynamo" for k in sys.modules)
