#!/usr/bin/env python3
"""Stage entrypoint using the native v3 config and checkpoint engine."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
from run_4d_pipeline import main

if __name__ == "__main__":
    raise SystemExit(main(["--stop-after", "alignment", *sys.argv[1:]]))
