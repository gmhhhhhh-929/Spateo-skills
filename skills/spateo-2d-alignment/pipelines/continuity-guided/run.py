#!/usr/bin/env python3
"""Validate annotation inputs, then run the preserved continuity-guided engine."""
from pathlib import Path
import subprocess
import sys


def main():
    import engine
    args = engine.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"Use a new output directory: {args.output_dir}")
    from validate_inputs import validate_inputs
    validate_inputs(args.slice_dir, "annotation-onehot", args.annotation_key,
                    args.spatial_key, args.representation_key, "continuity")
    subprocess.run([sys.executable, str(Path(__file__).resolve().with_name("engine.py")), *sys.argv[1:]], check=True)


if __name__ == "__main__":
    main()
