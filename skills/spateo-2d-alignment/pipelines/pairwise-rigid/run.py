#!/usr/bin/env python3
"""Portable full-data SN-S CLI around the unchanged pairwise release runners."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def clean_environment(runner, overrides):
    """Ignore ambient settings recognized by the preserved scientific runner."""
    configuration_keys = set()
    for node in ast.walk(ast.parse(runner.read_text())):
        candidate = None
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "get" and node.args:
            if ast.dump(node.func.value) == ast.dump(ast.parse("os.environ", mode="eval").body):
                candidate = node.args[0]
        elif isinstance(node, ast.Subscript) and ast.dump(node.value) == ast.dump(ast.parse("os.environ", mode="eval").body):
            candidate = node.slice
        if isinstance(candidate, ast.Constant) and isinstance(candidate.value, str):
            configuration_keys.add(candidate.value)
    environment = {key: value for key, value in os.environ.items() if key not in configuration_keys}
    environment.update(overrides)
    return environment, sorted(configuration_keys), sorted(configuration_keys.intersection(os.environ))


def plan(args):
    root = Path(__file__).resolve().parent
    spatial = args.representation == "spatial-only"
    runner = root / "runners" / ("spateo_rigid_sns_alignment_spatial_only.py" if spatial else "spateo_rigid_sns_alignment.py")
    environment = {
        "DATASET_ROOT": str(args.slice_dir.resolve()), "OUTDIR": str(args.output_dir.resolve() / "native"),
        "PAIR_START_ORDER": "1", "N_SLICES": "0", "MAX_CELLS_PER_SLICE": "0", "BATCH_SIZE": str(args.batch_size),
        "STAGE1_MODE": "SN-S", "STAGE1_MAX_ITER": str(args.max_iter), "STAGE2_MODE": "none",
        "USE_MORPHO_ALIGN_REF": "0", "USE_SPATEO_INTERNAL_DOWNSAMPLING": "0", "N_SAMPLING_REF": "40000",
        "INITIAL_PRETRANSFORM_RECIPE": "", "INITIAL_COORDINATES_CSV": "", "EXCLUDE_CELLTYPES_REGEX": "",
        "REQUIRE_CUDA": "0" if args.device == "cpu" else "1", "SPATEO_DEVICE": args.device,
        "RNG_SEED": str(args.seed), "Z_DISPLAY_SPACING": str(args.z_display_spacing),
        "EXPRESSION_MODE": "spatial_only" if spatial else "normal", "DUMMY_REP_DIM": "30",
        "SIGMA2_INIT_SCALE": "2", "SIGMA2_END": "", "SAVE_FULL_ASSIGNMENT": "0", "RUN_LABEL": args.representation,
    }
    bootstrap = "import os,random,runpy,numpy as np,torch; s=int(os.environ['RNG_SEED']); random.seed(s); np.random.seed(s); torch.manual_seed(s); torch.cuda.manual_seed_all(s); runpy.run_path(__import__('sys').argv[1],run_name='__main__')"
    return runner, environment, [sys.executable, "-c", bootstrap, str(runner)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slice-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--representation", choices=["expression-pca", "annotation-onehot", "spatial-only"], required=True)
    parser.add_argument("--annotation-key", default="anno")
    parser.add_argument("--pca-provenance", type=Path)
    parser.add_argument("--max-iter", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=800)
    parser.add_argument("--seed", type=int, default=20260817)
    parser.add_argument("--device", default="0", help="CUDA index such as 0, or cpu")
    parser.add_argument("--z-display-spacing", type=float, default=40, help="Viewer spacing only; never biological z")
    parser.add_argument("--dry-run", action="store_true", help="Print the explicit execution plan without reading inputs or running Spateo")
    args = parser.parse_args()
    if args.max_iter < 1 or args.batch_size < 1:
        parser.error("--max-iter and --batch-size must be positive")
    runner, parameters, command = plan(args)
    environment, configuration_keys, cleared_keys = clean_environment(runner, parameters)
    if args.dry_run:
        print(json.dumps({"command": command, "environment": parameters, "preset": "full-data-sns",
                          "runner_environment_keys_isolated": configuration_keys,
                          "algorithm_executed": False}, indent=2))
        return
    from validate_inputs import sha, validate_inputs
    import numpy as np
    import pandas as pd
    summary, records = validate_inputs(args.slice_dir, args.representation, args.annotation_key, pca_provenance=args.pca_provenance)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    invocation = {"preset": "full-data-sns", "command": command, "parameters": parameters,
                  "runner_sha256": sha(runner), "wrapper_sha256": sha(__file__), "input_audit": summary,
                  "runner_environment_keys_isolated": configuration_keys,
                  "discarded_ambient_runner_keys": cleared_keys,
                  "runtime_environment": {key: environment.get(key) for key in
                     ["CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "PYTHONPATH"]}}
    (args.output_dir / "invocation.json").write_text(json.dumps(invocation, indent=2) + "\n")
    mapping = pd.concat([pd.DataFrame({"workflow_cell_id": [f"{record['slice_id']}:{i}" for i in range(len(record["cell_ids"]))],
                        "cell_id": record["cell_ids"], "slice_id": record["slice_id"],
                        "raw_x": record["xy"][:, 0], "raw_y": record["xy"][:, 1]}) for record in records], ignore_index=True)
    mapping.to_csv(args.output_dir / "input_cell_id_map.csv.gz", index=False, compression="gzip")
    with (args.output_dir / "run.log").open("w") as log:
        subprocess.run(command, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
    effective_provenance = args.output_dir / "native/spateo_two_stage_provenance.json"
    invocation["effective_runner_provenance"] = json.loads(effective_provenance.read_text())
    (args.output_dir / "invocation.json").write_text(json.dumps(invocation, indent=2) + "\n")
    native = pd.read_csv(args.output_dir / "native/spateo_two_stage_aligned_points.csv", dtype={"cell_id": str, "slice_id": str})
    native = native.rename(columns={"cell_id": "workflow_cell_id"})
    merged = mapping.merge(native, on="workflow_cell_id", how="outer", validate="one_to_one", suffixes=("", "_native"), indicator=True)
    if not (merged["_merge"] == "both").all() or not (merged["slice_id"] == merged["slice_id_native"]).all():
        raise ValueError("Native output identity/slice mapping failed")
    if not np.allclose(merged[["raw_x", "raw_y"]], merged[["raw_x_native", "raw_y_native"]], atol=1e-4, rtol=1e-7):
        raise ValueError("Native raw XY differs from audited input")
    final = merged[["cell_id", "slice_id", "stage1_rigid_x", "stage1_rigid_y"]].rename(columns={"stage1_rigid_x": "aligned_x", "stage1_rigid_y": "aligned_y"})
    if not np.isfinite(final[["aligned_x", "aligned_y"]]).all().all() or not final["cell_id"].is_unique:
        raise ValueError("Invalid final coordinates")
    output = args.output_dir / "aligned_coordinates.csv.gz"
    final.to_csv(output, index=False, compression="gzip")
    frozen = {"status": "frozen", "pipeline": "pairwise-rigid", "n_cells": len(final),
              "coordinate_sha256": sha(output), "input_audit": summary, "runner_sha256": sha(runner),
              "wrapper_sha256": sha(__file__), "parameters": parameters, "ground_truth_used": False,
              "effective_runner_provenance_sha256": sha(effective_provenance),
              "cell_ids_sha256": hashlib.sha256(json.dumps(sorted(final["cell_id"]), separators=(",", ":")).encode()).hexdigest()}
    (args.output_dir / "frozen_manifest.json").write_text(json.dumps(frozen, indent=2) + "\n")
    print(json.dumps({"status": "frozen", "n_cells": len(final), "output_dir": str(args.output_dir)}))


if __name__ == "__main__":
    main()
