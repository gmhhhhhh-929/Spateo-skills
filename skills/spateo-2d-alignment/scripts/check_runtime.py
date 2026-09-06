#!/usr/bin/env python3
"""Check imports and alignment API keyword compatibility; no registration is run."""
import argparse
import importlib
import importlib.metadata
import inspect
import json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cpu", help="cpu or the intended CUDA device index")
    args = parser.parse_args()
    errors, versions = [], {}
    for distribution, module in [("numpy", "numpy"), ("pandas", "pandas"), ("scipy", "scipy"),
                                 ("anndata", "anndata"), ("h5py", "h5py"), ("matplotlib", "matplotlib"),
                                 ("torch", "torch"), ("POT", "ot"), ("spateo-release", "spateo")]:
        try:
            importlib.import_module(module)
            versions[distribution] = importlib.metadata.version(distribution)
        except Exception as exc:
            errors.append(f"{module}: {type(exc).__name__}: {exc}")
    try:
        from spateo.alignment.morpho_alignment import morpho_align, morpho_align_ref
        from spateo.alignment.methods import Morpho_pairwise
        required = {
            morpho_align: {"models", "rep_layer", "rep_field", "spatial_key", "key_added", "iter_key_added",
                          "vecfld_key_added", "mode", "dissimilarity", "max_iter", "dtype", "device", "verbose"},
            morpho_align_ref: {"models", "models_ref", "n_sampling", "sampling_method"},
            Morpho_pairwise.__init__: {"nn_init", "SVI_mode", "batch_size", "pre_compute_dist", "sparse_calculation_mode",
                    "sparse_top_k", "use_chunk", "chunk_capacity", "K", "beta", "lambdaVF", "partial_robust_level",
                    "sigma2_init_scale", "sigma2_end", "use_hvg", "init_transform", "init_layer", "init_field",
                    "allow_flip", "nonrigid_start_iter"},
        }
        for function, names in required.items():
            missing = names - set(inspect.signature(function).parameters)
            if missing:
                errors.append(f"{function.__qualname__} missing explicit parameters: {sorted(missing)}")
    except Exception as exc:
        errors.append(f"alignment API: {type(exc).__name__}: {exc}")
    if args.device != "cpu":
        try:
            import torch
            if not torch.cuda.is_available() or not 0 <= int(args.device) < torch.cuda.device_count():
                errors.append(f"Requested CUDA device is unavailable: {args.device}")
        except Exception as exc:
            errors.append(f"CUDA: {type(exc).__name__}: {exc}")
    print(json.dumps({"status": "failed" if errors else "passed", "versions": versions,
                      "errors": errors, "alignment_executed": False,
                      "scope": "Import and signature checks only; not a complete pipeline/GPU accuracy test"}, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
