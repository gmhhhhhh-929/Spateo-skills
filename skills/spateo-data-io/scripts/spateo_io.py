#!/usr/bin/env python3
"""Thin source-backed detect/convert CLI. No registration or silent correction."""
import argparse
import contextlib
import hashlib
import inspect
import json
import os
from pathlib import Path
import platform
import sys
import uuid


READERS = ["read_h5ad", "read_10x_h5", "read_10x_mtx", "read_auto_spatial",
           "read_visium", "read_visium_hd", "read_visium_hd_bin", "read_visium_hd_seg",
           "read_xenium", "read_atera", "read_merfish", "read_seqfish", "read_nanostring",
           "read_slideseq", "read_starmap_plus", "read_bgi", "read_bgi_agg", "read_seqscope"]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    detect = sub.add_parser("detect")
    detect.add_argument("input", type=Path)
    detect.add_argument("--technology")
    convert = sub.add_parser("convert")
    convert.add_argument("input", type=Path)
    convert.add_argument("output", type=Path)
    convert.add_argument("--reader", required=True, choices=READERS)
    convert.add_argument("--kwargs-json", type=Path, help="JSON object passed as reader keyword arguments")
    convert.add_argument("--matrix-semantics", choices=["unknown", "counts", "expression"], default="unknown")
    convert.add_argument("--coordinate-unit")
    convert.add_argument("--spatial-key", default="spatial")
    convert.add_argument("--require-spatial", action="store_true")
    convert.add_argument("--require-2d", action="store_true", help="Reject XYZ for an explicit 2D handoff")
    args = parser.parse_args()
    # Import only after argument handling; --help does not require scientific packages.
    import spateo
    import spateo.io as io
    from validate_anndata import ids_sha256, validate
    if args.command == "detect":
        matches = io.detect_spatial_technologies(str(args.input), technology=args.technology)
        print(json.dumps([{"technology": m.technology, "reader": m.reader, "path": str(m.path),
                           "kwargs": dict(m.kwargs), "confidence": m.confidence,
                           "evidence": list(m.evidence)} for m in matches], ensure_ascii=False, indent=2))
        return 0 if matches else 2
    output = args.output.resolve()
    sidecar = output.with_name(output.name + ".manifest.json")
    if output.exists() or sidecar.exists():
        raise FileExistsError("Refusing to overwrite output or manifest: " + str(output))
    if output.suffix != ".h5ad":
        raise ValueError("Output must have .h5ad suffix")
    kwargs = json.loads(args.kwargs_json.read_text()) if args.kwargs_json else {}
    if not isinstance(kwargs, dict):
        raise ValueError("--kwargs-json must contain an object")
    if "cache_file" in kwargs or kwargs.get("return_match", False):
        raise ValueError("CLI owns output persistence; omit cache_file and return_match")
    reader = getattr(io, args.reader)
    inspect.signature(reader).bind(str(args.input), **kwargs)
    with contextlib.redirect_stdout(sys.stderr):
        adata = reader(str(args.input), **kwargs)
    import anndata as ad
    if not isinstance(adata, ad.AnnData):
        raise TypeError("Selected reader did not return AnnData")
    if args.reader == "read_bgi_agg":
        # AGG is valid IO, but deliberately not accepted by the gene/spatial validator.
        audit = {"status": "not_applicable", "reason": "AGG XY grid; not a gene-matrix alignment input"}
    else:
        audit = validate(adata, matrix_semantics=args.matrix_semantics, spatial_key=args.spatial_key,
                         coordinate_unit=args.coordinate_unit, require_spatial=args.require_spatial,
                         require_2d=args.require_2d)
        if audit["status"] != "pass":
            raise ValueError(json.dumps(audit, ensure_ascii=False))
    manifest = {"reader": args.reader, "kwargs": kwargs, "source": str(args.input.resolve()),
                "input_file_sha256": sha256(args.input) if args.input.is_file() else None,
                "input_inventory": io.spatial_file_manifest(args.input),
                "n_obs": adata.n_obs, "n_vars": adata.n_vars,
                "cell_ids_sha256": ids_sha256(adata.obs_names),
                "spateo_version": str(spateo.__version__), "python": platform.python_version(),
                "reader_module": inspect.getsourcefile(inspect.unwrap(reader)),
                "source_commit_basis": "d6aa68addc475dd0b56f69cebe7823b1f79933a9",
                "source_commit_verified_by_cli": False, "validation": audit}
    if manifest["reader_module"]:
        manifest["reader_module_sha256"] = sha256(manifest["reader_module"])
    # Inventory paths are not hashes. Hash only the explicitly supplied file here;
    # large directory bundles require a separately chosen core-input hash manifest.
    manifest["directory_core_input_hashes_complete"] = args.input.is_file()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name("." + output.name + "." + uuid.uuid4().hex + ".tmp.h5ad")
    try:
        adata.write_h5ad(temporary)
        roundtrip = ad.read_h5ad(temporary)
        import numpy as np
        from scipy import sparse
        if list(roundtrip.obs_names) != list(adata.obs_names) or list(roundtrip.var_names) != list(adata.var_names):
            raise ValueError("H5AD roundtrip changed identity or order")
        difference = roundtrip.X - adata.X
        values = difference.data if sparse.issparse(difference) else np.asarray(difference)
        if np.any(values != 0):
            raise ValueError("H5AD roundtrip changed X")
        for key in adata.obsm:
            before, after = adata.obsm[key], roundtrip.obsm[key]
            if sparse.issparse(before):
                if (before != after).nnz:
                    raise ValueError("H5AD roundtrip changed obsm: " + key)
            elif not np.array_equal(np.asarray(before), np.asarray(after)):
                raise ValueError("H5AD roundtrip changed obsm: " + key)
        manifest["coordinate_and_matrix_roundtrip"] = "pass"
        manifest["output_sha256"] = sha256(temporary)
        # Hard link publishes atomically without replacing any existing output.
        os.link(temporary, output)
        with sidecar.open("x") as stream:
            json.dump(manifest, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    finally:
        temporary.unlink(missing_ok=True)
    print(json.dumps({"status": "saved", "output": str(output), "manifest": str(sidecar)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
