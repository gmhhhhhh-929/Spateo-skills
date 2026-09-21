#!/usr/bin/env python3
"""Contract-based spatial discovery/export and explicit single-reader conversion."""
import argparse
import contextlib
import hashlib
import inspect
import json
import os
from pathlib import Path
import sys
import uuid

SOURCE_COMMIT = "615644f88613bea8ceb2e2df1e2391d16de55ec1"
READERS = [
    "read_h5ad",
    "read_10x_h5",
    "read_10x_mtx",
    "read_visium",
    "read_visium_hd",
    "read_visium_hd_bin",
    "read_visium_hd_seg",
    "read_xenium",
    "read_atera",
    "read_merfish",
    "read_seqfish",
    "read_nanostring",
    "read_slideseq",
    "read_starmap_plus",
    "read_bgi",
    "read_bgi_agg",
    "read_stereoseq",
    "read_seqscope",
]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_equal(before, after, location="adata"):
    """Compare every persisted AnnData component, including images and layers."""
    import numpy as np
    import pandas as pd
    from scipy import sparse

    if (
        isinstance(before, dict)
        or hasattr(before, "keys")
        and not isinstance(before, pd.DataFrame)
    ):
        if set(before.keys()) != set(after.keys()):
            raise ValueError("Roundtrip changed keys: " + location)
        for k in before.keys():
            assert_equal(before[k], after[k], location + "/" + str(k))
    elif isinstance(before, pd.DataFrame):
        pd.testing.assert_frame_equal(
            before, after, check_dtype=False, check_categorical=False
        )
    elif sparse.issparse(before):
        if before.shape != after.shape or (before != after).nnz:
            raise ValueError("Roundtrip changed sparse values: " + location)
    elif before is None:
        if after is not None:
            raise ValueError("Roundtrip changed None: " + location)
    else:
        np.testing.assert_equal(np.asarray(before), np.asarray(after), err_msg=location)


def save_checked(adata, output, audit, provenance):
    import anndata as ad

    output = Path(output).resolve()
    sidecar = output.with_name(output.name + ".manifest.json")
    if output.exists() or sidecar.exists():
        raise FileExistsError(
            "Refusing to overwrite output or manifest: " + str(output)
        )
    if output.suffix != ".h5ad":
        raise ValueError("Output must use .h5ad")
    output.parent.mkdir(parents=True, exist_ok=True)
    temp = output.with_name("." + output.name + "." + uuid.uuid4().hex + ".h5ad")
    try:
        adata.write_h5ad(temp)
        loaded = ad.read_h5ad(temp)
        for attr in (
            "X",
            "obs",
            "var",
            "obsm",
            "varm",
            "obsp",
            "varp",
            "layers",
            "uns",
        ):
            assert_equal(getattr(adata, attr), getattr(loaded, attr), attr)
        if adata.raw is not None:
            if loaded.raw is None:
                raise ValueError("Roundtrip lost raw")
            for attr in ("X", "var", "varm"):
                assert_equal(
                    getattr(adata.raw, attr), getattr(loaded.raw, attr), "raw/" + attr
                )
        manifest = {
            **provenance,
            "validation": audit,
            "roundtrip": "pass",
            "output_sha256": sha256(temp),
        }
        os.link(temp, output)
        with sidecar.open("x") as stream:
            json.dump(manifest, stream, indent=2)
    finally:
        temp.unlink(missing_ok=True)
    return {"output": str(output), "manifest": str(sidecar)}


def select_feature_ids(data, column):
    if column is None:
        return data
    if column not in data.var:
        raise ValueError("Feature ID column is absent: " + column)
    ids = data.var[column]
    if (
        ids.isna().any()
        or not ids.astype(str).is_unique
        or any(not x.strip() for x in ids.astype(str))
    ):
        raise ValueError("Selected feature IDs must be unique and nonempty")
    if data.raw is not None:
        raise ValueError(
            "Feature-ID reassignment with raw requires an explicit matching raw-feature migration"
        )
    if "source_var_name" in data.var:
        raise ValueError(
            "source_var_name already exists; preserve the existing identity history"
        )
    data.var["source_var_name"] = data.var_names.to_numpy(dtype=str)
    data.var_names = ids.astype(str).to_numpy()
    return data


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("discover", "read"):
        p = sub.add_parser(name)
        p.add_argument("input", type=Path)
        p.add_argument("--technology")
        p.add_argument("--max-memory-bytes", type=int, default=1024**3)
        p.add_argument("--max-files", type=int, default=10000)
        p.add_argument("--max-depth", type=int, default=4)
        p.add_argument(
            "--load-images",
            action="store_true",
            help="Load optional raster assets; default is metadata-only",
        )
        p.add_argument("--stereoseq-bin-size", type=int)
        p.add_argument("--stereoseq-chemistry", choices=["V1", "V2"])
        if name == "read":
            p.add_argument("--output-dir", type=Path, required=True)
            p.add_argument(
                "--dataset-key",
                action="append",
                help="Explicit keys to export; full discovery status stays in report",
            )
            p.add_argument(
                "--matrix-semantics",
                choices=["unknown", "counts", "expression"],
                default="unknown",
            )
            p.add_argument("--coordinate-unit")
            p.add_argument(
                "--feature-id-column",
                help="Explicit stable ID column; preserve original var_names in source_var_name",
            )
    convert = sub.add_parser("convert")
    convert.add_argument("input", type=Path)
    convert.add_argument("output", type=Path)
    convert.add_argument("--reader", required=True, choices=READERS)
    convert.add_argument("--kwargs-json", type=Path)
    convert.add_argument(
        "--matrix-semantics",
        choices=["unknown", "counts", "expression"],
        default="unknown",
    )
    convert.add_argument("--coordinate-unit")
    convert.add_argument("--feature-id-column")
    convert.add_argument("--spatial-key", default="spatial")
    convert.add_argument("--require-spatial", action="store_true")
    convert.add_argument("--require-2d", action="store_true")
    args = parser.parse_args(argv)
    import spateo as st
    from validate_anndata import validate

    basis = {
        "skill_source_commit": SOURCE_COMMIT,
        "source_commit_verified_by_cli": False,
        "spateo_version": str(st.__version__),
        "source": str(args.input.resolve()),
        "input_file_sha256": sha256(args.input) if args.input.is_file() else None,
        "directory_content_hashes_complete": args.input.is_file(),
    }
    if args.command in ("discover", "read"):
        kwargs = {
            k: getattr(args, k)
            for k in (
                "technology",
                "max_memory_bytes",
                "max_files",
                "max_depth",
                "load_images",
                "stereoseq_bin_size",
                "stereoseq_chemistry",
            )
        }
        with contextlib.redirect_stdout(sys.stderr):
            result = st.io.read_spatial(
                args.input, load=args.command == "read", **kwargs
            )
        if args.command == "discover":
            print(json.dumps(result.report, indent=2))
            return (
                0
                if result.datasets
                and result.discovery.get("complete")
                and all(e.status == "deferred" for e in result.datasets.values())
                else 2
            )
        out = args.output_dir.resolve()
        out.mkdir(parents=True, exist_ok=False)
        selected = set(args.dataset_key or result.datasets)
        unknown = selected - set(result.datasets)
        exports, errors = {}, {}
        for key, entry in result.datasets.items():
            if key not in selected or entry.status != "ready":
                continue
            try:
                select_feature_ids(entry.adata, args.feature_id_column)
                audit = validate(
                    entry.adata,
                    matrix_semantics=args.matrix_semantics,
                    coordinate_unit=args.coordinate_unit,
                    require_spatial=True,
                )
                if audit["status"] != "pass":
                    raise ValueError(json.dumps(audit))
                filename = (
                    "dataset-" + hashlib.sha256(key.encode()).hexdigest()[:16] + ".h5ad"
                )
                exports[key] = save_checked(
                    entry.adata,
                    out / filename,
                    audit,
                    {
                        **basis,
                        "dataset_key": key,
                        "reader": "read_spatial",
                        "options": kwargs,
                        "feature_id_column": args.feature_id_column,
                    },
                )
            except Exception as exc:
                errors[key] = str(exc)
        report = {
            **result.report,
            "exports": exports,
            "export_errors": errors,
            "unknown_selected_keys": sorted(unknown),
        }
        report["export_status"] = (
            "complete"
            if selected and set(exports) == selected and not errors and not unknown
            else "incomplete"
        )
        (out / "read_report.json").write_text(json.dumps(report, indent=2))
        print(json.dumps(report, indent=2))
        # Partial scope is never reported as overall success, even with selected successful entries.
        return (
            0 if result.status == "ok" and report["export_status"] == "complete" else 2
        )
    kwargs = json.loads(args.kwargs_json.read_text()) if args.kwargs_json else {}
    if not isinstance(kwargs, dict) or "cache_file" in kwargs:
        raise ValueError(
            "Reader kwargs must be an object; CLI owns persistence, so omit cache_file"
        )
    # Preserve H5AD type/metadata exactly; Spateo direct H5AD reader sets UMI metadata.
    import anndata as ad

    reader = ad.read_h5ad if args.reader == "read_h5ad" else getattr(st.io, args.reader)
    inspect.signature(reader).bind(str(args.input), **kwargs)
    with contextlib.redirect_stdout(sys.stderr):
        data = reader(str(args.input), **kwargs)
    if not isinstance(data, ad.AnnData):
        raise TypeError("Explicit reader must return AnnData")
    select_feature_ids(data, args.feature_id_column)
    if args.reader == "read_bgi_agg":
        audit = {
            "status": "not_applicable",
            "reason": "AGG XY grid; not a cell-by-gene handoff",
        }
    else:
        audit = validate(
            data,
            matrix_semantics=args.matrix_semantics,
            coordinate_unit=args.coordinate_unit,
            spatial_key=args.spatial_key,
            require_spatial=args.require_spatial,
            require_2d=args.require_2d,
        )
        if audit["status"] != "pass":
            raise ValueError(json.dumps(audit))
    result = save_checked(
        data,
        args.output,
        audit,
        {
            **basis,
            "reader": args.reader,
            "options": kwargs,
            "feature_id_column": args.feature_id_column,
        },
    )
    print(json.dumps({"status": "saved", **result}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
