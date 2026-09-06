#!/usr/bin/env python3
"""Read-only AnnData structure/2D representation audit; no provenance inference."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse


def ids_sha256(values):
    payload = json.dumps(sorted(map(str, values)), ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def validate(adata, matrix_semantics="unknown", matrix_key="X", spatial_key="spatial",
             coordinate_unit=None, alignment_mode="none", representation_key="X_pca",
             annotation_key=None, categories=None, require_spatial=False, require_2d=False):
    errors, warnings = [], []
    report = {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars),
              "cell_ids_sha256": ids_sha256(adata.obs_names), "matrix_semantics": matrix_semantics,
              "matrix_key": matrix_key, "spatial_key": spatial_key,
              "coordinate_unit_declared": coordinate_unit, "alignment_mode": alignment_mode,
              "source_provenance_verified": False}
    if adata.n_obs == 0 or adata.n_vars == 0:
        errors.append("Empty observations or variables")
    for name, index in (("obs_names", adata.obs_names), ("var_names", adata.var_names)):
        if not index.is_unique or index.hasnans or any(not str(x).strip() for x in index):
            errors.append(name + " must be unique and nonempty; preserve an explicit identity mapping")
    if adata.uns.get("__type") == "AGG":
        errors.append("AGG XY grid is not an observation-by-gene matrix")
    if matrix_key != "X" and matrix_key not in adata.layers:
        errors.append("Missing matrix layer: " + matrix_key)
    else:
        matrix = adata.X if matrix_key == "X" else adata.layers[matrix_key]
        if matrix is None:
            errors.append("Selected matrix is absent")
        else:
            nonfinite = negative = noninteger = 0
            for first in range(0, adata.n_obs, 4096):
                block = matrix[first:first + 4096]
                values = block.data if sparse.issparse(block) else np.asarray(block).ravel()
                if values.dtype.kind not in "biuf":
                    errors.append("Matrix is not numeric")
                    break
                nonfinite += int(np.count_nonzero(~np.isfinite(values)))
                if matrix_semantics == "counts":
                    negative += int(np.count_nonzero(values < 0))
                    noninteger += int(np.count_nonzero(np.isfinite(values) & (values != np.floor(values))))
            report["matrix_nonfinite"] = nonfinite
            if nonfinite:
                errors.append("Matrix contains nonfinite values")
            if matrix_semantics == "counts" and (negative or noninteger):
                errors.append("Declared counts contain negative or fractional values")
    if matrix_semantics == "unknown":
        warnings.append("Matrix semantic origin is unknown; numeric shape cannot establish gene counts")
    need_xy = require_spatial or require_2d or alignment_mode != "none"
    if spatial_key not in adata.obsm:
        if need_xy:
            errors.append("Missing spatial key: " + spatial_key)
        else:
            warnings.append("No spatial key; expression-only AnnData has no inferred coordinates")
    else:
        xy = np.asarray(adata.obsm[spatial_key])
        report["spatial_shape"] = list(xy.shape)
        allowed_dims = (2,) if (require_2d or alignment_mode != "none") else (2, 3)
        if xy.ndim != 2 or xy.shape[0] != adata.n_obs or xy.shape[1] not in allowed_dims or xy.dtype.kind not in "biuf":
            errors.append("Spatial contract requires numeric (n_obs, %s) coordinates; retain native XYZ before deriving XY" % ("2" if allowed_dims == (2,) else "2 or 3"))
        elif not np.isfinite(xy).all():
            errors.append("Spatial coordinates contain nonfinite values")
        else:
            report["spatial_extent"] = np.ptp(xy, axis=0).tolist()
            if adata.n_obs > 1 and np.all(xy == xy[:1]):
                warnings.append("All coordinates coincide; inspect spatial origin and observation unit")
        if not coordinate_unit:
            warnings.append("Coordinate unit/frame is not declared; inspect platform metadata before alignment")
    if alignment_mode != "none":
        if not coordinate_unit:
            errors.append("Alignment contract requires an explicit coordinate unit")
        if representation_key not in adata.obsm:
            errors.append("Missing representation key: " + representation_key)
        else:
            rep = adata.obsm[representation_key]
            rep = rep.toarray() if sparse.issparse(rep) else np.asarray(rep)
            if rep.ndim != 2 or rep.shape[0] != adata.n_obs or rep.shape[1] < 1:
                errors.append("Invalid representation shape")
            elif rep.dtype.kind not in "biuf" or not np.isfinite(rep).all():
                errors.append("Representation must be numeric and finite")
            elif adata.n_obs:
                constant = bool(np.all(rep == rep[:1]))
                onehot = bool(np.all((rep == 0) | (rep == 1)) and np.all(rep.sum(axis=1) == 1))
                report.update({"representation_shape": list(rep.shape), "constant_rows": constant,
                               "strict_onehot": onehot})
                if alignment_mode == "expression-pca":
                    if constant or onehot:
                        errors.append("Expression PCA cannot be accepted as constant or strict annotation one-hot")
                    if np.any(np.all(rep == 0, axis=1)):
                        errors.append("Zero-norm PCA rows do not satisfy the Spateo 2D cosine-input contract")
                    warnings.append("PCA origin/shared basis not numerically provable; verify expression and fit manifest")
                elif alignment_mode == "spatial-only":
                    if not constant or not np.any(rep[0] != 0):
                        errors.append("Spatial-only requires identical finite nonzero vectors")
                elif alignment_mode == "annotation-onehot":
                    if not onehot:
                        errors.append("Annotation encoding must be exact one-hot")
                    if not annotation_key or annotation_key not in adata.obs:
                        errors.append("Annotation key is required and must exist")
                    if categories is None or not isinstance(categories, list) or not all(isinstance(c, str) for c in categories):
                        errors.append("Supply the shared ordered string category list via --categories-json")
                    elif len(categories) != rep.shape[1] or len(set(categories)) != len(categories):
                        errors.append("Category list must uniquely name every representation column")
                    elif annotation_key in adata.obs and onehot:
                        labels = adata.obs[annotation_key]
                        if labels.isna().any() or any(not str(v).strip() for v in labels):
                            errors.append("Missing/blank annotations cannot be silently encoded")
                        elif not np.array_equal(np.asarray(categories)[rep.argmax(axis=1)], labels.astype(str).to_numpy()):
                            errors.append("One-hot columns do not decode to the corresponding observation annotations")
    report.update({"status": "fail" if errors else "pass", "errors": errors, "warnings": warnings})
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h5ad", type=Path)
    parser.add_argument("--matrix-semantics", choices=["unknown", "counts", "expression"], default="unknown")
    parser.add_argument("--matrix-key", default="X", help="X or an existing layer name")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--require-spatial", action="store_true")
    parser.add_argument("--require-2d", action="store_true", help="Require XY; generic IO otherwise preserves XY or XYZ")
    parser.add_argument("--coordinate-unit")
    parser.add_argument("--alignment-mode", choices=["none", "expression-pca", "annotation-onehot", "spatial-only"], default="none")
    parser.add_argument("--representation-key", default="X_pca")
    parser.add_argument("--annotation-key")
    parser.add_argument("--categories-json", type=Path)
    args = parser.parse_args()
    import anndata as ad
    kwargs = vars(args).copy()
    path = kwargs.pop("h5ad")
    category_path = kwargs.pop("categories_json")
    kwargs["categories"] = json.loads(category_path.read_text()) if category_path else None
    # backed='r' keeps counts on disk; coordinates/metadata are still loaded by AnnData.
    adata = ad.read_h5ad(path, backed="r")
    try:
        report = validate(adata, **kwargs)
    finally:
        adata.file.close()
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
