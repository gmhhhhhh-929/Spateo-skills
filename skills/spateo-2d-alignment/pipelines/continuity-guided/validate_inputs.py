#!/usr/bin/env python3
"""Inspect blind slice identities and feature contracts without opening source data."""
import argparse
import hashlib
import json
import re
from pathlib import Path

import h5py
import numpy as np


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


def text_array(values):
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in values], dtype=str)


def column(node):
    if isinstance(node, h5py.Group) and "categories" in node:
        codes = node["codes"][()]
        if np.any(codes < 0):
            raise ValueError(f"Missing annotation: {node.name}")
        return text_array(node["categories"][()])[codes]
    return text_array(node[()])


def feature_sha(features):
    return hashlib.sha256(np.ascontiguousarray(features, dtype='<f4').tobytes(order='C')).hexdigest()


def cell_ids_sha(ids):
    return hashlib.sha256(json.dumps(sorted(map(str, ids)), ensure_ascii=False, separators=(',', ':')).encode()).hexdigest()


def validate_pca_manifest(path, records, dimension):
    path = Path(path)
    if not path.is_file():
        raise ValueError('Expression PCA requires expression_pca_manifest.json from a joint PCA preparation, or --pca-provenance')
    manifest = json.loads(path.read_text())
    if (manifest.get('schema_version') != 'expression-pca-v1'
            or manifest.get('feature_mode', manifest.get('mode')) != 'expression-pca'
            or manifest.get('shared_basis') is not True
            or manifest.get('fitted_on_spatial') is not False
            or manifest.get('fitted_on_annotation') is not False):
        raise ValueError('PCA provenance must declare one shared expression-only basis without coordinate/annotation fitting')
    basis = Path(manifest.get('basis_file', 'basis.npz'))
    if not basis.is_absolute():
        basis = path.parent / basis
    basis_hash = manifest.get('basis_sha256')
    if not basis.is_file() or sha(basis) != basis_hash:
        raise ValueError('PCA basis file hash mismatch or missing basis artifact')
    effective = manifest.get('n_pcs', manifest.get('effective_n_pcs'))
    if effective is not None and int(effective) != dimension:
        raise ValueError('PCA dimension differs from shared-basis provenance')
    entries = manifest.get('slices', [])
    by_name = {Path(e['path']).name: e for e in entries}
    if len(entries) != len(by_name) or set(by_name) != {r['path'].name for r in records}:
        raise ValueError('PCA provenance must cover exactly these slices once')
    for record in records:
        entry = by_name[record['path'].name]
        if entry.get('basis_sha256') != basis_hash:
            raise ValueError('Per-slice PCA bases differ from the declared joint basis')
        if entry.get('feature_key') != record['representation_key'] or entry.get('spatial_key') != record['spatial_key']:
            raise ValueError('PCA provenance feature/spatial keys differ')
        if entry.get('features_sha256') != record['features_sha256']:
            raise ValueError('PCA feature hash mismatch')
        if entry.get('output_h5ad_sha256') != record['sha256']:
            raise ValueError('PCA slice H5AD hash mismatch')
        if entry.get('cell_ids_sha256') != cell_ids_sha(record['cell_ids']) or int(entry.get('n_cells', -1)) != len(record['cell_ids']):
            raise ValueError('PCA provenance cell identity/count mismatch')
    return {'status': 'verified', 'path': str(path.resolve()), 'sha256': sha(path),
            'basis_sha256': basis_hash, 'slices_verified': len(records),
            'feature_hash_encoding': 'little-endian float32 C-order raw bytes'}


def validate_inputs(slice_dir, representation, annotation_key="anno", spatial_key="spatial",
                    representation_key="X_pca", filename_style="pairwise", annotation_qc=None,
                    pca_provenance=None):
    if representation not in {'annotation-onehot', 'expression-pca', 'spatial-only'}:
        raise ValueError('Unknown representation mode')
    annotation_qc = annotation_qc or ('provided' if representation == 'annotation-onehot' else 'off')
    if annotation_qc not in {'off', 'provided'} or (representation == 'annotation-onehot' and annotation_qc != 'provided'):
        raise ValueError('Annotation one-hot requires annotation-qc provided')
    paths = list(Path(slice_dir).glob("*.h5ad"))
    if filename_style == "pairwise" and set(Path(slice_dir).rglob("*.h5ad")) != set(paths):
        raise ValueError("Pairwise input must be a flat directory: nested H5ADs would be read recursively by the preserved runner")
    if len(paths) < 2:
        raise ValueError("At least two H5AD slices are required")
    ordered = []
    stages = set()
    for path in paths:
        pattern = r"CS(?P<stage>\d+)_SL(?P<sl>\d+)_Y\w+\.Spatial\.h5ad$" if filename_style == "pairwise" else r"SL(?P<sl>\d+)"
        match = re.search(pattern, path.name, re.IGNORECASE)
        if not match:
            raise ValueError(f"Filename does not meet {filename_style} ordering contract: {path.name}")
        if filename_style == "pairwise":
            stages.add(match.group("stage"))
        ordered.append((int(match.group("sl")), path))
    if len(stages) > 1 or len({number for number, _ in ordered}) != len(paths):
        raise ValueError("Use one specimen/stage with distinct slice numbers")
    records, ids_seen, label_vectors = [], set(), {}
    dimensions = set()
    forbidden_obs = {"align_x", "align_y", "aligned_x", "aligned_y", "truth_x", "truth_y", "manual_x", "manual_y"}
    for number, path in sorted(ordered):
        with h5py.File(path, "r") as handle:
            # Inspect key names before loading any feature/coordinate array.
            if set(handle["obsm"].keys()) != {spatial_key, representation_key}:
                raise ValueError(f"Use a sanitized blind copy with only requested spatial and representation obsm: {path}")
            if forbidden_obs.intersection(handle["obs"].keys()):
                raise ValueError(f"Reference/aligned coordinate columns present in obs: {path}")
            index_key = handle["obs"].attrs.get("_index", "_index")
            if isinstance(index_key, bytes):
                index_key = index_key.decode()
            ids = column(handle["obs"][index_key])
            if len(ids) == 0 or len(set(ids)) != len(ids) or ids_seen.intersection(ids) or np.any(np.char.strip(ids) == ""):
                raise ValueError(f"Cell IDs must be nonempty and globally unique: {path}")
            ids_seen.update(ids)
            xy = handle["obsm"][spatial_key][()]
            features = handle["obsm"][representation_key][()]
            if xy.shape != (len(ids), 2) or features.ndim != 2 or len(features) != len(ids) or features.shape[1] < 1:
                raise ValueError(f"Expected Nx2 spatial and NxD representation: {path}")
            if not np.isfinite(xy).all() or not np.isfinite(features).all() or np.any(np.linalg.norm(features, axis=1) == 0):
                raise ValueError(f"Nonfinite coordinates/features or zero-norm representation: {path}")
            dimensions.add(features.shape[1])
            features_f32 = np.asarray(features, dtype=np.float32)
            if (not np.isfinite(features_f32).all()
                    or np.any(np.linalg.norm(features_f32.astype(np.float64), axis=1) == 0)):
                raise ValueError('Representation cannot be represented as finite nonzero-norm float32')
            onehot = bool(np.isin(features, [0, 1]).all() and np.all(features.sum(axis=1) == 1))
            labels = None
            if annotation_qc == 'provided':
                if annotation_key not in handle['obs']:
                    raise ValueError(f'annotation-qc provided requires obs[{annotation_key!r}] in every slice')
                labels = column(handle['obs'][annotation_key])
                if len(labels) != len(ids) or np.isin(np.char.lower(np.char.strip(labels)), ['', 'nan', 'none', 'null']).any():
                    raise ValueError(f'Annotation labels must be real, nonmissing, nonblank values: {path}')
                if '__SPATEO_QC_UNLABELED__' in labels:
                    raise ValueError('Reserved internal non-biological QC sentinel cannot be supplied as an annotation')
            if representation == "spatial-only":
                if features.shape[1] != 30 or not np.all(features == 1):
                    raise ValueError(f"Spatial-only requires identical nonzero ones30 vectors: {path}")
            elif representation == "expression-pca":
                if onehot or (len(features) > 1 and np.all(features == features[0])):
                    raise ValueError(f"Expression PCA cannot be one-hot or constant: {path}")
            elif representation == "annotation-onehot":
                if not onehot or annotation_key not in handle["obs"]:
                    raise ValueError(f"Annotation mode requires exact one-hot and obs[{annotation_key!r}]: {path}")
                for label in np.unique(labels):
                    group = features[labels == label]
                    if not np.all(group == group[0]) or (label in label_vectors and not np.array_equal(group[0], label_vectors[label])):
                        raise ValueError(f"Annotation encoding differs within/across slices for {label!r}")
                    label_vectors[label] = group[0].copy()
            z_metadata = {}
            for key in ("z", "physical_z"):
                if key in handle["obs"]:
                    values = column(handle["obs"][key]).astype(float)
                    if len(values) != len(ids) or not np.isfinite(values).all() or len(np.unique(values)) != 1:
                        raise ValueError(f"obs[{key!r}] must be constant and finite within a slice: {path}")
                    z_metadata[key] = float(values[0])
            if len(z_metadata) == 2 and z_metadata["z"] != z_metadata["physical_z"]:
                raise ValueError(f"z and physical_z metadata disagree: {path}")
            records.append({"path": path, "slice_id": path.name.split(".Spatial")[0],
                            "slice_number": number, "cell_ids": ids, "xy": xy, "sha256": sha(path),
                            "features_sha256": feature_sha(features), "representation_key": representation_key,
                            "spatial_key": spatial_key,
                            "physical_z": next(iter(z_metadata.values()), None)})
    if len(dimensions) != 1:
        raise ValueError("All slices need the same representation dimension and shared feature basis")
    if label_vectors:
        columns = [int(np.argmax(value)) for value in label_vectors.values()]
        if len(set(columns)) != len(columns):
            raise ValueError("Distinct annotations share a one-hot column")
    z_values = [record["physical_z"] for record in records]
    if any(value is not None for value in z_values):
        if any(value is None for value in z_values) or not np.all(np.diff(z_values) > 0):
            raise ValueError("Slice filenames must encode ascending physical z consistently; create correctly named blind copies before alignment")
    pca_audit = None
    if representation == 'expression-pca':
        pca_audit = validate_pca_manifest(pca_provenance or Path(slice_dir).parent/'expression_pca_manifest.json', records, next(iter(dimensions)))
    summary = {"n_slices": len(records), "n_cells": len(ids_seen), "representation": representation,
               "annotation_qc": annotation_qc,
               "representation_dimension": next(iter(dimensions)), "global_ids_unique": True,
               "source_reference_files_read": False,
               "physical_z_order": "strictly_increasing" if all(value is not None for value in z_values) else "metadata_absent_SL_order_only",
               "inputs": [{"path": str(record["path"]), "sha256": record["sha256"],
                           "n_cells": len(record["cell_ids"])} for record in records],
               "shared_expression_basis_provenance": pca_audit,
               "annotation_source": f"obs[{annotation_key}]" if annotation_qc == 'provided' else "disabled; internal non-biological geometry sentinel only"}
    return summary, records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slice-dir", type=Path, required=True)
    parser.add_argument("--representation", choices=["expression-pca", "annotation-onehot", "spatial-only"], required=True)
    parser.add_argument("--annotation-key", default="anno")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--representation-key", default="X_pca")
    parser.add_argument("--filename-style", choices=["pairwise", "continuity"], default="pairwise")
    parser.add_argument('--annotation-qc', choices=['off', 'provided'])
    parser.add_argument('--pca-provenance', type=Path)
    args = parser.parse_args()
    summary, _ = validate_inputs(**vars(args))
    print(json.dumps(summary, indent=2))
