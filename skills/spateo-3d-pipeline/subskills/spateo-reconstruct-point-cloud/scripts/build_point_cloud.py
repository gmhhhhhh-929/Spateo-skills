#!/usr/bin/env python3
"""Build, preview, save, and round-trip a Spateo 3D point-cloud model."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Sequence

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

import anndata as ad
import matplotlib as mpl
import numpy as np
import pandas as pd
import pyvista as pv
import spateo as st

pv.OFF_SCREEN = True

SCHEMA_VERSION = "spateo-point-cloud-v1"
SUPPORTED_SOURCE_REVISION = "615644f88613bea8ceb2e2df1e2391d16de55ec1"
SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_value(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    return value


def _load_json(path: Optional[str], *, label: str) -> Any:
    if path is None:
        return None
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise ValueError(f"{label} JSON does not exist: {source}")
    try:
        return json.loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid {label} JSON: {source}: {exc}") from exc


def _prepare_output_dir(path: Path) -> None:
    if path.exists():
        if not path.is_dir():
            raise ValueError(f"Output path exists and is not a directory: {path}")
        if any(path.iterdir()):
            raise FileExistsError(
                f"Refusing to write into non-empty output directory: {path}"
            )
    else:
        path.mkdir(parents=True, exist_ok=False)


def _coordinates(
    adata: ad.AnnData, spatial_key: str, allow_planar: bool
) -> tuple[np.ndarray, dict[str, Any]]:
    if spatial_key not in adata.obsm:
        raise ValueError(f"Missing adata.obsm[{spatial_key!r}]")
    raw = adata.obsm[spatial_key]
    if hasattr(raw, "toarray"):
        raw = raw.toarray()
    try:
        coords = np.asarray(raw, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"adata.obsm[{spatial_key!r}] must be numeric") from exc
    expected = (adata.n_obs, 3)
    if coords.shape != expected:
        raise ValueError(
            f"adata.obsm[{spatial_key!r}] must have shape {expected}; observed {coords.shape}"
        )
    if adata.n_obs == 0:
        raise ValueError("Point-cloud input contains no observations")
    if not np.isfinite(coords).all():
        bad = int(coords.size - np.isfinite(coords).sum())
        raise ValueError(
            f"adata.obsm[{spatial_key!r}] contains {bad} non-finite coordinate values"
        )
    centered = coords - coords.mean(axis=0, keepdims=True)
    spans = np.ptp(coords, axis=0)
    scaled = centered / np.where(spans > 0, spans, 1.0)
    rank = int(np.linalg.matrix_rank(scaled))
    if rank < 3 and not allow_planar:
        raise ValueError(
            f"adata.obsm[{spatial_key!r}] has coordinate rank {rank}, not full 3D support; "
            "use --allow-planar only when this is intentional"
        )
    minimum = coords.min(axis=0)
    maximum = coords.max(axis=0)
    return coords, {
        "key": spatial_key,
        "shape": list(coords.shape),
        "coordinate_rank": rank,
        "allow_planar": bool(allow_planar),
        "bounds": {
            "x": [float(minimum[0]), float(maximum[0])],
            "y": [float(minimum[1]), float(maximum[1])],
            "z": [float(minimum[2]), float(maximum[2])],
        },
        "span": [float(v) for v in spans],
    }


def _read_external_labels(
    adata: ad.AnnData,
    filename: str,
    id_column: str,
    value_column: str,
) -> tuple[str, np.ndarray, dict[str, Any]]:
    path = Path(filename).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"External label table does not exist: {path}")
    separator = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    table = pd.read_csv(path, sep=separator)
    missing_columns = [
        column for column in (id_column, value_column) if column not in table.columns
    ]
    if missing_columns:
        raise ValueError(f"External label table is missing columns: {missing_columns}")
    ids = table[id_column].astype(str)
    if ids.duplicated().any():
        examples = ids[ids.duplicated(keep=False)].head(5).tolist()
        raise ValueError(f"External label IDs are not unique; examples: {examples}")
    expected = pd.Index(adata.obs_names.astype(str))
    observed = pd.Index(ids)
    missing = expected.difference(observed)
    extra = observed.difference(expected)
    if len(missing) or len(extra):
        raise ValueError(
            "External labels must cover exactly the AnnData observations; "
            f"missing={len(missing)}, extra={len(extra)}"
        )
    values = (
        pd.Series(table[value_column].to_numpy(), index=observed)
        .reindex(expected)
        .to_numpy()
    )
    if pd.isna(values).any():
        raise ValueError("External labels contain missing values")
    temporary_key = "__spateo_point_cloud_external__"
    while temporary_key in adata.obs:
        temporary_key += "_"
    adata.obs[temporary_key] = values
    return (
        temporary_key,
        values,
        {
            "path": str(path),
            "sha256": _sha256(path),
            "id_column": id_column,
            "value_column": value_column,
        },
    )


def _is_numeric(values: np.ndarray) -> bool:
    return bool(np.issubdtype(np.asarray(values).dtype, np.number))


def _validate_color_contract(
    values: np.ndarray,
    palette: Any,
    alphamap: Any,
    masks: Sequence[str],
) -> tuple[bool, list[str]]:
    values = np.asarray(values).reshape(-1)
    if pd.isna(values).any():
        raise ValueError("Color values contain missing values")
    numeric = _is_numeric(values)
    if numeric:
        if not np.isfinite(values.astype(float)).all():
            raise ValueError("Continuous color values contain non-finite values")
        if masks:
            raise ValueError("--mask is supported only for categorical labels")
        if isinstance(palette, dict):
            raise ValueError(
                "A continuous scalar cannot use a label-to-color palette dictionary"
            )
        if isinstance(palette, str) and palette not in mpl.colormaps:
            raise ValueError(
                f"Continuous values require a registered Matplotlib colormap; observed {palette!r}"
            )
        if isinstance(alphamap, dict):
            raise ValueError(
                "A continuous scalar cannot use a label-to-alpha dictionary"
            )
        if float(alphamap) != 1.0:
            raise ValueError(
                "Spateo does not encode alphamap for continuous scalars; control opacity when plotting instead"
            )
        return True, []

    labels = values.astype(str)
    categories = sorted(np.unique(labels).tolist())
    category_set = set(categories)
    missing_masks = sorted(set(masks) - category_set)
    if missing_masks:
        raise ValueError(f"Masked labels are absent from the data: {missing_masks}")
    if "mask" in category_set and "mask" not in masks:
        raise ValueError(
            "The literal categorical label 'mask' is reserved by Spateo for transparent points; "
            "rename it or pass --mask mask to confirm that intent"
        )
    visible = category_set - set(masks)
    if isinstance(palette, dict):
        missing = sorted(visible - {str(k) for k in palette})
        if missing:
            raise ValueError(
                f"Palette dictionary does not cover categorical labels: {missing[:20]}"
            )
    elif isinstance(palette, list) and len(palette) < len(categories):
        raise ValueError(
            f"Palette list has {len(palette)} colors for {len(categories)} categories"
        )
    if isinstance(alphamap, dict):
        missing = sorted(visible - {str(k) for k in alphamap})
        if missing:
            raise ValueError(
                f"Alpha dictionary does not cover categorical labels: {missing[:20]}"
            )
        if not all(
            isinstance(v, (int, float)) and 0 <= float(v) <= 1
            for v in alphamap.values()
        ):
            raise ValueError(
                "Alpha dictionary values must be numeric values between 0 and 1"
            )
    return False, categories


def _label_summary(values: np.ndarray) -> dict[str, Any]:
    values = np.asarray(values).reshape(-1)
    if _is_numeric(values):
        numeric = values.astype(float)
        quantiles = np.quantile(numeric, [0.0, 0.25, 0.5, 0.75, 1.0])
        return {
            "kind": "continuous",
            "minimum": float(quantiles[0]),
            "q25": float(quantiles[1]),
            "median": float(quantiles[2]),
            "q75": float(quantiles[3]),
            "maximum": float(quantiles[4]),
        }
    labels, counts = np.unique(values.astype(str), return_counts=True)
    order = np.argsort(-counts, kind="stable")
    limit = 100
    return {
        "kind": "categorical",
        "n_categories": int(len(labels)),
        "counts": {str(labels[i]): int(counts[i]) for i in order[:limit]},
        "counts_truncated": bool(len(labels) > limit),
    }


def _preview_model(
    pc: pv.PolyData, max_points: int, seed: int
) -> tuple[pv.PolyData, dict[str, Any]]:
    if max_points <= 0:
        raise ValueError("--preview-max-points must be positive")
    if pc.n_points <= max_points:
        return pc, {"sampling": "all", "n_points": int(pc.n_points), "seed": None}
    rng = np.random.default_rng(seed)
    indices = np.sort(rng.choice(pc.n_points, size=max_points, replace=False))
    preview = pv.PolyData(np.asarray(pc.points)[indices])
    for key in pc.point_data.keys():
        preview.point_data[key] = np.asarray(pc.point_data[key])[indices]
    return preview, {
        "sampling": "uniform_without_replacement",
        "n_points": int(max_points),
        "seed": int(seed),
    }


def _render_preview(
    pc: pv.PolyData,
    output: Path,
    key_added: str,
    plot_cmap: Any,
    max_points: int,
    seed: int,
    point_size: float,
) -> dict[str, Any]:
    if point_size <= 0:
        raise ValueError("--point-size must be positive")
    preview, sampling = _preview_model(pc, max_points=max_points, seed=seed)
    rgba_key = f"{key_added}_rgba"
    categorical = rgba_key in preview.point_data
    legend_entries: list[list[Any]] = []
    legend_summary: dict[str, Any] = {"included": False, "categories": []}
    if categorical:
        full_labels = np.asarray(pc.point_data[key_added]).astype(str)
        full_rgba = np.asarray(pc.point_data[rgba_key], dtype=float)
        for label in sorted(np.unique(full_labels)):
            indices = np.flatnonzero(full_labels == label)
            visible = indices[full_rgba[indices, 3] > 0]
            if len(visible):
                legend_entries.append(
                    [
                        label,
                        tuple(float(value) for value in full_rgba[int(visible[0]), :3]),
                    ]
                )
        if 0 < len(legend_entries) <= 20:
            legend_summary = {
                "included": True,
                "categories": [entry[0] for entry in legend_entries],
            }
        elif len(legend_entries) > 20:
            legend_summary = {
                "included": False,
                "categories": [entry[0] for entry in legend_entries],
                "reason": "more_than_20_visible_categories",
            }
    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(1400, 1200))
    plotter.set_background((247 / 255, 248 / 255, 248 / 255))
    views = [
        (0, 0, "XY", plotter.view_xy),
        (0, 1, "XZ", plotter.view_xz),
        (1, 0, "YZ", plotter.view_yz),
        (1, 1, "Isometric", plotter.view_isometric),
    ]
    for row, column, title, view_method in views:
        plotter.subplot(row, column)
        kwargs: dict[str, Any] = {
            "style": "points",
            "point_size": point_size,
            "render_points_as_spheres": False,
            "show_scalar_bar": False,
        }
        if categorical:
            kwargs.update({"scalars": rgba_key, "rgba": True})
        else:
            kwargs.update({"scalars": key_added, "cmap": plot_cmap})
            if title == "Isometric":
                kwargs["show_scalar_bar"] = True
        plotter.add_mesh(preview, **kwargs)
        plotter.add_text(title, position="upper_left", font_size=12, color="black")
        view_method()
        plotter.reset_camera()
        if title == "Isometric" and legend_summary["included"]:
            plotter.add_legend(
                legend_entries,
                bcolor=(247 / 255, 248 / 255, 248 / 255),
                loc="upper right",
            )
    plotter.screenshot(str(output), transparent_background=False)
    plotter.close()
    sampling.update(
        {
            "path": str(output),
            "sha256": _sha256(output),
            "views": ["xy", "xz", "yz", "isometric"],
            "legend": legend_summary,
        }
    )
    return sampling


def _assert_roundtrip(
    original: pv.PolyData, loaded: Any, key_added: str
) -> dict[str, Any]:
    if not isinstance(loaded, pv.PolyData):
        raise RuntimeError(
            f"Spateo reloaded {type(loaded).__name__}, expected pyvista.PolyData"
        )
    if loaded.n_points != original.n_points:
        raise RuntimeError(
            f"VTK round-trip changed point count: {original.n_points} -> {loaded.n_points}"
        )
    coordinate_error = np.abs(np.asarray(original.points) - np.asarray(loaded.points))
    maximum_error = float(coordinate_error.max(initial=0.0))
    if not np.allclose(original.points, loaded.points, rtol=1e-10, atol=1e-10):
        raise RuntimeError(
            f"VTK round-trip changed coordinates; max_abs_error={maximum_error}"
        )
    expected_arrays = set(original.point_data.keys())
    observed_arrays = set(loaded.point_data.keys())
    if expected_arrays != observed_arrays:
        raise RuntimeError(
            "VTK round-trip changed point-data arrays; "
            f"missing={sorted(expected_arrays - observed_arrays)}, extra={sorted(observed_arrays - expected_arrays)}"
        )
    for key in sorted(expected_arrays):
        before = np.asarray(original.point_data[key])
        after = np.asarray(loaded.point_data[key])
        if before.shape != after.shape:
            raise RuntimeError(f"VTK round-trip changed shape of point_data[{key!r}]")
        if _is_numeric(before):
            equal = np.allclose(before, after, rtol=1e-6, atol=1e-7, equal_nan=True)
        else:
            equal = np.array_equal(before.astype(str), after.astype(str))
        if not equal:
            raise RuntimeError(f"VTK round-trip changed point_data[{key!r}]")
    required = {"obs_index", key_added}
    missing_required = required - observed_arrays
    if missing_required:
        raise RuntimeError(
            f"Reloaded point cloud is missing required arrays: {sorted(missing_required)}"
        )
    return {
        "reader": "spateo.tdr.read_model",
        "dataset_type": type(loaded).__name__,
        "n_points_equal": True,
        "coordinates_allclose": True,
        "maximum_absolute_coordinate_error": maximum_error,
        "point_data_equal": True,
    }


def build_point_cloud(args: argparse.Namespace) -> dict[str, Any]:
    input_path = Path(args.input_h5ad).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    if not input_path.is_file():
        raise ValueError(f"Input H5AD does not exist: {input_path}")
    if not SAFE_NAME.fullmatch(args.name) or args.name.lower().endswith(".vtk"):
        raise ValueError("--name must be a safe filename stem without a .vtk extension")
    if args.key_added == "obs_index" or args.key_added.endswith("_rgba"):
        raise ValueError("--key-added conflicts with reserved point-data array names")
    _prepare_output_dir(output_dir)

    adata = ad.read_h5ad(input_path)
    if not adata.obs_names.is_unique:
        raise ValueError(
            "AnnData obs_names must be unique for point-to-observation traceability"
        )
    _, spatial_summary = _coordinates(adata, args.spatial_key, args.allow_planar)
    spatial_summary["unit"] = args.coordinate_unit
    spatial_summary["frame"] = args.coordinate_frame

    palette = _load_json(args.palette_json, label="palette")
    if palette is None:
        palette = args.colormap
    if not isinstance(palette, (str, list, dict)):
        raise ValueError(
            "Palette must be a color/cmap string, JSON list, or JSON dictionary"
        )
    alphamap = _load_json(args.alpha_json, label="alpha")
    if alphamap is None:
        alphamap = args.alpha
    if isinstance(alphamap, list):
        raise ValueError(
            "Alpha lists are ambiguous in the current Spateo API; use a scalar or full label dictionary"
        )
    if not isinstance(alphamap, (int, float, dict)):
        raise ValueError("Alpha must be a scalar or JSON dictionary")
    if isinstance(alphamap, (int, float)) and not 0 <= float(alphamap) <= 1:
        raise ValueError("--alpha must be between 0 and 1")

    groupby: Any = None
    source_values = np.asarray(["same"] * adata.n_obs)
    coloring: dict[str, Any] = {"mode": "uniform", "groupby": None}
    if args.obs is not None:
        if args.obs not in adata.obs:
            raise ValueError(f"Missing adata.obs[{args.obs!r}]")
        groupby = args.obs
        source_values = np.asarray(adata.obs[args.obs].to_numpy())
        if not _is_numeric(source_values):
            adata.obs[args.obs] = adata.obs[args.obs].astype(str)
            source_values = np.asarray(adata.obs[args.obs].to_numpy())
        coloring = {"mode": "obs", "groupby": args.obs}
    elif args.gene:
        if not adata.var_names.is_unique:
            raise ValueError("AnnData var_names must be unique for gene-based coloring")
        genes = list(args.gene)
        if len(set(genes)) != len(genes):
            raise ValueError("--gene values must be unique")
        missing_genes = sorted(set(genes) - set(adata.var_names.astype(str)))
        if missing_genes:
            raise ValueError(f"Genes are absent from adata.var_names: {missing_genes}")
        if args.layer != "X" and args.layer not in adata.layers:
            raise ValueError(f"Missing adata.layers[{args.layer!r}]")
        groupby = genes[0] if len(genes) == 1 else tuple(genes)
        matrix = (
            adata[:, genes].X
            if args.layer == "X"
            else adata[:, genes].layers[args.layer]
        )
        source_values = np.asarray(matrix.sum(axis=1)).reshape(-1)
        coloring = {
            "mode": "gene" if len(genes) == 1 else "gene_sum",
            "genes": genes,
            "layer": args.layer,
        }
    elif args.labels_file is not None:
        groupby, source_values, external = _read_external_labels(
            adata,
            filename=args.labels_file,
            id_column=args.labels_id_column,
            value_column=args.labels_value_column,
        )
        if not _is_numeric(source_values):
            adata.obs[groupby] = adata.obs[groupby].astype(str)
            source_values = np.asarray(adata.obs[groupby].to_numpy())
        coloring = {
            "mode": "external",
            "groupby": args.labels_value_column,
            "external_labels": external,
        }

    if coloring["mode"] == "uniform" and args.mask:
        raise ValueError("--mask requires categorical --obs or --labels-file values")

    numeric, raw_categories = _validate_color_contract(
        source_values, palette, alphamap, args.mask
    )
    construct_palette = palette.copy() if isinstance(palette, dict) else palette
    construct_alphamap = alphamap.copy() if isinstance(alphamap, dict) else alphamap
    if not numeric and args.mask:
        if isinstance(construct_palette, dict):
            for label in args.mask:
                construct_palette.setdefault(label, "gainsboro")
        if isinstance(construct_alphamap, dict):
            for label in args.mask:
                construct_alphamap.setdefault(label, 0.0)
    pc, plot_cmap = st.tdr.construct_pc(
        adata=adata,
        layer=args.layer,
        spatial_key=args.spatial_key,
        groupby=groupby,
        key_added=args.key_added,
        mask=None,
        colormap=construct_palette,
        alphamap=construct_alphamap,
    )
    if pc.n_points != adata.n_obs:
        raise RuntimeError(
            f"Spateo constructed {pc.n_points} points for {adata.n_obs} observations"
        )
    if "obs_index" not in pc.point_data or not np.array_equal(
        np.asarray(pc.point_data["obs_index"]).astype(str),
        adata.obs_names.astype(str).to_numpy(),
    ):
        raise RuntimeError(
            "Spateo point obs_index does not preserve the input observation order"
        )
    # Apply display masking without replacing the source biological label.
    # This also avoids the pinned source's scalar-alpha mask overwrite bug.
    rgba_key = f"{args.key_added}_rgba"
    if args.mask and rgba_key in pc.point_data:
        effective_labels = np.asarray(pc.point_data[args.key_added]).astype(str)
        rgba = np.asarray(pc.point_data[rgba_key], dtype=np.float32).copy()
        rgba[np.isin(effective_labels, np.asarray(args.mask, dtype=str)), 3] = 0.0
        pc.point_data[rgba_key] = rgba
    labels = np.asarray(pc.point_data[args.key_added])
    if len(labels) != adata.n_obs:
        raise RuntimeError(
            "Spateo point labels no longer match the input observation count"
        )
    if _is_numeric(labels) and not np.isfinite(labels.astype(float)).all():
        raise ValueError(
            "Constructed point-cloud scalar values contain non-finite values"
        )

    vtk_path = output_dir / f"{args.name}.vtk"
    preview_path = output_dir / f"{args.name}.preview.png"
    manifest_path = output_dir / f"{args.name}.manifest.json"
    st.tdr.save_model(pc, str(vtk_path), binary=not args.ascii)
    loaded = st.tdr.read_model(str(vtk_path))
    roundtrip = _assert_roundtrip(pc, loaded, key_added=args.key_added)

    preview: Optional[dict[str, Any]] = None
    if not args.skip_preview:
        preview = _render_preview(
            pc=loaded,
            output=preview_path,
            key_added=args.key_added,
            plot_cmap=plot_cmap,
            max_points=args.preview_max_points,
            seed=args.preview_seed,
            point_size=args.point_size,
        )

    coloring.update(
        {
            "key_added": args.key_added,
            "numeric": numeric,
            "raw_categories": raw_categories if not numeric else None,
            "mask": list(args.mask),
            "requested_colormap_or_palette": _json_value(palette),
            "effective_colormap_or_palette": _json_value(construct_palette),
            "requested_alphamap": _json_value(alphamap),
            "effective_alphamap": _json_value(construct_alphamap),
            "plot_cmap_returned_by_spateo": _json_value(plot_cmap),
            "stored_arrays": [args.key_added]
            + (
                [f"{args.key_added}_rgba"]
                if f"{args.key_added}_rgba" in pc.point_data
                else []
            ),
            "summary": _label_summary(labels),
        }
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "pass",
        "source_contract_revision": SUPPORTED_SOURCE_REVISION,
        "input": {
            "h5ad": str(input_path),
            "sha256": _sha256(input_path),
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "obs_names_unique": True,
        },
        "spatial": spatial_summary,
        "coloring": coloring,
        "model": {
            "type": type(pc).__name__,
            "n_points": int(pc.n_points),
            "n_cells": int(pc.n_cells),
            "point_data": sorted(pc.point_data.keys()),
        },
        "output": {
            "vtk": str(vtk_path),
            "vtk_sha256": _sha256(vtk_path),
            "vtk_binary": not args.ascii,
            "preview": preview,
            "manifest": str(manifest_path),
        },
        "roundtrip": roundtrip,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "anndata": ad.__version__,
            "pyvista": pv.__version__,
            "spateo": getattr(st, "__version__", "unknown"),
            "spateo_module": str(Path(st.__file__).resolve()),
        },
    }
    manifest_path.write_text(
        json.dumps(_json_value(manifest), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate an H5AD, construct a Spateo point cloud, save .vtk, and verify Spateo reload."
    )
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument(
        "--output-dir",
        required=True,
        help="A new or empty directory; existing contents are never overwritten.",
    )
    parser.add_argument(
        "--name", default="point_cloud", help="Output filename stem (without .vtk)."
    )
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument(
        "--coordinate-unit", help="Known coordinate unit to record; never inferred."
    )
    parser.add_argument(
        "--coordinate-frame",
        help="Known coordinate-frame description to record; never inferred.",
    )
    parser.add_argument("--key-added", default="groups")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--obs", help="Color by an adata.obs column.")
    source.add_argument(
        "--gene",
        action="append",
        help="Color by one gene or the sum of repeated --gene arguments.",
    )
    source.add_argument(
        "--labels-file",
        help="CSV/TSV with exact observation IDs and user-provided values.",
    )
    parser.add_argument(
        "--layer", default="X", help="X or an adata.layers key for gene coloring."
    )
    parser.add_argument("--labels-id-column", default="obs_index")
    parser.add_argument("--labels-value-column", default="value")
    parser.add_argument(
        "--colormap", default="rainbow", help="Matplotlib cmap name or a single color."
    )
    parser.add_argument(
        "--palette-json",
        help="JSON file containing a color list or complete label-to-color dictionary.",
    )
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument(
        "--alpha-json",
        help="JSON file containing a complete label-to-alpha dictionary.",
    )
    parser.add_argument(
        "--mask",
        action="append",
        default=[],
        help="Categorical value to make transparent; repeat as needed.",
    )
    parser.add_argument(
        "--allow-planar",
        action="store_true",
        help="Allow rank-deficient n_obs x 3 coordinates.",
    )
    parser.add_argument(
        "--ascii",
        action="store_true",
        help="Write ASCII VTK instead of the smaller default binary VTK.",
    )
    parser.add_argument("--skip-preview", action="store_true")
    parser.add_argument("--preview-max-points", type=int, default=200000)
    parser.add_argument("--preview-seed", type=int, default=0)
    parser.add_argument("--point-size", type=float, default=2.0)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = build_point_cloud(args)
    print(
        json.dumps(
            {
                "status": "pass",
                "vtk": manifest["output"]["vtk"],
                "manifest": manifest["output"]["manifest"],
            }
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
