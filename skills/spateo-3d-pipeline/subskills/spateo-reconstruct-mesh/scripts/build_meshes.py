#!/usr/bin/env python3
"""Build selectable surface meshes from a validated Spateo point-cloud VTK."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Sequence

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

import matplotlib as mpl
import numpy as np
import pyvista as pv
import spateo as st
from scipy import ndimage
from skimage import measure

pv.OFF_SCREEN = True

SCHEMA_VERSION = "spateo-surface-mesh-v1"
VALIDATED_SOURCE_REVISION = "615644f88613bea8ceb2e2df1e2391d16de55ec1"
VALIDATED_MESH_METHODS_SHA256 = (
    "1dd0ead7338df143903178c7357ddd81cdd1cda2d3c4591b7509c613a7fc76fa"
)
SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
DEFAULT_COLORS = [
    "#d98cb3",
    "#79a857",
    "#e9c46a",
    "#6aaed6",
    "#b89ad7",
    "#ef8354",
    "#62b6a7",
    "#f28482",
]


@dataclass(frozen=True)
class Target:
    name: str
    values: tuple[str, ...]
    indices: np.ndarray
    color: str
    opacity: float
    role: str


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
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _load_json_dict(path: Optional[str], label: str) -> dict[str, Any]:
    if path is None:
        return {}
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise ValueError(f"{label} JSON does not exist: {source}")
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid {label} JSON: {source}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} JSON must contain an object keyed by mesh name")
    return {str(key): item for key, item in value.items()}


def _prepare_output_dir(path: Path) -> Path:
    if path.exists():
        if not path.is_dir():
            raise ValueError(f"Output path exists and is not a directory: {path}")
        if any(path.iterdir()):
            raise FileExistsError(
                f"Refusing to write into non-empty output directory: {path}"
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=f".{path.name}.staging-", dir=str(path.parent)))


def _commit_output_dir(staging_dir: Path, output_dir: Path) -> None:
    if output_dir.exists():
        if not output_dir.is_dir() or any(output_dir.iterdir()):
            raise FileExistsError(
                f"Output destination changed while building: {output_dir}"
            )
        output_dir.rmdir()
    os.replace(staging_dir, output_dir)


def _parse_selection(text: str) -> tuple[str, tuple[str, ...]]:
    if "=" not in text:
        raise ValueError(
            f"Invalid --selection {text!r}; expected NAME=VALUE or NAME=VALUE1|VALUE2"
        )
    name, raw_values = text.split("=", 1)
    name = name.strip()
    values = tuple(value.strip() for value in raw_values.split("|") if value.strip())
    if not SAFE_NAME.fullmatch(name):
        raise ValueError(f"Invalid mesh name in --selection: {name!r}")
    if not values:
        raise ValueError(f"--selection {text!r} contains no label values")
    if len(set(values)) != len(values):
        raise ValueError(f"--selection {text!r} repeats a label value")
    return name, values


def _validate_color(value: Any, label: str) -> str:
    if not isinstance(value, str) or not mpl.colors.is_color_like(value):
        raise ValueError(f"Invalid color for {label}: {value!r}")
    return value


def _validate_opacity(value: Any, label: str) -> float:
    try:
        opacity = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid opacity for {label}: {value!r}") from exc
    if not 0 <= opacity <= 1:
        raise ValueError(f"Opacity for {label} must be between 0 and 1")
    return opacity


def _load_point_cloud(path: Path) -> tuple[pv.PolyData, np.ndarray, dict[str, Any]]:
    if not path.is_file():
        raise ValueError(f"Input point-cloud VTK does not exist: {path}")
    model = st.tdr.read_model(str(path))
    if not isinstance(model, pv.PolyData):
        raise ValueError(
            f"Input must reload as pyvista.PolyData; observed {type(model).__name__}"
        )
    if model.n_points == 0:
        raise ValueError("Input point cloud contains no points")
    if "obs_index" not in model.point_data:
        raise ValueError(
            "Input point cloud is missing point_data['obs_index']; build it with the point-cloud skill first"
        )
    coordinates = np.asarray(model.points, dtype=np.float64)
    if coordinates.shape != (model.n_points, 3) or not np.isfinite(coordinates).all():
        raise ValueError("Point-cloud coordinates must be finite n_points x 3 values")
    spans = np.ptp(coordinates, axis=0)
    scaled = (coordinates - coordinates.mean(axis=0)) / np.where(spans > 0, spans, 1)
    rank = int(np.linalg.matrix_rank(scaled))
    if rank < 3:
        raise ValueError(f"Point cloud has coordinate rank {rank}, not full 3D support")
    obs_ids = np.asarray(model.point_data["obs_index"]).astype(str)
    if len(np.unique(obs_ids)) != model.n_points:
        raise ValueError("point_data['obs_index'] must contain unique observation IDs")
    return (
        model,
        coordinates,
        {
            "path": str(path),
            "sha256": _sha256(path),
            "dataset_type": type(model).__name__,
            "n_points": int(model.n_points),
            "n_cells": int(model.n_cells),
            "coordinate_rank": rank,
            "bounds": [float(value) for value in model.bounds],
            "point_data": sorted(model.point_data.keys()),
        },
    )


def _make_targets(
    args: argparse.Namespace,
    point_cloud: pv.PolyData,
    palette: dict[str, Any],
    opacities: dict[str, Any],
) -> tuple[list[Target], Optional[np.ndarray], dict[str, int]]:
    parsed = [_parse_selection(text) for text in args.selection]
    names = [name for name, _ in parsed]
    folded_names = [name.casefold() for name in names]
    if len(folded_names) != len(set(folded_names)):
        raise ValueError(
            "Every --selection mesh name must be unique ignoring letter case"
        )
    if args.body_name.casefold() in set(folded_names):
        raise ValueError("--body-name conflicts with a --selection mesh name")

    labels: Optional[np.ndarray] = None
    label_counts: dict[str, int] = {}
    if parsed:
        if args.label_key is None:
            raise ValueError("--label-key is required when --selection is used")
        if args.label_key not in point_cloud.point_data:
            raise ValueError(
                f"Input point cloud is missing point_data[{args.label_key!r}]"
            )
        labels = np.asarray(point_cloud.point_data[args.label_key]).astype(str)
        unique, counts = np.unique(labels, return_counts=True)
        label_counts = {str(key): int(count) for key, count in zip(unique, counts)}

    targets: list[Target] = []
    if not args.no_body:
        color = _validate_color(
            palette.get(args.body_name, args.body_color), args.body_name
        )
        opacity = _validate_opacity(
            opacities.get(args.body_name, args.body_opacity), args.body_name
        )
        targets.append(
            Target(
                name=args.body_name,
                values=(),
                indices=np.arange(point_cloud.n_points, dtype=np.int64),
                color=color,
                opacity=opacity,
                role="body",
            )
        )

    for position, (name, values) in enumerate(parsed):
        assert labels is not None
        missing = [value for value in values if value not in label_counts]
        if missing:
            raise ValueError(
                f"Selection {name!r} references absent {args.label_key!r} values: {missing}"
            )
        indices = np.flatnonzero(np.isin(labels, np.asarray(values, dtype=str)))
        if len(indices) < args.minimum_selection_points:
            raise ValueError(
                f"Selection {name!r} has {len(indices)} points; require at least "
                f"{args.minimum_selection_points}"
            )
        color = _validate_color(
            palette.get(name, DEFAULT_COLORS[position % len(DEFAULT_COLORS)]), name
        )
        opacity = _validate_opacity(opacities.get(name, args.selection_opacity), name)
        targets.append(
            Target(
                name=name,
                values=values,
                indices=indices,
                color=color,
                opacity=opacity,
                role="selection",
            )
        )
    if not targets:
        raise ValueError("No meshes requested; remove --no-body or add --selection")
    return targets, labels, label_counts


def _shared_grid(
    coordinates: np.ndarray,
    target_long_axis: int,
    voxel_size: Optional[float],
    sigma: float,
    padding: Optional[int],
    max_voxels: int,
) -> dict[str, Any]:
    if sigma <= 0:
        raise ValueError("--sigma must be positive")
    spans = np.ptp(coordinates, axis=0)
    if voxel_size is None:
        if target_long_axis < 32:
            raise ValueError("--target-long-axis must be at least 32")
        voxel_size = float(spans.max() / target_long_axis)
    if not np.isfinite(voxel_size) or voxel_size <= 0:
        raise ValueError("--voxel-size must be positive and finite")
    if padding is None:
        padding = int(np.ceil(4 * sigma)) + 2
    if padding < int(np.ceil(3 * sigma)):
        raise ValueError("--padding must be at least ceil(3 * sigma)")
    origin = coordinates.min(axis=0) - padding * voxel_size
    shape = np.ceil(spans / voxel_size).astype(int) + 1 + 2 * padding
    # Python integers do not overflow; np.prod(..., dtype=int64) can wrap and
    # accidentally bypass the memory guard for very large requested grids.
    n_voxels = math.prod(int(value) for value in shape)
    if n_voxels > max_voxels:
        raise ValueError(
            f"Requested grid has {n_voxels:,} voxels, above --max-voxels={max_voxels:,}; "
            "increase --voxel-size or lower --target-long-axis"
        )
    # Store density samples at grid nodes and place every source point at its
    # nearest node.  Flooring would shift every splat by up to one full voxel
    # and measurably reduce final point containment near the surface.
    all_indices = np.rint((coordinates - origin) / voxel_size).astype(np.int64)
    if np.any(all_indices < 0) or np.any(all_indices >= shape):
        raise RuntimeError("Internal grid construction failed to contain all points")
    return {
        "origin": origin,
        "shape": tuple(int(value) for value in shape),
        "voxel_size": float(voxel_size),
        "padding": int(padding),
        "sigma": float(sigma),
        "n_voxels": n_voxels,
        "all_indices": all_indices,
    }


def _component_mask(
    density: np.ndarray,
    point_indices: np.ndarray,
    level: float,
    minimum_voxels: int,
    minimum_fraction: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    above = density >= level
    components, count = ndimage.label(
        above, structure=ndimage.generate_binary_structure(rank=3, connectivity=1)
    )
    sizes = np.bincount(components.ravel())
    foreground_sizes = sizes[1:]
    if count == 0 or len(foreground_sizes) == 0:
        return np.zeros_like(above), {
            "components_before_filter": 0,
            "components_kept": 0,
            "coverage": 0.0,
            "minimum_component_size": int(minimum_voxels),
            "voxel_connectivity": 6,
        }
    cutoff = max(
        int(minimum_voxels), int(np.ceil(foreground_sizes.max() * minimum_fraction))
    )
    keep_ids = np.flatnonzero(foreground_sizes >= cutoff) + 1
    kept = np.isin(components, keep_ids)
    coverage = float(
        kept[
            point_indices[:, 0],
            point_indices[:, 1],
            point_indices[:, 2],
        ].mean()
    )
    return kept, {
        "components_before_filter": int(count),
        "components_kept": int(len(keep_ids)),
        "component_voxels_before_filter": sorted(
            (int(value) for value in foreground_sizes), reverse=True
        )[:100],
        "component_voxels_kept": sorted(
            (int(sizes[index]) for index in keep_ids), reverse=True
        ),
        "minimum_component_size": int(cutoff),
        "voxel_connectivity": 6,
        "coverage": coverage,
    }


def _choose_level(
    density: np.ndarray,
    point_indices: np.ndarray,
    fixed_level: Optional[float],
    coverage_target: float,
    min_level: float,
    max_level: float,
    level_steps: int,
    minimum_voxels: int,
    minimum_fraction: float,
) -> tuple[float, np.ndarray, dict[str, Any]]:
    if not 0 < coverage_target <= 1:
        raise ValueError("Coverage targets must be in (0, 1]")
    if fixed_level is not None:
        if not 0 < fixed_level < 1:
            raise ValueError("--level must be between 0 and 1")
        levels = np.asarray([fixed_level], dtype=float)
        mode = "fixed"
    else:
        if not 0 < min_level < max_level < 1:
            raise ValueError("Require 0 < --min-level < --max-level < 1")
        if level_steps < 2:
            raise ValueError("--level-steps must be at least 2")
        levels = np.geomspace(max_level, min_level, level_steps)
        mode = "automatic_coverage"

    attempts: list[dict[str, Any]] = []
    chosen: Optional[tuple[float, np.ndarray, dict[str, Any]]] = None
    for candidate in levels:
        kept, summary = _component_mask(
            density=density,
            point_indices=point_indices,
            level=float(candidate),
            minimum_voxels=minimum_voxels,
            minimum_fraction=minimum_fraction,
        )
        attempts.append(
            {
                "level": float(candidate),
                "coverage": summary["coverage"],
                "components_kept": summary["components_kept"],
            }
        )
        chosen = (float(candidate), kept, summary)
        if fixed_level is not None or (
            summary["coverage"] >= coverage_target and summary["components_kept"] > 0
        ):
            break
    assert chosen is not None
    level, kept, summary = chosen
    summary.update(
        {
            "selection_mode": mode,
            "coverage_target": float(coverage_target),
            "target_met": bool(summary["coverage"] >= coverage_target),
            "attempts": attempts,
        }
    )
    return level, kept, summary


def _label_mesh(mesh: pv.PolyData, target: Target) -> pv.PolyData:
    mesh = mesh.extract_surface().triangulate().clean()
    # Some VTK filters retain auxiliary vertex/line cells alongside polygons.
    # Rebuild from the polygon connectivity so a delivered surface VTK contains
    # triangles only and mass/volume metrics do not emit mixed-cell warnings.
    mesh = pv.PolyData(np.asarray(mesh.points), np.asarray(mesh.faces)).clean()
    if mesh.n_points == 0 or mesh.n_cells == 0:
        raise ValueError(f"Reconstruction for {target.name!r} produced an empty mesh")
    mesh.cell_data["mesh_label"] = np.asarray([target.name] * mesh.n_cells)
    rgba = np.asarray(
        mpl.colors.to_rgba(target.color, alpha=target.opacity), dtype=float
    )
    mesh.cell_data["mesh_label_rgba"] = np.tile(rgba, (mesh.n_cells, 1))
    return mesh


def _repair_open_surface(
    mesh: pv.PolyData, enabled: bool, max_degenerate_cells: int = 16
) -> tuple[pv.PolyData, dict[str, Any]]:
    def component_count(surface: pv.PolyData) -> int:
        if surface.n_cells == 0:
            return 0
        connected = surface.connectivity()
        region_ids = np.asarray(connected.cell_data.get("RegionId", []))
        return int(len(np.unique(region_ids))) if len(region_ids) else 1

    before = int(mesh.n_open_edges)
    before_topology = {
        "n_points": int(mesh.n_points),
        "n_cells": int(mesh.n_cells),
        "components": component_count(mesh),
        "bounds": [float(value) for value in mesh.bounds],
    }
    applied = False
    component_repairs: list[dict[str, Any]] = []
    dropped_degenerate_components: list[dict[str, Any]] = []
    if enabled and before:
        try:
            import pymeshfix as mf
        except ImportError as exc:
            raise ImportError(
                "Open edges require pymeshfix; install the Spateo 3D dependencies "
                "or pass --skip-repair to retain them explicitly"
            ) from exc
        triangle_mesh = pv.PolyData(
            np.asarray(mesh.points), np.asarray(mesh.faces)
        ).triangulate()
        repaired_parts: list[pv.PolyData] = []
        for index, body in enumerate(triangle_mesh.split_bodies()):
            part = body.extract_surface().triangulate().clean()
            part = pv.PolyData(np.asarray(part.points), np.asarray(part.faces)).clean()
            part_open_before = int(part.n_open_edges)
            part_before = {
                "component": index,
                "n_points_before": int(part.n_points),
                "n_cells_before": int(part.n_cells),
                "open_edges_before": part_open_before,
                "surface_area_before": float(part.area),
            }
            repair_error = None
            if part_open_before:
                fixer = mf.MeshFix(part)
                # MeshFix otherwise defaults to remove_smallest_components=True.
                # Repair each connected body separately as an additional guard
                # against deleting legitimate bilateral/disconnected anatomy.
                try:
                    fixer.repair(joincomp=False, remove_smallest_components=False)
                    part = (
                        pv.PolyData(
                            np.asarray(fixer.mesh.points),
                            np.asarray(fixer.mesh.faces),
                        )
                        .triangulate()
                        .clean()
                    )
                except Exception as exc:  # pragma: no cover - backend detail
                    repair_error = f"{type(exc).__name__}: {exc}"
                    part = pv.PolyData()
            if part.n_points == 0 or part.n_cells == 0:
                if part_before["n_cells_before"] <= max_degenerate_cells:
                    dropped_degenerate_components.append(
                        {
                            **part_before,
                            "reason": "unrepairable_degenerate_fragment",
                            "repair_error": repair_error,
                        }
                    )
                    continue
                raise RuntimeError(
                    f"Mesh repair erased connected component {index} with "
                    f"{part_before['n_cells_before']} cells"
                )
            component_repairs.append(
                {
                    **part_before,
                    "open_edges_after": int(part.n_open_edges),
                    "n_points_after": int(part.n_points),
                    "n_cells_after": int(part.n_cells),
                }
            )
            repaired_parts.append(part)
        if not repaired_parts:
            raise RuntimeError("Mesh repair erased every connected component")
        mesh = repaired_parts[0]
        for part in repaired_parts[1:]:
            mesh = mesh.merge(part, merge_points=False)
        mesh = mesh.extract_surface().triangulate()
        applied = True
    after_topology = {
        "n_points": int(mesh.n_points),
        "n_cells": int(mesh.n_cells),
        "components": component_count(mesh),
        "bounds": [float(value) for value in mesh.bounds],
    }
    accounted_components = after_topology["components"] + len(
        dropped_degenerate_components
    )
    if applied and accounted_components < before_topology["components"]:
        raise RuntimeError(
            "Mesh repair reduced connected components from "
            f"{before_topology['components']} to {after_topology['components']}; "
            "refusing to discard or merge potentially valid disconnected anatomy"
        )
    return mesh, {
        "enabled": bool(enabled),
        "applied": applied,
        "open_edges_before": before,
        "open_edges_after": int(mesh.n_open_edges),
        "topology_before": before_topology,
        "topology_after": after_topology,
        "components_preserved": (
            after_topology["components"] >= before_topology["components"]
        ),
        "component_count_unchanged": (
            after_topology["components"] == before_topology["components"]
        ),
        "component_accounting_complete": (
            accounted_components >= before_topology["components"]
        ),
        "component_repairs": component_repairs,
        "dropped_degenerate_components": dropped_degenerate_components,
        "max_degenerate_cells": int(max_degenerate_cells),
    }


def _density_mesh(
    coordinates: np.ndarray,
    target: Target,
    grid: dict[str, Any],
    args: argparse.Namespace,
) -> tuple[pv.PolyData, dict[str, Any]]:
    target_points = coordinates[target.indices]
    point_indices = grid["all_indices"][target.indices]
    occupancy = np.zeros(grid["shape"], dtype=np.uint8)
    occupancy[point_indices[:, 0], point_indices[:, 1], point_indices[:, 2]] = 1
    density = ndimage.gaussian_filter(
        occupancy.astype(np.float32), sigma=grid["sigma"], mode="constant", cval=0
    )
    maximum = float(density.max(initial=0))
    if maximum <= 0:
        raise ValueError(f"Selection {target.name!r} produced an empty density field")
    density /= maximum
    coverage_target = (
        args.body_coverage_target
        if target.role == "body"
        else args.selection_coverage_target
    )
    level, kept, level_summary = _choose_level(
        density=density,
        point_indices=point_indices,
        fixed_level=args.level,
        coverage_target=coverage_target,
        min_level=args.min_level,
        max_level=args.max_level,
        level_steps=args.level_steps,
        minimum_voxels=args.min_component_voxels,
        minimum_fraction=args.min_component_fraction,
    )
    if args.level is None and not level_summary["target_met"]:
        raise ValueError(
            f"Selection {target.name!r} did not reach its density coverage target "
            f"{coverage_target:.3f}; best coverage was "
            f"{level_summary['coverage']:.3f}. Review component thresholds, sigma, "
            "or grid resolution instead of accepting an unqualified mesh."
        )
    if not kept.any():
        raise ValueError(
            f"Selection {target.name!r} has no supported density component at level {level}"
        )
    field = density.copy()
    field[~kept] = 0
    if (
        np.any(kept[[0, -1], :, :])
        or np.any(kept[:, [0, -1], :])
        or np.any(kept[:, :, [0, -1]])
    ):
        raise RuntimeError(
            f"Selection {target.name!r} reaches the grid boundary; increase --padding"
        )
    vertices, faces, _, _ = measure.marching_cubes(
        field,
        level=level,
        spacing=(grid["voxel_size"],) * 3,
        allow_degenerate=False,
        method="lewiner",
    )
    vertices += grid["origin"]
    vtk_faces = np.column_stack(
        [np.full(len(faces), 3, dtype=np.int64), faces.astype(np.int64)]
    )
    mesh = pv.PolyData(vertices, vtk_faces.ravel()).triangulate().clean()
    if args.smooth > 0:
        mesh = mesh.smooth_taubin(
            n_iter=args.smooth,
            pass_band=args.taubin_pass_band,
            boundary_smoothing=True,
            feature_smoothing=False,
        )
    if args.decimate > 0:
        mesh = mesh.decimate_pro(
            args.decimate,
            preserve_topology=True,
            boundary_vertex_deletion=False,
        )
    mesh, repair = _repair_open_surface(
        mesh,
        enabled=not args.skip_repair,
        max_degenerate_cells=args.repair_max_degenerate_cells,
    )
    mesh = _label_mesh(mesh, target)
    return mesh, {
        "algorithm": "occupancy_gaussian_marching_cubes",
        "source_points": int(len(target_points)),
        "occupied_voxels": int(occupancy.sum()),
        "density_maximum_before_normalization": maximum,
        "level": float(level),
        "level_selection": level_summary,
        "smooth_filter": "taubin",
        "smooth_iterations": int(args.smooth),
        "taubin_pass_band": float(args.taubin_pass_band),
        "decimate_reduction": float(args.decimate),
        "repair": repair,
    }


def _spateo_mesh(
    coordinates: np.ndarray,
    target: Target,
    args: argparse.Namespace,
) -> tuple[pv.PolyData, dict[str, Any]]:
    # construct_surface performs MeshFix and ACVD remeshing after the native
    # marching-cubes core.  On real serial-section data, split_bodies may hand
    # mixed VTK connectivity to MeshFix and fail before a surface is returned.
    # Use the pinned Spateo core directly, then normalize to triangles here.
    from spateo.tdr.models.models_individual.mesh_methods import marching_cube_mesh

    subset = pv.PolyData(coordinates[target.indices])
    subset.point_data["obs_index"] = np.asarray(
        [str(index) for index in target.indices]
    )
    # The pinned implementation samples with the legacy global NumPy RNG.
    # Derive a per-target seed so adding another selection cannot change an
    # existing target's geometry.
    target_seed = int.from_bytes(
        hashlib.sha256(f"{args.random_seed}:{target.name}".encode()).digest()[:4],
        byteorder="little",
    )
    random_state = np.random.get_state()
    try:
        np.random.seed(target_seed)
        surface = marching_cube_mesh(
            pc=subset,
            levelset=args.levelset,
            mc_scale_factor=args.mc_scale_factor,
            dist_sample_num=args.dist_sample_num,
        )
    finally:
        np.random.set_state(random_state)
    surface = (
        pv.PolyData(np.asarray(surface.points), np.asarray(surface.faces))
        .triangulate()
        .clean()
    )
    if args.smooth > 0:
        surface = surface.smooth(n_iter=args.smooth)
    surface, repair = _repair_open_surface(
        surface,
        enabled=not args.skip_repair,
        max_degenerate_cells=args.repair_max_degenerate_cells,
    )
    mesh = _label_mesh(surface, target)
    return mesh, {
        "algorithm": "spateo_marching_cube_mesh_core",
        "source_points": int(len(target.indices)),
        "levelset": float(args.levelset),
        "mc_scale_factor": float(args.mc_scale_factor),
        "dist_sample_num": args.dist_sample_num,
        "random_seed": int(args.random_seed),
        "derived_target_seed": target_seed,
        "smooth_filter": "laplacian",
        "smooth_iterations": int(args.smooth),
        "construct_surface_postprocessing": "bypassed_due_to_mixed_connectivity_meshfix_failure",
        "repair": repair,
    }


def _inside_fraction(
    points: np.ndarray, mesh: pv.PolyData, limit: int
) -> dict[str, Any]:
    if limit <= 0:
        return {"value": None, "n_points": 0, "sampling": "disabled"}
    if len(points) > limit:
        indices = np.linspace(0, len(points) - 1, limit, dtype=np.int64)
        points = points[indices]
        sampling = "even_index_sample"
    else:
        sampling = "all"
    try:
        selected = pv.PolyData(points).select_enclosed_points(
            surface=mesh, tolerance=0.0, check_surface=False
        )
        flags = np.asarray(selected.point_data["SelectedPoints"]).astype(bool)
        value: Optional[float] = float(flags.mean())
        error = None
    except Exception as exc:  # pragma: no cover - VTK-specific failure detail
        value = None
        error = f"{type(exc).__name__}: {exc}"
    return {
        "value": value,
        "n_points": int(len(points)),
        "sampling": sampling,
        "error": error,
    }


def _mesh_summary(
    mesh: pv.PolyData,
    source_points: np.ndarray,
    inside_check_points: int,
) -> dict[str, Any]:
    try:
        volume: Optional[float] = float(mesh.volume)
    except Exception:
        volume = None
    return {
        "dataset_type": type(mesh).__name__,
        "n_points": int(mesh.n_points),
        "n_cells": int(mesh.n_cells),
        "bounds": [float(value) for value in mesh.bounds],
        "surface_area": float(mesh.area),
        "volume": volume,
        "open_edges": int(mesh.n_open_edges),
        "cell_data": sorted(mesh.cell_data.keys()),
        "source_point_inside_fraction": _inside_fraction(
            source_points, mesh, inside_check_points
        ),
    }


def _save_and_verify(
    mesh: pv.PolyData, path: Path, binary: bool
) -> tuple[pv.PolyData, dict[str, Any]]:
    st.tdr.save_model(mesh, str(path), binary=binary)
    loaded = st.tdr.read_model(str(path))
    if not isinstance(loaded, pv.PolyData):
        raise RuntimeError(
            f"Spateo reloaded {type(loaded).__name__}, expected pyvista.PolyData"
        )
    if loaded.n_points != mesh.n_points or loaded.n_cells != mesh.n_cells:
        raise RuntimeError(
            f"VTK round-trip changed topology for {path.name}: "
            f"points {mesh.n_points}->{loaded.n_points}, cells {mesh.n_cells}->{loaded.n_cells}"
        )
    if not np.allclose(mesh.points, loaded.points, rtol=1e-8, atol=1e-8):
        raise RuntimeError(f"VTK round-trip changed coordinates for {path.name}")
    if not np.array_equal(np.asarray(mesh.faces), np.asarray(loaded.faces)):
        raise RuntimeError(f"VTK round-trip changed face connectivity for {path.name}")
    required = {"mesh_label", "mesh_label_rgba"}
    if not required.issubset(loaded.cell_data.keys()):
        raise RuntimeError(
            f"VTK round-trip lost required cell arrays for {path.name}: "
            f"{sorted(required - set(loaded.cell_data.keys()))}"
        )
    for key in sorted(required):
        expected = np.asarray(mesh.cell_data[key])
        observed = np.asarray(loaded.cell_data[key])
        if expected.shape != observed.shape:
            raise RuntimeError(
                f"VTK round-trip changed cell array shape for {path.name}:{key}"
            )
        if np.issubdtype(expected.dtype, np.number) and np.issubdtype(
            observed.dtype, np.number
        ):
            equal = np.allclose(expected, observed, rtol=1e-7, atol=1e-8)
        else:
            equal = np.array_equal(expected.astype(str), observed.astype(str))
        if not equal:
            raise RuntimeError(
                f"VTK round-trip changed cell array values for {path.name}:{key}"
            )
    return loaded, {
        "reader": "spateo.tdr.read_model",
        "dataset_type": type(loaded).__name__,
        "n_points_equal": True,
        "n_cells_equal": True,
        "coordinates_allclose": True,
        "face_connectivity_equal": True,
        "required_cell_data_present": True,
        "required_cell_data_values_equal": True,
        "open_edges": int(loaded.n_open_edges),
        "sha256": _sha256(path),
        "binary": bool(binary),
    }


def _render_preview(
    point_cloud: pv.PolyData,
    meshes: list[tuple[Target, pv.PolyData]],
    output: Path,
    label_key: Optional[str],
    max_points: int,
    point_size: float,
    point_opacity: float,
) -> dict[str, Any]:
    if max_points < 0:
        raise ValueError("--preview-max-points cannot be negative")
    if point_size <= 0:
        raise ValueError("--point-size must be positive")
    if not 0 <= point_opacity <= 1:
        raise ValueError("--point-opacity must be between 0 and 1")
    if max_points and point_cloud.n_points > max_points:
        sample_indices = np.linspace(
            0, point_cloud.n_points - 1, max_points, dtype=np.int64
        )
        preview_pc = pv.PolyData(np.asarray(point_cloud.points)[sample_indices])
        for key in point_cloud.point_data.keys():
            preview_pc.point_data[key] = np.asarray(point_cloud.point_data[key])[
                sample_indices
            ]
        sampling = "even_index_sample"
    else:
        preview_pc = point_cloud
        sample_indices = np.arange(point_cloud.n_points)
        sampling = "all" if max_points else "disabled"

    plotter = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(1600, 1200))
    plotter.set_background("#090b0f")
    views = [
        (0, 0, "Isometric", "iso"),
        (0, 1, "XY", "xy"),
        (1, 0, "XZ", "xz"),
        (1, 1, "YZ", "yz"),
    ]
    legend = [[target.name, target.color] for target, _ in meshes]
    for row, column, title, view in views:
        plotter.subplot(row, column)
        for target, mesh in meshes:
            plotter.add_mesh(
                mesh,
                color=target.color,
                opacity=target.opacity,
                smooth_shading=True,
                show_edges=False,
                ambient=0.25,
                diffuse=0.75,
                specular=0.35,
                specular_power=18,
                show_scalar_bar=False,
            )
        if max_points:
            point_kwargs: dict[str, Any] = {
                "style": "points_gaussian",
                "point_size": point_size,
                "opacity": point_opacity,
                "show_scalar_bar": False,
                "render_points_as_spheres": False,
            }
            rgba_key = f"{label_key}_rgba" if label_key else None
            if rgba_key and rgba_key in preview_pc.point_data:
                point_kwargs.update({"scalars": rgba_key, "rgba": True})
            else:
                point_kwargs["color"] = "#f5f5f5"
            plotter.add_points(preview_pc, **point_kwargs)
        plotter.add_text(title, position="upper_left", font_size=11, color="white")
        if view == "iso":
            plotter.view_isometric()
        elif view == "xy":
            plotter.view_xy()
        elif view == "xz":
            plotter.view_xz()
        else:
            plotter.view_yz()
        plotter.reset_camera()
        if view == "iso":
            plotter.add_legend(
                legend,
                loc="upper right",
                bcolor="#171b22",
                border=True,
                face="rectangle",
            )
    plotter.screenshot(str(output), transparent_background=False)
    plotter.close()
    return {
        "path": str(output),
        "sha256": _sha256(output),
        "views": ["isometric", "xy", "xz", "yz"],
        "point_sampling": sampling,
        "n_preview_points": int(len(sample_indices)) if max_points else 0,
        "mesh_names": [target.name for target, _ in meshes],
    }


def _validate_arguments(args: argparse.Namespace) -> None:
    if not SAFE_NAME.fullmatch(args.name):
        raise ValueError("--name must be a safe filename stem")
    if not SAFE_NAME.fullmatch(args.body_name):
        raise ValueError("--body-name must be a safe filename stem")
    if args.smooth < 0:
        raise ValueError("--smooth cannot be negative")
    if args.inside_check_points < 0:
        raise ValueError("--inside-check-points cannot be negative")
    if args.repair_max_degenerate_cells < 0:
        raise ValueError("--repair-max-degenerate-cells cannot be negative")
    if not 0 <= args.minimum_final_inside_fraction <= 1:
        raise ValueError("--minimum-final-inside-fraction must be in [0, 1]")
    if args.method == "density":
        if not 0 <= args.decimate < 1:
            raise ValueError("--decimate must be in [0, 1)")
        if not 0 < args.taubin_pass_band <= 2:
            raise ValueError("--taubin-pass-band must be in (0, 2]")
        if not 0 <= args.min_component_fraction <= 1:
            raise ValueError("--min-component-fraction must be in [0, 1]")
        if args.min_component_voxels < 1:
            raise ValueError("--min-component-voxels must be positive")
        if args.max_voxels < 1:
            raise ValueError("--max-voxels must be positive")
    else:
        if not 0 < args.levelset < 1:
            raise ValueError("--levelset must be between 0 and 1")
        if not np.isfinite(args.mc_scale_factor) or args.mc_scale_factor <= 0:
            raise ValueError("--mc-scale-factor must be positive and finite")
        if args.dist_sample_num is not None and args.dist_sample_num < 2:
            raise ValueError("--dist-sample-num must be at least 2")
        if not 0 <= args.random_seed <= np.iinfo(np.uint32).max:
            raise ValueError("--random-seed must be in [0, 2**32 - 1]")


def _build_meshes_in_staging(
    args: argparse.Namespace,
    point_cloud: pv.PolyData,
    coordinates: np.ndarray,
    input_summary: dict[str, Any],
    targets: list[Target],
    label_counts: dict[str, int],
    output_dir: Path,
    staging_dir: Path,
) -> dict[str, Any]:
    from spateo.tdr.models.models_individual import mesh_methods

    runtime_mesh_methods_path = Path(mesh_methods.__file__).resolve()
    runtime_mesh_methods_sha256 = _sha256(runtime_mesh_methods_path)
    grid: Optional[dict[str, Any]] = None
    if args.method == "density":
        grid = _shared_grid(
            coordinates=coordinates,
            target_long_axis=args.target_long_axis,
            voxel_size=args.voxel_size,
            sigma=args.sigma,
            padding=args.padding,
            max_voxels=args.max_voxels,
        )

    outputs: list[dict[str, Any]] = []
    warnings: list[str] = []
    if runtime_mesh_methods_sha256 != VALIDATED_MESH_METHODS_SHA256:
        warnings.append(
            "Runtime Spateo mesh_methods.py differs from the source revision "
            "validated by this skill"
        )
    preview_meshes: list[tuple[Target, pv.PolyData]] = []
    for target in targets:
        if args.method == "density":
            assert grid is not None
            mesh, reconstruction = _density_mesh(
                coordinates=coordinates,
                target=target,
                grid=grid,
                args=args,
            )
        else:
            mesh, reconstruction = _spateo_mesh(
                coordinates=coordinates,
                target=target,
                args=args,
            )

        filename = f"{args.name}.{target.name}.vtk"
        working_vtk_path = staging_dir / filename
        final_vtk_path = output_dir / filename
        loaded, roundtrip = _save_and_verify(
            mesh=mesh, path=working_vtk_path, binary=not args.ascii
        )
        if not args.skip_repair and loaded.n_open_edges:
            repaired, post_reload_repair = _repair_open_surface(
                loaded,
                enabled=True,
                max_degenerate_cells=args.repair_max_degenerate_cells,
            )
            repaired = _label_mesh(repaired, target)
            loaded, roundtrip = _save_and_verify(
                mesh=repaired, path=working_vtk_path, binary=not args.ascii
            )
            reconstruction["post_reload_repair"] = post_reload_repair
        if not args.skip_repair and loaded.n_open_edges:
            raise RuntimeError(
                f"Repaired VTK for {target.name!r} still has "
                f"{loaded.n_open_edges} open edges after Spateo reload"
            )
        dropped_fragments = sum(
            len(section.get("dropped_degenerate_components", []))
            for section in (
                reconstruction.get("repair", {}),
                reconstruction.get("post_reload_repair", {}),
            )
        )
        if dropped_fragments:
            warnings.append(
                f"{target.name}: repair removed {dropped_fragments} explicitly "
                "recorded unrepairable degenerate fragment(s)"
            )
        for stage_name, section in (
            ("pre-save", reconstruction.get("repair", {})),
            ("post-reload", reconstruction.get("post_reload_repair", {})),
        ):
            before_components = section.get("topology_before", {}).get("components")
            after_components = section.get("topology_after", {}).get("components")
            if (
                section.get("applied")
                and before_components is not None
                and after_components is not None
                and after_components > before_components
            ):
                warnings.append(
                    f"{target.name}: {stage_name} repair split connected "
                    f"components from {before_components} to {after_components}"
                )

        model_summary = _mesh_summary(
            mesh=loaded,
            source_points=coordinates[target.indices],
            inside_check_points=args.inside_check_points,
        )
        inside_value = model_summary["source_point_inside_fraction"]["value"]
        if inside_value is None:
            warning = f"{target.name}: final source-point containment was not available"
            if args.method == "density" and args.minimum_final_inside_fraction > 0:
                raise RuntimeError(warning)
            warnings.append(warning)
        elif inside_value < args.minimum_final_inside_fraction:
            message = (
                f"{target.name}: final source-point containment {inside_value:.3f} "
                f"is below {args.minimum_final_inside_fraction:.3f}"
            )
            if args.method == "density":
                raise RuntimeError(message)
            warnings.append(message)
        if args.skip_repair and loaded.n_open_edges:
            warnings.append(
                f"{target.name}: retained {loaded.n_open_edges} open edges by request"
            )
        if (
            args.method == "density"
            and args.level is not None
            and not reconstruction["level_selection"]["target_met"]
        ):
            warnings.append(
                f"{target.name}: fixed density level did not meet the requested "
                "voxel coverage target"
            )

        preview_meshes.append((target, loaded))
        outputs.append(
            {
                "name": target.name,
                "role": target.role,
                "selected_values": list(target.values),
                "color": target.color,
                "opacity": target.opacity,
                "vtk": str(final_vtk_path),
                "vtk_sha256": roundtrip["sha256"],
                "reconstruction": reconstruction,
                "model": model_summary,
                "roundtrip": roundtrip,
            }
        )

    preview = None
    if not args.skip_preview:
        preview_filename = f"{args.name}.preview.png"
        preview = _render_preview(
            point_cloud=point_cloud,
            meshes=preview_meshes,
            output=staging_dir / preview_filename,
            label_key=args.label_key,
            max_points=args.preview_max_points,
            point_size=args.point_size,
            point_opacity=args.point_opacity,
        )
        preview["path"] = str(output_dir / preview_filename)

    # Guard again at the file-contract boundary.  On a case-insensitive file
    # system, case-only name collisions otherwise overwrite an earlier VTK.
    final_paths = [Path(row["vtk"]) for row in outputs]
    folded_paths = [str(path).casefold() for path in final_paths]
    if len(folded_paths) != len(set(folded_paths)):
        raise RuntimeError("Mesh output paths are not unique ignoring letter case")
    for row, final_path in zip(outputs, final_paths):
        staged_path = staging_dir / final_path.name
        if not staged_path.is_file() or _sha256(staged_path) != row["vtk_sha256"]:
            raise RuntimeError(
                f"Staged mesh hash no longer matches manifest: {final_path.name}"
            )

    manifest_filename = f"{args.name}.manifest.json"
    final_manifest_path = output_dir / manifest_filename
    grid_summary = None
    if grid is not None:
        grid_summary = {
            "origin": grid["origin"],
            "shape": grid["shape"],
            "voxel_size": grid["voxel_size"],
            "padding": grid["padding"],
            "sigma_voxels": grid["sigma"],
            "n_voxels": grid["n_voxels"],
        }
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "pass_with_warnings" if warnings else "pass",
        "quality": {
            "minimum_final_inside_fraction": float(args.minimum_final_inside_fraction),
            "warnings": warnings,
        },
        "validated_against_revision": VALIDATED_SOURCE_REVISION,
        "input": input_summary,
        "selection": {
            "label_key": args.label_key,
            "available_label_counts": label_counts,
            "minimum_selection_points": int(args.minimum_selection_points),
        },
        "method": args.method,
        "parameters": {
            "active_groups": [
                "common",
                args.method,
                "repair",
                "preview",
                "serialization",
            ],
            "inactive_method_group": (
                "spateo-marching-cube" if args.method == "density" else "density"
            ),
            "common": {
                "smooth": int(args.smooth),
                "inside_check_points": int(args.inside_check_points),
                "minimum_final_inside_fraction": float(
                    args.minimum_final_inside_fraction
                ),
            },
            "density": {
                "target_long_axis": int(args.target_long_axis),
                "voxel_size": args.voxel_size,
                "voxel_size_mode": (
                    "explicit" if args.voxel_size is not None else "target_long_axis"
                ),
                "sigma": float(args.sigma),
                "padding": args.padding,
                "max_voxels": int(args.max_voxels),
                "level": args.level,
                "min_level": float(args.min_level),
                "max_level": float(args.max_level),
                "level_steps": int(args.level_steps),
                "body_coverage_target": float(args.body_coverage_target),
                "selection_coverage_target": float(args.selection_coverage_target),
                "min_component_voxels": int(args.min_component_voxels),
                "min_component_fraction": float(args.min_component_fraction),
                "voxel_connectivity": 6,
                "taubin_pass_band": float(args.taubin_pass_band),
                "decimate": float(args.decimate),
            },
            "spateo-marching-cube": {
                "levelset": float(args.levelset),
                "mc_scale_factor": float(args.mc_scale_factor),
                "dist_sample_num": args.dist_sample_num,
                "random_seed": int(args.random_seed),
            },
            "repair": {
                "enabled": not args.skip_repair,
                "max_degenerate_cells": int(args.repair_max_degenerate_cells),
            },
            "preview": {
                "enabled": not args.skip_preview,
                "max_points": int(args.preview_max_points),
                "point_size": float(args.point_size),
                "point_opacity": float(args.point_opacity),
            },
            "serialization": {"binary": not args.ascii},
        },
        "shared_grid": grid_summary,
        "meshes": outputs,
        "output": {
            "directory": str(output_dir),
            "preview": preview,
            "manifest": str(final_manifest_path),
        },
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "scikit_image": __import__("skimage").__version__,
            "pyvista": pv.__version__,
            "spateo": getattr(st, "__version__", "unknown"),
            "spateo_module": str(Path(st.__file__).resolve()),
            "runtime_mesh_methods_path": str(runtime_mesh_methods_path),
            "runtime_mesh_methods_sha256": runtime_mesh_methods_sha256,
            "validated_mesh_methods_sha256": VALIDATED_MESH_METHODS_SHA256,
            "runtime_mesh_methods_matches_validated_source": (
                runtime_mesh_methods_sha256 == VALIDATED_MESH_METHODS_SHA256
            ),
        },
    }
    (staging_dir / manifest_filename).write_text(
        json.dumps(_json_value(manifest), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


def build_meshes(args: argparse.Namespace) -> dict[str, Any]:
    input_path = Path(args.input_vtk).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    _validate_arguments(args)

    point_cloud, coordinates, input_summary = _load_point_cloud(input_path)
    palette = _load_json_dict(args.palette_json, "palette")
    opacities = _load_json_dict(args.opacity_json, "opacity")
    targets, _, label_counts = _make_targets(
        args=args,
        point_cloud=point_cloud,
        palette=palette,
        opacities=opacities,
    )
    staging_dir = _prepare_output_dir(output_dir)
    try:
        manifest = _build_meshes_in_staging(
            args=args,
            point_cloud=point_cloud,
            coordinates=coordinates,
            input_summary=input_summary,
            targets=targets,
            label_counts=label_counts,
            output_dir=output_dir,
            staging_dir=staging_dir,
        )
        _commit_output_dir(staging_dir=staging_dir, output_dir=output_dir)
    except Exception:
        shutil.rmtree(staging_dir, ignore_errors=True)
        raise
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a full-body surface and selected annotation surfaces from a "
            "round-trip-validated Spateo point-cloud VTK."
        )
    )
    parser.add_argument("--input-vtk", required=True)
    parser.add_argument(
        "--output-dir",
        required=True,
        help="A new or empty directory; existing contents are never overwritten.",
    )
    parser.add_argument("--name", default="surface_meshes")
    parser.add_argument(
        "--label-key",
        help="Point-data label array used by --selection (for example lineage or tissue).",
    )
    parser.add_argument(
        "--selection",
        action="append",
        default=[],
        help="Mesh selection as NAME=VALUE or NAME=VALUE1|VALUE2; repeat as needed.",
    )
    parser.add_argument("--no-body", action="store_true")
    parser.add_argument("--body-name", default="body")
    parser.add_argument("--body-color", default="#d8d1b8")
    parser.add_argument("--body-opacity", type=float, default=0.20)
    parser.add_argument("--selection-opacity", type=float, default=0.88)
    parser.add_argument("--palette-json")
    parser.add_argument("--opacity-json")
    parser.add_argument("--minimum-selection-points", type=int, default=50)
    parser.add_argument(
        "--method",
        choices=("density", "spateo-marching-cube"),
        default="density",
    )

    density = parser.add_argument_group("density marching-cubes method")
    density.add_argument("--target-long-axis", type=int, default=160)
    density.add_argument("--voxel-size", type=float)
    density.add_argument("--sigma", type=float, default=1.35)
    density.add_argument("--padding", type=int)
    density.add_argument("--max-voxels", type=int, default=32_000_000)
    density.add_argument("--level", type=float)
    density.add_argument("--min-level", type=float, default=0.015)
    density.add_argument("--max-level", type=float, default=0.35)
    density.add_argument("--level-steps", type=int, default=28)
    density.add_argument("--body-coverage-target", type=float, default=0.985)
    density.add_argument("--selection-coverage-target", type=float, default=0.97)
    density.add_argument("--min-component-voxels", type=int, default=12)
    density.add_argument("--min-component-fraction", type=float, default=0.005)
    density.add_argument("--taubin-pass-band", type=float, default=0.08)
    density.add_argument(
        "--decimate",
        type=float,
        default=0.0,
        help="Topology-preserving target reduction fraction in [0, 1).",
    )

    native = parser.add_argument_group("Spateo native marching-cube method")
    native.add_argument("--levelset", type=float, default=0.5)
    native.add_argument("--mc-scale-factor", type=float, default=0.8)
    native.add_argument("--dist-sample-num", type=int, default=100)
    native.add_argument(
        "--random-seed",
        type=int,
        default=0,
        help="Base seed for deterministic per-target native distance sampling.",
    )

    parser.add_argument(
        "--smooth",
        type=int,
        default=40,
        help=(
            "Taubin iterations for density mode or Spateo Laplacian iterations for "
            "native mode."
        ),
    )
    parser.add_argument("--inside-check-points", type=int, default=20000)
    parser.add_argument(
        "--minimum-final-inside-fraction",
        type=float,
        default=0.80,
        help=(
            "Minimum final source-point containment. Density mode fails below it; "
            "native comparison mode records pass_with_warnings."
        ),
    )
    parser.add_argument(
        "--skip-repair",
        action="store_true",
        help="Retain open edges instead of repairing a face-only mesh with MeshFix.",
    )
    parser.add_argument(
        "--repair-max-degenerate-cells",
        type=int,
        default=16,
        help=(
            "Allow repair to drop only an unrepairable connected fragment at or "
            "below this triangle count; every drop is recorded."
        ),
    )
    parser.add_argument("--skip-preview", action="store_true")
    parser.add_argument("--preview-max-points", type=int, default=30000)
    parser.add_argument("--point-size", type=float, default=2.0)
    parser.add_argument("--point-opacity", type=float, default=0.28)
    parser.add_argument("--ascii", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = build_meshes(args)
    print(
        json.dumps(
            {
                "status": manifest["status"],
                "method": manifest["method"],
                "meshes": [row["vtk"] for row in manifest["meshes"]],
                "preview": (
                    manifest["output"]["preview"]["path"]
                    if manifest["output"]["preview"]
                    else None
                ),
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
