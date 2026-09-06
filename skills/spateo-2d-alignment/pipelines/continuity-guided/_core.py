"""Private dependency closure for the packaged continuity-guided pipeline."""
from __future__ import annotations


import argparse


import csv


import gzip


import hashlib


import json


import math


import re


import sys


import traceback


from dataclasses import dataclass


from datetime import datetime


from pathlib import Path


from typing import Iterable


import anndata as ad


import matplotlib


matplotlib.use("Agg")


import matplotlib.pyplot as plt


import numpy as np


from scipy import sparse


from scipy.spatial import cKDTree


@dataclass
class SliceData:
    path: Path
    slice_id: int
    cell_ids: np.ndarray
    annotations: np.ndarray
    xy: np.ndarray
    sample_idx: np.ndarray
    content_sha256: str

    @property
    def sample_xy(self) -> np.ndarray:
        return self.xy[self.sample_idx]

    @property
    def sample_annotations(self) -> np.ndarray:
        return self.annotations[self.sample_idx]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_arrays(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.dtype).encode())
        digest.update(str(contiguous.shape).encode())
        digest.update(contiguous.view(np.uint8))
    return digest.hexdigest()


def slice_number(path: Path) -> int:
    match = re.search(r"SL(\d+)", path.name, flags=re.IGNORECASE)
    if match is None:
        raise ValueError(f"Cannot parse slice number from {path.name}")
    return int(match.group(1))


def deterministic_sample(labels: np.ndarray, cap: int, seed: int) -> np.ndarray:
    count = len(labels)
    if count <= cap:
        return np.arange(count, dtype=int)
    rng = np.random.default_rng(seed)
    selected: set[int] = set()
    unique, counts = np.unique(labels, return_counts=True)
    minimum_per_label = max(4, min(18, cap // max(1, 3 * len(unique))))
    for label, label_count in zip(unique, counts):
        indices = np.flatnonzero(labels == label)
        take = min(label_count, minimum_per_label)
        chosen = rng.choice(indices, size=take, replace=False)
        selected.update(int(index) for index in chosen)
    remaining = cap - len(selected)
    if remaining > 0:
        pool = np.array(sorted(set(range(count)) - selected), dtype=int)
        chosen = rng.choice(pool, size=min(remaining, len(pool)), replace=False)
        selected.update(int(index) for index in chosen)
    return np.array(sorted(selected), dtype=int)


def load_slices(args: argparse.Namespace) -> list[SliceData]:
    paths = sorted(args.slice_dir.glob("*.h5ad"), key=slice_number)
    if len(paths) < 2:
        raise RuntimeError(f"Need at least two slice H5AD files in {args.slice_dir}")
    slices: list[SliceData] = []
    all_ids: set[str] = set()
    for index, path in enumerate(paths):
        model = ad.read_h5ad(path, backed="r")
        try:
            if "spatial_3d" in model.obsm:
                raise RuntimeError(f"Blind input contains forbidden obsm['spatial_3d']: {path}")
            if args.spatial_key not in model.obsm:
                raise KeyError(f"Missing obsm[{args.spatial_key!r}] in {path}")
            if args.annotation_key not in model.obs:
                raise KeyError(f"Missing obs[{args.annotation_key!r}] in {path}")
            if not model.obs_names.is_unique:
                raise RuntimeError(f"Non-unique obs_names in {path}")
            cell_ids = model.obs_names.astype(str).to_numpy(copy=True)
            annotations = model.obs[args.annotation_key].astype(str).to_numpy(copy=True)
            xy = np.asarray(model.obsm[args.spatial_key], dtype=np.float64).copy()
        finally:
            model.file.close()
        if xy.ndim != 2 or xy.shape[1] < 2 or len(xy) != len(cell_ids):
            raise RuntimeError(f"Invalid XY array in {path}: {xy.shape}")
        xy = xy[:, :2]
        if not np.isfinite(xy).all():
            raise RuntimeError(f"Non-finite XY coordinate in {path}")
        overlap = all_ids.intersection(cell_ids.tolist())
        if overlap:
            raise RuntimeError(f"Cell ids occur in multiple slices; example={next(iter(overlap))}")
        all_ids.update(cell_ids.tolist())
        sample_idx = deterministic_sample(annotations, args.sample_cap, args.seed + index)
        content_hash = sha256_arrays(
            np.char.encode(cell_ids.astype("U"), "utf-8"),
            np.char.encode(annotations.astype("U"), "utf-8"),
            xy.astype("<f8"),
        )
        slices.append(
            SliceData(
                path=path,
                slice_id=slice_number(path),
                cell_ids=cell_ids,
                annotations=annotations,
                xy=xy,
                sample_idx=sample_idx,
                content_sha256=content_hash,
            )
        )
    return slices


def apply_matrix(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    return points @ matrix[:2, :2].T + matrix[:2, 2]


def rigid_matrix(angle_deg: float, translation: np.ndarray) -> np.ndarray:
    theta = math.radians(angle_deg)
    cosine, sine = math.cos(theta), math.sin(theta)
    matrix = np.eye(3, dtype=np.float64)
    matrix[:2, :2] = np.array([[cosine, -sine], [sine, cosine]])
    matrix[:2, 2] = np.asarray(translation, dtype=np.float64)
    return matrix


def robust_diagonal(points: np.ndarray) -> float:
    lower, upper = np.quantile(points, [0.02, 0.98], axis=0)
    return max(float(np.linalg.norm(upper - lower)), 1e-8)


def trimmed_mean(values: np.ndarray, quantile: float = 0.90) -> float:
    if len(values) == 0:
        return math.inf
    cutoff = np.quantile(values, quantile)
    retained = values[values <= cutoff]
    return float(np.mean(retained)) if len(retained) else float(cutoff)


def label_stats(points: np.ndarray, labels: np.ndarray, minimum: int = 8) -> dict[str, dict[str, object]]:
    stats: dict[str, dict[str, object]] = {}
    for label in np.unique(labels):
        mask = labels == label
        count = int(mask.sum())
        if count < minimum:
            continue
        subset = points[mask]
        center = np.median(subset, axis=0)
        spread = float(np.median(np.linalg.norm(subset - center, axis=1)))
        stats[str(label)] = {"count": count, "center": center, "spread": spread}
    return stats


def endpoint_landmarks(points: np.ndarray, reference: np.ndarray | None = None) -> list[np.ndarray]:
    """Return robust centers of the four ends of an oriented tissue outline."""
    orientation_source = points if reference is None else reference
    centered = orientation_source - np.median(orientation_source, axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    axes = [vt[0], -vt[0], vt[1], -vt[1]]
    landmarks: list[np.ndarray] = []
    for axis in axes:
        projections = points @ axis
        cutoff = np.quantile(projections, 0.94)
        end_points = points[projections >= cutoff]
        landmarks.append(np.median(end_points, axis=0))
    return landmarks


def score_alignment(
    moving: np.ndarray,
    moving_labels: np.ndarray,
    fixed: np.ndarray,
    fixed_labels: np.ndarray,
) -> dict[str, float]:
    scale = 0.5 * (robust_diagonal(moving) + robust_diagonal(fixed))
    fixed_tree = cKDTree(fixed)
    moving_tree = cKDTree(moving)
    m_to_f = fixed_tree.query(moving, k=1, workers=1)[0]
    f_to_m = moving_tree.query(fixed, k=1, workers=1)[0]
    geometry = 0.5 * (trimmed_mean(m_to_f) + trimmed_mean(f_to_m)) / scale

    moving_stats = label_stats(moving, moving_labels)
    fixed_stats = label_stats(fixed, fixed_labels)
    shared = sorted(set(moving_stats).intersection(fixed_stats))
    annotation_distances: list[float] = []
    annotation_weights: list[float] = []
    for label in shared:
        mstat, fstat = moving_stats[label], fixed_stats[label]
        distance = float(np.linalg.norm(np.asarray(mstat["center"]) - np.asarray(fstat["center"]))) / scale
        concentrated = 1.0 / (
            0.08 + (float(mstat["spread"]) + float(fstat["spread"])) / (2.0 * scale)
        )
        support = math.sqrt(min(int(mstat["count"]), int(fstat["count"])))
        annotation_distances.append(distance)
        annotation_weights.append(min(30.0, support * concentrated))
    if annotation_distances:
        annotation = float(np.average(annotation_distances, weights=annotation_weights))
    else:
        annotation = geometry * 1.5

    moving_ends = endpoint_landmarks(moving, reference=fixed)
    fixed_ends = endpoint_landmarks(fixed, reference=fixed)
    endpoint = float(
        np.mean(
            [np.linalg.norm(moving_end - fixed_end) for moving_end, fixed_end in zip(moving_ends, fixed_ends)]
        )
        / scale
    )

    moving_lower, moving_upper = np.quantile(moving, [0.02, 0.98], axis=0)
    fixed_lower, fixed_upper = np.quantile(fixed, [0.02, 0.98], axis=0)
    extent = float(np.linalg.norm((moving_upper - moving_lower) - (fixed_upper - fixed_lower))) / scale
    total = 0.42 * geometry + 0.30 * annotation + 0.18 * endpoint + 0.10 * extent
    return {
        "total": total,
        "geometry": geometry,
        "annotation": annotation,
        "endpoint": endpoint,
        "extent": extent,
        "shared_annotations": float(len(shared)),
        "scale": scale,
    }


def fit_proper_rigid(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    if source.shape != target.shape or source.ndim != 2 or source.shape[1] != 2:
        raise ValueError("Rigid fit requires paired N x 2 arrays")
    source_mean = source.mean(axis=0)
    target_mean = target.mean(axis=0)
    source_centered = source - source_mean
    target_centered = target - target_mean
    u, _, vt = np.linalg.svd(source_centered.T @ target_centered)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = vt.T @ u.T
    translation = target_mean - rotation @ source_mean
    matrix = np.eye(3)
    matrix[:2, :2] = rotation
    matrix[:2, 2] = translation
    return matrix


def matrix_rotation_deg(matrix: np.ndarray) -> float:
    return math.degrees(math.atan2(matrix[1, 0], matrix[0, 0]))


def aligned_sample(slice_data: SliceData, transform: np.ndarray) -> np.ndarray:
    return apply_matrix(slice_data.sample_xy, transform)


def compute_edge_scores(slices: list[SliceData], transforms: list[np.ndarray]) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for index in range(len(slices) - 1):
        left = aligned_sample(slices[index], transforms[index])
        right = aligned_sample(slices[index + 1], transforms[index + 1])
        score = score_alignment(
            right,
            slices[index + 1].sample_annotations,
            left,
            slices[index].sample_annotations,
        )
        rows.append({"boundary_index": float(index), **score})
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_coordinates(
    path: Path, slices: list[SliceData], transforms: list[np.ndarray]
) -> int:
    count = 0
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cell_id", "slice_id", "anno", "aligned_x", "aligned_y"])
        for slice_data, transform in zip(slices, transforms):
            aligned = apply_matrix(slice_data.xy, transform)
            for cell_id, annotation, xy in zip(slice_data.cell_ids, slice_data.annotations, aligned):
                writer.writerow(
                    [str(cell_id), slice_data.slice_id, str(annotation), f"{xy[0]:.10f}", f"{xy[1]:.10f}"]
                )
                count += 1
    return count


def plot_overview(
    path: Path,
    slices: list[SliceData],
    transforms: list[np.ndarray],
    seed: int,
) -> None:
    rng = np.random.default_rng(seed)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), dpi=180)
    color_values = plt.cm.turbo(np.linspace(0.02, 0.98, len(slices)))
    for index, (slice_data, transform) in enumerate(zip(slices, transforms)):
        count = min(1800, len(slice_data.xy))
        selected = rng.choice(len(slice_data.xy), size=count, replace=False)
        axes[0].scatter(
            slice_data.xy[selected, 0],
            slice_data.xy[selected, 1],
            s=0.55,
            alpha=0.45,
            linewidths=0,
            color=color_values[index],
        )
        aligned = apply_matrix(slice_data.xy[selected], transform)
        axes[1].scatter(
            aligned[:, 0],
            aligned[:, 1],
            s=0.55,
            alpha=0.45,
            linewidths=0,
            color=color_values[index],
        )
    axes[0].set_title("Prescribed-unaligned input")
    axes[1].set_title("Continuity-first frozen candidate")
    for axis in axes:
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_frame_on(False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_edge_scores(
    path: Path,
    slices: list[SliceData],
    before: list[dict[str, float]],
    after: list[dict[str, float]],
) -> None:
    labels = [f"{slices[i].slice_id}/{slices[i + 1].slice_id}" for i in range(len(slices) - 1)]
    x = np.arange(len(labels))
    fig, axis = plt.subplots(figsize=(max(7.0, len(labels) * 0.34), 3.2), dpi=180)
    axis.plot(x, [row["total"] for row in before], marker="o", ms=3, lw=1, label="before block rescue")
    axis.plot(x, [row["total"] for row in after], marker="o", ms=3, lw=1, label="frozen")
    axis.set_xticks(x)
    axis.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    axis.set_ylabel("Blind pair score (lower is better)")
    axis.legend(frameon=False, fontsize=7)
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
