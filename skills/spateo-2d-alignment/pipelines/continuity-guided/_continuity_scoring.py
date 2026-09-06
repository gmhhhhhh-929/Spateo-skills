"""Private dependency closure for the packaged continuity-guided pipeline."""
from __future__ import annotations


import argparse


import json


import math


import sys


import traceback


from dataclasses import dataclass


from datetime import datetime


from pathlib import Path


import matplotlib


matplotlib.use("Agg")


import matplotlib.pyplot as plt


import numpy as np


from scipy.optimize import least_squares


from scipy.spatial import cKDTree


import _core as core


def annotation_shape_score(
    moving: np.ndarray,
    moving_labels: np.ndarray,
    fixed: np.ndarray,
    fixed_labels: np.ndarray,
    scale: float,
) -> tuple[float, int]:
    """Score within-lineage outlines, not only lineage centroids."""
    records: list[float] = []
    weights: list[float] = []
    shared = sorted(set(moving_labels.astype(str)).intersection(fixed_labels.astype(str)))
    for label in shared:
        moving_group = moving[moving_labels.astype(str) == label]
        fixed_group = fixed[fixed_labels.astype(str) == label]
        if min(len(moving_group), len(fixed_group)) < 8:
            continue
        forward = cKDTree(fixed_group).query(moving_group, k=1, workers=1)[0]
        reverse = cKDTree(moving_group).query(fixed_group, k=1, workers=1)[0]
        distance = 0.5 * (core.trimmed_mean(forward, 0.85) + core.trimmed_mean(reverse, 0.85))
        moving_center = np.median(moving_group, axis=0)
        fixed_center = np.median(fixed_group, axis=0)
        spread = 0.5 * (
            np.median(np.linalg.norm(moving_group - moving_center, axis=1))
            + np.median(np.linalg.norm(fixed_group - fixed_center, axis=1))
        )
        concentration = 1.0 / (0.10 + float(spread) / scale)
        support = math.sqrt(min(len(moving_group), len(fixed_group)))
        records.append(float(distance / scale))
        # The cap prevents one abundant lineage from determining the pose.
        weights.append(min(18.0, support * concentration))
    if not records:
        return math.inf, 0
    return float(np.average(records, weights=weights)), len(records)


def longitudinal_endpoints(
    points: np.ndarray, reference: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return robust low/high endpoints along the reference longitudinal axis."""
    centered = reference - np.median(reference, axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    axis = vt[0]
    projection = points @ axis
    low = points[projection <= np.quantile(projection, 0.06)]
    high = points[projection >= np.quantile(projection, 0.94)]
    return np.median(low, axis=0), np.median(high, axis=0)


def balanced_endpoint_offset(moving: np.ndarray, fixed: np.ndarray) -> np.ndarray:
    moving_low, moving_high = longitudinal_endpoints(moving, fixed)
    fixed_low, fixed_high = longitudinal_endpoints(fixed, fixed)
    return 0.5 * ((fixed_low - moving_low) + (fixed_high - moving_high))


def longitudinal_endpoint_score(
    moving: np.ndarray, fixed: np.ndarray, scale: float
) -> tuple[float, float, float]:
    moving_low, moving_high = longitudinal_endpoints(moving, fixed)
    fixed_low, fixed_high = longitudinal_endpoints(fixed, fixed)
    low = float(np.linalg.norm(moving_low - fixed_low) / scale)
    high = float(np.linalg.norm(moving_high - fixed_high) / scale)
    return low, high, max(low, high)


def score_alignment_v2(
    moving: np.ndarray,
    moving_labels: np.ndarray,
    fixed: np.ndarray,
    fixed_labels: np.ndarray,
) -> dict[str, float]:
    base = core.score_alignment(moving, moving_labels, fixed, fixed_labels)
    semantic_shape, semantic_labels = annotation_shape_score(
        moving, moving_labels, fixed, fixed_labels, base["scale"]
    )
    if not np.isfinite(semantic_shape):
        semantic_shape = base["annotation"] * 1.25
    endpoint_low, endpoint_high, endpoint_max = longitudinal_endpoint_score(
        moving, fixed, base["scale"]
    )
    total = (
        0.22 * base["geometry"]
        + 0.18 * base["annotation"]
        + 0.18 * semantic_shape
        + 0.34 * endpoint_max
        + 0.08 * base["extent"]
    )
    return {
        **base,
        "total": float(total),
        "semantic_shape": float(semantic_shape),
        "semantic_labels": float(semantic_labels),
        "endpoint_low": endpoint_low,
        "endpoint_high": endpoint_high,
        "endpoint_max": endpoint_max,
    }
