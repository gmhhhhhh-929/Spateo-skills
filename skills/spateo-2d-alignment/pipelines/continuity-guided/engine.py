#!/usr/bin/env python3
"""Blind Spateo-first serial alignment with autonomous continuity repair."""

from __future__ import annotations

import argparse
import json
import math
import sys
import traceback
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import numpy as np

import _core as core


PIPELINE_VERSION = "spateo-continuity-v3.7.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True)
    parser.add_argument("--slice-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--annotation-key", default="anno")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--representation-key", default="X_pca")
    parser.add_argument("--sample-cap", type=int, default=1800)
    parser.add_argument("--seed", type=int, default=20260817)
    parser.add_argument("--max-iter", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=800)
    parser.add_argument("--device", default="0")
    parser.add_argument("--priority-annotation", default="ROD")
    parser.add_argument("--terminal-passes", type=int, default=2)
    parser.add_argument("--internal-block-passes", type=int, default=2)
    parser.add_argument("--internal-block-min-side-slices", type=int, default=3)
    parser.add_argument("--internal-block-score-ratio", type=float, default=1.75)
    parser.add_argument("--internal-block-mad-multiplier", type=float, default=3.0)
    parser.add_argument("--internal-block-endpoint-ratio", type=float, default=1.50)
    parser.add_argument("--internal-block-endpoint-floor", type=float, default=0.06)
    parser.add_argument("--internal-block-min-improvement", type=float, default=0.05)
    parser.add_argument("--internal-block-min-endpoint-mean-improvement", type=float, default=0.25)
    parser.add_argument(
        "--repair-policy",
        choices=["auto", "learned_probe"],
        default="auto",
        help="learned_probe freezes the previously successful three-step recipe for validation",
    )
    parser.add_argument(
        "--postprocess",
        choices=["none", "auto"],
        default="auto",
        help="apply blind interface and annotation-continuity repairs after the Spateo baseline",
    )
    return parser.parse_args()


def load_models(
    slices: list[core.SliceData], args: argparse.Namespace
) -> list[ad.AnnData]:
    models: list[ad.AnnData] = []
    annotation_vectors: dict[str, np.ndarray] = {}
    for slice_data in slices:
        model = ad.read_h5ad(slice_data.path)
        if "spatial_3d" in model.obsm:
            raise RuntimeError(f"Blind input contains forbidden obsm['spatial_3d']: {slice_data.path}")
        if args.representation_key not in model.obsm:
            raise KeyError(
                f"Missing obsm[{args.representation_key!r}] required for Spateo: {slice_data.path}"
            )
        representation = np.asarray(model.obsm[args.representation_key], dtype=np.float32)
        if representation.ndim != 2 or representation.shape[0] != model.n_obs:
            raise RuntimeError(f"Invalid representation in {slice_data.path}: {representation.shape}")
        if not np.isfinite(representation).all():
            raise RuntimeError(f"Non-finite representation in {slice_data.path}")
        labels = model.obs[args.annotation_key].astype(str).to_numpy()
        for label in np.unique(labels):
            group = representation[labels == label]
            representative = group[0]
            if np.max(np.abs(group - representative)) > 1e-6:
                raise RuntimeError(
                    f"obsm[{args.representation_key!r}] is not a deterministic annotation encoding "
                    f"for {label!r} in {slice_data.path}"
                )
            if np.linalg.norm(representative) <= 1e-8:
                raise RuntimeError(f"Zero annotation vector for {label!r} in {slice_data.path}")
            if label in annotation_vectors and not np.allclose(
                annotation_vectors[label], representative, atol=1e-6, rtol=0.0
            ):
                raise RuntimeError(f"Inconsistent annotation vector for {label!r} across slices")
            annotation_vectors[label] = representative.copy()
        model.obsm[args.representation_key] = representation
        model.obsm[args.spatial_key] = np.asarray(
            model.obsm[args.spatial_key], dtype=np.float32
        )[:, :2]
        models.append(model)
    labels = sorted(annotation_vectors)
    for index, label in enumerate(labels):
        for other in labels[index + 1 :]:
            if np.allclose(annotation_vectors[label], annotation_vectors[other], atol=1e-6, rtol=0.0):
                raise RuntimeError(
                    f"Annotations {label!r} and {other!r} share the same representation vector"
                )
    return models


def run_spateo_serial(
    models: list[ad.AnnData], args: argparse.Namespace
) -> tuple[list[ad.AnnData], list[np.ndarray]]:
    import torch
    import spateo as st

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
    aligned, assignments = st.align.morpho_align(
        models=models,
        rep_layer=args.representation_key,
        rep_field="obsm",
        spatial_key=args.spatial_key,
        key_added="align_spatial",
        iter_key_added=None,
        vecfld_key_added="VecFld_morpho",
        mode="SN-S",
        dissimilarity="cos",
        max_iter=args.max_iter,
        dtype="float32",
        device=args.device,
        verbose=False,
        batch_size=args.batch_size,
        beta=0.05,
        use_hvg=False,
        SVI_mode=True,
        pre_compute_dist=True,
        nn_init=True,
        init_transform=True,
        init_layer=args.representation_key,
        init_field="obsm",
        allow_flip=False,
        nonrigid_start_iter=80,
    )
    return aligned, assignments


def fractional_rigid(matrix: np.ndarray, fraction: float) -> np.ndarray:
    angle = fraction * core.matrix_rotation_deg(matrix)
    return core.rigid_matrix(angle, fraction * matrix[:2, 2])


def make_work_model(
    model: ad.AnnData,
    coordinates: np.ndarray,
    args: argparse.Namespace,
    mask: np.ndarray | None = None,
) -> ad.AnnData:
    if mask is None:
        work = model.copy()
        work.obsm[args.spatial_key] = np.asarray(coordinates, dtype=np.float32)
    else:
        work = model[mask].copy()
        work.obsm[args.spatial_key] = np.asarray(coordinates[mask], dtype=np.float32)
    return work


def run_spateo_candidate_models(
    models: list[ad.AnnData], args: argparse.Namespace
) -> list[ad.AnnData]:
    import torch
    import spateo as st

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
    aligned, _ = st.align.morpho_align(
        models=models,
        rep_layer=args.representation_key,
        rep_field="obsm",
        spatial_key=args.spatial_key,
        key_added="align_spatial",
        iter_key_added=None,
        vecfld_key_added="VecFld_candidate",
        mode="SN-S",
        dissimilarity="cos",
        max_iter=args.max_iter,
        dtype="float32",
        device=args.device,
        verbose=False,
        batch_size=args.batch_size,
        beta=0.05,
        use_hvg=False,
        SVI_mode=True,
        pre_compute_dist=True,
        # Candidate coordinates already live in the baseline aligned frame.
        # Re-running expression NN initialization is both unnecessary and
        # undefined for a single-lineage subset whose cosine distances are all
        # zero.
        nn_init=False,
        init_transform=False,
        init_layer=args.representation_key,
        init_field="obsm",
        allow_flip=False,
        nonrigid_start_iter=80,
    )
    return aligned


def aligned_points(slices: list[core.SliceData], transforms: list[np.ndarray], index: int) -> np.ndarray:
    return core.apply_matrix(slices[index].xy, transforms[index])


def sampled_pair_score(
    slices: list[core.SliceData], transforms: list[np.ndarray], left: int, right: int
) -> dict[str, float]:
    import _continuity_scoring as pose

    return pose.score_alignment_v2(
        core.aligned_sample(slices[right], transforms[right]),
        slices[right].sample_annotations,
        core.aligned_sample(slices[left], transforms[left]),
        slices[left].sample_annotations,
    )


def local_objective(
    slices: list[core.SliceData], transforms: list[np.ndarray], boundaries: list[int]
) -> float:
    scores = [sampled_pair_score(slices, transforms, index, index + 1) for index in boundaries]
    return float(np.mean([score["total"] for score in scores]))


def reverse_boundary_fractional_rescue(
    slices: list[core.SliceData],
    models: list[ad.AnnData],
    transforms: list[np.ndarray],
    boundary: int,
    args: argparse.Namespace,
) -> tuple[list[np.ndarray], dict[str, object]]:
    fixed = boundary + 1
    moving = boundary
    fixed_points = aligned_points(slices, transforms, fixed)
    moving_points = aligned_points(slices, transforms, moving)
    work = [
        make_work_model(models[fixed], fixed_points, args),
        make_work_model(models[moving], moving_points, args),
    ]
    aligned = run_spateo_candidate_models(work, args)
    delta = core.fit_proper_rigid(
        moving_points, np.asarray(aligned[1].obsm["align_spatial"], dtype=np.float64)
    )
    boundaries = sorted(set([max(0, boundary - 1), boundary]))
    baseline_objective = local_objective(slices, transforms, boundaries)
    baseline_interface = sampled_pair_score(slices, transforms, boundary, fixed)
    candidates: list[tuple[float, float, list[np.ndarray], dict[str, float], np.ndarray]] = []
    for fraction in [0.0, 0.25, 0.50, 0.75, 1.0]:
        trial = [matrix.copy() for matrix in transforms]
        fraction_delta = fractional_rigid(delta, fraction)
        trial[moving] = fraction_delta @ trial[moving]
        interface = sampled_pair_score(slices, trial, boundary, fixed)
        objective = local_objective(slices, trial, boundaries)
        candidates.append((objective, fraction, trial, interface, fraction_delta))
    candidates.sort(key=lambda item: item[0])
    if args.repair_policy == "learned_probe":
        objective, fraction, trial, interface, selected_delta = next(
            item for item in candidates if item[1] == 0.25
        )
    else:
        eligible = []
        for item in candidates:
            item_objective, item_fraction, _, item_interface, _ = item
            item_improvement = (baseline_objective - item_objective) / max(
                baseline_objective, 1e-12
            )
            item_endpoint_safe = (
                item_interface["endpoint_low"] <= baseline_interface["endpoint_low"] * 1.08
                and item_interface["endpoint_high"] <= baseline_interface["endpoint_high"] * 1.08
            )
            item_annotation_safe = (
                item_interface["annotation"] <= baseline_interface["annotation"] * 1.15
            )
            if (
                item_fraction > 0
                and item_improvement >= 0.015
                and item_endpoint_safe
                and item_annotation_safe
            ):
                eligible.append(item)
        objective, fraction, trial, interface, selected_delta = (
            # Use the smallest correction that clears all blind gates.  The
            # boundary slice belongs to an internally coherent terminal block;
            # chasing the lowest single-interface score over-corrects it.
            min(eligible, key=lambda item: item[1])
            if eligible
            else next(item for item in candidates if item[1] == 0.0)
        )
    endpoint_safe = (
        interface["endpoint_low"] <= baseline_interface["endpoint_low"] * 1.08
        and interface["endpoint_high"] <= baseline_interface["endpoint_high"] * 1.08
    )
    annotation_safe = interface["annotation"] <= baseline_interface["annotation"] * 1.15
    improvement = (baseline_objective - objective) / max(baseline_objective, 1e-12)
    accepted = (
        fraction > 0 and improvement >= 0.015 and endpoint_safe and annotation_safe
    ) or (args.repair_policy == "learned_probe" and fraction == 0.25)
    proposal_fraction = fraction
    proposal_objective = objective
    proposal_interface = interface
    proposal_delta = selected_delta.copy()
    if not accepted:
        trial = transforms
        fraction = 0.0
        objective = baseline_objective
        selected_delta = np.eye(3)
    return trial, {
        "operation": "boundary_slice_fractional_reverse_spateo",
        "left_slice": slices[boundary].slice_id,
        "right_slice": slices[fixed].slice_id,
        "affected_slices": str(slices[moving].slice_id),
        "accepted": accepted,
        "selected_fraction": fraction,
        "score_before": baseline_objective,
        "score_after": objective,
        "improvement_fraction": improvement if accepted else 0.0,
        "proposal_fraction": proposal_fraction,
        "proposal_score": proposal_objective,
        "proposal_endpoint_max": proposal_interface["endpoint_max"],
        "proposal_annotation": proposal_interface["annotation"],
        "proposal_delta_rotation_deg": core.matrix_rotation_deg(proposal_delta),
        "proposal_delta_tx": float(proposal_delta[0, 2]),
        "proposal_delta_ty": float(proposal_delta[1, 2]),
        "delta_rotation_deg": core.matrix_rotation_deg(selected_delta),
        "delta_tx": float(selected_delta[0, 2]),
        "delta_ty": float(selected_delta[1, 2]),
    }


def lineage_pair_distance(
    points_left: np.ndarray,
    labels_left: np.ndarray,
    points_right: np.ndarray,
    labels_right: np.ndarray,
    label: str,
) -> float:
    from scipy.spatial import cKDTree

    left = points_left[labels_left.astype(str) == label]
    right = points_right[labels_right.astype(str) == label]
    if min(len(left), len(right)) < 8:
        return math.inf
    scale = 0.5 * (core.robust_diagonal(points_left) + core.robust_diagonal(points_right))
    forward = cKDTree(left).query(right, k=1, workers=1)[0]
    reverse = cKDTree(right).query(left, k=1, workers=1)[0]
    return 0.5 * (core.trimmed_mean(forward, 0.85) + core.trimmed_mean(reverse, 0.85)) / scale


def terminal_candidate_score(
    slices: list[core.SliceData],
    transforms: list[np.ndarray],
    terminal_indices: list[int],
    priority_label: str,
) -> tuple[float, float, float]:
    full_scores: list[float] = []
    lineage_scores: list[float] = []
    for left, right in zip(terminal_indices[:-1], terminal_indices[1:]):
        score = sampled_pair_score(slices, transforms, left, right)
        full_scores.append(score["total"])
        lineage_scores.append(
            lineage_pair_distance(
                core.aligned_sample(slices[left], transforms[left]),
                slices[left].sample_annotations,
                core.aligned_sample(slices[right], transforms[right]),
                slices[right].sample_annotations,
                priority_label,
            )
        )
    full = float(np.mean(full_scores))
    finite = [value for value in lineage_scores if np.isfinite(value)]
    lineage = float(np.mean(finite)) if finite else full
    return 0.52 * full + 0.48 * lineage, full, lineage


def priority_lineage_terminal_rescue(
    slices: list[core.SliceData],
    models: list[ad.AnnData],
    transforms: list[np.ndarray],
    terminal_indices: list[int],
    args: argparse.Namespace,
) -> tuple[list[np.ndarray], dict[str, object]]:
    shared = set(slices[terminal_indices[0]].annotations.astype(str))
    for index in terminal_indices[1:]:
        shared.intersection_update(slices[index].annotations.astype(str))
    eligible = [
        label
        for label in sorted(shared)
        if min(np.sum(slices[index].annotations.astype(str) == label) for index in terminal_indices) >= 40
    ]
    priority = args.priority_annotation if args.priority_annotation in eligible else ""
    if not priority and eligible:
        # Select a well-supported, spatially concentrated lineage without using
        # any reference coordinate.
        ranked: list[tuple[float, str]] = []
        for label in eligible:
            supports = []
            relative_spreads = []
            for index in terminal_indices:
                mask = slices[index].sample_annotations.astype(str) == label
                group = slices[index].sample_xy[mask]
                if len(group) < 8:
                    continue
                supports.append(len(group))
                relative_spreads.append(core.robust_diagonal(group) / core.robust_diagonal(slices[index].sample_xy))
            if supports:
                ranked.append((math.sqrt(min(supports)) / (0.15 + np.median(relative_spreads)), label))
        priority = max(ranked)[1] if ranked else ""
    if not priority:
        return transforms, {
            "operation": "priority_lineage_terminal_spateo",
            "accepted": False,
            "reason": "no_eligible_shared_annotation",
        }

    baseline_score, baseline_full, baseline_lineage = terminal_candidate_score(
        slices, transforms, terminal_indices, priority
    )
    ordered = list(reversed(terminal_indices))
    work: list[ad.AnnData] = []
    for index in ordered:
        mask = models[index].obs[args.annotation_key].astype(str).to_numpy() == priority
        work.append(
            make_work_model(models[index], aligned_points(slices, transforms, index), args, mask)
        )
    aligned = run_spateo_candidate_models(work, args)
    trial = [matrix.copy() for matrix in transforms]
    for position, index in enumerate(ordered):
        if position == 0:
            continue
        mask = models[index].obs[args.annotation_key].astype(str).to_numpy() == priority
        current_subset = aligned_points(slices, transforms, index)[mask]
        output_subset = np.asarray(aligned[position].obsm["align_spatial"], dtype=np.float64)
        delta = core.fit_proper_rigid(current_subset, output_subset)
        trial[index] = delta @ trial[index]
    score, full, lineage = terminal_candidate_score(slices, trial, terminal_indices, priority)
    lineage_improvement = (baseline_lineage - lineage) / max(baseline_lineage, 1e-12)
    delta_rotations: list[float] = []
    delta_translation_ratios: list[float] = []
    for index in terminal_indices[:-1]:
        delta = trial[index] @ np.linalg.inv(transforms[index])
        delta_rotations.append(abs(core.matrix_rotation_deg(delta)))
        before_points = aligned_points(slices, transforms, index)
        after_points = core.apply_matrix(before_points, delta)
        center_displacement = np.linalg.norm(
            np.median(after_points, axis=0) - np.median(before_points, axis=0)
        )
        delta_translation_ratios.append(
            float(center_displacement / core.robust_diagonal(before_points))
        )
    bounded = (
        max(delta_rotations, default=0.0) <= 12.0
        and max(delta_translation_ratios, default=0.0) <= 0.18
    )
    accepted = (
        lineage_improvement >= 0.15 and bounded
    ) or args.repair_policy == "learned_probe"
    proposal_score, proposal_full, proposal_lineage = score, full, lineage
    if not accepted:
        trial = transforms
        score, full, lineage = baseline_score, baseline_full, baseline_lineage
    changed = terminal_indices[:-1] if accepted else []
    return trial, {
        "operation": "priority_lineage_terminal_spateo",
        "priority_annotation": priority,
        "terminal_slices": ",".join(str(slices[index].slice_id) for index in terminal_indices),
        "affected_slices": ",".join(str(slices[index].slice_id) for index in changed),
        "accepted": accepted,
        "score_before": baseline_score,
        "score_after": score,
        "full_score_before": baseline_full,
        "full_score_after": full,
        "lineage_score_before": baseline_lineage,
        "lineage_score_after": lineage,
        "lineage_improvement_fraction": lineage_improvement if accepted else 0.0,
        "max_delta_rotation_deg": max(delta_rotations, default=0.0),
        "max_delta_translation_ratio": max(delta_translation_ratios, default=0.0),
        "proposal_score": proposal_score,
        "proposal_full_score": proposal_full,
        "proposal_lineage_score": proposal_lineage,
    }


def endpoint_balanced_block_translation(
    slices: list[core.SliceData],
    transforms: list[np.ndarray],
    boundary: int,
    args: argparse.Namespace,
) -> tuple[list[np.ndarray], dict[str, object]]:
    import _continuity_scoring as pose

    fixed = boundary + 1
    moving_points = core.aligned_sample(slices[boundary], transforms[boundary])
    fixed_points = core.aligned_sample(slices[fixed], transforms[fixed])
    translation = pose.balanced_endpoint_offset(moving_points, fixed_points)
    baseline = sampled_pair_score(slices, transforms, boundary, fixed)
    candidates: list[tuple[float, float, list[np.ndarray], dict[str, float], np.ndarray]] = []
    for multiplier in [0.0, 0.50, 0.75, 1.0, 1.25]:
        trial = [matrix.copy() for matrix in transforms]
        delta = core.rigid_matrix(0.0, multiplier * translation)
        for index in range(0, boundary + 1):
            trial[index] = delta @ trial[index]
        score = sampled_pair_score(slices, trial, boundary, fixed)
        candidates.append((score["total"], multiplier, trial, score, delta))
    candidates.sort(key=lambda item: item[0])
    if args.repair_policy == "learned_probe":
        _, multiplier, trial, score, delta = next(item for item in candidates if item[1] == 1.0)
    else:
        eligible = []
        for item in candidates:
            _, item_multiplier, _, item_score, _ = item
            item_endpoint_improvement = (
                baseline["endpoint_max"] - item_score["endpoint_max"]
            ) / max(baseline["endpoint_max"], 1e-12)
            item_endpoint_mean_improvement = (
                0.5 * (baseline["endpoint_low"] + baseline["endpoint_high"])
                - 0.5 * (item_score["endpoint_low"] + item_score["endpoint_high"])
            ) / max(0.5 * (baseline["endpoint_low"] + baseline["endpoint_high"]), 1e-12)
            item_annotation_safe = item_score["annotation"] <= baseline["annotation"] * 1.12
            item_objective_improvement = (
                baseline["total"] - item_score["total"]
            ) / max(baseline["total"], 1e-12)
            if (
                item_multiplier > 0
                and item_endpoint_improvement >= 0.15
                and item_endpoint_mean_improvement >= 0.05
                and item_annotation_safe
                and item_objective_improvement >= 0.02
            ):
                eligible.append(item)
        _, multiplier, trial, score, delta = (
            min(eligible, key=lambda item: item[0])
            if eligible
            else next(item for item in candidates if item[1] == 0.0)
        )
    endpoint_improvement = (baseline["endpoint_max"] - score["endpoint_max"]) / max(
        baseline["endpoint_max"], 1e-12
    )
    endpoint_mean_improvement = (
        0.5 * (baseline["endpoint_low"] + baseline["endpoint_high"])
        - 0.5 * (score["endpoint_low"] + score["endpoint_high"])
    ) / max(0.5 * (baseline["endpoint_low"] + baseline["endpoint_high"]), 1e-12)
    annotation_safe = score["annotation"] <= baseline["annotation"] * 1.12
    objective_improvement = (baseline["total"] - score["total"]) / max(baseline["total"], 1e-12)
    accepted = (
        multiplier > 0
        and endpoint_improvement >= 0.15
        and endpoint_mean_improvement >= 0.05
        and annotation_safe
        and objective_improvement >= 0.02
    ) or (args.repair_policy == "learned_probe" and multiplier == 1.0)
    proposal_multiplier = multiplier
    proposal_score = score.copy()
    proposal_delta = delta.copy()
    if not accepted:
        trial = transforms
        multiplier = 0.0
        score = baseline
        delta = np.eye(3)
    return trial, {
        "operation": "endpoint_balanced_block_translation",
        "left_slice": slices[boundary].slice_id,
        "right_slice": slices[fixed].slice_id,
        "affected_slices": (
            ",".join(str(slices[index].slice_id) for index in range(boundary + 1))
            if accepted
            else ""
        ),
        "accepted": accepted,
        "selected_multiplier": multiplier,
        "score_before": baseline["total"],
        "score_after": score["total"],
        "endpoint_max_before": baseline["endpoint_max"],
        "endpoint_max_after": score["endpoint_max"],
        "endpoint_improvement_fraction": endpoint_improvement if accepted else 0.0,
        "endpoint_mean_improvement_fraction": endpoint_mean_improvement if accepted else 0.0,
        "proposal_multiplier": proposal_multiplier,
        "proposal_score": proposal_score["total"],
        "proposal_endpoint_max": proposal_score["endpoint_max"],
        "proposal_annotation": proposal_score["annotation"],
        "proposal_delta_tx": float(proposal_delta[0, 2]),
        "proposal_delta_ty": float(proposal_delta[1, 2]),
        "delta_tx": float(delta[0, 2]),
        "delta_ty": float(delta[1, 2]),
    }


def sparse_terminal_singleton_rescue(
    slices: list[core.SliceData],
    models: list[ad.AnnData],
    transforms: list[np.ndarray],
    moving: int,
    fixed: int,
    median_edge_score: float,
    args: argparse.Namespace,
) -> tuple[list[np.ndarray], dict[str, object]]:
    left, right = sorted([moving, fixed])
    baseline = sampled_pair_score(slices, transforms, left, right)
    count_ratio = len(slices[moving].xy) / max(len(slices[fixed].xy), 1)
    triggered = count_ratio < 0.50 and baseline["total"] > 1.45 * median_edge_score
    if not triggered:
        return transforms, {
            "operation": "sparse_terminal_singleton_spateo",
            "moving_slice": slices[moving].slice_id,
            "fixed_slice": slices[fixed].slice_id,
            "affected_slices": "",
            "triggered": False,
            "accepted": False,
            "count_ratio": count_ratio,
            "score_before": baseline["total"],
            "median_edge_score": median_edge_score,
            "reason": "trigger_not_met",
        }

    moving_points = aligned_points(slices, transforms, moving)
    fixed_points = aligned_points(slices, transforms, fixed)
    work = [
        make_work_model(models[fixed], fixed_points, args),
        make_work_model(models[moving], moving_points, args),
    ]
    aligned = run_spateo_candidate_models(work, args)
    full_delta = core.fit_proper_rigid(
        moving_points, np.asarray(aligned[1].obsm["align_spatial"], dtype=np.float64)
    )
    candidates: list[tuple[float, float, list[np.ndarray], dict[str, float], np.ndarray]] = []
    for fraction in [0.0, 0.25, 0.50, 0.75, 1.0]:
        trial = [matrix.copy() for matrix in transforms]
        delta = fractional_rigid(full_delta, fraction)
        trial[moving] = delta @ trial[moving]
        score = sampled_pair_score(slices, trial, left, right)
        candidates.append((score["total"], fraction, trial, score, delta))
    candidates.sort(key=lambda item: item[0])
    _, fraction, trial, score, delta = candidates[0]
    improvement = (baseline["total"] - score["total"]) / max(baseline["total"], 1e-12)
    after_points = core.apply_matrix(moving_points, delta)
    center_displacement_ratio = float(
        np.linalg.norm(np.median(after_points, axis=0) - np.median(moving_points, axis=0))
        / core.robust_diagonal(moving_points)
    )
    bounded = (
        abs(core.matrix_rotation_deg(delta)) <= 20.0
        and center_displacement_ratio <= 0.25
    )
    accepted = fraction > 0 and improvement >= 0.08 and bounded
    if not accepted:
        trial = transforms
        fraction = 0.0
        score = baseline
        delta = np.eye(3)
    return trial, {
        "operation": "sparse_terminal_singleton_spateo",
        "moving_slice": slices[moving].slice_id,
        "fixed_slice": slices[fixed].slice_id,
        "affected_slices": str(slices[moving].slice_id) if accepted else "",
        "triggered": True,
        "accepted": accepted,
        "count_ratio": count_ratio,
        "score_before": baseline["total"],
        "score_after": score["total"],
        "improvement_fraction": improvement if accepted else 0.0,
        "selected_fraction": fraction,
        "median_edge_score": median_edge_score,
        "center_displacement_ratio": center_displacement_ratio,
        "delta_rotation_deg": core.matrix_rotation_deg(delta),
        "delta_tx": float(delta[0, 2]),
        "delta_ty": float(delta[1, 2]),
    }


def robust_high_threshold(
    values: list[float], ratio: float, mad_multiplier: float
) -> tuple[float, float, float]:
    array = np.asarray(values, dtype=np.float64)
    median = float(np.median(array))
    mad = float(np.median(np.abs(array - median)))
    robust_sigma = 1.4826 * mad
    threshold = max(float(ratio) * median, median + float(mad_multiplier) * robust_sigma)
    return threshold, median, mad


def internal_block_translation_repairs(
    slices: list[core.SliceData],
    transforms: list[np.ndarray],
    args: argparse.Namespace,
) -> tuple[list[np.ndarray], list[dict[str, object]]]:
    """Repair a translation discontinuity between two internally coherent blocks.

    A prefix translation is gauge-equivalent to the opposite suffix translation,
    while preserving every within-block relative pose. Candidate construction and
    acceptance use only the frozen blind coordinates and annotation-aware scores.
    """
    operations: list[dict[str, object]] = []
    minimum = max(2, int(args.internal_block_min_side_slices))
    if len(slices) < 2 * minimum:
        return transforms, operations

    for pass_index in range(max(0, int(args.internal_block_passes))):
        score_rows = [
            sampled_pair_score(slices, transforms, index, index + 1)
            for index in range(len(slices) - 1)
        ]
        totals = [float(row["total"]) for row in score_rows]
        endpoints = [float(row["endpoint_max"]) for row in score_rows]
        score_threshold, score_median, score_mad = robust_high_threshold(
            totals,
            args.internal_block_score_ratio,
            args.internal_block_mad_multiplier,
        )
        endpoint_median = float(np.median(endpoints))
        endpoint_threshold = max(
            float(args.internal_block_endpoint_floor),
            float(args.internal_block_endpoint_ratio) * endpoint_median,
        )
        candidates = [
            index
            for index, row in enumerate(score_rows)
            if index + 1 >= minimum
            and len(slices) - index - 1 >= minimum
            and float(row["total"]) >= score_threshold
            and float(row["endpoint_max"]) >= endpoint_threshold
        ]
        candidates.sort(
            key=lambda index: (
                float(score_rows[index]["total"]) / max(score_threshold, 1e-12)
                + float(score_rows[index]["endpoint_max"]) / max(endpoint_threshold, 1e-12)
            ),
            reverse=True,
        )
        if not candidates:
            break

        accepted_this_pass = False
        for boundary in candidates:
            before = [matrix.copy() for matrix in transforms]
            trial, operation = endpoint_balanced_block_translation(
                slices, before, boundary, args
            )
            operation["operation"] = "internal_block_endpoint_translation"
            operation["triggered"] = True
            operation["internal_pass"] = pass_index + 1
            operation["score_threshold"] = score_threshold
            operation["score_median"] = score_median
            operation["score_mad"] = score_mad
            operation["endpoint_threshold"] = endpoint_threshold
            operation["endpoint_median"] = endpoint_median
            operation["left_block_n_slices"] = boundary + 1
            operation["right_block_n_slices"] = len(slices) - boundary - 1
            local_scale = 0.5 * (
                core.robust_diagonal(core.aligned_sample(slices[boundary], before[boundary]))
                + core.robust_diagonal(core.aligned_sample(slices[boundary + 1], before[boundary + 1]))
            )
            displacement = float(
                np.hypot(
                    float(operation.get("proposal_delta_tx", 0.0)),
                    float(operation.get("proposal_delta_ty", 0.0)),
                )
            )
            displacement_ratio = displacement / max(local_scale, 1e-12)
            operation["proposal_center_displacement_ratio"] = displacement_ratio
            strict_improvement = (
                float(operation["score_before"]) - float(operation["proposal_score"])
            ) / max(float(operation["score_before"]), 1e-12)
            operation["proposal_objective_improvement_fraction"] = strict_improvement
            strict_accepted = (
                bool(operation.get("accepted", False))
                and strict_improvement >= float(args.internal_block_min_improvement)
                and float(operation.get("endpoint_mean_improvement_fraction", 0.0))
                >= float(args.internal_block_min_endpoint_mean_improvement)
                and displacement_ratio <= 0.25
            )
            if not strict_accepted:
                if bool(operation.get("accepted", False)):
                    operation["reason"] = "internal_block_strict_gate_not_met"
                operation["accepted"] = False
                operation["affected_slices"] = ""
                operation["selected_multiplier"] = 0.0
                operation["score_after"] = operation["score_before"]
                operation["objective_improvement_fraction"] = 0.0
                operation["delta_tx"] = 0.0
                operation["delta_ty"] = 0.0
                trial = before
            else:
                operation["objective_improvement_fraction"] = strict_improvement
            operations.append(operation)
            if strict_accepted:
                transforms = trial
                accepted_this_pass = True
                break
        if not accepted_this_pass:
            break
    return transforms, operations


def automatic_continuity_repairs(
    slices: list[core.SliceData],
    models: list[ad.AnnData],
    transforms: list[np.ndarray],
    args: argparse.Namespace,
) -> tuple[list[np.ndarray], list[dict[str, object]], int]:
    anchor = int(np.argmax([len(item.xy) for item in slices]))
    operations: list[dict[str, object]] = []
    # This repair chain is specifically for a short sparse terminal block next
    # to the largest high-support slice.  A longer side requires its own pose
    # graph; treating it as a three-slice end would create a new discontinuity.
    if 2 <= anchor <= 3:
        boundary = anchor - 1
        try:
            transforms, operation = reverse_boundary_fractional_rescue(
                slices, models, transforms, boundary, args
            )
        except Exception as error:
            operation = {
                "operation": "boundary_slice_fractional_reverse_spateo",
                "accepted": False,
                "error": f"{type(error).__name__}: {error}",
                "fallback": "retain_pre_operation_coordinates",
            }
        operations.append(operation)
        terminal = list(range(0, anchor))
        try:
            transforms, operation = priority_lineage_terminal_rescue(
                slices, models, transforms, terminal, args
            )
        except Exception as error:
            operation = {
                "operation": "priority_lineage_terminal_spateo",
                "accepted": False,
                "error": f"{type(error).__name__}: {error}",
                "fallback": "retain_pre_operation_coordinates",
            }
        operations.append(operation)
        try:
            transforms, operation = endpoint_balanced_block_translation(
                slices, transforms, boundary, args
            )
        except Exception as error:
            operation = {
                "operation": "endpoint_balanced_block_translation",
                "accepted": False,
                "error": f"{type(error).__name__}: {error}",
                "fallback": "retain_pre_operation_coordinates",
            }
        operations.append(operation)
    try:
        transforms, internal_operations = internal_block_translation_repairs(
            slices, transforms, args
        )
    except Exception as error:
        internal_operations = [
            {
                "operation": "internal_block_endpoint_translation",
                "accepted": False,
                "error": f"{type(error).__name__}: {error}",
                "fallback": "retain_pre_operation_coordinates",
            }
        ]
    operations.extend(internal_operations)
    current_scores = [
        sampled_pair_score(slices, transforms, index, index + 1)["total"]
        for index in range(len(slices) - 1)
    ]
    median_edge_score = float(np.median(current_scores))
    for moving, fixed in [(0, 1), (len(slices) - 1, len(slices) - 2)]:
        for pass_index in range(args.terminal_passes):
            try:
                transforms, operation = sparse_terminal_singleton_rescue(
                    slices, models, transforms, moving, fixed, median_edge_score, args
                )
            except Exception as error:
                operation = {
                    "operation": "sparse_terminal_singleton_spateo",
                    "moving_slice": slices[moving].slice_id,
                    "fixed_slice": slices[fixed].slice_id,
                    "accepted": False,
                    "error": f"{type(error).__name__}: {error}",
                    "fallback": "retain_pre_operation_coordinates",
                }
            operation["terminal_pass"] = pass_index + 1
            operations.append(operation)
            if not bool(operation.get("accepted", False)):
                break
    return transforms, operations, anchor


def normalize_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    return [{key: row.get(key, "") for key in fields} for row in rows]


def extract_transforms(
    slices: list[core.SliceData], aligned_models: list[ad.AnnData]
) -> tuple[list[np.ndarray], list[dict[str, object]]]:
    transforms: list[np.ndarray] = []
    rows: list[dict[str, object]] = []
    for index, (slice_data, model) in enumerate(zip(slices, aligned_models)):
        aligned = np.asarray(model.obsm["align_spatial"], dtype=np.float64)
        matrix = core.fit_proper_rigid(slice_data.xy, aligned)
        replay = core.apply_matrix(slice_data.xy, matrix)
        residual = np.linalg.norm(replay - aligned, axis=1)
        transforms.append(matrix)
        vecfld = model.uns.get("VecFld_morpho", {})
        rows.append(
            {
                "slice_id": slice_data.slice_id,
                "role": "fixed_anchor" if index == 0 else "moving",
                "rotation_deg": core.matrix_rotation_deg(matrix),
                "tx": float(matrix[0, 2]),
                "ty": float(matrix[1, 2]),
                "rigid_replay_rmse": float(np.sqrt(np.mean(residual**2))),
                "rigid_replay_max": float(np.max(residual)),
                "sigma2": vecfld.get("sigma2", "") if isinstance(vecfld, dict) else "",
            }
        )
    return transforms, rows


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise FileExistsError(
            f"Output directory already exists; use a new immutable path: {args.output_dir}"
        )
    args.output_dir.mkdir(parents=True)
    slices = core.load_slices(args)
    models = load_models(slices, args)
    aligned_models, assignments = run_spateo_serial(models, args)
    transforms, pair_rows = extract_transforms(slices, aligned_models)
    baseline_transforms = [matrix.copy() for matrix in transforms]
    operations: list[dict[str, object]] = []
    continuity_anchor = 0
    if args.postprocess == "auto":
        transforms, operations, continuity_anchor = automatic_continuity_repairs(
            slices, models, transforms, args
        )

    coordinate_path = args.output_dir / "aligned_coordinates.csv.gz"
    transform_path = args.output_dir / "slice_transforms.csv"
    pair_path = args.output_dir / "spateo_pair_summary.csv"
    operation_path = args.output_dir / "interface_operations.csv"
    qc_path = args.output_dir / "blind_qc_summary.json"
    frozen_path = args.output_dir / "frozen_manifest.json"
    overview_path = args.output_dir / "blind_alignment_overview.png"
    edge_plot_path = args.output_dir / "interface_scores.png"

    n_cells = core.write_coordinates(coordinate_path, slices, transforms)
    transform_rows: list[dict[str, object]] = []
    for slice_data, matrix in zip(slices, transforms):
        transform_rows.append(
            {
                "slice_id": slice_data.slice_id,
                "source_file": str(slice_data.path),
                "n_cells": len(slice_data.xy),
                "m00": matrix[0, 0],
                "m01": matrix[0, 1],
                "m02": matrix[0, 2],
                "m10": matrix[1, 0],
                "m11": matrix[1, 1],
                "m12": matrix[1, 2],
                "rotation_deg": core.matrix_rotation_deg(matrix),
            }
        )
    core.write_csv(transform_path, transform_rows)
    core.write_csv(pair_path, pair_rows)
    core.write_csv(operation_path, normalize_rows(operations))
    core.plot_overview(overview_path, slices, transforms, args.seed)
    baseline_edges = core.compute_edge_scores(slices, baseline_transforms)
    edges = core.compute_edge_scores(slices, transforms)
    core.plot_edge_scores(edge_plot_path, slices, baseline_edges, edges)

    parameters = {key: value for key, value in vars(args).items() if key != "output_dir"}
    parameters["slice_dir"] = str(parameters["slice_dir"])
    assignment_shapes = [list(item.shape) for item in assignments]
    qc = {
        "pipeline_version": PIPELINE_VERSION,
        "stage": args.stage,
        "status": "pass",
        "blind_contract": {
            "spatial_3d_available_to_runner": False,
            "ground_truth_used_for_generation": False,
            "ground_truth_used_for_selection": False,
            "identity_source": "obs_names",
            "annotation_source": f"obsm[{args.representation_key}] derived from obs[{args.annotation_key}]",
            "coordinate_source": f"obsm[{args.spatial_key}]",
        },
        "n_slices": len(slices),
        "n_cells": n_cells,
        "slice_ids": [item.slice_id for item in slices],
        "anchor_slice": slices[0].slice_id,
        "continuity_anchor_slice": slices[continuity_anchor].slice_id,
        "parameters": parameters,
        "assignment_shapes": assignment_shapes,
        "baseline_edge_score_mean": float(np.mean([row["total"] for row in baseline_edges])),
        "edge_score_mean": float(np.mean([row["total"] for row in edges])),
        "edge_score_max": float(np.max([row["total"] for row in edges])),
        "accepted_operations": sum(bool(row.get("accepted", False)) for row in operations),
        "internal_block_candidates": sum(
            row.get("operation") == "internal_block_endpoint_translation"
            for row in operations
        ),
        "accepted_internal_block_repairs": sum(
            row.get("operation") == "internal_block_endpoint_translation"
            and bool(row.get("accepted", False))
            for row in operations
        ),
        "operations": operations,
        "baseline_edges": baseline_edges,
        "edges": edges,
    }
    qc_path.write_text(json.dumps(qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    frozen = {
        "pipeline_version": PIPELINE_VERSION,
        "stage": args.stage,
        "status": "frozen_before_reference_evaluation",
        "frozen_at": datetime.now().astimezone().isoformat(),
        "ground_truth_used_before_freeze": False,
        "spatial_3d_read_before_freeze": False,
        "output": {"path": str(coordinate_path), "sha256": core.sha256_file(coordinate_path), "n_cells": n_cells},
        "recipe": {"path": str(transform_path), "sha256": core.sha256_file(transform_path)},
        "runner": {"path": str(Path(__file__).resolve()), "sha256": core.sha256_file(Path(__file__).resolve()), "python": sys.executable},
        "inputs": [
            {
                "slice_id": item.slice_id,
                "path": str(item.path),
                "n_cells": len(item.xy),
                "blind_content_sha256": item.content_sha256,
            }
            for item in slices
        ],
        "parameters": parameters,
        "artifacts": {
            "spateo_pair_summary": str(pair_path),
            "interface_operations": str(operation_path),
            "blind_qc_summary": str(qc_path),
            "blind_alignment_overview": str(overview_path),
            "interface_scores": str(edge_plot_path),
        },
    }
    frozen_path.write_text(json.dumps(frozen, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": "frozen",
                "stage": args.stage,
                "output_dir": str(args.output_dir),
                "n_cells": n_cells,
                "anchor_slice": slices[0].slice_id,
                "continuity_anchor_slice": slices[continuity_anchor].slice_id,
                "accepted_operations": qc["accepted_operations"],
                "baseline_edge_score_mean": qc["baseline_edge_score_mean"],
                "edge_score_mean": qc["edge_score_mean"],
            },
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise
