#!/usr/bin/env python
from __future__ import annotations

import json
import os
import re
import datetime as dt
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


if "DATASET_ROOT" not in os.environ:
    raise SystemExit("DATASET_ROOT is required and must point to the input h5ad directory.")
if "OUTDIR" not in os.environ:
    raise SystemExit("OUTDIR is required and must point to a writable result directory.")

ROOT = Path(os.environ["DATASET_ROOT"])
OUTDIR = Path(os.environ["OUTDIR"])
PAIR_START_ORDER = int(os.environ.get("PAIR_START_ORDER", "1"))
N_SLICES = int(os.environ.get("N_SLICES", "0"))
MAX_CELLS_PER_SLICE = int(os.environ.get("MAX_CELLS_PER_SLICE", "0"))
BATCH_SIZE = int(os.environ.get("BATCH_SIZE", "800"))
STAGE1_MAX_ITER = int(os.environ.get("STAGE1_MAX_ITER", os.environ.get("MAX_ITER", "300")))
STAGE2_MAX_ITER = int(os.environ.get("STAGE2_MAX_ITER", os.environ.get("MAX_ITER", "300")))
STAGE1_MODE = os.environ.get("STAGE1_MODE", "SN-S")
STAGE2_MODE = os.environ.get("STAGE2_MODE", "none")
N_SAMPLING_REF = int(os.environ.get("N_SAMPLING_REF", "40000"))
USE_MORPHO_ALIGN_REF = os.environ.get("USE_MORPHO_ALIGN_REF", "1") == "1"
SAVE_FULL_ASSIGNMENT = os.environ.get("SAVE_FULL_ASSIGNMENT", "0") == "1"
Z_DISPLAY_SPACING = float(os.environ.get("Z_DISPLAY_SPACING", "40"))
EXCLUDE_CELLTYPES_REGEX = os.environ.get("EXCLUDE_CELLTYPES_REGEX", "")
RNG_SEED = int(os.environ.get("RNG_SEED", "20260612"))
SPATEO_DEVICE = os.environ.get("SPATEO_DEVICE", "").strip()
REQUIRE_CUDA = os.environ.get("REQUIRE_CUDA", "1") == "1"
USE_SPATEO_INTERNAL_DOWNSAMPLING = os.environ.get("USE_SPATEO_INTERNAL_DOWNSAMPLING", "1") == "1"
INITIAL_PRETRANSFORM_RECIPE = os.environ.get("INITIAL_PRETRANSFORM_RECIPE", "").strip()
INITIAL_COORDINATES_CSV = os.environ.get("INITIAL_COORDINATES_CSV", "").strip()
INITIAL_COORDINATE_CELL_ID_COL = os.environ.get("INITIAL_COORDINATE_CELL_ID_COL", "cell_id").strip()
INITIAL_COORDINATE_X_COL = os.environ.get("INITIAL_COORDINATE_X_COL", "manual_x").strip()
INITIAL_COORDINATE_Y_COL = os.environ.get("INITIAL_COORDINATE_Y_COL", "manual_y").strip()
ALLOW_PARTIAL_INITIAL_COORDINATES = os.environ.get("ALLOW_PARTIAL_INITIAL_COORDINATES", "0") == "1"


def stage2_is_enabled() -> bool:
    return STAGE2_MODE.strip().lower() not in {"", "none", "skip", "off", "false", "0"}


def decode_array(arr):
    arr = np.asarray(arr)
    if arr.dtype.kind == "S":
        return arr.astype(str)
    if arr.dtype.kind == "O":
        return np.array(
            [x.decode("utf-8", "replace") if isinstance(x, (bytes, bytearray)) else str(x) for x in arr],
            dtype=object,
        )
    return arr


def read_obs_column(obs, name: str, idx: np.ndarray) -> np.ndarray:
    obj = obs[name]
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        cats = decode_array(obj["categories"][()])
        codes = obj["codes"][()][idx]
        values = np.array(["NA"] * len(idx), dtype=object)
        valid = codes >= 0
        values[valid] = cats[codes[valid].astype(int)]
        return values
    values = decode_array(obj[()])
    return values[idx].astype(str)


def collect_global_order() -> pd.DataFrame:
    rows = []
    for fp in sorted(ROOT.rglob("*.h5ad")):
        match = re.search(r"CS(?P<stage>\d+)_SL(?P<sl>\d+)_(?P<chip>Y\w+)\.Spatial", fp.name)
        if not match:
            continue
        rows.append(
            {
                "stage": f"CS{match.group('stage')}",
                "sl_number": int(match.group("sl")),
                "chip_id": match.group("chip"),
                "slice_id": fp.name.split(".Spatial")[0],
                "file": str(fp),
            }
        )
    if not rows:
        raise SystemExit(f"No h5ad files discovered under {ROOT}")
    df = pd.DataFrame(rows).sort_values(["stage", "sl_number"]).reset_index(drop=True)
    df["global_order_by_SL"] = df.index + 1
    df["relative_z_by_SL"] = df["sl_number"] - df["sl_number"].min()
    df["z_display"] = (df["global_order_by_SL"] - 1) * Z_DISPLAY_SPACING
    return df


def load_slice(row: pd.Series, rng: np.random.Generator) -> ad.AnnData:
    with h5py.File(row["file"], "r") as handle:
        spatial = handle["obsm"]["spatial"][()]
        x_pca = handle["obsm"]["X_pca"][()]
        n_obs = spatial.shape[0]

        candidate_idx = np.arange(n_obs)
        if EXCLUDE_CELLTYPES_REGEX and "celltype" in handle["obs"]:
            pattern = re.compile(EXCLUDE_CELLTYPES_REGEX)
            all_celltypes = read_obs_column(handle["obs"], "celltype", candidate_idx)
            keep = np.array([pattern.search(str(x)) is None for x in all_celltypes], dtype=bool)
            candidate_idx = candidate_idx[keep]

        n_candidates = len(candidate_idx)
        if MAX_CELLS_PER_SLICE > 0 and n_candidates > MAX_CELLS_PER_SLICE:
            idx = np.sort(rng.choice(candidate_idx, size=MAX_CELLS_PER_SLICE, replace=False))
            sampling_mode = f"sampled_{MAX_CELLS_PER_SLICE}"
        else:
            idx = candidate_idx
            sampling_mode = "full"

        obs = pd.DataFrame(index=[f"{row['slice_id']}:{i}" for i in idx])
        obs["slice_id"] = row["slice_id"]
        obs["chip_id"] = row["chip_id"]
        obs["stage"] = row["stage"]
        obs["sl_number"] = int(row["sl_number"])
        obs["relative_z_by_SL"] = float(row["relative_z_by_SL"])
        obs["global_order_by_SL"] = int(row["global_order_by_SL"])
        obs["z_display"] = float(row["z_display"])
        obs["sampling_mode"] = sampling_mode
        if "celltype" in handle["obs"]:
            obs["celltype"] = read_obs_column(handle["obs"], "celltype", idx)

        adata = ad.AnnData(X=sparse.csr_matrix((len(idx), 1), dtype=np.float32), obs=obs)
        adata.var_names = ["dummy"]
        adata.obsm["spatial"] = spatial[idx].astype(np.float32)
        adata.obsm["X_pca"] = x_pca[idx].astype(np.float32)
        return adata


def to_jsonable(value: Any):
    if hasattr(value, "detach") and hasattr(value, "cpu"):
        return value.detach().cpu().numpy().tolist()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, dict):
        return {k: to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(v) for v in value]
    return value


def norm_slice_token(value: Any) -> str:
    text = str(value)
    match = re.search(r"(\d+)", text)
    return match.group(1) if match else text


def step_matches_slice(step: dict[str, Any], adata: ad.AnnData) -> bool:
    scope = step.get("scope") or {}
    labels = scope.get("slice_labels") or []
    start = scope.get("slice_start")
    end = scope.get("slice_end")
    slice_id = str(adata.obs["slice_id"].iloc[0])
    sl_number = int(adata.obs["sl_number"].iloc[0])
    if labels:
        wanted = {norm_slice_token(x) for x in labels}
        return norm_slice_token(slice_id) in wanted or str(sl_number) in wanted
    if start is not None or end is not None:
        s = int(norm_slice_token(start if start is not None else end))
        e = int(norm_slice_token(end if end is not None else start))
        lo, hi = sorted([s, e])
        return lo <= sl_number <= hi
    return False


def pretransform_center(step: dict[str, Any], affected_xy: np.ndarray) -> np.ndarray:
    center = step.get("center") or {}
    if center.get("mode") == "manual":
        return np.array([float(center.get("x", 0)), float(center.get("y", 0))], dtype=np.float32)
    return affected_xy.mean(axis=0).astype(np.float32)


def apply_initial_pretransform(raw_slices: list[ad.AnnData], recipe_path: str) -> dict[str, Any]:
    path = Path(recipe_path).expanduser()
    if not path.exists():
        raise SystemExit(f"INITIAL_PRETRANSFORM_RECIPE does not exist: {path}")
    with path.open() as handle:
        recipe = json.load(handle)

    for adata in raw_slices:
        adata.obsm["raw_input_spatial"] = adata.obsm["spatial"].copy()

    applied_steps = []
    for step_index, step in enumerate(recipe.get("steps", []), start=1):
        op = step.get("operation")
        if op not in {"translate", "rotate_z", "scale", "mirror_x", "mirror_y"}:
            continue
        matched = [adata for adata in raw_slices if step_matches_slice(step, adata)]
        if not matched:
            applied_steps.append({"step_index": step_index, "operation": op, "matched_slices": [], "affected_cells": 0})
            continue
        affected_xy = np.concatenate([adata.obsm["spatial"][:, :2] for adata in matched], axis=0)
        center = pretransform_center(step, affected_xy)
        affected_cells = 0
        params = step.get("parameters") or {}
        for adata in matched:
            xy = adata.obsm["spatial"][:, :2].astype(np.float32, copy=True)
            if op == "translate":
                xy[:, 0] += float(params.get("dx", 0))
                xy[:, 1] += float(params.get("dy", 0))
            elif op == "rotate_z":
                theta = np.deg2rad(float(params.get("degrees", 0)))
                rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]], dtype=np.float32)
                xy = (xy - center) @ rot.T + center
            elif op == "scale":
                scale = float(params.get("scale", 1))
                xy = (xy - center) * scale + center
            elif op == "mirror_x":
                xy[:, 0] = 2 * center[0] - xy[:, 0]
            elif op == "mirror_y":
                xy[:, 1] = 2 * center[1] - xy[:, 1]
            adata.obsm["spatial"][:, :2] = xy
            affected_cells += adata.n_obs
        applied_steps.append(
            {
                "step_index": step_index,
                "operation": op,
                "matched_slices": [str(adata.obs["slice_id"].iloc[0]) for adata in matched],
                "affected_cells": int(affected_cells),
                "center_xy": center.tolist(),
                "parameters": params,
            }
        )

    for adata in raw_slices:
        adata.obsm["initial_pretransform_spatial"] = adata.obsm["spatial"].copy()

    summary = {
        "recipe_path": str(path),
        "recipe_version": recipe.get("version", ""),
        "created_at": recipe.get("created_at", ""),
        "applied_at": dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds"),
        "steps": applied_steps,
    }
    with (OUTDIR / "initial_pretransform_applied.json").open("w") as handle:
        json.dump(to_jsonable(summary), handle, indent=2)
    return summary


def apply_initial_coordinates_csv(raw_slices: list[ad.AnnData], csv_path: str) -> dict[str, Any]:
    path = Path(csv_path).expanduser()
    if not path.exists():
        raise SystemExit(f"INITIAL_COORDINATES_CSV does not exist: {path}")

    usecols = [
        INITIAL_COORDINATE_CELL_ID_COL,
        INITIAL_COORDINATE_X_COL,
        INITIAL_COORDINATE_Y_COL,
    ]
    try:
        coords = pd.read_csv(path, usecols=usecols)
    except ValueError as exc:
        raise SystemExit(
            "INITIAL_COORDINATES_CSV is missing required columns: "
            f"{', '.join(usecols)}"
        ) from exc
    if coords[INITIAL_COORDINATE_CELL_ID_COL].duplicated().any():
        dup = coords.loc[coords[INITIAL_COORDINATE_CELL_ID_COL].duplicated(), INITIAL_COORDINATE_CELL_ID_COL].iloc[0]
        raise SystemExit(f"INITIAL_COORDINATES_CSV has duplicated cell id: {dup}")

    coords = coords.set_index(INITIAL_COORDINATE_CELL_ID_COL)
    applied = []
    total_missing = 0
    for adata in raw_slices:
        if "raw_input_spatial" not in adata.obsm:
            adata.obsm["raw_input_spatial"] = adata.obsm["spatial"].copy()
        aligned = coords.reindex(adata.obs_names)
        missing = aligned[INITIAL_COORDINATE_X_COL].isna() | aligned[INITIAL_COORDINATE_Y_COL].isna()
        n_missing = int(missing.sum())
        total_missing += n_missing
        if n_missing and not ALLOW_PARTIAL_INITIAL_COORDINATES:
            sid = str(adata.obs["slice_id"].iloc[0])
            raise SystemExit(
                "INITIAL_COORDINATES_CSV does not cover all loaded cells. "
                f"First incomplete slice={sid}, missing={n_missing}/{adata.n_obs}. "
                "Use a full coordinate table, or set ALLOW_PARTIAL_INITIAL_COORDINATES=1 "
                "only for explicit downsampled tests."
            )
        xy = adata.obsm["spatial"][:, :2].astype(np.float32, copy=True)
        present = ~missing.to_numpy()
        if present.any():
            xy[present, 0] = aligned.loc[present, INITIAL_COORDINATE_X_COL].to_numpy(dtype=np.float32)
            xy[present, 1] = aligned.loc[present, INITIAL_COORDINATE_Y_COL].to_numpy(dtype=np.float32)
            adata.obsm["spatial"][:, :2] = xy
        adata.obsm["initial_coordinate_spatial"] = adata.obsm["spatial"].copy()
        applied.append(
            {
                "slice_id": str(adata.obs["slice_id"].iloc[0]),
                "sl_number": int(adata.obs["sl_number"].iloc[0]),
                "n_cells": int(adata.n_obs),
                "matched_cells": int(present.sum()),
                "missing_cells": n_missing,
            }
        )

    summary = {
        "csv_path": str(path),
        "cell_id_column": INITIAL_COORDINATE_CELL_ID_COL,
        "x_column": INITIAL_COORDINATE_X_COL,
        "y_column": INITIAL_COORDINATE_Y_COL,
        "allow_partial": ALLOW_PARTIAL_INITIAL_COORDINATES,
        "total_missing_cells": int(total_missing),
        "applied_at": dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds"),
        "slices": applied,
    }
    with (OUTDIR / "initial_coordinates_applied.json").open("w") as handle:
        json.dump(to_jsonable(summary), handle, indent=2)
    return summary


def save_stage_audit(
    aligned_slices: list[ad.AnnData],
    pis: list[np.ndarray],
    stage_name: str,
    vecfld_key: str,
) -> None:
    audit_dir = OUTDIR / "spateo_pair_audit" / stage_name
    audit_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, pi in enumerate(pis):
        fixed = aligned_slices[i]
        moving = aligned_slices[i + 1]
        fixed_id = fixed.obs["slice_id"].iloc[0]
        moving_id = moving.obs["slice_id"].iloc[0]
        vecfld = moving.uns.get(vecfld_key, {})
        pair_name = f"pair_{i + 1:03d}_{fixed_id}__{moving_id}"
        small_keys = [
            "beta",
            "normalize_c",
            "dissimilarity",
            "sigma2",
            "gamma",
            "NA",
            "sigma2_variance",
            "method",
            "kernel_type",
        ]
        summary = {
            "stage": stage_name,
            "pair_index": i + 1,
            "fixed_slice": fixed_id,
            "moving_slice": moving_id,
            "assignment_shape": list(pi.shape),
            "saved_full_assignment": SAVE_FULL_ASSIGNMENT,
            "vecfld_summary": {key: to_jsonable(vecfld.get(key)) for key in small_keys if key in vecfld},
        }
        matrix_payload = {}
        for key in [
            "R",
            "t",
            "optimal_R",
            "optimal_t",
            "init_R",
            "init_t",
            "Coff",
            "inducing_variables",
            "normalize_scales",
            "normalize_means",
        ]:
            if key in vecfld and vecfld[key] is not None:
                matrix_payload[key] = np.asarray(vecfld[key])
        if SAVE_FULL_ASSIGNMENT:
            matrix_payload["P"] = np.asarray(pi)
        elif pi.size:
            row_cap = min(64, pi.shape[0])
            col_cap = min(512, pi.shape[1])
            matrix_payload["P_sample"] = np.asarray(pi[:row_cap, :col_cap])
            summary["assignment_sample_shape"] = [row_cap, col_cap]
        np.savez_compressed(audit_dir / f"{pair_name}.npz", **matrix_payload)
        with (audit_dir / f"{pair_name}.json").open("w") as handle:
            json.dump(to_jsonable(summary), handle, indent=2)
        rows.append(summary)
    with (audit_dir / "spateo_pair_audit_summary.json").open("w") as handle:
        json.dump(to_jsonable(rows), handle, indent=2)


def run_morpho_stage(
    slices: list[ad.AnnData],
    spatial_key: str,
    key_added: str,
    vecfld_key: str,
    mode: str,
    max_iter: int,
    device: str,
):
    from spateo.alignment.morpho_alignment import morpho_align, morpho_align_ref

    def make_reference_models(models: list[ad.AnnData]) -> list[ad.AnnData]:
        rng = np.random.default_rng(RNG_SEED)
        refs = []
        for model in models:
            n_obs = model.n_obs
            n = min(N_SAMPLING_REF, n_obs)
            if n < n_obs:
                idx = np.sort(rng.choice(np.arange(n_obs), size=n, replace=False))
                ref = model[idx, :].copy()
            else:
                ref = model.copy()
            ref.obs["reference_sampling_mode"] = f"random_{n}" if n < n_obs else "full"
            refs.append(ref)
        return refs

    common_kwargs = dict(
        spatial_key=spatial_key,
        key_added=key_added,
        iter_key_added=None,
        vecfld_key_added=vecfld_key,
        mode=mode,
        device=device,
        verbose=True,
        rep_layer="X_pca",
        rep_field="obsm",
        dissimilarity="cos",
        max_iter=max_iter,
        nn_init=False,
        SVI_mode=True,
        batch_size=BATCH_SIZE,
        pre_compute_dist=False,
        sparse_calculation_mode=False,
        sparse_top_k=256,
        use_chunk=True,
        chunk_capacity=1,
        K=50,
        beta=0.05,
        lambdaVF=1,
        partial_robust_level=1,
        sigma2_init_scale=2,
    )
    if USE_MORPHO_ALIGN_REF:
        refs = None if USE_SPATEO_INTERNAL_DOWNSAMPLING else make_reference_models(slices)
        aligned, aligned_ref, pis, pis_ref = morpho_align_ref(
            models=slices,
            models_ref=refs,
            n_sampling=N_SAMPLING_REF if refs is None else None,
            sampling_method="random",
            **common_kwargs,
        )
        return aligned, pis, aligned_ref, pis_ref
    aligned, pis = morpho_align(models=slices, **common_kwargs)
    return aligned, pis, None, None


def copy_stage_coordinates(
    stage1: list[ad.AnnData],
    stage2_input: list[ad.AnnData],
    stage1_key: str,
) -> None:
    for rigid_model, target_model in zip(stage1, stage2_input):
        rigid_xy = rigid_model.obsm[stage1_key].astype(np.float32)
        target_model.obsm["stage1_rigid_spatial"] = rigid_xy.copy()
        target_model.obsm["spatial"] = target_model.obsm["spatial"].astype(np.float32)
        target_model.obsm["raw_spatial"] = target_model.obsm["spatial"].copy()


def make_3d(xy: np.ndarray, obs: pd.DataFrame) -> np.ndarray:
    z = obs["z_display"].to_numpy(dtype=np.float32)
    return np.column_stack([xy[:, 0], xy[:, 1], z])


def write_outputs(
    raw_slices: list[ad.AnnData],
    stage1_slices: list[ad.AnnData],
    stage2_slices: list[ad.AnnData],
    stage1_key: str,
    stage2_key: str,
) -> None:
    rows = []
    point_frames = []
    combined = []
    for raw, s1, s2 in zip(raw_slices, stage1_slices, stage2_slices):
        sid = raw.obs["slice_id"].iloc[0]
        raw_xy = raw.obsm.get("raw_input_spatial", raw.obsm["spatial"])
        initial_xy = raw.obsm.get("initial_pretransform_spatial", raw.obsm["spatial"])
        coordinate_init_xy = raw.obsm.get("initial_coordinate_spatial", raw.obsm["spatial"])
        s1_xy = s1.obsm[stage1_key]
        s2_xy = s2.obsm[stage2_key]
        s2_rigid_xy = s2.obsm.get(f"{stage2_key}_rigid", s2_xy)
        s2_nonrigid_xy = s2.obsm.get(f"{stage2_key}_nonrigid", s2_xy)
        stage2_delta = s2_nonrigid_xy - s1_xy
        delta_norm = np.sqrt(np.sum(stage2_delta**2, axis=1))
        rows.append(
            {
                "slice_id": sid,
                "n_cells": raw.n_obs,
                "sampling_mode": raw.obs["sampling_mode"].iloc[0],
                "sl_number": int(raw.obs["sl_number"].iloc[0]),
                "z_relative": float(raw.obs["relative_z_by_SL"].iloc[0]),
                "z_display": float(raw.obs["z_display"].iloc[0]),
                "stage1_x_min": float(s1_xy[:, 0].min()),
                "stage1_x_max": float(s1_xy[:, 0].max()),
                "stage1_y_min": float(s1_xy[:, 1].min()),
                "stage1_y_max": float(s1_xy[:, 1].max()),
                "stage2_x_min": float(s2_xy[:, 0].min()),
                "stage2_x_max": float(s2_xy[:, 0].max()),
                "stage2_y_min": float(s2_xy[:, 1].min()),
                "stage2_y_max": float(s2_xy[:, 1].max()),
                "stage2_nonrigid_delta_mean": float(delta_norm.mean()),
                "stage2_nonrigid_delta_max": float(delta_norm.max()),
            }
        )
        point_frame = pd.DataFrame(
            {
                "cell_id": raw.obs_names,
                "slice_id": raw.obs["slice_id"].to_numpy(),
                "stage": raw.obs["stage"].to_numpy(),
                "chip_id": raw.obs["chip_id"].to_numpy(),
                "sl_number": raw.obs["sl_number"].to_numpy(dtype=int),
                "z": raw.obs["relative_z_by_SL"].to_numpy(dtype=np.float32),
                "z_display": raw.obs["z_display"].to_numpy(dtype=np.float32),
                "raw_x": raw_xy[:, 0],
                "raw_y": raw_xy[:, 1],
                "initial_pretransform_x": initial_xy[:, 0],
                "initial_pretransform_y": initial_xy[:, 1],
                "initial_coordinate_x": coordinate_init_xy[:, 0],
                "initial_coordinate_y": coordinate_init_xy[:, 1],
                "stage1_rigid_x": s1_xy[:, 0],
                "stage1_rigid_y": s1_xy[:, 1],
                "stage2_rigid_residual_x": s2_rigid_xy[:, 0],
                "stage2_rigid_residual_y": s2_rigid_xy[:, 1],
                "stage2_nonrigid_x": s2_xy[:, 0],
                "stage2_nonrigid_y": s2_xy[:, 1],
                "stage2_nonrigid_delta": delta_norm,
                "sampling_mode": raw.obs["sampling_mode"].to_numpy(),
            }
        )
        if "celltype" in raw.obs:
            point_frame["celltype"] = raw.obs["celltype"].astype(str).to_numpy()
        point_frames.append(point_frame)

        out = raw.copy()
        out.obsm["raw_spatial"] = raw_xy.astype(np.float32)
        out.obsm["initial_pretransform_spatial"] = initial_xy.astype(np.float32)
        out.obsm["initial_coordinate_spatial"] = coordinate_init_xy.astype(np.float32)
        out.obsm["spateo_stage1_rigid"] = s1_xy.astype(np.float32)
        out.obsm["spateo_stage2_rigid_residual"] = s2_rigid_xy.astype(np.float32)
        out.obsm["spateo_stage2_nonrigid"] = s2_xy.astype(np.float32)
        out.obsm["spateo_stage1_rigid_3D"] = make_3d(s1_xy, raw.obs).astype(np.float32)
        out.obsm["spateo_stage2_nonrigid_3D"] = make_3d(s2_xy, raw.obs).astype(np.float32)
        out.obsm["aligned_spatial_3D"] = (
            out.obsm["spateo_stage2_nonrigid_3D"]
            if stage2_is_enabled()
            else out.obsm["spateo_stage1_rigid_3D"]
        )
        combined.append(out)

    pd.DataFrame(rows).to_csv(OUTDIR / "spateo_two_stage_coordinate_summary.csv", index=False)
    pd.concat(point_frames, ignore_index=True).to_csv(
        OUTDIR / "spateo_two_stage_aligned_points.csv",
        index=False,
    )
    aligned_adata = ad.concat(
        combined,
        join="outer",
        label="slice_key",
        keys=[x.obs["slice_id"].iloc[0] for x in combined],
        index_unique=None,
    )
    aligned_adata.uns["spateo_two_stage"] = {
        "stage1_mode": STAGE1_MODE,
        "stage2_mode": STAGE2_MODE,
        "stage2_skipped": not stage2_is_enabled(),
        "stage1_key": stage1_key,
        "stage2_key": stage2_key,
        "z_display_spacing": Z_DISPLAY_SPACING,
        "use_morpho_align_ref": USE_MORPHO_ALIGN_REF,
        "n_sampling_ref": N_SAMPLING_REF,
    }
    aligned_adata.write_h5ad(OUTDIR / "spateo_two_stage_aligned_lightweight.h5ad")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    if STAGE1_MODE != "SN-S":
        raise SystemExit("STAGE1_MODE should be SN-S for rigid output.")
    if stage2_is_enabled() and STAGE2_MODE != "SN-N":
        raise SystemExit("STAGE2_MODE should be SN-N, none, skip, off, false, or 0.")
    if INITIAL_PRETRANSFORM_RECIPE and INITIAL_COORDINATES_CSV:
        raise SystemExit("Use either INITIAL_PRETRANSFORM_RECIPE or INITIAL_COORDINATES_CSV, not both.")

    from importlib.metadata import version

    import torch
    import ot

    print("spateo_version", version("spateo-release"), flush=True)
    print("POT_version", ot.__version__, flush=True)
    print("torch_version", torch.__version__, flush=True)
    print("torch_cuda_version", torch.version.cuda, flush=True)
    print("torch_cuda_available", torch.cuda.is_available(), flush=True)
    print("torch_cuda_device_count", torch.cuda.device_count(), flush=True)
    if REQUIRE_CUDA and not torch.cuda.is_available():
        raise SystemExit("CUDA is required for this run, but torch.cuda.is_available() is False.")
    device = SPATEO_DEVICE or ("0" if torch.cuda.is_available() else "cpu")
    print("device", device, flush=True)

    order = collect_global_order()
    total_slices = len(order)
    n_slices = N_SLICES if N_SLICES > 0 else total_slices
    start_idx = PAIR_START_ORDER - 1
    selected_order = order.iloc[start_idx : start_idx + n_slices].copy()
    if len(selected_order) != n_slices:
        raise SystemExit(f"Need {n_slices} slices from PAIR_START_ORDER={PAIR_START_ORDER}.")
    selected_order.to_csv(OUTDIR / "spateo_two_stage_slices.csv", index=False)
    print("selected_slices", " -> ".join(selected_order["slice_id"].tolist()), flush=True)

    provenance = {
        "dataset_root": str(ROOT),
        "outdir": str(OUTDIR),
        "pair_start_order": PAIR_START_ORDER,
        "n_slices": n_slices,
        "max_cells_per_slice": MAX_CELLS_PER_SLICE,
        "batch_size": BATCH_SIZE,
        "stage1_mode": STAGE1_MODE,
        "stage1_max_iter": STAGE1_MAX_ITER,
        "stage2_mode": STAGE2_MODE,
        "stage2_max_iter": STAGE2_MAX_ITER,
        "stage2_skipped": not stage2_is_enabled(),
        "use_morpho_align_ref": USE_MORPHO_ALIGN_REF,
        "use_spateo_internal_downsampling": USE_SPATEO_INTERNAL_DOWNSAMPLING,
        "n_sampling_ref": N_SAMPLING_REF,
        "initial_pretransform_recipe": INITIAL_PRETRANSFORM_RECIPE,
        "initial_coordinates_csv": INITIAL_COORDINATES_CSV,
        "initial_coordinate_cell_id_col": INITIAL_COORDINATE_CELL_ID_COL,
        "initial_coordinate_x_col": INITIAL_COORDINATE_X_COL,
        "initial_coordinate_y_col": INITIAL_COORDINATE_Y_COL,
        "allow_partial_initial_coordinates": ALLOW_PARTIAL_INITIAL_COORDINATES,
        "exclude_celltypes_regex": EXCLUDE_CELLTYPES_REGEX,
        "save_full_assignment": SAVE_FULL_ASSIGNMENT,
        "z_display_spacing": Z_DISPLAY_SPACING,
        "slice_order": "filename-derived SL number ascending",
        "method": (
            "Spateo SN-S rigid alignment only."
            if not stage2_is_enabled()
            else "Two-stage Spateo morpho alignment: SN-S rigid coordinates, then SN-N non-rigid from rigid coordinates."
        ),
        "note": "MAX_CELLS_PER_SLICE=0 means complete slices, no subsampling.",
    }
    with (OUTDIR / "spateo_two_stage_provenance.json").open("w") as handle:
        json.dump(provenance, handle, indent=2)

    rng = np.random.default_rng(RNG_SEED)
    raw_slices = [load_slice(row, rng) for _, row in selected_order.iterrows()]
    for adata in raw_slices:
        print(
            "loaded",
            adata.obs["slice_id"].iloc[0],
            "n_cells",
            adata.n_obs,
            "sampling_mode",
            adata.obs["sampling_mode"].iloc[0],
            "spatial_shape",
            adata.obsm["spatial"].shape,
            "x_pca_shape",
            adata.obsm["X_pca"].shape,
            flush=True,
        )
    if INITIAL_PRETRANSFORM_RECIPE:
        pretransform_summary = apply_initial_pretransform(raw_slices, INITIAL_PRETRANSFORM_RECIPE)
        print("initial_pretransform_applied", json.dumps(to_jsonable(pretransform_summary)), flush=True)
    if INITIAL_COORDINATES_CSV:
        coordinate_summary = apply_initial_coordinates_csv(raw_slices, INITIAL_COORDINATES_CSV)
        print("initial_coordinates_applied", json.dumps(to_jsonable(coordinate_summary)), flush=True)

    stage1_key = "spateo_stage1_rigid"
    stage2_key = "spateo_stage2_nonrigid"

    print("stage1_start", STAGE1_MODE, "spatial_key=spatial", flush=True)
    stage1_slices, stage1_pis, stage1_ref, stage1_ref_pis = run_morpho_stage(
        slices=raw_slices,
        spatial_key="spatial",
        key_added=stage1_key,
        vecfld_key="VecFld_stage1_rigid",
        mode=STAGE1_MODE,
        max_iter=STAGE1_MAX_ITER,
        device=device,
    )
    save_stage_audit(stage1_slices, stage1_pis, "stage1_SN-S_rigid", "VecFld_stage1_rigid")
    if stage1_ref is not None and stage1_ref_pis is not None:
        save_stage_audit(stage1_ref, stage1_ref_pis, "stage1_SN-S_rigid_reference", "VecFld_stage1_rigid")

    if stage2_is_enabled():
        stage2_input = [adata.copy() for adata in raw_slices]
        copy_stage_coordinates(stage1_slices, stage2_input, stage1_key)
        print("stage2_start", STAGE2_MODE, "spatial_key=stage1_rigid_spatial", flush=True)
        stage2_slices, stage2_pis, stage2_ref, stage2_ref_pis = run_morpho_stage(
            slices=stage2_input,
            spatial_key="stage1_rigid_spatial",
            key_added=stage2_key,
            vecfld_key="VecFld_stage2_nonrigid",
            mode=STAGE2_MODE,
            max_iter=STAGE2_MAX_ITER,
            device=device,
        )
        save_stage_audit(stage2_slices, stage2_pis, "stage2_SN-N_nonrigid", "VecFld_stage2_nonrigid")
        if stage2_ref is not None and stage2_ref_pis is not None:
            save_stage_audit(stage2_ref, stage2_ref_pis, "stage2_SN-N_nonrigid_reference", "VecFld_stage2_nonrigid")
    else:
        print("stage2_skipped", STAGE2_MODE, "rigid_only_output=stage1", flush=True)
        stage2_slices = []
        for model in stage1_slices:
            copied = model.copy()
            xy = copied.obsm[stage1_key].astype(np.float32)
            copied.obsm[stage2_key] = xy.copy()
            copied.obsm[f"{stage2_key}_rigid"] = xy.copy()
            copied.obsm[f"{stage2_key}_nonrigid"] = xy.copy()
            stage2_slices.append(copied)

    write_outputs(raw_slices, stage1_slices, stage2_slices, stage1_key, stage2_key)
    print("wrote_results", OUTDIR, flush=True)


if __name__ == "__main__":
    main()
