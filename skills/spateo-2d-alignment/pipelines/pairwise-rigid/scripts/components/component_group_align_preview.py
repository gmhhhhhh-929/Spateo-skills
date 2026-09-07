#!/usr/bin/env python3
"""Build non-destructive component-group alignment previews."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import numpy as np
    import pandas as pd
    from PIL import Image, ImageDraw, ImageFont
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "component_group_align_preview.py requires numpy, pandas, and Pillow. "
        "Use the Codex bundled Python runtime if the system Python lacks them."
    ) from exc


def natural_slice(value: Any) -> str:
    text = str(value)
    return text[:-2] if text.endswith(".0") else text


def read_pair(args: argparse.Namespace, fixed_slice: str, moving_slice: str) -> pd.DataFrame:
    wanted = {natural_slice(fixed_slice), natural_slice(moving_slice)}
    labels = pd.read_csv(args.components)
    labels[args.slice_col] = labels[args.slice_col].map(natural_slice)
    labels = labels[labels[args.slice_col].isin(wanted)].copy()
    if args.cell_id_col not in labels.columns:
        raise SystemExit(f"components file lacks {args.cell_id_col}")

    usecols = {args.cell_id_col, args.slice_col, args.x_col, args.y_col}
    if args.z_col:
        usecols.add(args.z_col)
    if args.celltype_col:
        usecols.add(args.celltype_col)
    chunks = []
    for chunk in pd.read_csv(args.coordinates, usecols=lambda c: c in usecols, chunksize=args.chunksize):
        chunk[args.slice_col] = chunk[args.slice_col].map(natural_slice)
        chunk = chunk[chunk[args.slice_col].isin(wanted)]
        if not chunk.empty:
            chunks.append(chunk)
    if not chunks:
        raise SystemExit("No coordinate rows matched fixed/moving slices.")
    coords = pd.concat(chunks, ignore_index=True)
    merged = coords.merge(
        labels.drop(columns=[c for c in [args.x_col, args.y_col, args.z_col] if c in labels.columns], errors="ignore"),
        on=[args.cell_id_col, args.slice_col],
        how="left",
        suffixes=("", "_component"),
    )
    if "component_rank" not in merged.columns and "component_id" not in merged.columns:
        raise SystemExit("components file must contain component_rank or component_id")
    return merged


def nearest_neighbors(query: np.ndarray, ref: np.ndarray, batch_size: int = 384) -> tuple[np.ndarray, np.ndarray]:
    dists = np.empty(len(query), dtype=float)
    idxs = np.empty(len(query), dtype=int)
    for start in range(0, len(query), batch_size):
        q = query[start : start + batch_size]
        diff = q[:, None, :] - ref[None, :, :]
        sq = np.einsum("qri,qri->qr", diff, diff, optimize=True)
        idx = np.argmin(sq, axis=1)
        idxs[start : start + len(q)] = idx
        dists[start : start + len(q)] = np.sqrt(sq[np.arange(len(q)), idx])
    return dists, idxs


def rigid_fit(src: np.ndarray, dst: np.ndarray, allow_reflection: bool = False) -> tuple[np.ndarray, np.ndarray]:
    src_center = src.mean(axis=0)
    dst_center = dst.mean(axis=0)
    x = src - src_center
    y = dst - dst_center
    u, _, vt = np.linalg.svd(x.T @ y)
    r = u @ vt
    if not allow_reflection and np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = u @ vt
    t = dst_center - src_center @ r
    return r, t


def apply_transform(points: np.ndarray, r: np.ndarray, t: np.ndarray) -> np.ndarray:
    return points @ r + t


def transform_summary(r: np.ndarray, t: np.ndarray, residual: np.ndarray | None = None) -> dict[str, Any]:
    out = {
        "rotation_degrees": float(math.degrees(math.atan2(float(r[0, 1]), float(r[0, 0])))),
        "determinant": float(np.linalg.det(r)),
        "r00": float(r[0, 0]),
        "r01": float(r[0, 1]),
        "r10": float(r[1, 0]),
        "r11": float(r[1, 1]),
        "tx": float(t[0]),
        "ty": float(t[1]),
    }
    if residual is not None and len(residual):
        out["fit_residual_median"] = float(np.median(residual))
        out["fit_residual_p95"] = float(np.percentile(residual, 95))
        out["fit_residual_max"] = float(np.max(residual))
    return out


def estimate_icp(
    src: np.ndarray,
    dst: np.ndarray,
    iterations: int,
    trim_fraction: float,
    allow_reflection: bool,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    if len(src) < 3 or len(dst) < 3:
        raise SystemExit("ICP groups need at least 3 points in source and target.")
    r, t = np.eye(2), dst.mean(axis=0) - src.mean(axis=0)
    for _ in range(max(iterations, 1)):
        moved = apply_transform(src, r, t)
        dist, idx = nearest_neighbors(moved, dst)
        cutoff = np.quantile(dist, min(max(trim_fraction, 0.05), 1.0))
        keep = dist <= cutoff
        if keep.sum() < 3:
            keep = np.ones(len(src), dtype=bool)
        r, t = rigid_fit(src[keep], dst[idx[keep]], allow_reflection)
    residual, _ = nearest_neighbors(apply_transform(src, r, t), dst)
    return r, t, transform_summary(r, t, residual)


def estimate_external(points_csv: Path, slice_col: str, moving_slice: str) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    df = pd.read_csv(points_csv)
    if slice_col not in df.columns:
        raise SystemExit(f"external aligned points missing slice column: {slice_col}")
    moving = df[df[slice_col].map(natural_slice).eq(natural_slice(moving_slice))].copy()
    if moving.empty:
        raise SystemExit(f"external aligned points contain no rows for moving slice {moving_slice}")
    required = ["initial_coordinate_x", "initial_coordinate_y", "stage1_rigid_x", "stage1_rigid_y"]
    missing = [col for col in required if col not in moving.columns]
    if missing:
        raise SystemExit(f"external aligned points missing columns: {missing}")
    src = moving[["initial_coordinate_x", "initial_coordinate_y"]].to_numpy(float)
    dst = moving[["stage1_rigid_x", "stage1_rigid_y"]].to_numpy(float)
    r, t = rigid_fit(src, dst, allow_reflection=False)
    residual = np.sqrt(np.sum((apply_transform(src, r, t) - dst) ** 2, axis=1))
    summary = transform_summary(r, t, residual)
    summary["external_aligned_points"] = str(points_csv)
    summary["external_points"] = int(len(src))
    return r, t, summary


def nn_summary(query: np.ndarray, ref: np.ndarray) -> dict[str, float]:
    if len(query) == 0 or len(ref) == 0:
        return {"median": math.nan, "mean": math.nan, "p90": math.nan, "p95": math.nan, "max": math.nan}
    dist, _ = nearest_neighbors(query, ref)
    return {
        "median": float(np.median(dist)),
        "mean": float(np.mean(dist)),
        "p90": float(np.percentile(dist, 90)),
        "p95": float(np.percentile(dist, 95)),
        "max": float(np.max(dist)),
    }


def component_mask(df: pd.DataFrame, slice_col: str, slice_value: str, component_key: str, components: list[Any]) -> pd.Series:
    wanted = {str(x) for x in components}
    return df[slice_col].map(natural_slice).eq(natural_slice(slice_value)) & df[component_key].astype("Int64").astype(str).isin(wanted)


def load_font(size: int) -> ImageFont.ImageFont:
    for candidate in ["/System/Library/Fonts/Supplemental/Arial.ttf", "/System/Library/Fonts/Helvetica.ttc"]:
        try:
            return ImageFont.truetype(candidate, size)
        except Exception:
            pass
    return ImageFont.load_default()


def to_screen(points: np.ndarray, bounds: tuple[float, float, float, float], panel: tuple[int, int, int, int]) -> np.ndarray:
    xmin, xmax, ymin, ymax = bounds
    x0, y0, x1, y1 = panel
    pad = 30
    scale = min((x1 - x0 - 2 * pad) / max(xmax - xmin, 1e-6), (y1 - y0 - 2 * pad) / max(ymax - ymin, 1e-6))
    used_w = (xmax - xmin) * scale
    used_h = (ymax - ymin) * scale
    ox = x0 + (x1 - x0 - used_w) / 2
    oy = y0 + (y1 - y0 - used_h) / 2
    return np.column_stack([ox + (points[:, 0] - xmin) * scale, oy + (ymax - points[:, 1]) * scale])


def draw_points(img: Image.Image, pts: np.ndarray, color: tuple[int, int, int, int], radius: int = 1) -> None:
    layer = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(layer)
    for x, y in pts:
        draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=color)
    img.alpha_composite(layer)


def make_preview(out: pd.DataFrame, args: argparse.Namespace, plan: dict[str, Any], metrics: dict[str, Any]) -> None:
    fixed = out[out[args.slice_col].map(natural_slice).eq(natural_slice(plan["fixed_slice"]))][[args.x_col, args.y_col]].to_numpy(float)
    moving_before = out[out[args.slice_col].map(natural_slice).eq(natural_slice(plan["moving_slice"]))][[args.x_col, args.y_col]].to_numpy(float)
    moving_after = out[out[args.slice_col].map(natural_slice).eq(natural_slice(plan["moving_slice"]))][["aligned_x", "aligned_y"]].to_numpy(float)
    group_labels = out[out[args.slice_col].map(natural_slice).eq(natural_slice(plan["moving_slice"]))]["alignment_group"].astype(str).to_numpy()
    allpts = np.vstack([fixed, moving_before, moving_after])
    mx = (allpts[:, 0].max() - allpts[:, 0].min()) * 0.08
    my = (allpts[:, 1].max() - allpts[:, 1].min()) * 0.08
    bounds = (float(allpts[:, 0].min() - mx), float(allpts[:, 0].max() + mx), float(allpts[:, 1].min() - my), float(allpts[:, 1].max() + my))
    img = Image.new("RGBA", (2400, 1360), (246, 248, 250, 255))
    draw = ImageDraw.Draw(img)
    draw.text((52, 34), f"Component-group preview {plan['moving_slice']} -> {plan['fixed_slice']}", fill=(20, 28, 38), font=load_font(34))
    draw.text((52, 78), "Non-destructive candidate. Fixed slice is blue; moving groups are colored after transform.", fill=(74, 84, 96), font=load_font(21))
    panels = [
        ((42, 126, 1168, 1288), "Before", moving_before, None),
        ((1232, 126, 2358, 1288), "After group-wise transforms", moving_after, group_labels),
    ]
    palette = [(218, 91, 44, 95), (38, 166, 91, 120), (170, 80, 200, 120), (230, 170, 20, 120), (60, 180, 190, 120)]
    for panel, title, moving, labels in panels:
        draw.rectangle(panel, fill=(255, 255, 255, 255), outline=(210, 215, 220, 255), width=2)
        draw.text((panel[0] + 26, panel[1] + 22), title, fill=(20, 28, 38), font=load_font(29))
        plot_box = (panel[0] + 8, panel[1] + 76, panel[2] - 8, panel[3] - 16)
        draw_points(img, to_screen(fixed, bounds, plot_box), (38, 103, 207, 80), 1)
        if labels is None:
            draw_points(img, to_screen(moving, bounds, plot_box), (218, 91, 44, 85), 1)
        else:
            for i, group_id in enumerate(sorted(set(labels))):
                pts = moving[labels == group_id]
                color = (130, 130, 130, 45) if group_id == "ignored" else palette[i % len(palette)]
                draw_points(img, to_screen(pts, bounds, plot_box), color, 1 if group_id != "ignored" else 1)
    lines = []
    for group in metrics["groups"]:
        tr = group["transform"]
        after = group["moving_to_fixed_after"]
        lines.append(f"{group['group_id']}: rot {tr['rotation_degrees']:.2f} deg, t=({tr['tx']:.1f}, {tr['ty']:.1f}); after median {after['median']:.1f}, p95 {after['p95']:.1f}")
    y = 1220
    for line in lines[:4]:
        draw.text((72, y), line, fill=(35, 42, 52), font=load_font(20))
        y += 28
    img.convert("RGB").save(args.output_dir / "component_group_preview.png", quality=95)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coordinates", required=True, type=Path)
    ap.add_argument("--components", required=True, type=Path)
    ap.add_argument("--group-plan", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--slice-col", default="sl_number")
    ap.add_argument("--cell-id-col", default="cell_id")
    ap.add_argument("--x-col", default="manual_x")
    ap.add_argument("--y-col", default="manual_y")
    ap.add_argument("--z-col", default="")
    ap.add_argument("--celltype-col", default="")
    ap.add_argument("--chunksize", type=int, default=250_000)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.overwrite:
        raise SystemExit(f"Output directory exists and is not empty: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plan = json.loads(args.group_plan.read_text())
    fixed_slice = natural_slice(plan["fixed_slice"])
    moving_slice = natural_slice(plan["moving_slice"])
    df = read_pair(args, fixed_slice, moving_slice)
    df["aligned_x"] = df[args.x_col].to_numpy(float)
    df["aligned_y"] = df[args.y_col].to_numpy(float)
    df["alignment_group"] = "ignored"
    metrics: dict[str, Any] = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "coordinates": str(args.coordinates),
        "components": str(args.components),
        "group_plan": str(args.group_plan),
        "fixed_slice": fixed_slice,
        "moving_slice": moving_slice,
        "groups": [],
    }
    recipe_steps = []

    for group in plan["groups"]:
        component_key = group.get("component_key", plan.get("component_key", "component_rank"))
        if component_key not in df.columns:
            raise SystemExit(f"component key not found: {component_key}")
        fixed_est_mask = component_mask(df, args.slice_col, fixed_slice, component_key, group["estimate_from_fixed"])
        moving_est_mask = component_mask(df, args.slice_col, moving_slice, component_key, group["estimate_from_moving"])
        fixed_apply_mask = component_mask(df, args.slice_col, fixed_slice, component_key, group["fixed_components"])
        moving_apply_mask = component_mask(df, args.slice_col, moving_slice, component_key, group["apply_to_moving"])
        method = group.get("method", "icp")
        if method == "icp":
            src = df.loc[moving_est_mask, [args.x_col, args.y_col]].to_numpy(float)
            dst = df.loc[fixed_est_mask, [args.x_col, args.y_col]].to_numpy(float)
            r, t, tr = estimate_icp(
                src,
                dst,
                int(group.get("max_iterations", 60)),
                float(group.get("trim_fraction", 0.9)),
                bool(group.get("allow_reflection", False)),
            )
        elif method == "external":
            if "external_aligned_points" not in group:
                raise SystemExit(f"group {group['group_id']} method external requires external_aligned_points")
            r, t, tr = estimate_external(Path(group["external_aligned_points"]), args.slice_col, moving_slice)
        else:
            raise SystemExit(f"Unsupported method for group {group['group_id']}: {method}")

        pts = df.loc[moving_apply_mask, [args.x_col, args.y_col]].to_numpy(float)
        moved = apply_transform(pts, r, t)
        df.loc[moving_apply_mask, "aligned_x"] = moved[:, 0]
        df.loc[moving_apply_mask, "aligned_y"] = moved[:, 1]
        df.loc[moving_apply_mask, "alignment_group"] = str(group["group_id"])

        fixed_eval = df.loc[fixed_apply_mask, [args.x_col, args.y_col]].to_numpy(float)
        moving_before = df.loc[moving_apply_mask, [args.x_col, args.y_col]].to_numpy(float)
        moving_after = df.loc[moving_apply_mask, ["aligned_x", "aligned_y"]].to_numpy(float)
        group_metrics = {
            "group_id": group["group_id"],
            "method": method,
            "component_key": component_key,
            "fixed_components": group["fixed_components"],
            "moving_components": group["moving_components"],
            "estimate_from_fixed": group["estimate_from_fixed"],
            "estimate_from_moving": group["estimate_from_moving"],
            "apply_to_moving": group["apply_to_moving"],
            "fixed_eval_points": int(len(fixed_eval)),
            "moving_apply_points": int(len(moving_after)),
            "transform": tr,
            "moving_to_fixed_before": nn_summary(moving_before, fixed_eval),
            "moving_to_fixed_after": nn_summary(moving_after, fixed_eval),
        }
        metrics["groups"].append(group_metrics)
        recipe_steps.append(
            {
                "group_id": group["group_id"],
                "operation": "rigid_2d",
                "scope": {"slice": moving_slice, "component_key": component_key, "components": group["apply_to_moving"]},
                "transform": tr,
                "method": method,
            }
        )

    ignored = df[df[args.slice_col].map(natural_slice).eq(moving_slice) & df["alignment_group"].eq("ignored")]
    metrics["ignored_moving_points"] = int(len(ignored))
    df.to_csv(args.output_dir / "component_group_preview_points.csv", index=False)
    (args.output_dir / "component_group_metrics.json").write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    recipe = {
        "operation_type": "component_group_alignment_preview",
        "created_at": metrics["created_at"],
        "fixed_slice": fixed_slice,
        "moving_slice": moving_slice,
        "source_coordinates": str(args.coordinates),
        "source_components": str(args.components),
        "steps": recipe_steps,
        "note": "Candidate recipe only. Confirm visually before composing into final coordinates.",
    }
    (args.output_dir / "component_group_recipe.json").write_text(json.dumps(recipe, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    make_preview(df, args, plan, metrics)
    print(json.dumps({
        "preview_points": str(args.output_dir / "component_group_preview_points.csv"),
        "metrics": str(args.output_dir / "component_group_metrics.json"),
        "recipe": str(args.output_dir / "component_group_recipe.json"),
        "figure": str(args.output_dir / "component_group_preview.png"),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
