#!/usr/bin/env python3
"""Detect per-slice connected tissue components from spatial coordinates."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from collections import Counter, deque
from datetime import datetime, timezone
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    from PIL import Image, ImageDraw, ImageFont
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "detect_components.py requires numpy, pandas, and Pillow. "
        "Use the Codex bundled Python runtime if the system Python lacks them."
    ) from exc


def parse_slices(value: str | None):
    if not value:
        return None
    out = []
    for item in value.split(","):
        item = item.strip()
        if item:
            out.append(item)
    return set(out)


def estimate_bin_size(df: pd.DataFrame, x_col: str, y_col: str, sample_n: int = 12000) -> float:
    pts = df[[x_col, y_col]].dropna().to_numpy(dtype=float)
    if len(pts) < 2:
        return 1.0
    if len(pts) > sample_n:
        rng = np.random.default_rng(20260624)
        pts = pts[rng.choice(len(pts), size=sample_n, replace=False)]
    scale = np.ptp(pts, axis=0)
    nonzero = scale[scale > 0]
    if len(nonzero) == 0:
        return 1.0
    # Grid hashing gives a robust local spacing estimate without scipy.
    coarse = max(float(np.median(nonzero)) / 160.0, 1.0)
    gx = np.floor(pts[:, 0] / coarse).astype(np.int64)
    gy = np.floor(pts[:, 1] / coarse).astype(np.int64)
    buckets: dict[tuple[int, int], list[int]] = {}
    for i, key in enumerate(zip(gx, gy)):
        buckets.setdefault(key, []).append(i)
    dists = []
    for i, (cx, cy) in enumerate(zip(gx, gy)):
        best = math.inf
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in buckets.get((cx + dx, cy + dy), []):
                    if i == j:
                        continue
                    ddx = pts[i, 0] - pts[j, 0]
                    ddy = pts[i, 1] - pts[j, 1]
                    dist = ddx * ddx + ddy * ddy
                    if dist < best:
                        best = dist
        if math.isfinite(best) and best > 0:
            dists.append(math.sqrt(best))
    if not dists:
        return max(float(np.median(nonzero)) / 200.0, 1.0)
    return max(float(np.median(dists)) * 1.5, 1.0)


def dilate(mask: np.ndarray, radius: int) -> np.ndarray:
    if radius <= 0:
        return mask.copy()
    out = mask.copy()
    coords = [(dx, dy) for dx in range(-radius, radius + 1) for dy in range(-radius, radius + 1) if dx * dx + dy * dy <= radius * radius]
    h, w = mask.shape
    ys, xs = np.nonzero(mask)
    for dx, dy in coords:
        yy = ys + dy
        xx = xs + dx
        valid = (0 <= yy) & (yy < h) & (0 <= xx) & (xx < w)
        out[yy[valid], xx[valid]] = True
    return out


def erode(mask: np.ndarray, radius: int) -> np.ndarray:
    if radius <= 0:
        return mask.copy()
    return ~dilate(~mask, radius)


def close_mask(mask: np.ndarray, radius: int) -> np.ndarray:
    if radius <= 0:
        return mask.copy()
    return erode(dilate(mask, radius), radius)


def label_components(mask: np.ndarray) -> tuple[np.ndarray, int, Counter]:
    labels = np.zeros(mask.shape, dtype=np.int32)
    current = 0
    sizes: Counter[int] = Counter()
    h, w = mask.shape
    for y0, x0 in zip(*np.nonzero(mask)):
        if labels[y0, x0] != 0:
            continue
        current += 1
        q = deque([(int(y0), int(x0))])
        labels[y0, x0] = current
        while q:
            y, x = q.popleft()
            sizes[current] += 1
            for yy in (y - 1, y, y + 1):
                for xx in (x - 1, x, x + 1):
                    if yy == y and xx == x:
                        continue
                    if 0 <= yy < h and 0 <= xx < w and mask[yy, xx] and labels[yy, xx] == 0:
                        labels[yy, xx] = current
                        q.append((yy, xx))
    return labels, current, sizes


def nearest_component_for_empty(labels: np.ndarray, gx: np.ndarray, gy: np.ndarray, max_r: int = 4) -> np.ndarray:
    assigned = labels[gy, gx].astype(np.int32)
    missing_idx = np.where(assigned == 0)[0]
    h, w = labels.shape
    for idx in missing_idx:
        x = int(gx[idx])
        y = int(gy[idx])
        found = 0
        for r in range(1, max_r + 1):
            y0 = max(0, y - r)
            y1 = min(h, y + r + 1)
            x0 = max(0, x - r)
            x1 = min(w, x + r + 1)
            vals = labels[y0:y1, x0:x1]
            nz = vals[vals > 0]
            if len(nz):
                found = int(Counter(nz.tolist()).most_common(1)[0][0])
                break
        assigned[idx] = found
    return assigned


def color_for(idx: int) -> tuple[int, int, int]:
    palette = [
        (230, 25, 75), (60, 180, 75), (0, 130, 200), (245, 130, 48),
        (145, 30, 180), (70, 240, 240), (240, 50, 230), (210, 245, 60),
        (250, 190, 190), (0, 128, 128), (230, 190, 255), (170, 110, 40),
    ]
    return palette[(idx - 1) % len(palette)]


def draw_slice_png(points: pd.DataFrame, x_col: str, y_col: str, label_col: str, path: Path, title: str) -> None:
    w, h = 1600, 1400
    margin = 80
    img = Image.new("RGB", (w, h), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    xmin, xmax = float(points[x_col].min()), float(points[x_col].max())
    ymin, ymax = float(points[y_col].min()), float(points[y_col].max())
    if xmax == xmin:
        xmax += 1
    if ymax == ymin:
        ymax += 1
    sx = (w - 2 * margin) / (xmax - xmin)
    sy = (h - 2 * margin) / (ymax - ymin)
    s = min(sx, sy)
    ox = margin + ((w - 2 * margin) - (xmax - xmin) * s) / 2
    oy = margin + ((h - 2 * margin) - (ymax - ymin) * s) / 2
    draw.text((margin, 25), title, fill=(0, 0, 0), font=font)
    draw.rectangle((margin, margin, w - margin, h - margin), outline=(220, 220, 220))
    xs = (ox + (points[x_col].to_numpy(float) - xmin) * s).astype(int)
    ys = (h - (oy + (points[y_col].to_numpy(float) - ymin) * s)).astype(int)
    labs = points[label_col].to_numpy(int)
    for x, y, lab in zip(xs, ys, labs):
        color = (180, 180, 180) if lab <= 0 else color_for(int(lab))
        draw.point((int(x), int(y)), fill=color)
    legend_y = 50
    rank_lookup = {}
    if "component_rank" in points.columns:
        for lab, rank in zip(points[label_col].to_numpy(int), points["component_rank"].to_numpy(int)):
            rank_lookup.setdefault(int(lab), int(rank))
    for lab, count in Counter(labs.tolist()).most_common(12):
        color = (180, 180, 180) if lab <= 0 else color_for(int(lab))
        draw.rectangle((w - margin - 180, legend_y, w - margin - 165, legend_y + 10), fill=color)
        rank = rank_lookup.get(int(lab), 0)
        label = "noise" if lab <= 0 else f"rank {rank} id {lab}"
        draw.text((w - margin - 158, legend_y - 2), f"{label}: {count}", fill=(0, 0, 0), font=font)
        legend_y += 16
    img.save(path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--slice-col", default="sl_number")
    ap.add_argument("--x-col", default="manual_x")
    ap.add_argument("--y-col", default="manual_y")
    ap.add_argument("--cell-id-col", default="cell_id")
    ap.add_argument("--celltype-col", default="")
    ap.add_argument("--slices", help="Comma-separated slice labels to process.")
    ap.add_argument("--bin-size", default="auto", help="'auto' or numeric coordinate units.")
    ap.add_argument("--close-radius-bins", type=int, default=2)
    ap.add_argument("--min-points", type=int, default=200)
    ap.add_argument("--min-fraction", type=float, default=0.0)
    ap.add_argument(
        "--max-components",
        type=int,
        default=0,
        help="Keep at most this many largest non-noise components per slice. Use 0 to keep all components passing thresholds.",
    )
    ap.add_argument("--chunksize", type=int, default=250_000)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.overwrite:
        raise SystemExit(f"Output directory exists and is not empty: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = args.output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)

    wanted = parse_slices(args.slices)
    usecols = [args.slice_col, args.x_col, args.y_col]
    if args.cell_id_col:
        usecols.append(args.cell_id_col)
    if args.celltype_col:
        usecols.append(args.celltype_col)
    chunks = []
    for chunk in pd.read_csv(args.input, usecols=lambda c: c in set(usecols), chunksize=args.chunksize):
        chunk[args.slice_col] = chunk[args.slice_col].astype(str)
        if wanted is not None:
            chunk = chunk[chunk[args.slice_col].isin(wanted)]
        if not chunk.empty:
            chunks.append(chunk)
    if not chunks:
        raise SystemExit("No rows matched the requested slices.")
    df = pd.concat(chunks, ignore_index=True)
    df["_row_id"] = np.arange(len(df), dtype=np.int64)

    all_labels = []
    summaries = []
    manifest = {
        "version": "detect_spatial_components_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input": str(args.input),
        "parameters": vars(args) | {"input": str(args.input), "output_dir": str(args.output_dir)},
        "slices": [],
    }

    for slice_label, sdf in df.groupby(args.slice_col, sort=True):
        sdf = sdf.reset_index(drop=True)
        bin_size = estimate_bin_size(sdf, args.x_col, args.y_col) if args.bin_size == "auto" else float(args.bin_size)
        xmin = float(sdf[args.x_col].min())
        ymin = float(sdf[args.y_col].min())
        gx = np.floor((sdf[args.x_col].to_numpy(float) - xmin) / bin_size).astype(np.int32)
        gy = np.floor((sdf[args.y_col].to_numpy(float) - ymin) / bin_size).astype(np.int32)
        width = int(gx.max()) + 3
        height = int(gy.max()) + 3
        mask = np.zeros((height, width), dtype=bool)
        mask[gy, gx] = True
        closed = close_mask(mask, args.close_radius_bins)
        labels_grid, _, grid_sizes = label_components(closed)
        point_component = nearest_component_for_empty(labels_grid, gx, gy)
        point_counts = Counter(point_component.tolist())
        min_count_by_fraction = int(math.ceil(len(sdf) * args.min_fraction))
        min_count = max(args.min_points, min_count_by_fraction)
        thresholded = [lab for lab, n in point_counts.most_common() if lab > 0 and n >= min_count]
        if args.max_components and args.max_components > 0:
            thresholded = thresholded[: args.max_components]
        kept = set(thresholded)
        point_component = np.array([lab if lab in kept else 0 for lab in point_component], dtype=np.int32)
        kept_counts = Counter(point_component.tolist())
        ranked = [lab for lab, _ in kept_counts.most_common() if lab > 0]
        rank_map = {lab: i + 1 for i, lab in enumerate(ranked)}
        component_rank = np.array([rank_map.get(int(lab), 0) for lab in point_component], dtype=np.int32)

        out = sdf.copy()
        out["component_id"] = point_component
        out["component_rank"] = component_rank
        out["component_points"] = [kept_counts[int(lab)] for lab in point_component]
        centroid_x = {}
        centroid_y = {}
        for lab in sorted(set(point_component.tolist())):
            if lab <= 0:
                continue
            part = out[out["component_id"] == lab]
            centroid_x[lab] = float(part[args.x_col].mean())
            centroid_y[lab] = float(part[args.y_col].mean())
        out["component_centroid_x"] = [centroid_x.get(int(lab), math.nan) for lab in point_component]
        out["component_centroid_y"] = [centroid_y.get(int(lab), math.nan) for lab in point_component]
        all_labels.append(out)

        for lab in sorted(set(point_component.tolist())):
            part = out[out["component_id"] == lab]
            if part.empty:
                continue
            rec = {
                "slice": slice_label,
                "component_id": int(lab),
                "component_rank": int(rank_map.get(int(lab), 0)),
                "points": int(len(part)),
                "fraction": float(len(part) / len(out)),
                "centroid_x": float(part[args.x_col].mean()),
                "centroid_y": float(part[args.y_col].mean()),
                "min_x": float(part[args.x_col].min()),
                "max_x": float(part[args.x_col].max()),
                "min_y": float(part[args.y_col].min()),
                "max_y": float(part[args.y_col].max()),
                "bin_size": float(bin_size),
                "close_radius_bins": int(args.close_radius_bins),
            }
            if args.celltype_col and args.celltype_col in part.columns and lab > 0:
                rec["top_celltypes"] = ";".join(f"{k}:{v}" for k, v in Counter(part[args.celltype_col].astype(str)).most_common(8))
            summaries.append(rec)

        fig_path = fig_dir / f"slice_{slice_label}_components.png"
        draw_slice_png(out, args.x_col, args.y_col, "component_id", fig_path, f"Slice {slice_label}: detected components")
        manifest["slices"].append(
            {
                "slice": slice_label,
                "points": int(len(sdf)),
                "bin_size": float(bin_size),
                "grid_shape": [int(height), int(width)],
                "components_kept": int(len(ranked)),
                "components_detected_before_filter": int(sum(1 for lab in point_counts if lab > 0)),
                "min_count_threshold": int(min_count),
                "min_count_by_points": int(args.min_points),
                "min_count_by_fraction": int(min_count_by_fraction),
                "max_components_cap": int(args.max_components),
                "noise_points": int(kept_counts.get(0, 0)),
                "figure": str(fig_path),
            }
        )

    labels = pd.concat(all_labels, ignore_index=True).sort_values("_row_id")
    labels = labels.drop(columns=["_row_id"])
    labels.to_csv(args.output_dir / "component_labels.csv", index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / "component_summary.csv", index=False)
    manifest["outputs"] = {
        "component_labels": str(args.output_dir / "component_labels.csv"),
        "component_summary": str(args.output_dir / "component_summary.csv"),
        "figures": str(fig_dir),
    }
    (args.output_dir / "component_detection_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest["outputs"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
