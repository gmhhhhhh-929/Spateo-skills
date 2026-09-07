#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from pairwise_core import edge_masks, json_dump, load_json_or_none, natural_key, read_csv_rows, read_point_table, sanitize_token


def load_edges(run_dir: Path, edge_ids: set[str] | None) -> list[dict[str, object]]:
    out = []
    edges_root = run_dir / "edges"
    for edge_dir in sorted(edges_root.iterdir(), key=lambda p: natural_key(p.name)) if edges_root.exists() else []:
        if not edge_dir.is_dir():
            continue
        if edge_ids and edge_dir.name not in edge_ids:
            continue
        obj = load_json_or_none(edge_dir / "edge_transform.json")
        if isinstance(obj, dict):
            out.append(obj)
    return out


def sample(mask: np.ndarray, max_points: int, rng: np.random.Generator) -> np.ndarray:
    ids = np.where(mask)[0]
    if max_points > 0 and len(ids) > max_points:
        ids = np.sort(rng.choice(ids, size=max_points, replace=False))
    return ids


def render_pair_png(path: Path, fixed_xy: np.ndarray, moving_xy: np.ndarray, title: str) -> None:
    width, height = 1100, 1000
    margin = 80
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img, "RGBA")
    all_xy = np.vstack([arr for arr in [fixed_xy, moving_xy] if len(arr)]) if len(fixed_xy) or len(moving_xy) else np.zeros((0, 2))
    if len(all_xy):
        lo = np.nanmin(all_xy, axis=0)
        hi = np.nanmax(all_xy, axis=0)
        span = np.maximum(hi - lo, 1.0)
        scale = min((width - 2 * margin) / span[0], (height - 2 * margin) / span[1])

        def project(xy: np.ndarray) -> np.ndarray:
            out = np.zeros_like(xy)
            out[:, 0] = margin + (xy[:, 0] - lo[0]) * scale
            out[:, 1] = height - margin - (xy[:, 1] - lo[1]) * scale
            return out

        for pts, color in [(project(fixed_xy), (31, 119, 180, 95)), (project(moving_xy), (214, 39, 40, 95))]:
            for x, y in pts:
                draw.ellipse((x - 2, y - 2, x + 2, y + 2), fill=color)
    draw.text((24, 18), title, fill=(20, 20, 20, 255), font=ImageFont.load_default())
    draw.rectangle((24, 52, 44, 72), fill=(31, 119, 180, 160))
    draw.text((52, 55), "fixed", fill=(20, 20, 20, 255), font=ImageFont.load_default())
    draw.rectangle((120, 52, 140, 72), fill=(214, 39, 40, 160))
    draw.text((148, 55), "moving", fill=(20, 20, 20, 255), font=ImageFont.load_default())
    img.save(path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Render visual QC overlays for suspicious pairwise edges/cell types.")
    parser.add_argument("--pairwise-run-dir", required=True)
    parser.add_argument("--points-csv", required=True)
    parser.add_argument("--qc-csv", help="suspicious_celltypes.csv; if omitted, render edge all-cell views.")
    parser.add_argument("--edge-id", action="append")
    parser.add_argument("--celltype", action="append")
    parser.add_argument("--x-col", default="stage1_rigid_x")
    parser.add_argument("--y-col", default="stage1_rigid_y")
    parser.add_argument("--z-col", default="z_display")
    parser.add_argument("--slice-col", default="sl_number")
    parser.add_argument("--celltype-col", default="celltype")
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--max-points-per-layer", type=int, default=6000)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--rng-seed", type=int, default=20260616)
    args = parser.parse_args()

    run_dir = Path(args.pairwise_run_dir).expanduser().resolve()
    outdir = Path(args.output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    table = read_point_table(
        Path(args.points_csv).expanduser().resolve(),
        args.x_col,
        args.y_col,
        args.z_col,
        args.slice_col,
        args.celltype_col,
        args.cell_id_col,
    )
    selected_edges = set(args.edge_id or [])
    selected_celltypes = set(args.celltype or [])
    qc_rows = read_csv_rows(Path(args.qc_csv).expanduser().resolve()) if args.qc_csv else []
    if qc_rows:
        qc_rows = qc_rows[: args.top_k]
        selected_edges.update(str(row.get("edge_id", "")) for row in qc_rows if row.get("edge_id"))
    edges = load_edges(run_dir, selected_edges or None)
    if args.edge_id:
        edges = [edge for edge in edges if str(edge.get("edge_id")) in selected_edges]
    if not edges:
        raise SystemExit("No edges selected for visual QC.")
    rng = np.random.default_rng(args.rng_seed)
    outputs = []

    for edge in edges[: args.top_k]:
        edge_id = str(edge.get("edge_id"))
        fixed_mask, moving_mask = edge_masks(table, edge)
        edge_dir = outdir / sanitize_token(edge_id)
        edge_dir.mkdir(parents=True, exist_ok=True)
        views = [("all_cells", None)]
        for row in qc_rows:
            if str(row.get("edge_id")) == edge_id and row.get("celltype"):
                ct = str(row["celltype"])
                if selected_celltypes and ct not in selected_celltypes:
                    continue
                if ("celltype", ct) not in views:
                    views.append(("celltype", ct))
        if selected_celltypes:
            for ct in selected_celltypes:
                if ("celltype", ct) not in views:
                    views.append(("celltype", ct))

        image_paths = []
        for kind, celltype in views[:8]:
            ct_mask = np.ones(table.n_obs, dtype=bool) if celltype is None else table.celltypes == celltype
            f_ids = sample(fixed_mask & ct_mask, args.max_points_per_layer, rng)
            m_ids = sample(moving_mask & ct_mask, args.max_points_per_layer, rng)
            title = f"{edge_id} · {celltype or 'all cells'}"
            name = f"{sanitize_token(celltype or 'all_cells')}_pair_overlay.png"
            path = edge_dir / name
            render_pair_png(path, table.coords[f_ids, :2], table.coords[m_ids, :2], title)
            image_paths.append(str(path))

        outputs.append({"edge_id": edge_id, "images": image_paths})
        json_dump(
            edge_dir / "visual_assessment.template.json",
            {
                "edge_id": edge_id,
                "images": image_paths,
                "agent_visual_assessment": {
                    "looks_aligned": None,
                    "possible_rotation": None,
                    "possible_translation": None,
                    "split_or_discontinuous": None,
                    "confidence": None,
                    "notes": "",
                },
            },
        )
    json_dump(outdir / "visual_qc_index.json", {"outputs": outputs})
    print(outdir / "visual_qc_index.json")


if __name__ == "__main__":
    main()
