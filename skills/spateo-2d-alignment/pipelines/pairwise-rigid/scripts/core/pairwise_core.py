#!/usr/bin/env python3
from __future__ import annotations

import csv
import datetime as dt
import json
import math
import re
import shlex
from pathlib import Path
from typing import Any, Iterable

import numpy as np


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def natural_key(value: Any) -> list[Any]:
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", str(value))]


def parse_sl_number(value: Any) -> int | None:
    match = re.search(r"SL\s*0*(\d+)|\b0*(\d+)\b", str(value), flags=re.IGNORECASE)
    if not match:
        return None
    return int(match.group(1) or match.group(2))


def slice_token(value: Any) -> str:
    number = parse_sl_number(value)
    return f"SL{number}" if number is not None else sanitize_token(value)


def sanitize_token(value: Any) -> str:
    text = str(value).strip()
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return text.strip("_") or "item"


def shell_quote(value: Any) -> str:
    return shlex.quote(str(value))


def maybe_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def json_dump(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(to_jsonable(obj), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_json_or_none(path: Path) -> Any | None:
    if not path.exists():
        return None
    try:
        return load_json(path)
    except Exception:
        return None


def to_jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(v) for v in value]
    return value


def require_columns(header: Iterable[str], required: Iterable[str]) -> None:
    missing = [col for col in required if col and col not in header]
    if missing:
        raise SystemExit(f"Missing required columns: {', '.join(missing)}")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def write_csv_rows(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = []
        seen: set[str] = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    fieldnames.append(key)
                    seen.add(key)
    if not fieldnames:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def discover_slice_order_from_root(dataset_root: Path) -> list[dict[str, Any]]:
    rows = []
    for fp in sorted(dataset_root.rglob("*.h5ad"), key=lambda p: natural_key(p.name)):
        match = re.search(r"(?P<stage>CS\d+)_SL(?P<sl>\d+)_(?P<chip>Y\w+)\.Spatial", fp.name)
        if not match:
            continue
        rows.append(
            {
                "stage": match.group("stage"),
                "sl_number": int(match.group("sl")),
                "chip_id": match.group("chip"),
                "slice_id": fp.name.split(".Spatial")[0],
                "file": str(fp),
            }
        )
    rows.sort(key=lambda row: (natural_key(row["stage"]), int(row["sl_number"])))
    for i, row in enumerate(rows, start=1):
        row["global_order_by_SL"] = i
    return rows


def read_slice_order_csv(path: Path) -> list[dict[str, Any]]:
    rows = read_csv_rows(path)
    if not rows:
        return []
    require_columns(rows[0].keys(), ["slice_id", "sl_number"])
    out = []
    for i, row in enumerate(rows, start=1):
        sl = parse_sl_number(row.get("sl_number") or row.get("slice_id"))
        if sl is None:
            raise SystemExit(f"Cannot parse SL number from row {i}: {row}")
        item = dict(row)
        item["sl_number"] = sl
        item["stage"] = item.get("stage") or infer_stage(item.get("slice_id", ""))
        item["global_order_by_SL"] = int(item.get("global_order_by_SL") or item.get("order") or i)
        out.append(item)
    out.sort(key=lambda row: (natural_key(row.get("stage", "")), int(row["global_order_by_SL"])))
    return out


def infer_stage(value: Any) -> str:
    match = re.search(r"CS\d+", str(value), flags=re.IGNORECASE)
    return match.group(0).upper() if match else ""


def select_order_rows(
    dataset_root: Path,
    slice_order_csv: Path | None = None,
    stage: str | None = None,
    window_start: str | None = None,
    window_end: str | None = None,
) -> list[dict[str, Any]]:
    rows = read_slice_order_csv(slice_order_csv) if slice_order_csv else discover_slice_order_from_root(dataset_root)
    if stage:
        rows = [row for row in rows if str(row.get("stage", "")).lower() == str(stage).lower()]
    if window_start or window_end:
        start = parse_sl_number(window_start or window_end)
        end = parse_sl_number(window_end or window_start)
        if start is None or end is None:
            raise SystemExit("--window-start/--window-end must include an SL number.")
        lo, hi = sorted([start, end])
        rows = [row for row in rows if lo <= int(row["sl_number"]) <= hi]
    rows.sort(key=lambda row: (natural_key(row.get("stage", "")), int(row["global_order_by_SL"])))
    return rows


def build_adjacent_edges(order_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in order_rows:
        grouped.setdefault(str(row.get("stage", "")), []).append(row)
    edges: list[dict[str, Any]] = []
    used: set[str] = set()
    for stage in sorted(grouped, key=natural_key):
        rows = sorted(grouped[stage], key=lambda row: int(row["global_order_by_SL"]))
        for fixed, moving in zip(rows, rows[1:]):
            base_id = f"{slice_token(fixed['sl_number'])}__{slice_token(moving['sl_number'])}"
            edge_id = base_id if base_id not in used else f"{sanitize_token(stage)}__{base_id}"
            used.add(edge_id)
            edges.append(
                {
                    "edge_id": edge_id,
                    "stage": stage,
                    "fixed_slice_id": fixed.get("slice_id", ""),
                    "moving_slice_id": moving.get("slice_id", ""),
                    "fixed_sl_number": int(fixed["sl_number"]),
                    "moving_sl_number": int(moving["sl_number"]),
                    "fixed_global_order_by_SL": int(fixed["global_order_by_SL"]),
                    "moving_global_order_by_SL": int(moving["global_order_by_SL"]),
                    "pair_start_order": int(fixed["global_order_by_SL"]),
                    "n_slices": 2,
                    "fixed_file": fixed.get("file", ""),
                    "moving_file": moving.get("file", ""),
                }
            )
    return edges


def find_edge_dirs(run_dir: Path | None = None, edges_root: Path | None = None, edge_dir: Path | None = None) -> list[Path]:
    if edge_dir:
        return [edge_dir.expanduser().resolve()]
    roots = []
    if edges_root:
        roots.append(edges_root.expanduser().resolve())
    if run_dir:
        roots.append(run_dir.expanduser().resolve() / "edges")
    edge_dirs: list[Path] = []
    for root in roots:
        if root.exists():
            edge_dirs.extend([p for p in root.iterdir() if p.is_dir()])
    return sorted(edge_dirs, key=lambda p: natural_key(p.name))


def summary_pairs_from_json(obj: Any) -> list[dict[str, Any]]:
    if isinstance(obj, list):
        return [x for x in obj if isinstance(x, dict)]
    if isinstance(obj, dict):
        for key in ["pairs", "edges", "pair_summaries"]:
            if isinstance(obj.get(key), list):
                return [x for x in obj[key] if isinstance(x, dict)]
    return []


def infer_dsub_job_id(edge_dir: Path) -> str:
    ids: set[str] = set()
    for fp in (edge_dir / "logs").glob("*"):
        match = re.search(r"\.(\d+)\.(?:out|err)$", fp.name)
        if match:
            ids.add(match.group(1))
    return ";".join(sorted(ids, key=natural_key))


def find_pair_audit_artifacts(
    edge_dir: Path,
    stage: str,
    provenance: dict[str, Any],
) -> tuple[Path | None, dict[str, Any] | None, Path | None, Path | None]:
    audit_dir = edge_dir / "spateo_pair_audit" / stage
    summary_path = audit_dir / "spateo_pair_audit_summary.json"
    pairs = summary_pairs_from_json(load_json_or_none(summary_path))
    fixed = provenance.get("fixed_slice_id")
    moving = provenance.get("moving_slice_id")
    pair_summary = None
    for item in pairs:
        if (not fixed or item.get("fixed_slice") == fixed) and (not moving or item.get("moving_slice") == moving):
            pair_summary = item
            break
    if pair_summary is None and pairs:
        pair_summary = pairs[0]
    pair_json_path = None
    if audit_dir.exists():
        for path in sorted(audit_dir.glob("pair_*.json"), key=lambda p: natural_key(p.name)):
            if path.name == "spateo_pair_audit_summary.json":
                continue
            obj = load_json_or_none(path)
            if not isinstance(obj, dict):
                continue
            if (not fixed or obj.get("fixed_slice") == fixed) and (not moving or obj.get("moving_slice") == moving):
                pair_json_path = path
                if pair_summary is None:
                    pair_summary = obj
                break
        if pair_json_path is None:
            jsons = [p for p in sorted(audit_dir.glob("pair_*.json"), key=lambda p: natural_key(p.name)) if p.name != "spateo_pair_audit_summary.json"]
            if jsons:
                pair_json_path = jsons[0]
                if pair_summary is None:
                    obj = load_json_or_none(pair_json_path)
                    pair_summary = obj if isinstance(obj, dict) else None
    npz_path = None
    if pair_json_path and pair_json_path.with_suffix(".npz").exists():
        npz_path = pair_json_path.with_suffix(".npz")
    elif audit_dir.exists():
        npzs = sorted(audit_dir.glob("pair_*.npz"), key=lambda p: natural_key(p.name))
        if npzs:
            npz_path = npzs[0]
    return (summary_path if summary_path.exists() else None), pair_summary, pair_json_path, npz_path


def read_transform_npz(path: Path | None) -> tuple[dict[str, Any], list[str]]:
    if path is None:
        return {}, ["missing_npz"]
    try:
        data = np.load(path, allow_pickle=True)
    except Exception as exc:
        return {}, [f"npz_read_error:{exc}"]
    out: dict[str, Any] = {}
    warnings = []
    with data:
        for key in ["R", "t", "optimal_R", "optimal_t", "init_R", "init_t", "normalize_scales", "normalize_means"]:
            if key in data:
                out[key] = np.asarray(data[key]).tolist()
        for key in ["P", "P_sample", "Coff", "inducing_variables"]:
            if key in data:
                out[f"{key}_shape"] = list(np.asarray(data[key]).shape)
    if "R" not in out:
        warnings.append("missing_R")
    if "t" not in out:
        warnings.append("missing_t")
    return out, warnings


def write_edge_transform(edge_dir: Path, stage: str) -> dict[str, Any]:
    provenance_obj = load_json_or_none(edge_dir / "provenance.json")
    provenance = provenance_obj if isinstance(provenance_obj, dict) else {}
    summary_path, pair_summary, pair_json_path, npz_path = find_pair_audit_artifacts(edge_dir, stage, provenance)
    transform, warnings = read_transform_npz(npz_path)
    vecfld = (pair_summary or {}).get("vecfld_summary") or {}
    fixed_slice = (pair_summary or {}).get("fixed_slice") or provenance.get("fixed_slice_id", "")
    moving_slice = (pair_summary or {}).get("moving_slice") or provenance.get("moving_slice_id", "")
    status = "ok"
    if summary_path is None and pair_json_path is None:
        status = "missing_audit"
        warnings.append("missing_audit")
    elif "missing_npz" in warnings:
        status = "partial_missing_npz"
    elif "missing_R" in warnings or "missing_t" in warnings:
        status = "partial_missing_transform"
    dsub = provenance.get("dsub") if isinstance(provenance.get("dsub"), dict) else {}
    job_id = dsub.get("job_id") or infer_dsub_job_id(edge_dir)
    edge_transform = {
        "version": "spateo_pairwise_edge_transform_v1",
        "created_at": now_iso(),
        "status": status,
        "warnings": sorted(set(warnings)),
        "edge_id": provenance.get("edge_id") or edge_dir.name,
        "stage": stage,
        "fixed_slice_id": fixed_slice,
        "moving_slice_id": moving_slice,
        "fixed_sl_number": parse_sl_number(fixed_slice) or provenance.get("fixed_sl_number", ""),
        "moving_sl_number": parse_sl_number(moving_slice) or provenance.get("moving_sl_number", ""),
        "pair_start_order": provenance.get("pair_start_order", ""),
        "n_slices": provenance.get("n_slices", 2),
        "R": transform.get("R"),
        "t": transform.get("t"),
        "optimal_R": transform.get("optimal_R"),
        "optimal_t": transform.get("optimal_t"),
        "init_R": transform.get("init_R"),
        "init_t": transform.get("init_t"),
        "normalize_scales": transform.get("normalize_scales"),
        "normalize_means": transform.get("normalize_means"),
        "sigma2": vecfld.get("sigma2"),
        "gamma": vecfld.get("gamma"),
        "sigma2_variance": vecfld.get("sigma2_variance"),
        "assignment_shape": (pair_summary or {}).get("assignment_shape"),
        "assignment_sample_shape": (pair_summary or {}).get("assignment_sample_shape"),
        "saved_full_assignment": (pair_summary or {}).get("saved_full_assignment"),
        "matrix_shapes": {k: transform.get(k) for k in ["P_shape", "P_sample_shape", "Coff_shape", "inducing_variables_shape"] if k in transform},
        "input_h5ad_path": {"fixed": provenance.get("fixed_file", ""), "moving": provenance.get("moving_file", "")},
        "dataset_root": provenance.get("dataset_root", ""),
        "coordinate_key": provenance.get("coordinate_key", ""),
        "runner_script": provenance.get("runner_script", ""),
        "runner_parameters": provenance.get("runner_parameters", {}),
        "dsub": {
            "project_account": dsub.get("project_account", ""),
            "resource": dsub.get("resource", ""),
            "job_name": dsub.get("job_name", ""),
            "job_id": job_id,
        },
        "paths": {
            "edge_dir": str(edge_dir),
            "audit_summary": str(summary_path) if summary_path else "",
            "pair_audit_json": str(pair_json_path) if pair_json_path else "",
            "pair_audit_npz": str(npz_path) if npz_path else "",
            "runner_points_csv": str(edge_dir / "spateo_two_stage_aligned_points.csv") if (edge_dir / "spateo_two_stage_aligned_points.csv").exists() else "",
            "runner_provenance_json": str(edge_dir / "spateo_two_stage_provenance.json") if (edge_dir / "spateo_two_stage_provenance.json").exists() else "",
            "planner_provenance_json": str(edge_dir / "provenance.json") if (edge_dir / "provenance.json").exists() else "",
        },
    }
    json_dump(edge_dir / "edge_transform.json", edge_transform)
    return edge_transform


def edge_transform_row(item: dict[str, Any]) -> dict[str, Any]:
    shape = item.get("assignment_shape") or []
    return {
        "edge_id": item.get("edge_id", ""),
        "status": item.get("status", ""),
        "warnings": ";".join(item.get("warnings") or []),
        "stage": item.get("stage", ""),
        "fixed_slice_id": item.get("fixed_slice_id", ""),
        "moving_slice_id": item.get("moving_slice_id", ""),
        "fixed_sl_number": item.get("fixed_sl_number", ""),
        "moving_sl_number": item.get("moving_sl_number", ""),
        "pair_start_order": item.get("pair_start_order", ""),
        "sigma2": item.get("sigma2", ""),
        "gamma": item.get("gamma", ""),
        "sigma2_variance": json.dumps(item.get("sigma2_variance", ""), ensure_ascii=False),
        "assignment_rows": shape[0] if len(shape) > 0 else "",
        "assignment_cols": shape[1] if len(shape) > 1 else "",
        "saved_full_assignment": item.get("saved_full_assignment", ""),
        "has_R": item.get("R") is not None,
        "has_t": item.get("t") is not None,
        "has_optimal_R": item.get("optimal_R") is not None,
        "has_optimal_t": item.get("optimal_t") is not None,
        "dsub_job_id": (item.get("dsub") or {}).get("job_id", ""),
        "edge_dir": (item.get("paths") or {}).get("edge_dir", ""),
        "audit_summary": (item.get("paths") or {}).get("audit_summary", ""),
        "pair_audit_npz": (item.get("paths") or {}).get("pair_audit_npz", ""),
        "runner_points_csv": (item.get("paths") or {}).get("runner_points_csv", ""),
    }


class PointTable:
    def __init__(
        self,
        rows: list[dict[str, str]],
        x_col: str,
        y_col: str,
        z_col: str,
        slice_col: str,
        celltype_col: str | None,
        cell_id_col: str | None,
    ) -> None:
        if rows:
            require_columns(rows[0].keys(), [x_col, y_col, z_col, slice_col])
            if celltype_col:
                require_columns(rows[0].keys(), [celltype_col])
            if cell_id_col:
                require_columns(rows[0].keys(), [cell_id_col])
        self.rows = rows
        self.x_col = x_col
        self.y_col = y_col
        self.z_col = z_col
        self.slice_col = slice_col
        self.celltype_col = celltype_col
        self.cell_id_col = cell_id_col
        self.coords = np.array([[float(r[x_col]), float(r[y_col]), float(r[z_col])] for r in rows], dtype=float)
        self.slices = np.array([str(r[slice_col]) for r in rows], dtype=object)
        self.celltypes = np.array([str(r[celltype_col]) for r in rows], dtype=object) if celltype_col else np.array(["all"] * len(rows), dtype=object)
        if cell_id_col:
            self.cell_ids = np.array([str(r[cell_id_col]) for r in rows], dtype=object)
        elif rows and "cell_id" in rows[0]:
            self.cell_ids = np.array([str(r["cell_id"]) for r in rows], dtype=object)
        else:
            self.cell_ids = np.array([str(i) for i in range(len(rows))], dtype=object)

    @property
    def n_obs(self) -> int:
        return len(self.rows)

    @property
    def slice_labels(self) -> list[str]:
        return sorted({str(x) for x in self.slices}, key=natural_key)

    @property
    def celltype_labels(self) -> list[str]:
        return sorted({str(x) for x in self.celltypes}, key=natural_key)


def read_point_table(
    path: Path,
    x_col: str,
    y_col: str,
    z_col: str,
    slice_col: str,
    celltype_col: str | None = None,
    cell_id_col: str | None = None,
) -> PointTable:
    return PointTable(read_csv_rows(path), x_col, y_col, z_col, slice_col, celltype_col, cell_id_col)


def slice_mask(table: PointTable, slice_value: Any) -> np.ndarray:
    wanted = {str(slice_value)}
    number = parse_sl_number(slice_value)
    if number is not None:
        wanted.add(str(number))
        wanted.add(f"SL{number}")
    out = []
    for value in table.slices:
        sl = parse_sl_number(value)
        out.append(str(value) in wanted or (number is not None and sl == number))
    return np.array(out, dtype=bool)


def edge_masks(table: PointTable, edge: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    fixed = edge.get("fixed_sl_number") or edge.get("fixed_slice_id") or edge.get("fixed_slice")
    moving = edge.get("moving_sl_number") or edge.get("moving_slice_id") or edge.get("moving_slice")
    return slice_mask(table, fixed), slice_mask(table, moving)


def centroid(points: np.ndarray) -> np.ndarray:
    return points.mean(axis=0) if len(points) else np.array([np.nan, np.nan, np.nan], dtype=float)


def radius(points: np.ndarray) -> float:
    if len(points) < 2:
        return float("nan")
    c = centroid(points)
    return float(np.median(np.linalg.norm(points[:, :2] - c[:2], axis=1)))


def principal_angle(points: np.ndarray) -> float:
    if len(points) < 3:
        return float("nan")
    xy = points[:, :2] - points[:, :2].mean(axis=0)
    cov = np.cov(xy.T)
    vals, vecs = np.linalg.eigh(cov)
    vec = vecs[:, int(np.argmax(vals))]
    return float(math.degrees(math.atan2(vec[1], vec[0])))


def angle_delta(a: float, b: float) -> float:
    if not np.isfinite(a) or not np.isfinite(b):
        return float("nan")
    delta = abs((a - b + 90.0) % 180.0 - 90.0)
    return float(delta)


def pair_metrics(fixed_points: np.ndarray, moving_points: np.ndarray) -> dict[str, Any]:
    cf = centroid(fixed_points)
    cm = centroid(moving_points)
    rf = radius(fixed_points)
    rm = radius(moving_points)
    distance = float(np.linalg.norm(cm[:2] - cf[:2])) if len(fixed_points) and len(moving_points) else float("nan")
    finite_radii = [x for x in [rf, rm] if np.isfinite(x)]
    denom = float(np.median(finite_radii)) if finite_radii else float("nan")
    normalized = float(distance / max(denom, 1e-6)) if np.isfinite(distance) and np.isfinite(denom) else float("nan")
    af = principal_angle(fixed_points)
    am = principal_angle(moving_points)
    ratio = float(rm / rf) if np.isfinite(rf) and rf > 1e-6 and np.isfinite(rm) else float("nan")
    return {
        "fixed_centroid_x": cf[0],
        "fixed_centroid_y": cf[1],
        "moving_centroid_x": cm[0],
        "moving_centroid_y": cm[1],
        "centroid_dx": cm[0] - cf[0],
        "centroid_dy": cm[1] - cf[1],
        "centroid_distance": distance,
        "fixed_radius": rf,
        "moving_radius": rm,
        "normalized_centroid_jump": normalized,
        "fixed_axis_angle_deg": af,
        "moving_axis_angle_deg": am,
        "axis_angle_delta_deg": angle_delta(af, am),
        "radius_ratio": ratio,
    }


def rigid_transform_about_center(points: np.ndarray, degrees: float, dx: float, dy: float, center: np.ndarray) -> np.ndarray:
    out = points.copy()
    theta = math.radians(degrees)
    rot = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]], dtype=float)
    out[:, :2] = (out[:, :2] - center[:2]) @ rot.T + center[:2]
    out[:, 0] += dx
    out[:, 1] += dy
    return out


def estimate_rigid_2d(source: np.ndarray, target: np.ndarray) -> tuple[float, float, float]:
    if len(source) < 2 or len(target) < 2:
        return 0.0, 0.0, 0.0
    n = min(len(source), len(target))
    src = source[:n, :2]
    dst = target[:n, :2]
    cs = src.mean(axis=0)
    cd = dst.mean(axis=0)
    xs = src - cs
    xd = dst - cd
    h = xs.T @ xd
    u, _, vt = np.linalg.svd(h)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    theta = math.degrees(math.atan2(r[1, 0], r[0, 0]))
    t = cd - cs @ r.T
    return float(theta), float(t[0]), float(t[1])


def estimate_cloud_axis_rigid_2d(source: np.ndarray, target: np.ndarray) -> tuple[float, float, float]:
    if len(source) < 2 or len(target) < 2:
        return 0.0, 0.0, 0.0
    cs = source[:, :2].mean(axis=0)
    ct = target[:, :2].mean(axis=0)
    a_source = principal_angle(source)
    a_target = principal_angle(target)
    if np.isfinite(a_source) and np.isfinite(a_target):
        degrees = (a_target - a_source + 90.0) % 180.0 - 90.0
    else:
        degrees = 0.0
    delta = ct - cs
    return float(degrees), float(delta[0]), float(delta[1])


def sample_points(points: np.ndarray, max_points: int, rng: np.random.Generator) -> np.ndarray:
    if max_points <= 0 or len(points) <= max_points:
        return points
    idx = np.sort(rng.choice(np.arange(len(points)), size=max_points, replace=False))
    return points[idx]
