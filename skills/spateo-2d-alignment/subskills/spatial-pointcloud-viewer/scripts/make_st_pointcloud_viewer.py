#!/usr/bin/env python3
import argparse
import base64
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def hex_color(value: str):
    value = value.lstrip("#")
    return [int(value[i : i + 2], 16) / 255.0 for i in (0, 2, 4)]


def hsv_color(h: float, sat: float, val: float):
    h = h % 1.0
    i = int(h * 6.0)
    f = h * 6.0 - i
    p = val * (1.0 - sat)
    q = val * (1.0 - f * sat)
    t = val * (1.0 - (1.0 - f) * sat)
    i %= 6
    if i == 0:
        return [val, t, p]
    if i == 1:
        return [q, val, p]
    if i == 2:
        return [p, val, t]
    if i == 3:
        return [p, q, val]
    if i == 4:
        return [t, p, val]
    return [val, p, q]


def srgb_to_linear(x: float) -> float:
    return x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4


def rgb_to_lab(rgb):
    r, g, b = [srgb_to_linear(float(x)) for x in rgb]
    x = r * 0.4124564 + g * 0.3575761 + b * 0.1804375
    y = r * 0.2126729 + g * 0.7151522 + b * 0.0721750
    z = r * 0.0193339 + g * 0.1191920 + b * 0.9503041
    x /= 0.95047
    z /= 1.08883

    def f(t: float) -> float:
        return t ** (1.0 / 3.0) if t > 0.008856 else 7.787 * t + 16.0 / 116.0

    fx, fy, fz = f(x), f(y), f(z)
    return np.asarray([116.0 * fy - 16.0, 500.0 * (fx - fy), 200.0 * (fy - fz)], dtype=np.float32)


def relative_luminance(rgb) -> float:
    r, g, b = [srgb_to_linear(float(x)) for x in rgb]
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def distinct_palette(count: int):
    if count <= 0:
        return []
    seed = hex_color("#00A7FF")
    candidates = []
    for hue in range(0, 360, 3):
        for sat in (0.55, 0.70, 0.85, 1.0):
            for val in (0.62, 0.76, 0.90, 1.0):
                rgb = hsv_color(hue / 360.0, sat, val)
                if relative_luminance(rgb) >= 0.08:
                    candidates.append(rgb)
    candidate_rgb = np.asarray(candidates, dtype=np.float32)
    candidate_lab = np.asarray([rgb_to_lab(c) for c in candidate_rgb], dtype=np.float32)
    selected = [seed]
    min_dist = np.linalg.norm(candidate_lab - rgb_to_lab(seed), axis=1)
    for _ in range(1, count):
        idx = int(np.argmax(min_dist))
        selected.append(candidate_rgb[idx].astype(float).tolist())
        dist = np.linalg.norm(candidate_lab - candidate_lab[idx], axis=1)
        min_dist = np.minimum(min_dist, dist)
    return selected


def optional_column(fieldnames, requested, fallbacks):
    if requested and str(requested).lower() not in {"auto", "none", "off", ""}:
        return requested if requested in fieldnames else None
    if requested and str(requested).lower() in {"none", "off", ""}:
        return None
    for col in fallbacks:
        if col in fieldnames:
            return col
    return None


def read_csv_points(path, x_col, y_col, z_col, celltype_col, slice_col, component_col, component_rank_col):
    points = []
    celltypes = []
    slices = []
    components = []
    component_ranks = []
    with Path(path).open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        missing = [c for c in [x_col, y_col, z_col, celltype_col] if c not in fieldnames]
        if missing:
            raise SystemExit(f"Missing CSV columns: {', '.join(missing)}")
        slice_enabled = bool(slice_col and slice_col in fieldnames)
        resolved_component_col = optional_column(fieldnames, component_col, ["sample_component_id", "component_id"])
        resolved_component_rank_col = optional_column(fieldnames, component_rank_col, ["sample_component_rank", "component_rank"])
        for row in reader:
            points.append((float(row[x_col]), float(row[y_col]), float(row[z_col])))
            celltypes.append(str(row[celltype_col]))
            if slice_enabled:
                slices.append(str(row[slice_col]))
            if resolved_component_col:
                components.append(str(row[resolved_component_col]))
                component_ranks.append(str(row[resolved_component_rank_col]) if resolved_component_rank_col else "")
    slice_values = np.asarray(slices, dtype=object) if slices else None
    component_values = np.asarray(components, dtype=object) if components else None
    component_rank_values = np.asarray(component_ranks, dtype=object) if component_ranks else None
    resolved = {
        "component_col": resolved_component_col,
        "component_rank_col": resolved_component_rank_col,
    }
    return np.asarray(points, dtype=np.float32), np.asarray(celltypes, dtype=object), slice_values, component_values, component_rank_values, resolved


def read_h5ad_points(path, obsm_key, celltype_col, slice_col, component_col, component_rank_col):
    try:
        import anndata as ad  # type: ignore
    except Exception as exc:
        raise SystemExit("Reading h5ad requires anndata in the active Python environment.") from exc
    adata = ad.read_h5ad(path)
    if obsm_key not in adata.obsm:
        raise SystemExit(f"Missing obsm key: {obsm_key}")
    if celltype_col not in adata.obs:
        raise SystemExit(f"Missing obs column: {celltype_col}")
    coords = np.asarray(adata.obsm[obsm_key], dtype=np.float32)
    if coords.ndim != 2 or coords.shape[1] < 3:
        raise SystemExit(f"obsm[{obsm_key!r}] must have at least three columns")
    celltypes = adata.obs[celltype_col].astype(str).to_numpy(dtype=object)
    slice_values = None
    if slice_col and slice_col in adata.obs:
        slice_values = adata.obs[slice_col].astype(str).to_numpy(dtype=object)
    obs_cols = set(map(str, adata.obs.columns))
    resolved_component_col = optional_column(obs_cols, component_col, ["sample_component_id", "component_id"])
    resolved_component_rank_col = optional_column(obs_cols, component_rank_col, ["sample_component_rank", "component_rank"])
    component_values = adata.obs[resolved_component_col].astype(str).to_numpy(dtype=object) if resolved_component_col else None
    component_rank_values = adata.obs[resolved_component_rank_col].astype(str).to_numpy(dtype=object) if resolved_component_rank_col else None
    resolved = {
        "component_col": resolved_component_col,
        "component_rank_col": resolved_component_rank_col,
    }
    return coords[:, :3].astype(np.float32, copy=False), celltypes, slice_values, component_values, component_rank_values, resolved


def downsample_indices(celltypes, max_points, seed):
    if max_points <= 0 or celltypes.size <= max_points:
        return np.arange(celltypes.size, dtype=np.int64)
    rng = np.random.default_rng(seed)
    order = ordered_unique(celltypes)
    counts = np.asarray([np.sum(celltypes == ct) for ct in order], dtype=float)
    weights = np.sqrt(np.maximum(counts, 1.0))
    targets = np.maximum(1, np.floor(max_points * weights / weights.sum()).astype(int))
    diff = max_points - int(targets.sum())
    if diff > 0:
        for i in np.argsort(-counts)[:diff]:
            targets[i] += 1
    elif diff < 0:
        for i in np.argsort(-targets):
            if diff == 0:
                break
            take = min(targets[i] - 1, -diff)
            targets[i] -= take
            diff += take
    selected = []
    for ct, target in zip(order, targets):
        idx = np.where(celltypes == ct)[0]
        if idx.size > target:
            idx = rng.choice(idx, size=int(target), replace=False)
        selected.append(idx)
    return np.sort(np.concatenate(selected).astype(np.int64))


def ordered_unique(values):
    seen = set()
    out = []
    for value in values.tolist():
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


_SLICE_SORT_RE = re.compile(r"-?\d+")


def slice_sort_key(label):
    text = str(label)
    match = _SLICE_SORT_RE.search(text)
    if match:
        try:
            return (0, int(match.group(0)), text)
        except ValueError:
            pass
    return (1, text)


def numeric_text(value):
    text = str(value).strip()
    try:
        number = float(text)
    except ValueError:
        return text
    if number.is_integer():
        return str(int(number))
    return text


def slice_label_text(value):
    text = numeric_text(value)
    return text if str(text).upper().startswith("SL") else f"SL{text}"


def is_noise_component(component, rank) -> bool:
    for value in (component, rank):
        text = str(value).strip().lower()
        if text in {"", "0", "0.0", "nan", "none", "noise"}:
            return True
    return False


def component_display_labels(components, component_ranks, slices):
    if components is None:
        return None
    labels = []
    ranks = component_ranks if component_ranks is not None else np.asarray([""] * len(components), dtype=object)
    slice_values = slices if slices is not None else np.asarray([""] * len(components), dtype=object)
    for comp, rank, sl in zip(components.tolist(), ranks.tolist(), slice_values.tolist()):
        prefix = f"{slice_label_text(sl)} " if str(sl).strip() else ""
        comp_text = numeric_text(comp)
        rank_text = numeric_text(rank)
        if is_noise_component(comp_text, rank_text):
            labels.append(f"{prefix}noise".strip())
        elif rank_text:
            labels.append(f"{prefix}rank{rank_text} id{comp_text}".strip())
        else:
            labels.append(f"{prefix}component id{comp_text}".strip())
    return np.asarray(labels, dtype=object)


_COMPONENT_SORT_RE = re.compile(r"(?:SL)?(-?\d+).*?rank(-?\d+).*?id(.+)$", re.IGNORECASE)


def component_sort_key(label):
    text = str(label)
    if "noise" in text.lower():
        slice_match = _SLICE_SORT_RE.search(text)
        return (0, int(slice_match.group(0)) if slice_match else 0, 10**9, text)
    match = _COMPONENT_SORT_RE.search(text)
    if match:
        try:
            return (0, int(match.group(1)), int(match.group(2)), match.group(3))
        except ValueError:
            return (0, match.group(1), match.group(2), match.group(3))
    return (1, text)


def display_policy_for(args) -> str:
    if args.component_balanced:
        return "balanced300k"
    if args.max_points != 0:
        return "debug_random_sample"
    kind = str(args.viewer_kind).lower()
    if "review_full_points" in kind or "candidate_full_points" in kind or "manual_full_points" in kind:
        return "manual_full_points"
    return "debug_all_points"


def slices_from_z(points):
    z = points[:, 2]
    unique = sorted(float(v) for v in np.unique(z))
    mapping = {value: i + 1 for i, value in enumerate(unique)}
    labels = [f"z{i:03d}" for i in range(1, len(unique) + 1)]
    label_by_z = {value: label for value, label in zip(unique, labels)}
    return np.asarray([label_by_z[float(v)] for v in z], dtype=object), mapping


def bounds_meta(points):
    pmin = points.min(axis=0)
    pmax = points.max(axis=0)
    center = (pmin + pmax) / 2.0
    span = np.maximum(pmax - pmin, 1e-6)
    return {
        "min": [float(x) for x in pmin],
        "max": [float(x) for x in pmax],
        "center": [float(x) for x in center],
        "span": [float(x) for x in span],
        "max_span": float(span.max()),
    }


def encode_f32(points):
    arr = points.astype("<f4", copy=False)
    return base64.b64encode(arr.tobytes()).decode("ascii")


def encode_index(values):
    return base64.b64encode(values.tobytes()).decode("ascii")


def label_ids(values, groups):
    name_to_id = {g["name"]: i for i, g in enumerate(groups)}
    dtype = np.uint32 if len(groups) > np.iinfo(np.uint16).max else np.uint16
    ids = np.asarray([name_to_id[str(v)] for v in values.tolist()], dtype=dtype)
    return ids, "uint32" if dtype == np.uint32 else "uint16"


def build_mode_groups(values, *, sort_key=None, points=None):
    if values is None:
        return None, None, None
    order = ordered_unique(values)
    if sort_key is not None:
        order = sorted(order, key=sort_key)
    if not order:
        return np.asarray([], dtype=np.uint16), [], "uint16"
    index_dtype = np.uint32 if values.shape[0] > np.iinfo(np.uint16).max else np.uint16
    palette = distinct_palette(len(order))
    parts = []
    groups = []
    offset = 0
    for i, label in enumerate(order):
        idx = np.where(values == label)[0].astype(index_dtype, copy=False)
        parts.append(idx)
        color = palette[i]
        groups.append(
            {
                "name": str(label),
                "count": int(idx.size),
                "offset": int(offset),
                "color": [round(float(c), 6) for c in color],
            }
        )
        if points is not None and idx.size:
            centroid = np.mean(points[idx], axis=0)
            groups[-1]["centroid"] = [round(float(v), 6) for v in centroid.tolist()]
        offset += int(idx.size)
    indices = np.concatenate(parts).astype(index_dtype, copy=False) if parts else np.asarray([], dtype=index_dtype)
    return indices, groups, "uint32" if index_dtype == np.uint32 else "uint16"


def render_html(
    meta,
    points_b64,
    celltype_index_b64,
    slice_index_b64,
    component_index_b64,
    point_celltype_id_b64,
    point_slice_id_b64,
    point_component_id_b64,
):
    meta_json = json.dumps(meta, ensure_ascii=False, separators=(",", ":"))
    doc = r"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
:root { color-scheme: dark; --bg:#080b10; --panel:#111823; --line:rgba(255,255,255,.13); --text:#edf3ff; --muted:#9ca8ba; }
* { box-sizing:border-box; }
html,body { margin:0;height:100%;overflow:hidden;background:var(--bg);color:var(--text);font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }
#app { height:100vh;display:grid;grid-template-columns:minmax(320px,380px) 1fr;min-height:0; }
aside { min-height:0;overflow:hidden;background:var(--panel);border-right:1px solid var(--line);display:flex;flex-direction:column; }
header { padding:14px;border-bottom:1px solid var(--line); }
h1 { margin:0;font-size:17px;letter-spacing:0; }
.toolbar { display:grid;grid-template-columns:1fr 1fr;gap:8px;padding:10px 14px;border-bottom:1px solid var(--line); }
.wide { grid-column:1 / -1; }
button { border:1px solid var(--line);background:#151e2a;color:var(--text);border-radius:7px;min-height:32px;padding:0 9px;font:inherit;font-size:12px;cursor:pointer; }
button:hover { border-color:rgba(255,255,255,.36);background:#1c2a3b; }
button.active { border-color:rgba(80,175,255,.82);background:#18334a; }
button:disabled { opacity:.42;cursor:not-allowed; }
select { min-width:0;width:100%;height:32px;border:1px solid var(--line);border-radius:7px;background:#0d1420;color:var(--text);padding:0 8px;font:inherit;font-size:12px; }
.pairTools { display:grid;grid-template-columns:54px minmax(0,1fr) minmax(0,1fr) 54px;gap:8px;align-items:center; }
.pairTools[hidden] { display:none; }
.pairMeta { grid-column:1 / -1;color:var(--muted);font-size:11px;min-height:15px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap; }
.slider { display:grid;grid-template-columns:auto 1fr auto;align-items:center;gap:8px;padding:3px 2px;color:var(--muted);font-size:12px; }
input[type=range] { width:100%; }
.section { flex:1;min-height:0;border-bottom:1px solid var(--line);display:flex;flex-direction:column; }
.sectionTitle { display:flex;align-items:center;justify-content:space-between;gap:10px;padding:8px 14px 2px;color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:0; }
.panel { flex:1;min-height:0;overflow:auto;padding:8px 10px 14px; }
.row { display:grid;grid-template-columns:18px 14px 1fr auto;gap:7px;align-items:center;min-height:30px;padding:3px 5px;border-radius:6px;font-size:12px; }
.row:hover { background:rgba(255,255,255,.06); }
.row.disabled { opacity:.38; }
.sw { width:12px;height:12px;border-radius:3px;border:1px solid rgba(255,255,255,.38); }
.name { overflow:hidden;text-overflow:ellipsis;white-space:nowrap; }
.num { color:var(--muted);font-variant-numeric:tabular-nums; }
main { position:relative;min-width:0;min-height:0; }
canvas { width:100%;height:100%;display:block;background:radial-gradient(circle at 50% 42%,#142033 0,#080b10 58%,#040609 100%);cursor:grab; }
canvas.dragging { cursor:grabbing; }
#hud,#status { position:absolute;background:rgba(5,8,12,.68);border:1px solid rgba(255,255,255,.12);border-radius:8px;color:#dce7f7;font-size:12px;line-height:1.45;backdrop-filter:blur(8px); }
#hud { left:14px;bottom:14px;max-width:min(740px,calc(100% - 28px));padding:9px 10px; }
#status { right:14px;top:14px;border-radius:999px;padding:7px 10px; }
#labelLayer { position:absolute;inset:0;pointer-events:none;overflow:hidden; }
.componentLabel { position:absolute;transform:translate(-50%,-50%);max-width:150px;padding:3px 6px;border:1px solid rgba(255,255,255,.45);border-radius:6px;background:rgba(5,8,12,.74);color:#f5f8ff;font-size:11px;line-height:1.2;white-space:nowrap;text-overflow:ellipsis;overflow:hidden;box-shadow:0 2px 10px rgba(0,0,0,.25); }
.componentLabel::before { content:"";display:inline-block;width:8px;height:8px;margin-right:5px;border-radius:2px;background:var(--label-color);vertical-align:0; }
@media (max-width:760px) {
  #app { grid-template-columns:1fr;grid-template-rows:46vh 54vh; }
  aside { order:2;border-right:0;border-top:1px solid var(--line); }
  main { order:1; }
}
</style>
</head>
<body>
<div id="app">
<aside>
<header><h1>__TITLE__</h1></header>
<div class="toolbar">
<button id="overallViewBtn" class="active">Overall</button>
<button id="pairViewBtn">Pair</button>
<div id="pairTools" class="pairTools wide" hidden>
<button id="prevPairBtn">Prev</button>
<select id="pairASelect"></select>
<select id="pairBSelect"></select>
<button id="nextPairBtn">Next</button>
<div id="pairMeta" class="pairMeta"></div>
</div>
<button id="cellModeBtn" class="active">Cell type</button>
<button id="sliceModeBtn">Slice</button>
<button id="componentModeBtn">Component</button>
<button id="labelBtn">Labels</button>
<button id="allBtn">All on</button>
<button id="noneBtn">All off</button>
<button id="resetBtn">Reset view</button>
<button id="invertBtn">Invert</button>
<div class="slider wide"><span>Point</span><input id="pointSlider" type="range" min="1" max="8" value="3" step="0.5"><span id="pointVal">3.0</span></div>
</div>
<div class="section">
<div class="sectionTitle"><span>Cell type</span><span id="cellCount"></span></div>
<div id="cellPanel" class="panel"></div>
</div>
<div class="section" id="sliceSection">
<div class="sectionTitle"><span>Slice</span><span id="sliceCount"></span></div>
<div id="slicePanel" class="panel"></div>
</div>
<div class="section" id="componentSection">
<div class="sectionTitle"><span>Component</span><span id="componentCount"></span></div>
<div id="componentPanel" class="panel"></div>
</div>
</aside>
<main>
<canvas id="gl"></canvas>
<div id="labelLayer"></div>
<div id="status">Loading...</div>
<div id="hud"></div>
</main>
</div>
<script>
const META = __META_JSON__;
const POINTS_B64 = __POINTS_B64__;
const CELLTYPE_INDEX_B64 = __CELLTYPE_INDEX_B64__;
const SLICE_INDEX_B64 = __SLICE_INDEX_B64__;
const COMPONENT_INDEX_B64 = __COMPONENT_INDEX_B64__;
const POINT_CELLTYPE_ID_B64 = __POINT_CELLTYPE_ID_B64__;
const POINT_SLICE_ID_B64 = __POINT_SLICE_ID_B64__;
const POINT_COMPONENT_ID_B64 = __POINT_COMPONENT_ID_B64__;
function decodeBytes(b64) {
  if (!b64) return new ArrayBuffer(0);
  const bin = atob(b64);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return bytes.buffer;
}
function decodeF32(b64) {
  return new Float32Array(decodeBytes(b64));
}
function decodeIndex(b64, dtype) {
  if (!b64) return null;
  const buf = decodeBytes(b64);
  return dtype === "uint32" ? new Uint32Array(buf) : new Uint16Array(buf);
}
function css(c) { return `rgb(${Math.round(c[0]*255)},${Math.round(c[1]*255)},${Math.round(c[2]*255)})`; }
function esc(s) { return String(s).replaceAll('&','&amp;').replaceAll('"','&quot;').replaceAll('<','&lt;').replaceAll('>','&gt;'); }
const points = decodeF32(POINTS_B64);
const modes = {
  celltype: {
    label: "Cell type",
    groups: META.celltypes.map(x => ({...x})),
    index: decodeIndex(CELLTYPE_INDEX_B64, META.celltype_index_dtype),
    pointIds: decodeIndex(POINT_CELLTYPE_ID_B64, META.point_celltype_id_dtype),
    visible: new Map(META.celltypes.map(c => [c.name, true])),
  },
  slice: {
    label: "Slice",
    groups: (META.slices || []).map(x => ({...x})),
    index: decodeIndex(SLICE_INDEX_B64, META.slice_index_dtype),
    pointIds: decodeIndex(POINT_SLICE_ID_B64, META.point_slice_id_dtype),
    visible: new Map((META.slices || []).map(c => [c.name, true])),
  },
  component: {
    label: "Component",
    groups: (META.components || []).map(x => ({...x})),
    index: decodeIndex(COMPONENT_INDEX_B64, META.component_index_dtype),
    pointIds: decodeIndex(POINT_COMPONENT_ID_B64, META.point_component_id_dtype),
    visible: new Map((META.components || []).map(c => [c.name, true])),
  },
};
let colorMode = "celltype";
let viewMode = "overall";
let pairA = 0;
let pairB = modes.slice.groups.length > 1 ? 1 : 0;
let yaw = -0.7, pitch = -0.22, zoom = 3.65, pointSize = 3.0 * (devicePixelRatio || 1);
let showComponentLabels = true;
const canvas = document.getElementById("gl");
const statusEl = document.getElementById("status");
const hud = document.getElementById("hud");
const labelLayer = document.getElementById("labelLayer");
const cellPanel = document.getElementById("cellPanel");
const slicePanel = document.getElementById("slicePanel");
const componentPanel = document.getElementById("componentPanel");
const cellCount = document.getElementById("cellCount");
const sliceCount = document.getElementById("sliceCount");
const componentCount = document.getElementById("componentCount");
const cellModeBtn = document.getElementById("cellModeBtn");
const sliceModeBtn = document.getElementById("sliceModeBtn");
const componentModeBtn = document.getElementById("componentModeBtn");
const labelBtn = document.getElementById("labelBtn");
const overallViewBtn = document.getElementById("overallViewBtn");
const pairViewBtn = document.getElementById("pairViewBtn");
const pairTools = document.getElementById("pairTools");
const pairASelect = document.getElementById("pairASelect");
const pairBSelect = document.getElementById("pairBSelect");
const pairMeta = document.getElementById("pairMeta");
const prevPairBtn = document.getElementById("prevPairBtn");
const nextPairBtn = document.getElementById("nextPairBtn");
const gl = canvas.getContext("webgl", {antialias:true, alpha:false});
if (!gl) throw new Error("WebGL unavailable");
const needUint32 = META.celltype_index_dtype === "uint32" || META.slice_index_dtype === "uint32" || META.component_index_dtype === "uint32";
const uint32Ext = needUint32 ? gl.getExtension("OES_element_index_uint") : true;
if (!uint32Ext) throw new Error("This viewer requires OES_element_index_uint for more than 65,535 points.");
const vs = `
attribute vec3 aPos;
uniform vec3 uCenter;
uniform float uScale;
uniform float uYaw;
uniform float uPitch;
uniform float uZoom;
uniform float uAspect;
uniform float uPointSize;
void main() {
  vec3 centered = aPos - uCenter;
  vec3 p = vec3(centered.x, centered.z, -centered.y) * uScale;
  float cy = cos(uYaw), sy = sin(uYaw);
  p = vec3(cy*p.x + sy*p.z, p.y, -sy*p.x + cy*p.z);
  float cx = cos(uPitch), sx = sin(uPitch);
  p = vec3(p.x, cx*p.y - sx*p.z, sx*p.y + cx*p.z);
  float z = p.z + uZoom;
  float f = 1.72;
  float near = .1;
  float far = 20.;
  float clipZ = ((far + near) / (far - near)) * z - (2. * far * near) / (far - near);
  gl_Position = vec4((f / uAspect) * p.x, f * p.y, clipZ, z);
  gl_PointSize = uPointSize;
}`;
const fs = `
precision mediump float;
uniform vec4 uColor;
void main() {
  vec2 d = gl_PointCoord - vec2(.5);
  if (dot(d,d) > .25) discard;
  gl_FragColor = uColor;
}`;
function shader(type, src) {
  const s = gl.createShader(type);
  gl.shaderSource(s, src);
  gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
  return s;
}
function program() {
  const p = gl.createProgram();
  gl.attachShader(p, shader(gl.VERTEX_SHADER, vs));
  gl.attachShader(p, shader(gl.FRAGMENT_SHADER, fs));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  return p;
}
const pr = program();
const aPos = gl.getAttribLocation(pr, "aPos");
const uni = Object.fromEntries(["uCenter","uScale","uYaw","uPitch","uZoom","uAspect","uPointSize","uColor"].map(n => [n, gl.getUniformLocation(pr, n)]));
const pointBuffer = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, pointBuffer);
gl.bufferData(gl.ARRAY_BUFFER, points, gl.STATIC_DRAW);
const celltypeIndexBuffer = gl.createBuffer();
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, celltypeIndexBuffer);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, modes.celltype.index, gl.STATIC_DRAW);
const sliceIndexBuffer = gl.createBuffer();
if (modes.slice.index) {
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, sliceIndexBuffer);
  gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, modes.slice.index, gl.STATIC_DRAW);
}
const componentIndexBuffer = gl.createBuffer();
if (modes.component.index) {
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, componentIndexBuffer);
  gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, modes.component.index, gl.STATIC_DRAW);
}
const panels = {celltype: cellPanel, slice: slicePanel, component: componentPanel};
const counts = {celltype: cellCount, slice: sliceCount, component: componentCount};
const modeButtons = {celltype: cellModeBtn, slice: sliceModeBtn, component: componentModeBtn};
const indexBuffers = {celltype: celltypeIndexBuffer, slice: sliceIndexBuffer, component: componentIndexBuffer};
const indexDtypes = {celltype: META.celltype_index_dtype, slice: META.slice_index_dtype, component: META.component_index_dtype};
function renderGroupPanel(panel, modeKey) {
  const mode = modes[modeKey];
  panel.innerHTML = "";
  const pairSliceNames = currentPairSliceNames();
  for (const c of mode.groups) {
    const lockedOut = viewMode === "pair" && modeKey === "slice" && !pairSliceNames.has(c.name);
    const row = document.createElement("label");
    row.className = `row${lockedOut ? " disabled" : ""}`;
    row.innerHTML = `<input type="checkbox" data-mode="${modeKey}" data-name="${esc(c.name)}" ${mode.visible.get(c.name) ? "checked" : ""} ${viewMode === "pair" && modeKey === "slice" ? "disabled" : ""}><span class="sw" style="background:${css(c.color)}"></span><span class="name">${esc(c.name)}</span><span class="num">${c.count.toLocaleString()}</span>`;
    panel.appendChild(row);
  }
  const shown = mode.groups.filter(g => mode.visible.get(g.name)).length;
  const target = counts[modeKey];
  target.textContent = `${shown}/${mode.groups.length}`;
}
function renderPanels() {
  renderGroupPanel(cellPanel, "celltype");
  renderGroupPanel(slicePanel, "slice");
  renderGroupPanel(componentPanel, "component");
  sliceModeBtn.disabled = !modes.slice.groups.length;
  componentModeBtn.disabled = !modes.component.groups.length;
  labelBtn.disabled = !modes.component.groups.length;
  labelBtn.classList.toggle("active", showComponentLabels);
  pairViewBtn.disabled = modes.slice.groups.length < 2;
  document.getElementById("sliceSection").style.display = modes.slice.groups.length ? "" : "none";
  document.getElementById("componentSection").style.display = modes.component.groups.length ? "" : "none";
}
function handlePanelChange(e) {
  const modeKey = e.target.dataset.mode;
  const name = e.target.dataset.name;
  if (!modeKey || !name) return;
  modes[modeKey].visible.set(name, e.target.checked);
  renderGroupPanel(panels[modeKey], modeKey);
  draw();
}
cellPanel.addEventListener("change", handlePanelChange);
slicePanel.addEventListener("change", handlePanelChange);
componentPanel.addEventListener("change", handlePanelChange);
function activeFilterMode() {
  return modes[colorMode];
}
function currentPairSliceNames() {
  const out = new Set();
  if (viewMode !== "pair" || !modes.slice.groups.length) return out;
  const a = modes.slice.groups[pairA];
  const b = modes.slice.groups[pairB];
  if (a) out.add(a.name);
  if (b) out.add(b.name);
  return out;
}
function applyPairSliceVisibility() {
  if (!modes.slice.groups.length) return;
  const selected = currentPairSliceNames();
  if (viewMode === "pair") {
    for (const g of modes.slice.groups) modes.slice.visible.set(g.name, selected.has(g.name));
  }
}
function renderPairControls() {
  const slices = modes.slice.groups;
  pairTools.hidden = viewMode !== "pair";
  if (!slices.length) return;
  pairA = Math.max(0, Math.min(pairA, slices.length - 1));
  pairB = Math.max(0, Math.min(pairB, slices.length - 1));
  if (slices.length > 1 && pairA === pairB) pairB = Math.min(pairA + 1, slices.length - 1);
  const opts = slices.map((s, i) => `<option value="${i}">${esc(s.name)} · ${s.count.toLocaleString()}</option>`).join("");
  pairASelect.innerHTML = opts;
  pairBSelect.innerHTML = opts;
  pairASelect.value = String(pairA);
  pairBSelect.value = String(pairB);
  const lo = Math.min(pairA, pairB);
  const hi = Math.max(pairA, pairB);
  prevPairBtn.disabled = lo <= 0;
  nextPairBtn.disabled = hi >= slices.length - 1;
  const a = slices[pairA];
  const b = slices[pairB];
  const adjacent = Math.abs(pairA - pairB) === 1;
  pairMeta.textContent = a && b ? `${a.name} vs ${b.name} · ${adjacent ? "adjacent" : "custom pair"}` : "";
}
function setViewMode(nextMode) {
  if (nextMode === "pair" && modes.slice.groups.length < 2) return;
  viewMode = nextMode;
  overallViewBtn.classList.toggle("active", viewMode === "overall");
  pairViewBtn.classList.toggle("active", viewMode === "pair");
  if (viewMode === "pair") {
    applyPairSliceVisibility();
  } else {
    for (const g of modes.slice.groups) modes.slice.visible.set(g.name, true);
  }
  renderPairControls();
  renderPanels();
  draw();
}
function setPairIndices(a, b) {
  if (!modes.slice.groups.length) return;
  pairA = Math.max(0, Math.min(Number(a), modes.slice.groups.length - 1));
  pairB = Math.max(0, Math.min(Number(b), modes.slice.groups.length - 1));
  if (modes.slice.groups.length > 1 && pairA === pairB) {
    pairB = pairA < modes.slice.groups.length - 1 ? pairA + 1 : pairA - 1;
  }
  applyPairSliceVisibility();
  renderPairControls();
  renderPanels();
  draw();
}
function stepPair(delta) {
  if (modes.slice.groups.length < 2) return;
  const lo = Math.min(pairA, pairB);
  const next = Math.max(0, Math.min(lo + delta, modes.slice.groups.length - 2));
  setPairIndices(next, next + 1);
}
function setModeVisible(modeKey, value) {
  for (const c of modes[modeKey].groups) modes[modeKey].visible.set(c.name, value);
  if (viewMode === "pair" && modeKey === "slice") applyPairSliceVisibility();
  renderGroupPanel(panels[modeKey], modeKey);
  draw();
}
function invertModeVisible(modeKey) {
  for (const c of modes[modeKey].groups) modes[modeKey].visible.set(c.name, !modes[modeKey].visible.get(c.name));
  if (viewMode === "pair" && modeKey === "slice") applyPairSliceVisibility();
  renderGroupPanel(panels[modeKey], modeKey);
  draw();
}
function setColorMode(modeKey) {
  if (!modes[modeKey] || !modes[modeKey].groups.length) return;
  colorMode = modeKey;
  for (const key of Object.keys(modeButtons)) modeButtons[key].classList.toggle("active", colorMode === key);
  draw();
}
cellModeBtn.onclick = () => setColorMode("celltype");
sliceModeBtn.onclick = () => setColorMode("slice");
componentModeBtn.onclick = () => setColorMode("component");
labelBtn.onclick = () => {
  showComponentLabels = !showComponentLabels;
  labelBtn.classList.toggle("active", showComponentLabels);
  draw();
};
overallViewBtn.onclick = () => setViewMode("overall");
pairViewBtn.onclick = () => setViewMode("pair");
pairASelect.onchange = () => setPairIndices(pairASelect.value, pairBSelect.value);
pairBSelect.onchange = () => setPairIndices(pairASelect.value, pairBSelect.value);
prevPairBtn.onclick = () => stepPair(-1);
nextPairBtn.onclick = () => stepPair(1);
document.getElementById("allBtn").onclick = () => setModeVisible(colorMode, true);
document.getElementById("noneBtn").onclick = () => setModeVisible(colorMode, false);
document.getElementById("invertBtn").onclick = () => invertModeVisible(colorMode);
document.getElementById("resetBtn").onclick = () => { yaw = -0.7; pitch = -0.22; zoom = 3.65; draw(); };
const pointSlider = document.getElementById("pointSlider");
pointSlider.oninput = () => {
  pointSize = Number(pointSlider.value) * (devicePixelRatio || 1);
  document.getElementById("pointVal").textContent = Number(pointSlider.value).toFixed(1);
  draw();
};
let dragging = false, lx = 0, ly = 0;
canvas.onpointerdown = e => {
  dragging = true;
  lx = e.clientX; ly = e.clientY;
  canvas.classList.add("dragging");
  canvas.setPointerCapture(e.pointerId);
};
canvas.onpointermove = e => {
  if (!dragging) return;
  yaw += (e.clientX - lx) * .008;
  pitch = Math.max(-1.45, Math.min(1.45, pitch + (e.clientY - ly) * .008));
  lx = e.clientX; ly = e.clientY;
  draw();
};
canvas.onpointerup = () => { dragging = false; canvas.classList.remove("dragging"); };
canvas.addEventListener("wheel", e => {
  e.preventDefault();
  zoom = Math.max(2.0, Math.min(9.0, zoom * Math.exp(e.deltaY * .001)));
  draw();
}, {passive:false});
function resize() {
  const dpr = devicePixelRatio || 1;
  const w = Math.max(1, Math.floor(canvas.clientWidth * dpr));
  const h = Math.max(1, Math.floor(canvas.clientHeight * dpr));
  if (canvas.width !== w || canvas.height !== h) {
    canvas.width = w;
    canvas.height = h;
    gl.viewport(0, 0, w, h);
  }
}
function draw() {
  resize();
  gl.clearColor(.015, .019, .027, 1);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  gl.enable(gl.DEPTH_TEST);
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  gl.useProgram(pr);
  const m = META.bounds;
  gl.uniform3f(uni.uCenter, m.center[0], m.center[1], m.center[2]);
  gl.uniform1f(uni.uScale, 2.22 / m.max_span);
  gl.uniform1f(uni.uYaw, yaw);
  gl.uniform1f(uni.uPitch, pitch);
  gl.uniform1f(uni.uZoom, zoom);
  gl.uniform1f(uni.uAspect, canvas.width / canvas.height);
  gl.uniform1f(uni.uPointSize, pointSize);
  gl.bindBuffer(gl.ARRAY_BUFFER, pointBuffer);
  gl.enableVertexAttribArray(aPos);
  gl.vertexAttribPointer(aPos, 3, gl.FLOAT, false, 0, 0);
  const active = activeFilterMode();
  const filterModes = Object.keys(modes)
    .filter(k => k !== colorMode && modes[k].groups.length)
    .map(k => ({
      key: k,
      pointIds: modes[k].pointIds,
      visibleIds: new Set(modes[k].groups.map((g, i) => modes[k].visible.get(g.name) ? i : -1).filter(i => i >= 0)),
      total: modes[k].groups.length,
    }))
    .filter(f => f.pointIds && f.visibleIds.size < f.total);
  const useOtherFilter = filterModes.length > 0;
  const indexBuffer = indexBuffers[colorMode];
  const indexType = indexDtypes[colorMode] === "uint32" ? gl.UNSIGNED_INT : gl.UNSIGNED_SHORT;
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer);
  let shown = 0;
  for (const c of active.groups) {
    if (!active.visible.get(c.name)) continue;
    gl.uniform4f(uni.uColor, c.color[0], c.color[1], c.color[2], .93);
    if (!useOtherFilter) {
      gl.drawElements(gl.POINTS, c.count, indexType, c.offset * (indexType === gl.UNSIGNED_INT ? 4 : 2));
      shown += c.count;
      continue;
    }
    const index = active.index;
    let runStart = -1;
    let runCount = 0;
    for (let i = c.offset; i < c.offset + c.count; i++) {
      const pointIdx = index[i];
      let ok = true;
      for (const f of filterModes) {
        if (!f.visibleIds.has(f.pointIds[pointIdx])) {
          ok = false;
          break;
        }
      }
      if (ok) {
        if (runStart < 0) runStart = i;
        runCount += 1;
      } else if (runStart >= 0) {
        gl.drawElements(gl.POINTS, runCount, indexType, runStart * (indexType === gl.UNSIGNED_INT ? 4 : 2));
        shown += runCount;
        runStart = -1;
        runCount = 0;
      }
    }
    if (runStart >= 0) {
      gl.drawElements(gl.POINTS, runCount, indexType, runStart * (indexType === gl.UNSIGNED_INT ? 4 : 2));
      shown += runCount;
    }
  }
  const pairLabel = viewMode === "pair" && modes.slice.groups[pairA] && modes.slice.groups[pairB]
    ? ` · ${modes.slice.groups[pairA].name} vs ${modes.slice.groups[pairB].name}`
    : "";
  statusEl.textContent = `${viewMode === "pair" ? "Pair" : "Overall"} · ${modes[colorMode].label} · ${shown.toLocaleString()} / ${META.total_points.toLocaleString()} cells`;
  hud.innerHTML = `${META.title}<br>${viewMode === "pair" ? "pair" : "overall"}${pairLabel} · color ${modes[colorMode].label.toLowerCase()} · x ${m.span[0].toFixed(1)}, y ${m.span[1].toFixed(1)}, z ${m.span[2].toFixed(1)}`;
  renderComponentLabels();
}
window.onresize = draw;
function normSliceName(value) {
  const text = String(value || "").trim();
  const m = text.match(/-?\d+/);
  return m ? String(Number(m[0])) : text;
}
function componentSliceName(name) {
  return normSliceName(String(name || "").split(" ")[0]);
}
function projectWorldPoint(world) {
  if (!world || world.length < 3) return null;
  const m = META.bounds;
  const scale = 2.22 / m.max_span;
  let px = (world[0] - m.center[0]) * scale;
  let py = (world[2] - m.center[2]) * scale;
  let pz = -(world[1] - m.center[1]) * scale;
  const cy = Math.cos(yaw), sy = Math.sin(yaw);
  const rx = cy * px + sy * pz;
  const rz = -sy * px + cy * pz;
  px = rx; pz = rz;
  const cx = Math.cos(pitch), sx = Math.sin(pitch);
  const ry = cx * py - sx * pz;
  const rz2 = sx * py + cx * pz;
  py = ry; pz = rz2;
  const z = pz + zoom;
  if (z <= 0.05) return null;
  const aspect = Math.max(1e-6, canvas.clientWidth / Math.max(1, canvas.clientHeight));
  const ndcX = (1.72 / aspect) * px / z;
  const ndcY = 1.72 * py / z;
  if (ndcX < -1.12 || ndcX > 1.12 || ndcY < -1.12 || ndcY > 1.12) return null;
  return {
    x: (ndcX * 0.5 + 0.5) * canvas.clientWidth,
    y: (0.5 - ndcY * 0.5) * canvas.clientHeight,
    z,
  };
}
function renderComponentLabels() {
  labelLayer.innerHTML = "";
  if (!showComponentLabels || colorMode !== "component" || !modes.component.groups.length) return;
  const visibleSlices = new Set(modes.slice.groups.filter(g => modes.slice.visible.get(g.name)).map(g => normSliceName(g.name)));
  const candidates = modes.component.groups
    .filter(g => modes.component.visible.get(g.name))
    .filter(g => !String(g.name).toLowerCase().includes("noise"))
    .filter(g => !visibleSlices.size || visibleSlices.has(componentSliceName(g.name)))
    .filter(g => g.centroid);
  const maxLabels = viewMode === "pair" ? 28 : 42;
  for (const g of candidates.slice(0, maxLabels)) {
    const pos = projectWorldPoint(g.centroid);
    if (!pos) continue;
    const el = document.createElement("div");
    el.className = "componentLabel";
    el.textContent = `${g.name} · ${g.count.toLocaleString()}`;
    el.style.left = `${pos.x}px`;
    el.style.top = `${pos.y}px`;
    el.style.setProperty("--label-color", css(g.color));
    labelLayer.appendChild(el);
  }
}
renderPairControls();
renderPanels();
draw();
</script>
</body>
</html>
"""
    return (
        doc.replace("__TITLE__", str(meta["title"]))
        .replace("__META_JSON__", meta_json)
        .replace("__POINTS_B64__", json.dumps(points_b64))
        .replace("__CELLTYPE_INDEX_B64__", json.dumps(celltype_index_b64))
        .replace("__SLICE_INDEX_B64__", json.dumps(slice_index_b64))
        .replace("__COMPONENT_INDEX_B64__", json.dumps(component_index_b64))
        .replace("__POINT_CELLTYPE_ID_B64__", json.dumps(point_celltype_id_b64))
        .replace("__POINT_SLICE_ID_B64__", json.dumps(point_slice_id_b64))
        .replace("__POINT_COMPONENT_ID_B64__", json.dumps(point_component_id_b64))
    )


def build_viewer(args):
    if "full_preview" in str(args.viewer_kind).lower():
        raise SystemExit(
            "viewer_kind must not contain 'full_preview'; use review_full_points, "
            "candidate_full_points, review_balanced300k, candidate_balanced300k, "
            "or debug_all_points."
        )
    if args.component_balanced and args.max_points:
        raise SystemExit("component-balanced viewers must render a precomputed balanced sample without --max-points.")
    input_path = Path(args.input)
    if args.format == "csv":
        points, celltypes, slices, components, component_ranks, resolved_optional_cols = read_csv_points(
            input_path,
            args.x_col,
            args.y_col,
            args.z_col,
            args.celltype_col,
            args.slice_col,
            args.component_col,
            args.component_rank_col,
        )
    else:
        points, celltypes, slices, components, component_ranks, resolved_optional_cols = read_h5ad_points(
            input_path,
            args.obsm_key,
            args.celltype_col,
            args.slice_col,
            args.component_col,
            args.component_rank_col,
        )
    source_count = int(points.shape[0])
    input_sha256 = args.input_sha256 or sha256_file(input_path)
    selected = downsample_indices(celltypes, args.max_points, args.seed)
    points = points[selected]
    celltypes = celltypes[selected]
    if slices is not None:
        slices = slices[selected]
    if components is not None:
        components = components[selected]
    if component_ranks is not None:
        component_ranks = component_ranks[selected]
    z_slice_mapping = None
    if slices is None and args.slice_from_z:
        slices, z_slice_mapping = slices_from_z(points)
    component_labels = component_display_labels(components, component_ranks, slices)
    celltype_indices, cell_meta, celltype_index_dtype = build_mode_groups(celltypes)
    slice_indices, slice_meta, slice_index_dtype = build_mode_groups(slices, sort_key=slice_sort_key)
    component_indices, component_meta, component_index_dtype = build_mode_groups(
        component_labels,
        sort_key=component_sort_key,
        points=points,
    )
    point_celltype_ids, point_celltype_id_dtype = label_ids(celltypes, cell_meta)
    if slices is not None and slice_meta:
        point_slice_ids, point_slice_id_dtype = label_ids(slices, slice_meta)
    else:
        slice_indices = np.asarray([], dtype=np.uint16)
        slice_meta = []
        slice_index_dtype = "uint16"
        point_slice_ids = np.asarray([], dtype=np.uint16)
        point_slice_id_dtype = "uint16"
    if component_labels is not None and component_meta:
        point_component_ids, point_component_id_dtype = label_ids(component_labels, component_meta)
    else:
        component_indices = np.asarray([], dtype=np.uint16)
        component_meta = []
        component_index_dtype = "uint16"
        point_component_ids = np.asarray([], dtype=np.uint16)
        point_component_id_dtype = "uint16"
    source_full_csv = args.source_full_csv or str(input_path)
    source_full_sha256 = args.source_full_sha256
    if not source_full_sha256 and Path(source_full_csv).expanduser() == input_path.expanduser():
        source_full_sha256 = input_sha256
    full_points = int(args.source_full_row_count) if args.source_full_row_count else source_count
    meta = {
        "title": args.title,
        "viewer_kind": args.viewer_kind,
        "display_policy": display_policy_for(args),
        "source": str(input_path),
        "input_sha256": input_sha256,
        "source_full_csv": source_full_csv,
        "source_full_sha256": source_full_sha256,
        "full_points": full_points,
        "format": args.format,
        "total_points": int(points.shape[0]),
        "displayed_points": int(points.shape[0]),
        "source_points": source_count,
        "sampled": int(points.shape[0]) < source_count,
        "component_balanced": bool(args.component_balanced),
        "sampling_method": args.sampling_method,
        "sampling_seed": int(args.sampling_seed if args.sampling_seed is not None else args.seed),
        "max_points": int(args.max_points),
        "celltypes": cell_meta,
        "slices": slice_meta,
        "components": component_meta,
        "bounds": bounds_meta(points),
        "coordinate_columns": [args.x_col, args.y_col, args.z_col] if args.format == "csv" else None,
        "obsm_key": args.obsm_key if args.format == "h5ad" else None,
        "celltype_col": args.celltype_col,
        "slice_col": args.slice_col,
        "component_col": resolved_optional_cols["component_col"],
        "component_rank_col": resolved_optional_cols["component_rank_col"],
        "component_label_format": "SL<slice> rank<component_rank> id<component_id>",
        "component_label_overlay": bool(component_meta),
        "slice_from_z": bool(args.slice_from_z and slices is not None),
        "z_slice_mapping": z_slice_mapping,
        "celltype_index_dtype": celltype_index_dtype,
        "slice_index_dtype": slice_index_dtype,
        "component_index_dtype": component_index_dtype,
        "point_celltype_id_dtype": point_celltype_id_dtype,
        "point_slice_id_dtype": point_slice_id_dtype,
        "point_component_id_dtype": point_component_id_dtype,
    }
    html = render_html(
        meta,
        encode_f32(points),
        encode_index(celltype_indices),
        encode_index(slice_indices),
        encode_index(component_indices),
        encode_index(point_celltype_ids),
        encode_index(point_slice_ids),
        encode_index(point_component_ids),
    )
    out_html = Path(args.output_html)
    out_summary = Path(args.output_summary)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")
    out_summary.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"wrote {out_html}")
    print(f"wrote {out_summary}")
    print(f"points {points.shape[0]:,} / source {source_count:,}")


def main():
    parser = argparse.ArgumentParser(description="Create a self-contained WebGL 3D ST point cloud viewer.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--format", choices=["csv", "h5ad"], required=True)
    parser.add_argument("--x-col", default="x")
    parser.add_argument("--y-col", default="y")
    parser.add_argument("--z-col", default="z")
    parser.add_argument("--obsm-key", default="spatial")
    parser.add_argument("--celltype-col", default="celltype")
    parser.add_argument("--slice-col", default="slice")
    parser.add_argument("--slice-from-z", action="store_true", help="Use discrete z values as slice labels when --slice-col is absent.")
    parser.add_argument("--component-col", default="auto", help="Component id column. 'auto' uses sample_component_id or component_id when present; 'none' disables component mode.")
    parser.add_argument("--component-rank-col", default="auto", help="Component rank column. 'auto' uses sample_component_rank or component_rank when present.")
    parser.add_argument("--title", default="CS13 ST")
    parser.add_argument("--max-points", type=int, default=0, help="0 means all points.")
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--viewer-kind", default="review_balanced300k")
    parser.add_argument("--component-balanced", action="store_true")
    parser.add_argument("--sampling-method", default="")
    parser.add_argument("--sampling-seed", type=int)
    parser.add_argument("--input-sha256", default="")
    parser.add_argument("--source-full-csv", default="")
    parser.add_argument("--source-full-sha256", default="")
    parser.add_argument("--source-full-row-count", type=int, default=0)
    parser.add_argument("--output-html", required=True)
    parser.add_argument("--output-summary", required=True)
    args = parser.parse_args()
    build_viewer(args)


if __name__ == "__main__":
    main()
