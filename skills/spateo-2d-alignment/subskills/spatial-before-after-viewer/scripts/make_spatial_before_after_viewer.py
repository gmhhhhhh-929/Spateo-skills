#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import csv
import json
import math
import re
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


def encode_f32(values: np.ndarray) -> str:
    return base64.b64encode(values.astype("<f4", copy=False).tobytes()).decode("ascii")


def encode_index(values: np.ndarray) -> str:
    return base64.b64encode(values.tobytes()).decode("ascii")


def ordered_unique(values: np.ndarray) -> list[str]:
    out: list[str] = []
    seen = set()
    for value in values.tolist():
        text = str(value)
        if text not in seen:
            out.append(text)
            seen.add(text)
    return out


def hsv_color(h: float, sat: float, val: float) -> list[float]:
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


def distinct_palette(count: int) -> list[list[float]]:
    if count <= 0:
        return []
    return [hsv_color(i / max(count, 1), 0.72, 0.95) for i in range(count)]


_SLICE_SORT_RE = re.compile(r"-?\d+")


def slice_sort_key(label: str) -> tuple[int, int | str, str]:
    match = _SLICE_SORT_RE.search(str(label))
    if match:
        return (0, int(match.group(0)), str(label))
    return (1, str(label), str(label))


def build_groups(values: np.ndarray, *, sort_slices: bool = False) -> tuple[np.ndarray, list[dict], str, np.ndarray, str]:
    order = ordered_unique(values)
    if sort_slices:
        order = sorted(order, key=slice_sort_key)
    palette = distinct_palette(len(order))
    index_dtype = np.uint32 if values.shape[0] > np.iinfo(np.uint16).max else np.uint16
    point_dtype = np.uint32 if len(order) > np.iinfo(np.uint16).max else np.uint16
    index_parts = []
    point_ids = np.empty(values.shape[0], dtype=point_dtype)
    groups = []
    offset = 0
    for i, label in enumerate(order):
        idx = np.where(values.astype(str) == label)[0].astype(index_dtype, copy=False)
        index_parts.append(idx)
        point_ids[idx] = i
        groups.append(
            {
                "name": label,
                "count": int(idx.size),
                "offset": int(offset),
                "color": [round(float(c), 6) for c in palette[i]],
            }
        )
        offset += int(idx.size)
    if index_parts:
        index = np.concatenate(index_parts)
    else:
        index = np.asarray([], dtype=index_dtype)
    return index, groups, "uint32" if index_dtype == np.uint32 else "uint16", point_ids, "uint32" if point_dtype == np.uint32 else "uint16"


def downsample_indices(n: int, max_points: int, seed: int) -> np.ndarray:
    if max_points <= 0 or n <= max_points:
        return np.arange(n, dtype=np.int64)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, size=max_points, replace=False)).astype(np.int64)


def bounds(points: np.ndarray) -> dict:
    pmin = points.min(axis=0)
    pmax = points.max(axis=0)
    center = (pmin + pmax) / 2
    span = np.maximum(pmax - pmin, 1e-6)
    return {
        "min": [float(x) for x in pmin],
        "max": [float(x) for x in pmax],
        "center": [float(x) for x in center],
        "span": [float(x) for x in span],
        "max_span": float(span.max()),
    }


def load_font(size: int) -> ImageFont.ImageFont:
    for path in ["/System/Library/Fonts/Supplemental/Arial.ttf", "/System/Library/Fonts/Helvetica.ttc"]:
        try:
            return ImageFont.truetype(path, size)
        except Exception:
            pass
    return ImageFont.load_default()


def sample_points(points: np.ndarray, max_n: int, seed: int) -> np.ndarray:
    if len(points) <= max_n:
        return points.astype(float, copy=True)
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(len(points), size=max_n, replace=False))
    return points[idx].astype(float, copy=True)


def screen(points: np.ndarray, plot_bounds: tuple[float, float, float, float], panel: tuple[int, int, int, int]) -> np.ndarray:
    xmin, xmax, ymin, ymax = plot_bounds
    x0, y0, x1, y1 = panel
    pad = 28
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


def make_contact_sheet(before: np.ndarray, after: np.ndarray, slices: np.ndarray, slice_labels: list[str], out_path: Path, title: str) -> None:
    pairs = list(zip(slice_labels[:-1], slice_labels[1:]))
    if not pairs:
        return
    thumbs = []
    for a, b in pairs:
        mask_a = slices.astype(str) == str(a)
        mask_b = slices.astype(str) == str(b)
        if not mask_a.any() or not mask_b.any():
            continue
        bfa = before[mask_a][:, :2]
        bfb = before[mask_b][:, :2]
        afa = after[mask_a][:, :2]
        afb = after[mask_b][:, :2]
        allpts = np.vstack([bfa, bfb, afa, afb])
        mx = (allpts[:, 0].max() - allpts[:, 0].min()) * 0.08
        my = (allpts[:, 1].max() - allpts[:, 1].min()) * 0.08
        plot_bounds = (float(allpts[:, 0].min() - mx), float(allpts[:, 0].max() + mx), float(allpts[:, 1].min() - my), float(allpts[:, 1].max() + my))
        panel_img = Image.new("RGBA", (980, 500), (255, 255, 255, 255))
        draw = ImageDraw.Draw(panel_img)
        draw.text((18, 12), f"SL{a}/SL{b}", fill=(25, 32, 42), font=load_font(20))
        panels = [((18, 50, 472, 482), "before", bfa, bfb), ((508, 50, 962, 482), "after", afa, afb)]
        for panel, label, p_a, p_b in panels:
            draw.rectangle(panel, fill=(255, 255, 255, 255), outline=(210, 215, 220, 255), width=1)
            draw.text((panel[0] + 12, panel[1] + 10), label, fill=(55, 64, 76), font=load_font(14))
            box = (panel[0] + 5, panel[1] + 34, panel[2] - 5, panel[3] - 5)
            draw_points(panel_img, screen(sample_points(p_a, 15000, int(re.sub(r"[^0-9-]", "", str(a)) or 0)), plot_bounds, box), (38, 103, 207, 72), 1)
            draw_points(panel_img, screen(sample_points(p_b, 15000, int(re.sub(r"[^0-9-]", "", str(b)) or 1)), plot_bounds, box), (218, 91, 44, 88), 1)
        thumbs.append((f"SL{a}_SL{b}", panel_img.convert("RGB")))
    cols = 2
    rows = math.ceil(len(thumbs) / cols)
    out = Image.new("RGB", (1080, 72 + rows * 300), (246, 248, 250))
    draw = ImageDraw.Draw(out)
    draw.text((30, 24), title, fill=(20, 28, 38), font=load_font(26))
    for i, (label, im) in enumerate(thumbs):
        im.thumbnail((510, 260), Image.LANCZOS)
        x = 30 + (i % cols) * 525
        y = 74 + (i // cols) * 300
        draw.text((x, y), label, fill=(35, 42, 52), font=load_font(15))
        out.paste(im, (x, y + 22))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.save(out_path, quality=94)


def read_points(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    required = [
        args.before_x,
        args.before_y,
        args.before_z,
        args.after_x,
        args.after_y,
        args.after_z,
        args.slice_col,
        args.color_col,
    ]
    before = []
    after = []
    slices = []
    colors = []
    with Path(args.input).open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fields = set(reader.fieldnames or [])
        missing = [col for col in required if col not in fields]
        if missing:
            raise SystemExit(f"Missing columns: {', '.join(missing)}")
        for row in reader:
            before.append([float(row[args.before_x]), float(row[args.before_y]), float(row[args.before_z])])
            after.append([float(row[args.after_x]), float(row[args.after_y]), float(row[args.after_z])])
            slices.append(str(row[args.slice_col]))
            colors.append(str(row[args.color_col]))
    return (
        np.asarray(before, dtype=np.float32),
        np.asarray(after, dtype=np.float32),
        np.asarray(slices, dtype=object),
        np.asarray(colors, dtype=object),
    )


def apply_display_z(points: np.ndarray, slices: np.ndarray, args: argparse.Namespace) -> tuple[np.ndarray, dict[str, float]]:
    out = points.copy()
    order = sorted(ordered_unique(slices), key=slice_sort_key)
    if args.physical_z_step > 0:
        mapping = {label: i * float(args.physical_z_step) for i, label in enumerate(order)}
        z = np.asarray([mapping[str(v)] for v in slices], dtype=np.float32)
    else:
        z = out[:, 2].astype(np.float32, copy=True)
    display = z * float(args.z_display_multiplier)
    out[:, 2] = display
    return out, {label: float(value) for label, value in (mapping.items() if args.physical_z_step > 0 else [])}


def render_html(meta: dict, before_b64: str, after_b64: str, slice_index_b64: str, color_index_b64: str, point_slice_ids_b64: str, point_color_ids_b64: str) -> str:
    meta_json = json.dumps(meta, ensure_ascii=False)
    doc = r"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
:root { color-scheme:dark; --bg:#070b10; --panel:#101722; --line:#263241; --text:#e7edf6; --muted:#93a4b8; }
* { box-sizing:border-box; }
body { margin:0;height:100vh;display:grid;grid-template-columns:290px 1fr;background:var(--bg);color:var(--text);font:13px/1.4 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }
aside { min-width:0;border-right:1px solid var(--line);background:var(--panel);display:flex;flex-direction:column;overflow:hidden; }
header { padding:13px 14px 9px;border-bottom:1px solid var(--line);display:grid;gap:9px; }
h1 { margin:0;font-size:16px;letter-spacing:0; }
.toolbar { display:grid;grid-template-columns:repeat(3,1fr);gap:6px; }
button { border:1px solid var(--line);background:#151e2a;color:var(--text);border-radius:7px;min-height:32px;padding:0 8px;font:inherit;font-size:12px;cursor:pointer; }
button.active { background:#2b6fd6;border-color:#4d8eff; }
button:disabled { opacity:.45;cursor:not-allowed; }
select { min-width:0;width:100%;height:32px;border:1px solid var(--line);border-radius:7px;background:#0d1420;color:var(--text);padding:0 8px;font:inherit;font-size:12px; }
.pairControls { display:grid;grid-template-columns:1fr 1fr;gap:6px; }
.pairButtons { display:grid;grid-template-columns:1fr 1fr;gap:6px; }
.hint { color:var(--muted);font-size:11px; }
.sectionTitle { display:flex;align-items:center;justify-content:space-between;gap:10px;padding:8px 14px 2px;color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:0; }
.panel { overflow:auto;padding:4px 8px 10px; }
.row { display:grid;grid-template-columns:18px 14px 1fr auto;gap:7px;align-items:center;min-height:30px;padding:3px 5px;border-radius:6px;font-size:12px; }
.sw { width:10px;height:10px;border-radius:999px;border:1px solid rgba(255,255,255,.25); }
.name { overflow:hidden;text-overflow:ellipsis;white-space:nowrap; }
.num { color:var(--muted);font-variant-numeric:tabular-nums; }
main { position:relative;min-width:0;min-height:0; }
canvas { width:100%;height:100%;display:block;cursor:grab; }
canvas.dragging { cursor:grabbing; }
#hud,#status { position:absolute;background:rgba(5,8,12,.68);border:1px solid rgba(255,255,255,.12);border-radius:8px;color:#dce7f7;font-size:12px;line-height:1.45;backdrop-filter:blur(8px); }
#hud { left:14px;top:14px;padding:9px 11px;max-width:min(540px,calc(100% - 28px)); }
#status { left:14px;bottom:14px;padding:7px 10px; }
</style>
</head>
<body>
<aside>
<header>
<h1>__TITLE__</h1>
<div class="toolbar">
<button id="beforeBtn" class="active">Before</button>
<button id="afterBtn">After</button>
<button id="deltaBtn">Delta</button>
</div>
<div class="toolbar">
<button id="overallBtn" class="active">Overall</button>
<button id="pairBtn">Pair</button>
<button id="resetBtn">Reset</button>
</div>
<div class="toolbar">
<button id="sliceColorBtn" class="active">Slice</button>
<button id="groupColorBtn">Group</button>
<button id="smallPointBtn">Point</button>
</div>
<div class="pairControls">
<select id="pairA"></select>
<select id="pairB"></select>
</div>
<div class="pairButtons">
<button id="prevPair">Prev</button>
<button id="nextPair">Next</button>
</div>
<div class="hint" id="pairMeta"></div>
</header>
<div class="sectionTitle"><span>Slice</span><span id="sliceCount"></span></div>
<div id="slicePanel" class="panel"></div>
<div class="sectionTitle"><span>Color group</span><span id="colorCount"></span></div>
<div id="colorPanel" class="panel"></div>
</aside>
<main>
<canvas id="gl"></canvas>
<div id="hud"></div>
<div id="status"></div>
</main>
<script>
const META = __META_JSON__;
const BEFORE_B64 = __BEFORE_B64__;
const AFTER_B64 = __AFTER_B64__;
const SLICE_INDEX_B64 = __SLICE_INDEX_B64__;
const COLOR_INDEX_B64 = __COLOR_INDEX_B64__;
const POINT_SLICE_IDS_B64 = __POINT_SLICE_IDS_B64__;
const POINT_COLOR_IDS_B64 = __POINT_COLOR_IDS_B64__;
function decodeF32(b64) {
  const bin = atob(b64); const bytes = new Uint8Array(bin.length);
  for (let i=0;i<bin.length;i++) bytes[i] = bin.charCodeAt(i);
  return new Float32Array(bytes.buffer);
}
function decodeIndex(b64, dtype) {
  const bin = atob(b64); const bytes = new Uint8Array(bin.length);
  for (let i=0;i<bin.length;i++) bytes[i] = bin.charCodeAt(i);
  return dtype === "uint32" ? new Uint32Array(bytes.buffer) : new Uint16Array(bytes.buffer);
}
const beforePoints = decodeF32(BEFORE_B64);
const afterPoints = decodeF32(AFTER_B64);
const sliceIndex = decodeIndex(SLICE_INDEX_B64, META.slice_index_dtype);
const colorIndex = decodeIndex(COLOR_INDEX_B64, META.color_index_dtype);
const pointSliceIds = decodeIndex(POINT_SLICE_IDS_B64, META.point_slice_id_dtype);
const pointColorIds = decodeIndex(POINT_COLOR_IDS_B64, META.point_color_id_dtype);
const canvas = document.getElementById("gl");
const gl = canvas.getContext("webgl", {antialias:true});
if (!gl) throw new Error("WebGL unavailable");
const uint32Ext = gl.getExtension("OES_element_index_uint");
if ((META.slice_index_dtype === "uint32" || META.color_index_dtype === "uint32") && !uint32Ext) throw new Error("OES_element_index_uint unavailable");
const vs = `
attribute vec3 aPos;
uniform vec3 uCenter;
uniform float uScale,uYaw,uPitch,uZoom,uAspect,uPointSize;
void main() {
  vec3 centered = aPos - uCenter;
  vec3 p = vec3(centered.x, centered.z, -centered.y) * uScale;
  float cy=cos(uYaw), sy=sin(uYaw), cx=cos(uPitch), sx=sin(uPitch);
  p = vec3(cy*p.x + sy*p.z, p.y, -sy*p.x + cy*p.z);
  p = vec3(p.x, cx*p.y - sx*p.z, sx*p.y + cx*p.z);
  float z = p.z + uZoom;
  float f = 1.35, near=.08, far=30.;
  float clipZ = ((far + near) / (far - near)) * z - (2. * far * near) / (far - near);
  gl_Position = vec4((f/uAspect)*p.x, f*p.y, clipZ, z);
  gl_PointSize = uPointSize;
}`;
const fs = `precision mediump float; uniform vec4 uColor; void main(){ vec2 d=gl_PointCoord-vec2(.5); if(dot(d,d)>.25) discard; gl_FragColor=uColor; }`;
function shader(type, src){ const s=gl.createShader(type); gl.shaderSource(s,src); gl.compileShader(s); if(!gl.getShaderParameter(s,gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s)); return s; }
const pr=gl.createProgram(); gl.attachShader(pr, shader(gl.VERTEX_SHADER, vs)); gl.attachShader(pr, shader(gl.FRAGMENT_SHADER, fs)); gl.linkProgram(pr);
if(!gl.getProgramParameter(pr, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(pr));
const aPos=gl.getAttribLocation(pr,"aPos");
const uni=Object.fromEntries(["uCenter","uScale","uYaw","uPitch","uZoom","uAspect","uPointSize","uColor"].map(n=>[n,gl.getUniformLocation(pr,n)]));
const pointBuffer=gl.createBuffer();
const sliceIndexBuffer=gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, sliceIndexBuffer); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, sliceIndex, gl.STATIC_DRAW);
const colorIndexBuffer=gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, colorIndexBuffer); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, colorIndex, gl.STATIC_DRAW);
let coordMode="before", viewMode="overall";
let colorMode="slice";
let yaw=-0.72, pitch=-0.34, zoom=3.7, pointSize=3*(devicePixelRatio||1);
let pairA=0, pairB=META.slices.length>1?1:0;
const visibleSlices = new Map(META.slices.map(s=>[s.name,true]));
const visibleColors = new Map(META.colors.map(s=>[s.name,true]));
function css(c){ return `rgb(${Math.round(c[0]*255)},${Math.round(c[1]*255)},${Math.round(c[2]*255)})`; }
function esc(s){ return String(s).replace(/[&<>"']/g, m=>({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#039;"}[m])); }
function currentPoints(){
  if (coordMode === "before") return beforePoints;
  if (coordMode === "after") return afterPoints;
  const arr = new Float32Array(beforePoints.length);
  for (let i=0;i<beforePoints.length;i+=3) {
    arr[i] = beforePoints[i] + (afterPoints[i]-beforePoints[i]) * META.delta_scale;
    arr[i+1] = beforePoints[i+1] + (afterPoints[i+1]-beforePoints[i+1]) * META.delta_scale;
    arr[i+2] = beforePoints[i+2] + (afterPoints[i+2]-beforePoints[i+2]);
  }
  return arr;
}
function renderPanel(panel, groups, visible, kind){
  panel.innerHTML="";
  for (const g of groups) {
    const row=document.createElement("label"); row.className="row";
    row.innerHTML=`<input type="checkbox" data-kind="${kind}" data-name="${esc(g.name)}" ${visible.get(g.name)?"checked":""}><span class="sw" style="background:${css(g.color)}"></span><span class="name">${esc(g.name)}</span><span class="num">${g.count.toLocaleString()}</span>`;
    panel.appendChild(row);
  }
}
function renderPanels(){
  renderPanel(document.getElementById("slicePanel"), META.slices, visibleSlices, "slice");
  renderPanel(document.getElementById("colorPanel"), META.colors, visibleColors, "color");
  document.getElementById("sliceCount").textContent = `${[...visibleSlices.values()].filter(Boolean).length}/${META.slices.length}`;
  document.getElementById("colorCount").textContent = `${[...visibleColors.values()].filter(Boolean).length}/${META.colors.length}`;
}
function setPair(a,b){
  pairA=Math.max(0,Math.min(Number(a),META.slices.length-1));
  pairB=Math.max(0,Math.min(Number(b),META.slices.length-1));
  if (pairA===pairB && META.slices.length>1) pairB = pairA < META.slices.length-1 ? pairA+1 : pairA-1;
  for (const s of META.slices) visibleSlices.set(s.name, viewMode==="pair" ? false : true);
  if (viewMode==="pair") { visibleSlices.set(META.slices[pairA].name,true); visibleSlices.set(META.slices[pairB].name,true); }
  renderPairControls(); renderPanels(); draw();
}
function renderPairControls(){
  const opts=META.slices.map((s,i)=>`<option value="${i}">${esc(s.name)} · ${s.count.toLocaleString()}</option>`).join("");
  const a=document.getElementById("pairA"), b=document.getElementById("pairB");
  a.innerHTML=opts; b.innerHTML=opts; a.value=pairA; b.value=pairB;
  document.getElementById("pairMeta").textContent = META.slices[pairA] && META.slices[pairB] ? `${META.slices[pairA].name} vs ${META.slices[pairB].name}` : "";
}
document.getElementById("slicePanel").onchange=e=>{ if(!e.target.dataset.kind)return; visibleSlices.set(e.target.dataset.name,e.target.checked); renderPanels(); draw(); };
document.getElementById("colorPanel").onchange=e=>{ if(!e.target.dataset.kind)return; visibleColors.set(e.target.dataset.name,e.target.checked); renderPanels(); draw(); };
document.getElementById("beforeBtn").onclick=()=>setCoord("before");
document.getElementById("afterBtn").onclick=()=>setCoord("after");
document.getElementById("deltaBtn").onclick=()=>setCoord("delta");
function setCoord(m){ coordMode=m; for (const id of ["beforeBtn","afterBtn","deltaBtn"]) document.getElementById(id).classList.toggle("active", id.startsWith(m)); draw(); }
document.getElementById("overallBtn").onclick=()=>{ viewMode="overall"; document.getElementById("overallBtn").classList.add("active"); document.getElementById("pairBtn").classList.remove("active"); for(const s of META.slices) visibleSlices.set(s.name,true); renderPanels(); draw(); };
document.getElementById("pairBtn").onclick=()=>{ viewMode="pair"; document.getElementById("pairBtn").classList.add("active"); document.getElementById("overallBtn").classList.remove("active"); setPair(pairA,pairB); };
document.getElementById("pairA").onchange=e=>setPair(e.target.value,pairB);
document.getElementById("pairB").onchange=e=>setPair(pairA,e.target.value);
document.getElementById("prevPair").onclick=()=>{ const lo=Math.max(0,Math.min(pairA,pairB)-1); setPair(lo,lo+1); };
document.getElementById("nextPair").onclick=()=>{ const lo=Math.min(META.slices.length-2,Math.min(pairA,pairB)+1); setPair(lo,lo+1); };
document.getElementById("resetBtn").onclick=()=>{ yaw=-0.72; pitch=-0.34; zoom=3.7; draw(); };
document.getElementById("sliceColorBtn").onclick=()=>setColorMode("slice");
document.getElementById("groupColorBtn").onclick=()=>setColorMode("group");
document.getElementById("smallPointBtn").onclick=()=>{ pointSize = pointSize > 2*(devicePixelRatio||1) ? 1.5*(devicePixelRatio||1) : 3*(devicePixelRatio||1); draw(); };
function setColorMode(m){ colorMode=m; document.getElementById("sliceColorBtn").classList.toggle("active",m==="slice"); document.getElementById("groupColorBtn").classList.toggle("active",m==="group"); draw(); }
let dragging=false,lastX=0,lastY=0;
canvas.onpointerdown=e=>{dragging=true;lastX=e.clientX;lastY=e.clientY;canvas.classList.add("dragging");canvas.setPointerCapture(e.pointerId);};
canvas.onpointermove=e=>{if(!dragging)return;const dx=e.clientX-lastX,dy=e.clientY-lastY;lastX=e.clientX;lastY=e.clientY;yaw+=dx*.006;pitch=Math.max(-1.45,Math.min(1.45,pitch+dy*.006));draw();};
canvas.onpointerup=()=>{dragging=false;canvas.classList.remove("dragging");};
canvas.onwheel=e=>{e.preventDefault();zoom=Math.max(2,Math.min(10,zoom*Math.exp(e.deltaY*.001)));draw();},{passive:false};
function resize(){ const dpr=devicePixelRatio||1,w=Math.max(1,Math.floor(canvas.clientWidth*dpr)),h=Math.max(1,Math.floor(canvas.clientHeight*dpr)); if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;gl.viewport(0,0,w,h);} }
function draw(){
  resize(); gl.clearColor(.015,.019,.027,1); gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT); gl.enable(gl.DEPTH_TEST); gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA);
  const pts=currentPoints(); gl.bindBuffer(gl.ARRAY_BUFFER, pointBuffer); gl.bufferData(gl.ARRAY_BUFFER, pts, gl.DYNAMIC_DRAW);
  gl.useProgram(pr); const m=META.bounds;
  gl.uniform3f(uni.uCenter,m.center[0],m.center[1],m.center[2]); gl.uniform1f(uni.uScale,2.22/m.max_span); gl.uniform1f(uni.uYaw,yaw); gl.uniform1f(uni.uPitch,pitch); gl.uniform1f(uni.uZoom,zoom); gl.uniform1f(uni.uAspect,canvas.width/canvas.height); gl.uniform1f(uni.uPointSize,pointSize);
  gl.enableVertexAttribArray(aPos); gl.vertexAttribPointer(aPos,3,gl.FLOAT,false,0,0);
  const activeGroups = colorMode==="slice" ? META.slices : META.colors;
  const activeIndex = colorMode==="slice" ? sliceIndex : colorIndex;
  const indexBuffer = colorMode==="slice" ? sliceIndexBuffer : colorIndexBuffer;
  const indexType = (colorMode==="slice" ? META.slice_index_dtype : META.color_index_dtype)==="uint32"?gl.UNSIGNED_INT:gl.UNSIGNED_SHORT;
  const otherGroups = colorMode==="slice" ? META.colors : META.slices;
  const otherPointIds = colorMode==="slice" ? pointColorIds : pointSliceIds;
  const activeVisible = colorMode==="slice" ? visibleSlices : visibleColors;
  const otherVisible = colorMode==="slice" ? visibleColors : visibleSlices;
  const otherVisibleIds = new Set(otherGroups.map((g,i)=>otherVisible.get(g.name)?i:-1).filter(i=>i>=0));
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,indexBuffer);
  let shown=0;
  for (let gi=0; gi<activeGroups.length; gi++) {
    const g=activeGroups[gi]; if(!activeVisible.get(g.name)) continue;
    gl.uniform4f(uni.uColor,g.color[0],g.color[1],g.color[2],.9);
    let run=-1,count=0;
    for(let i=g.offset;i<g.offset+g.count;i++){
      const pointIdx=activeIndex[i], otherId=otherPointIds[pointIdx];
      if(otherVisibleIds.has(otherId)){ if(run<0)run=i; count++; }
      else if(run>=0){ gl.drawElements(gl.POINTS,count,indexType,run*(indexType===gl.UNSIGNED_INT?4:2)); shown+=count; run=-1; count=0; }
    }
    if(run>=0){ gl.drawElements(gl.POINTS,count,indexType,run*(indexType===gl.UNSIGNED_INT?4:2)); shown+=count; }
  }
  document.getElementById("status").textContent=`${coordMode} · ${viewMode} · color ${colorMode} · ${shown.toLocaleString()} / ${META.total_points.toLocaleString()} cells`;
  document.getElementById("hud").innerHTML=`${META.title}<br>${coordMode} · color ${colorMode} · x ${m.span[0].toFixed(1)}, y ${m.span[1].toFixed(1)}, z ${m.span[2].toFixed(1)} · z display ${META.z_display_multiplier}x`;
}
window.onresize=draw;
renderPairControls(); renderPanels(); draw();
</script>
</body>
</html>"""
    return (
        doc.replace("__TITLE__", str(meta["title"]))
        .replace("__META_JSON__", meta_json)
        .replace("__BEFORE_B64__", json.dumps(before_b64))
        .replace("__AFTER_B64__", json.dumps(after_b64))
        .replace("__SLICE_INDEX_B64__", json.dumps(slice_index_b64))
        .replace("__COLOR_INDEX_B64__", json.dumps(color_index_b64))
        .replace("__POINT_SLICE_IDS_B64__", json.dumps(point_slice_ids_b64))
        .replace("__POINT_COLOR_IDS_B64__", json.dumps(point_color_ids_b64))
    )


def build(args: argparse.Namespace) -> None:
    before_raw, after_raw, slices, colors = read_points(args)
    before, z_mapping = apply_display_z(before_raw, slices, args)
    after, _ = apply_display_z(after_raw, slices, args)
    source_count = int(before.shape[0])
    slice_labels_full = sorted(ordered_unique(slices), key=slice_sort_key)
    if args.output_contact_sheet:
        make_contact_sheet(before, after, slices, slice_labels_full, Path(args.output_contact_sheet), args.title)
    selected = downsample_indices(source_count, args.max_points, args.seed)
    before = before[selected]
    after = after[selected]
    slices = slices[selected]
    colors = colors[selected]
    slice_index, slice_groups, slice_index_dtype, point_slice_ids, point_slice_id_dtype = build_groups(slices, sort_slices=True)
    color_index, color_groups, color_index_dtype, point_color_ids, point_color_id_dtype = build_groups(colors)
    all_points = np.vstack([before, after])
    meta = {
        "title": args.title,
        "source": str(Path(args.input)),
        "total_points": int(before.shape[0]),
        "source_points": source_count,
        "sampled": int(before.shape[0]) < source_count,
        "max_points": int(args.max_points),
        "slices": slice_groups,
        "colors": color_groups,
        "bounds": bounds(all_points),
        "before_columns": [args.before_x, args.before_y, args.before_z],
        "after_columns": [args.after_x, args.after_y, args.after_z],
        "slice_col": args.slice_col,
        "color_col": args.color_col,
        "physical_z_step": float(args.physical_z_step),
        "z_display_multiplier": float(args.z_display_multiplier),
        "z_slice_mapping": z_mapping,
        "delta_scale": float(args.delta_scale),
        "slice_index_dtype": slice_index_dtype,
        "color_index_dtype": color_index_dtype,
        "point_slice_id_dtype": point_slice_id_dtype,
        "point_color_id_dtype": point_color_id_dtype,
    }
    html = render_html(
        meta,
        encode_f32(before),
        encode_f32(after),
        encode_index(slice_index),
        encode_index(color_index),
        encode_index(point_slice_ids),
        encode_index(point_color_ids),
    )
    out_html = Path(args.output_html)
    out_summary = Path(args.output_summary)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")
    out_summary.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"wrote {out_html}")
    print(f"wrote {out_summary}")
    print(f"points {before.shape[0]:,} / source {source_count:,}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a self-contained before/after spatial point cloud viewer.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--slice-col", default="sl_number")
    parser.add_argument("--color-col", default="celltype")
    parser.add_argument("--before-x", required=True)
    parser.add_argument("--before-y", required=True)
    parser.add_argument("--before-z", required=True)
    parser.add_argument("--after-x", required=True)
    parser.add_argument("--after-y", required=True)
    parser.add_argument("--after-z", required=True)
    parser.add_argument("--physical-z-step", type=float, default=0.0)
    parser.add_argument("--z-display-multiplier", type=float, default=1.0)
    parser.add_argument("--delta-scale", type=float, default=6.0)
    parser.add_argument("--title", default="Spatial Before/After")
    parser.add_argument("--max-points", type=int, default=0)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--output-html", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-contact-sheet", default="")
    args = parser.parse_args()
    build(args)


if __name__ == "__main__":
    main()
