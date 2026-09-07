#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


def slice_key(value: str) -> int:
    return int(float(str(value).replace("SL", "")))


def read_scan(path: Path, include_statuses: set[str], top_n: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if include_statuses and row.get("status") not in include_statuses:
                continue
            row["triad_outlier_score"] = float(row["triad_outlier_score"])
            row["component_rank"] = int(float(row["component_rank"]))
            row["mid_points"] = int(float(row["mid_points"]))
            for key in ["mid_expected_x", "mid_expected_y", "mid_centroid_x", "mid_centroid_y"]:
                row[key] = float(row[key])
            rows.append(row)
    rows.sort(key=lambda r: (-float(r["triad_outlier_score"]), r["triad_id"], int(r["component_rank"])))
    if top_n > 0:
        rows = rows[:top_n]
    return rows


def triad_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return str(row["prev_slice"]), str(row["mid_slice"]), str(row["next_slice"])


def read_component_labels(path: Path, wanted_slices: set[str]) -> dict[str, dict[str, str]]:
    labels: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"cell_id", "sl_number", "component_rank", "component_id"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Component labels missing columns: {sorted(missing)}")
        for row in reader:
            sl = str(int(float(row["sl_number"])))
            if sl in wanted_slices:
                labels[row["cell_id"]] = {
                    "component_rank": str(int(float(row["component_rank"]))),
                    "component_id": str(row["component_id"]),
                }
    return labels


def choose_points(points: list[dict[str, Any]], max_points: int, max_highlight_points: int) -> list[dict[str, Any]]:
    if max_points <= 0 or len(points) <= max_points:
        return points
    highlighted = [p for p in points if p["highlight"]]
    other = [p for p in points if not p["highlight"]]
    if len(highlighted) > max_highlight_points:
        step = max(1, len(highlighted) // max_highlight_points)
        highlighted = highlighted[::step][:max_highlight_points]
    keep_other = max(0, max_points - len(highlighted))
    if len(other) <= keep_other:
        return highlighted + other
    step = max(1, len(other) // keep_other)
    sampled = other[::step][:keep_other]
    return highlighted + sampled


def load_points(
    coordinates: Path,
    components: Path,
    triads: list[dict[str, Any]],
    args: argparse.Namespace,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    wanted_slices = {str(row[k]) for row in triads for k in ["prev_slice", "mid_slice", "next_slice"]}
    labels = read_component_labels(components, wanted_slices)
    highlight_by_mid: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in triads:
        if row.get("status") in {"flagged", "watch"}:
            highlight_by_mid[(str(row["mid_slice"]), str(row["component_rank"]))].add(str(row["triad_id"]))
    triad_slices = {str(row["triad_id"]): {str(row["prev_slice"]), str(row["mid_slice"]), str(row["next_slice"])} for row in triads}
    slice_to_triads: dict[str, list[str]] = defaultdict(list)
    for triad_id, slices in triad_slices.items():
        for sl in slices:
            slice_to_triads[sl].append(triad_id)
    triad_points: dict[str, list[dict[str, Any]]] = {str(row["triad_id"]): [] for row in triads}
    with coordinates.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"cell_id", args.slice_col, args.x_col, args.y_col, args.z_col, args.celltype_col}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Coordinates missing columns: {sorted(missing)}")
        for row in reader:
            sl = str(int(float(row[args.slice_col])))
            if sl not in wanted_slices:
                continue
            label = labels.get(row["cell_id"])
            if label is None:
                continue
            rank = label["component_rank"]
            for triad_id in slice_to_triads.get(sl, []):
                highlight = triad_id in highlight_by_mid.get((sl, rank), set())
                triad_points[triad_id].append(
                    {
                        "x": float(row[args.x_col]),
                        "y": float(row[args.y_col]),
                        "z": float(row[args.z_col]),
                        "slice": sl,
                        "component_rank": rank,
                        "component_id": label["component_id"],
                        "celltype": row.get(args.celltype_col, ""),
                        "highlight": highlight,
                    }
                )
    payload_triads: list[dict[str, Any]] = []
    total_points = 0
    for row in triads:
        triad_id = str(row["triad_id"])
        pts = choose_points(triad_points.get(triad_id, []), args.max_points_per_triad, args.max_highlight_points_per_triad)
        total_points += len(pts)
        payload_triads.append(
            {
                "triad_id": triad_id,
                "prev_slice": str(row["prev_slice"]),
                "mid_slice": str(row["mid_slice"]),
                "next_slice": str(row["next_slice"]),
                "status": str(row["status"]),
                "issue_type": str(row["issue_type"]),
                "component_rank": int(row["component_rank"]),
                "mid_component_id": str(row["mid_component_id"]),
                "mid_points": int(row["mid_points"]),
                "score": round(float(row["triad_outlier_score"]), 3),
                "reason": str(row["reason"]),
                "mid_expected_x": float(row["mid_expected_x"]),
                "mid_expected_y": float(row["mid_expected_y"]),
                "mid_centroid_x": float(row["mid_centroid_x"]),
                "mid_centroid_y": float(row["mid_centroid_y"]),
                "points": pts,
            }
        )
    summary = {
        "triads": len(payload_triads),
        "total_rendered_points": total_points,
        "max_points_per_triad": args.max_points_per_triad,
        "wanted_slices": sorted(wanted_slices, key=slice_key),
    }
    return payload_triads, summary


def render_html(payload: dict[str, Any]) -> str:
    data = json.dumps(payload, ensure_ascii=False)
    title = html.escape(payload["title"])
    return f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
<style>
:root {{ color-scheme:dark; --bg:#080b0f; --panel:#101820; --panel2:#16212b; --line:#2b3a46; --text:#e8eef5; --muted:#96a7b8; --flag:#ffb000; --watch:#58a6ff; --normal:#7d8590; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; height:100vh; display:grid; grid-template-columns:360px 1fr; background:var(--bg); color:var(--text); font:13px/1.4 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }}
aside {{ min-width:0; border-right:1px solid var(--line); background:var(--panel); display:flex; flex-direction:column; overflow:hidden; }}
header {{ padding:12px 14px; border-bottom:1px solid var(--line); display:grid; gap:9px; }}
h1 {{ margin:0; font-size:16px; letter-spacing:0; }}
.toolbar {{ display:grid; grid-template-columns:1fr 1fr; gap:6px; }}
.toolbar.three {{ grid-template-columns:1fr 1fr 1fr; }}
button {{ border:1px solid var(--line); background:#15202b; color:var(--text); border-radius:7px; min-height:32px; padding:0 8px; font:inherit; font-size:12px; cursor:pointer; }}
button.active {{ background:#2b6fd6; border-color:#4d8eff; }}
button:disabled {{ opacity:.45; cursor:not-allowed; }}
.hint {{ color:var(--muted); font-size:11px; }}
#list {{ overflow:auto; padding:8px; display:grid; gap:6px; }}
.triad {{ text-align:left; min-height:64px; border:1px solid var(--line); background:var(--panel2); border-radius:8px; padding:8px; display:grid; gap:3px; cursor:pointer; }}
.triad.active {{ outline:2px solid #78a9ff; }}
.topline {{ display:flex; align-items:center; justify-content:space-between; gap:8px; }}
.id {{ font-weight:700; font-size:13px; }}
.status {{ font-size:11px; padding:2px 6px; border-radius:999px; color:#071018; font-weight:700; }}
.status.flagged {{ background:var(--flag); }}
.status.watch {{ background:var(--watch); }}
.status.normal {{ background:var(--normal); color:#fff; }}
.meta {{ color:var(--muted); font-size:12px; }}
main {{ display:grid; grid-template-rows:auto 1fr; min-width:0; min-height:0; }}
#info {{ min-width:0; padding:10px 14px; border-bottom:1px solid var(--line); background:#0d131b; display:flex; justify-content:space-between; gap:16px; align-items:center; }}
#info strong {{ font-size:14px; }}
#canvasWrap {{ position:relative; min-width:0; min-height:0; overflow:hidden; }}
canvas {{ width:100%; height:100%; display:block; background:#070b10; }}
#legend {{ position:absolute; left:14px; bottom:14px; background:rgba(5,8,12,.72); border:1px solid rgba(255,255,255,.12); border-radius:8px; padding:8px 10px; display:grid; gap:5px; backdrop-filter:blur(8px); }}
.leg {{ display:grid; grid-template-columns:12px auto; gap:7px; align-items:center; }}
.sw {{ width:12px; height:12px; border-radius:3px; }}
#hud {{ position:absolute; left:14px; top:14px; max-width:min(640px,calc(100% - 28px)); background:rgba(5,8,12,.72); border:1px solid rgba(255,255,255,.12); border-radius:8px; color:#dce7f7; padding:8px 10px; font-size:12px; line-height:1.45; }}
</style>
</head>
<body>
<aside>
<header>
<h1>{title}</h1>
<div class="toolbar">
<button id="prevFlag">Prev flagged</button>
<button id="nextFlag">Next flagged</button>
</div>
<div class="toolbar">
<button id="flagOnly" class="active">Flagged only</button>
<button id="showAll">Show all</button>
</div>
<div class="toolbar">
<button id="sliceColor" class="active">Color by slice</button>
<button id="componentColor">Color by component</button>
</div>
<div class="toolbar">
<button id="focusBtn" class="active">Focus suspicious</button>
<button id="viewBtn">3D triad view</button>
</div>
<div class="hint" id="summary"></div>
</header>
<div id="list"></div>
</aside>
<main>
<div id="info"><div><strong id="triadTitle"></strong><div class="meta" id="triadMeta"></div></div><div class="meta" id="triadReason"></div></div>
<div id="canvasWrap"><canvas id="canvas"></canvas><div id="hud"></div><div id="legend"></div></div>
</main>
<script>
const DATA = {data};
let filteredOnly = true;
let colorMode = "slice";
let focus = true;
let mode3d = false;
let active = 0;
const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
const sliceColors = {{prev:"#4ea1ff", mid:"#ffb000", next:"#2fd17c"}};
const componentPalette = ["#58a6ff","#ff7b72","#a371f7","#56d364","#e3b341","#f778ba","#39c5cf","#d2a8ff","#ffa657","#79c0ff"];
function triads() {{ return filteredOnly ? DATA.triads.filter(t => t.status !== "normal") : DATA.triads; }}
function resize() {{ const dpr = devicePixelRatio || 1; const w = Math.max(1, Math.floor(canvas.clientWidth*dpr)); const h = Math.max(1, Math.floor(canvas.clientHeight*dpr)); if (canvas.width !== w || canvas.height !== h) {{ canvas.width=w; canvas.height=h; }} }}
function esc(s) {{
  const map = {{"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#039;"}};
  return String(s).replace(/[&<>"']/g, m=>map[m]);
}}
function renderList() {{
  const list = document.getElementById("list"); const rows = triads();
  list.innerHTML = rows.map((t,i)=>`<div class="triad ${{i===active?'active':''}}" data-i="${{i}}"><div class="topline"><span class="id">${{esc(t.triad_id.replaceAll('_','-'))}}</span><span class="status ${{t.status}}">${{t.status}}</span></div><div class="meta">mid=SL${{t.mid_slice}} · rank${{t.component_rank}} · score=${{t.score}}</div><div class="meta">${{esc(t.reason)}}</div></div>`).join("");
  for (const el of list.querySelectorAll(".triad")) el.onclick = () => {{ active = Number(el.dataset.i); render(); }};
  document.getElementById("summary").textContent = `${{rows.length}} shown / ${{DATA.triads.length}} triad-component rows · rendered points are capped per triad`;
}}
function bounds(points) {{
  let xs=[], ys=[];
  for (const p of points) {{
    let x=p.x, y=p.y;
    if (mode3d) {{ x = p.x + (Number(p.slice)-Number(points[0].slice))*90; y = p.y - (Number(p.slice)-Number(points[0].slice))*28; }}
    xs.push(x); ys.push(y);
  }}
  const minX=Math.min(...xs), maxX=Math.max(...xs), minY=Math.min(...ys), maxY=Math.max(...ys);
  const mx=(maxX-minX)*0.08+1, my=(maxY-minY)*0.08+1;
  return [minX-mx,maxX+mx,minY-my,maxY+my];
}}
function project(p,b,w,h) {{
  let x=p.x, y=p.y;
  if (mode3d) {{ x = p.x + (Number(p.slice)-Number(current().prev_slice))*90; y = p.y - (Number(p.slice)-Number(current().prev_slice))*28; }}
  const [minX,maxX,minY,maxY]=b; const pad=36; const scale=Math.min((w-2*pad)/(maxX-minX), (h-2*pad)/(maxY-minY));
  const ox=(w-(maxX-minX)*scale)/2; const oy=(h-(maxY-minY)*scale)/2;
  return [ox+(x-minX)*scale, oy+(maxY-y)*scale];
}}
function current() {{ const rows=triads(); if (!rows.length) return null; active=Math.max(0,Math.min(active,rows.length-1)); return rows[active]; }}
function pointColor(t,p) {{
  if (colorMode === "component") return componentPalette[Math.abs(Number(p.component_rank)||0) % componentPalette.length];
  if (String(p.slice) === String(t.prev_slice)) return sliceColors.prev;
  if (String(p.slice) === String(t.mid_slice)) return sliceColors.mid;
  return sliceColors.next;
}}
function drawTriad(t) {{
  resize(); const w=canvas.width,h=canvas.height; ctx.clearRect(0,0,w,h);
  ctx.fillStyle="#070b10"; ctx.fillRect(0,0,w,h);
  if (!t || !t.points.length) return;
  const pts = focus ? t.points.filter(p => p.highlight || String(p.slice) !== String(t.mid_slice) || t.status === "normal") : t.points;
  const b = bounds(pts.length ? pts : t.points);
  for (const p of t.points) {{
    const [x,y]=project(p,b,w,h);
    const hi = p.highlight;
    ctx.globalAlpha = hi ? 0.98 : (focus ? 0.16 : 0.34);
    ctx.fillStyle = hi ? "#fff3b0" : pointColor(t,p);
    const r = hi ? 2.4*(devicePixelRatio||1) : 1.4*(devicePixelRatio||1);
    ctx.beginPath(); ctx.arc(x,y,r,0,Math.PI*2); ctx.fill();
  }}
  const c = project({{x:t.mid_centroid_x,y:t.mid_centroid_y,slice:t.mid_slice}}, b, w, h);
  const e = project({{x:t.mid_expected_x,y:t.mid_expected_y,slice:t.mid_slice}}, b, w, h);
  ctx.globalAlpha = 1; ctx.strokeStyle="#ffdf5d"; ctx.lineWidth=2*(devicePixelRatio||1);
  ctx.beginPath(); ctx.moveTo(c[0],c[1]); ctx.lineTo(e[0],e[1]); ctx.stroke();
  ctx.fillStyle="#ffdf5d"; ctx.beginPath(); ctx.arc(c[0],c[1],5*(devicePixelRatio||1),0,Math.PI*2); ctx.fill();
  ctx.fillStyle="#ffffff"; ctx.beginPath(); ctx.rect(e[0]-4*(devicePixelRatio||1),e[1]-4*(devicePixelRatio||1),8*(devicePixelRatio||1),8*(devicePixelRatio||1)); ctx.fill();
}}
function renderInfo(t) {{
  document.getElementById("triadTitle").textContent = t ? `${{t.triad_id.replaceAll('_','-')}} · ${{t.status}}` : "No triads";
  document.getElementById("triadMeta").textContent = t ? `mid=SL${{t.mid_slice}} · component rank${{t.component_rank}} · mid component ${{t.mid_component_id}} · ${{t.mid_points.toLocaleString()}} cells · score ${{t.score}}` : "";
  document.getElementById("triadReason").textContent = t ? t.reason : "";
  document.getElementById("hud").innerHTML = t ? `Highlight: SL${{t.mid_slice}} rank${{t.component_rank}}<br>Line: current mid centroid to prev/next consensus<br>Mode: ${{mode3d ? '3D triad view' : '2D overlay'}} · Color: ${{colorMode}}` : "";
  document.getElementById("legend").innerHTML = `<div class="leg"><span class="sw" style="background:${{sliceColors.prev}}"></span><span>prev slice</span></div><div class="leg"><span class="sw" style="background:${{sliceColors.mid}}"></span><span>mid slice</span></div><div class="leg"><span class="sw" style="background:${{sliceColors.next}}"></span><span>next slice</span></div><div class="leg"><span class="sw" style="background:#fff3b0"></span><span>suspicious mid component</span></div>`;
}}
function render() {{ renderList(); const t=current(); renderInfo(t); drawTriad(t); }}
function moveFlag(delta) {{ const rows=triads(); if (!rows.length) return; let i=active; for (let n=0;n<rows.length;n++) {{ i=(i+delta+rows.length)%rows.length; if (rows[i].status !== "normal") {{ active=i; break; }} }} render(); }}
document.getElementById("prevFlag").onclick=()=>moveFlag(-1);
document.getElementById("nextFlag").onclick=()=>moveFlag(1);
document.getElementById("flagOnly").onclick=()=>{{ filteredOnly=true; active=0; document.getElementById("flagOnly").classList.add("active"); document.getElementById("showAll").classList.remove("active"); render(); }};
document.getElementById("showAll").onclick=()=>{{ filteredOnly=false; active=0; document.getElementById("showAll").classList.add("active"); document.getElementById("flagOnly").classList.remove("active"); render(); }};
document.getElementById("sliceColor").onclick=()=>{{ colorMode="slice"; document.getElementById("sliceColor").classList.add("active"); document.getElementById("componentColor").classList.remove("active"); render(); }};
document.getElementById("componentColor").onclick=()=>{{ colorMode="component"; document.getElementById("componentColor").classList.add("active"); document.getElementById("sliceColor").classList.remove("active"); render(); }};
document.getElementById("focusBtn").onclick=()=>{{ focus=!focus; document.getElementById("focusBtn").classList.toggle("active", focus); render(); }};
document.getElementById("viewBtn").onclick=()=>{{ mode3d=!mode3d; document.getElementById("viewBtn").textContent = mode3d ? "2D overlay" : "3D triad view"; render(); }};
window.onresize=()=>drawTriad(current());
render();
</script>
</body>
</html>"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--components", required=True, type=Path)
    parser.add_argument("--scan", required=True, type=Path)
    parser.add_argument("--output-html", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--title", default="Triad component review")
    parser.add_argument("--x-col", default="manual_x")
    parser.add_argument("--y-col", default="manual_y")
    parser.add_argument("--z-col", default="manual_z")
    parser.add_argument("--slice-col", default="sl_number")
    parser.add_argument("--celltype-col", default="celltype")
    parser.add_argument("--max-triads", type=int, default=120)
    parser.add_argument("--max-points-per-triad", type=int, default=12000)
    parser.add_argument("--max-highlight-points-per-triad", type=int, default=8000)
    parser.add_argument("--include-status", default="flagged,watch,normal")
    args = parser.parse_args()

    statuses = {s.strip() for s in args.include_status.split(",") if s.strip()}
    scan_rows = read_scan(args.scan, statuses, args.max_triads)
    payload_triads, summary = load_points(args.coordinates, args.components, scan_rows, args)
    payload = {
        "title": args.title,
        "coordinates": str(args.coordinates),
        "components": str(args.components),
        "scan": str(args.scan),
        "triads": payload_triads,
        "summary": summary,
        "z_spacing": 40,
    }
    args.output_html.parent.mkdir(parents=True, exist_ok=True)
    args.output_html.write_text(render_html(payload), encoding="utf-8")
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.write_text(json.dumps({**summary, "html": str(args.output_html)}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
