#!/usr/bin/env python3
"""Build a portable interactive dashboard for a tracked Spateo 4D run."""

from __future__ import annotations

import argparse
import html
import json
import sys
from pathlib import Path
from typing import Any


FEATURE_CANDIDATES = (
    "displacement",
    "velocity",
    "acceleration",
    "curl",
    "divergence",
    "torsion",
    "curvature",
)


def json_for_script(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":")).replace(
        "</", "<\\/"
    )


def stratified_indices(labels: list[str], maximum: int, seed: int) -> list[int]:
    import numpy as np

    count = len(labels)
    if maximum <= 0 or count <= maximum:
        return list(range(count))
    rng = np.random.default_rng(seed)
    array = np.asarray(labels)
    groups = {label: np.flatnonzero(array == label) for label in np.unique(array)}
    chosen: list[int] = []
    for indices in groups.values():
        size = max(1, int(round(maximum * len(indices) / count)))
        chosen.extend(
            rng.choice(indices, size=min(size, len(indices)), replace=False).tolist()
        )
    if len(chosen) > maximum:
        mandatory = [int(rng.choice(indices)) for indices in groups.values()]
        if len(mandatory) >= maximum:
            return sorted(mandatory[:maximum])
        remaining = np.setdiff1d(
            np.asarray(chosen), np.asarray(mandatory), assume_unique=False
        )
        chosen = (
            mandatory
            + rng.choice(
                remaining, size=maximum - len(mandatory), replace=False
            ).tolist()
        )
    return sorted(set(chosen))


def finite_rows(values: Any) -> list[Any]:
    import numpy as np

    array = np.asarray(values, dtype=float)
    if array.ndim == 1:
        return [None if not np.isfinite(value) else float(value) for value in array]
    result: list[Any] = []
    for row in array:
        result.append(
            [None if not np.isfinite(value) else float(value) for value in row]
        )
    return result


def scalar_feature(
    adata: Any,
    key: str,
    indices: list[int],
    *,
    vector_key: str,
    mapped_key: str,
    spatial_key: str,
) -> list[float | None] | None:
    import numpy as np

    if key == "displacement":
        if vector_key in adata.obsm:
            values = np.linalg.norm(np.asarray(adata.obsm[vector_key])[indices], axis=1)
        elif mapped_key in adata.obsm:
            mapped = np.asarray(adata.obsm[mapped_key])[indices]
            coords = np.asarray(adata.obsm[spatial_key])[indices]
            values = np.linalg.norm(mapped - coords, axis=1)
        else:
            return None
    elif key in adata.obs:
        try:
            values = np.asarray(adata.obs[key].iloc[indices], dtype=float)
        except (TypeError, ValueError):
            return None
    elif key in adata.obsm:
        raw = np.asarray(adata.obsm[key])[indices]
        values = np.abs(raw) if raw.ndim == 1 else np.linalg.norm(raw, axis=1)
    else:
        return None
    values = np.asarray(values, dtype=float).ravel()
    return [None if not np.isfinite(value) else float(value) for value in values]


def read_dataset(
    path: Path,
    *,
    spatial_key: str,
    groupby: str,
    maximum: int,
    seed: int,
    features: list[str] | None = None,
    vector_key: str | None = None,
    mapped_key: str | None = None,
) -> dict[str, Any]:
    import anndata as ad
    import numpy as np

    adata = ad.read_h5ad(path, backed="r")
    try:
        if spatial_key not in adata.obsm:
            raise KeyError(
                f"Missing obsm[{spatial_key!r}] in {path}; available={list(adata.obsm.keys())}"
            )
        if groupby not in adata.obs:
            raise KeyError(
                f"Missing obs[{groupby!r}] in {path}; available={list(adata.obs.columns)}"
            )
        all_coords = np.asarray(adata.obsm[spatial_key], dtype=float)
        if all_coords.ndim != 2 or all_coords.shape[1] != 3:
            raise ValueError(
                f"Expected 3-D coordinates in {path}, got {all_coords.shape}"
            )
        all_labels = adata.obs[groupby].astype(str).tolist()
        indices = stratified_indices(all_labels, maximum, seed)
        coords = all_coords[indices]
        labels = [all_labels[index] for index in indices]
        payload: dict[str, Any] = {
            "source": str(path.resolve()),
            "coords": finite_rows(coords),
            "labels": labels,
            "obsNames": [str(adata.obs_names[index]) for index in indices],
            "displayed": len(indices),
            "total": len(all_labels),
        }
        if vector_key and vector_key in adata.obsm:
            vectors = np.asarray(adata.obsm[vector_key])[indices]
            if vectors.shape == coords.shape:
                payload["vectors"] = finite_rows(vectors)
        if mapped_key and mapped_key in adata.obsm:
            mapped = np.asarray(adata.obsm[mapped_key])[indices]
            if mapped.shape == coords.shape:
                payload["mapped"] = finite_rows(mapped)
        if features is not None and vector_key and mapped_key:
            payload["features"] = {}
            for key in dict.fromkeys(features):
                values = scalar_feature(
                    adata,
                    key,
                    indices,
                    vector_key=vector_key,
                    mapped_key=mapped_key,
                    spatial_key=spatial_key,
                )
                if values is not None:
                    payload["features"][key] = values
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()
    return payload


def read_mapping_summary(path: Path | None) -> list[dict[str, Any]]:
    if path is None or not path.is_file():
        return []
    import pandas as pd

    frame = pd.read_csv(path)
    if frame.empty:
        return []
    group_column = next(
        (
            column
            for column in frame.columns
            if column.lower() in {"lineage", "anno", "group", "groupby"}
        ),
        frame.columns[0],
    )
    value_column = next(
        (column for column in ("median", "mean", "displacement") if column in frame),
        None,
    )
    if value_column is None:
        return []
    rows: list[dict[str, Any]] = []
    for _, row in frame.iterrows():
        rows.append(
            {
                "group": str(row[group_column]),
                "value": float(row[value_column]),
                "n": int(row["n"]) if "n" in frame and row["n"] == row["n"] else None,
                "p05": (
                    float(row["p05"])
                    if "p05" in frame and row["p05"] == row["p05"]
                    else None
                ),
                "p95": (
                    float(row["p95"])
                    if "p95" in frame and row["p95"] == row["p95"]
                    else None
                ),
            }
        )
    return rows


def read_metric_summary(path: Path | None) -> list[dict[str, Any]]:
    if path is None or not path.is_file():
        return []
    import pandas as pd

    frame = pd.read_csv(path)
    keep = [
        column
        for column in (
            "metric",
            "anno",
            "lineage",
            "group",
            "n",
            "mean",
            "median",
            "p05",
            "p95",
        )
        if column in frame
    ]
    return frame[keep].where(frame[keep].notna(), None).to_dict(orient="records")


def read_glm_tables(directory: Path | None, limit: int = 100) -> dict[str, Any]:
    if directory is None or not directory.is_dir():
        return {}
    import pandas as pd

    tables: dict[str, Any] = {}
    for path in sorted(directory.glob("glm_degs_*.csv")):
        frame = pd.read_csv(path)
        if frame.empty:
            continue
        feature = path.stem.removeprefix("glm_degs_")
        gene_column = next(
            (
                column
                for column in frame.columns
                if column.lower() in {"gene", "genes", "gene_id", "index"}
            ),
            frame.columns[0],
        )
        q_column = next(
            (
                column
                for column in frame.columns
                if column.lower()
                in {"qval", "q_value", "qvalue", "padj", "adjusted_p_value"}
            ),
            None,
        )
        score_column = next(
            (
                column
                for column in frame.columns
                if column.lower()
                in {"llf", "log-likelihood", "score", "wald", "deviance", "coefficient"}
            ),
            None,
        )
        if q_column:
            frame = frame.sort_values(q_column, ascending=True, na_position="last")
        rows = []
        for _, row in frame.head(limit).iterrows():
            rows.append(
                {
                    "gene": str(row[gene_column]),
                    "q": (
                        None
                        if not q_column or row[q_column] != row[q_column]
                        else float(row[q_column])
                    ),
                    "score": (
                        None
                        if not score_column or row[score_column] != row[score_column]
                        else float(row[score_column])
                    ),
                }
            )
        tables[feature] = rows
    return tables


def plotly_script(use_cdn: bool) -> str:
    if not use_cdn:
        try:
            from plotly.offline.offline import get_plotlyjs

            return f"<script>{get_plotlyjs()}</script>"
        except ImportError as exc:
            raise RuntimeError(
                "Offline dashboard requires plotly; install it or explicitly select CDN mode"
            ) from exc
    return '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'


def build_payload(args: argparse.Namespace) -> dict[str, Any]:
    source = read_dataset(
        args.adata,
        spatial_key=args.spatial_key,
        groupby=args.groupby,
        maximum=args.max_points,
        seed=args.seed,
        features=args.features,
        vector_key=args.vector_key,
        mapped_key=args.mapped_key,
    )
    target = None
    if args.target_adata:
        target = read_dataset(
            args.target_adata,
            spatial_key=args.target_spatial_key or args.spatial_key,
            groupby=args.groupby,
            maximum=args.max_target_points,
            seed=args.seed + 1,
        )
    manifest = read_object(args.manifest) if args.manifest else {}
    return {
        "runId": manifest.get("run_id", args.adata.stem),
        "runStatus": manifest.get("status", "untracked"),
        "groupby": args.groupby,
        "source": source,
        "target": target,
        "mappingSummary": read_mapping_summary(args.mapping_summary),
        "metricSummary": read_metric_summary(args.metric_summary),
        "glm": read_glm_tables(args.glm_dir),
        "stages": [
            {
                "name": name,
                "status": stage.get("status", "unknown"),
                "reusedFrom": stage.get("reused_from"),
            }
            for name, stage in manifest.get("stages", {}).items()
        ],
        "vectorLimit": args.max_vectors,
        "defaultFeature": args.default_feature,
    }


def read_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    return value if isinstance(value, dict) else {}


def render_html(payload: dict[str, Any], use_cdn: bool, title: str | None) -> str:
    page_title = html.escape(title or f"Spateo 4D run {payload['runId']}")
    data = json_for_script(payload)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{page_title}</title>
{plotly_script(use_cdn)}
<style>
:root {{ color-scheme:light dark; --bg:#f7f8fa; --panel:#fff; --fg:#182230; --muted:#667085; --line:#d0d5dd; --accent:#2563eb; --source:#7c3aed; --target:#0f9d8b; --mapped:#f59e0b; --ok:#15803d; --warn:#b45309; --bad:#b42318; }}
@media (prefers-color-scheme:dark) {{ :root {{ --bg:#0f172a; --panel:#172033; --fg:#f2f4f7; --muted:#cbd5e1; --line:#475467; --accent:#60a5fa; --source:#a78bfa; --target:#2dd4bf; --mapped:#fbbf24; --ok:#4ade80; --warn:#fbbf24; --bad:#fb7185; }} }}
* {{ box-sizing:border-box; }}
body {{ margin:0; overflow-x:hidden; background:var(--bg); color:var(--fg); font:14px/1.45 system-ui,-apple-system,sans-serif; }}
body > #js-plotly-tester {{ position:fixed!important; left:0!important; top:0!important; width:1px!important; height:1px!important; overflow:hidden!important; }}
.shell {{ max-width:1440px; margin:auto; padding:16px; }}
.controls {{ display:flex; flex-wrap:wrap; gap:10px; align-items:end; margin-bottom:10px; }}
.field {{ display:grid; gap:4px; min-width:145px; }} label,.muted {{ color:var(--muted); font-size:12px; }}
select,input,button {{ font:inherit; color:var(--fg); background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:7px 9px; }}
button {{ cursor:pointer; }} button.primary {{ color:white; background:var(--accent); border-color:var(--accent); }}
.stats {{ display:flex; flex-wrap:wrap; gap:16px; padding:8px 0 12px; color:var(--muted); }} .stats strong {{ color:var(--fg); font-weight:500; }}
.layout {{ display:grid; grid-template-columns:minmax(0,2.2fr) minmax(300px,1fr); gap:12px; align-items:start; }}
.panel {{ min-width:0; background:var(--panel); border:1px solid var(--line); border-radius:12px; padding:10px; }}
#spatial {{ height:680px; }} #mapping {{ height:250px; }}
.side {{ display:grid; gap:12px; }}
.stage-list {{ display:flex; flex-wrap:wrap; gap:6px; }}
.stage {{ border:1px solid var(--line); border-radius:999px; padding:3px 8px; font-size:12px; }}
.stage.completed,.stage.reused {{ color:var(--ok); }} .stage.running,.stage.pending {{ color:var(--warn); }} .stage.failed {{ color:var(--bad); }}
table {{ width:100%; border-collapse:collapse; }} th,td {{ text-align:left; padding:6px; border-bottom:1px solid var(--line); }} th {{ color:var(--muted); font-weight:500; }}
#status {{ min-height:20px; color:var(--muted); }}
@media(max-width:900px) {{ .layout {{ grid-template-columns:1fr; }} #spatial {{ height:520px; }} }}
@media(max-width:520px) {{ .shell {{ padding:10px; }} .field {{ min-width:calc(50% - 6px); flex:1; }} #spatial {{ height:430px; }} }}
</style>
</head>
<body>
<main class="shell">
  <div class="controls">
    <div class="field"><label for="view">View</label><select id="view"><option value="flow">Developmental flow</option><option value="overlay">Stage overlay</option><option value="source">Source stage</option><option value="target">Target stage</option></select></div>
    <div class="field"><label for="group">Group</label><select id="group"></select></div>
    <div class="field"><label for="feature">Source color</label><select id="feature"></select></div>
    <div class="field"><label for="vectorScale">Vector scale <span id="vectorValue">1.0</span></label><input id="vectorScale" type="range" min="0" max="5" step="0.1" value="1"></div>
    <div class="field"><label for="pointSize">Point size <span id="pointValue">3</span></label><input id="pointSize" type="range" min="1" max="9" step="1" value="3"></div>
    <button id="ask" class="primary" type="button">Ask AI about selection</button>
  </div>
  <div class="stats"><span>Run <strong id="run"></strong></span><span>Status <strong id="runStatus"></strong></span><span>Shown <strong id="shown"></strong></span><span>Source median displacement <strong id="median">—</strong></span></div>
  <div class="layout">
    <section class="panel"><div id="spatial" role="img" aria-label="Interactive 3D source, target, mapped coordinates, and vector field"></div></section>
    <aside class="side">
      <section class="panel"><div class="stage-list" id="stages" aria-label="Pipeline stage status"></div></section>
      <section class="panel"><div id="mapping" role="img" aria-label="Mapped displacement summary by group"></div></section>
      <section class="panel">
        <div class="field"><label for="glmFeature">GLM feature</label><select id="glmFeature"></select></div>
        <table><thead><tr><th>Gene</th><th>q-value</th><th>Score</th></tr></thead><tbody id="glmRows"></tbody></table>
      </section>
      <div id="status" aria-live="polite"></div>
    </aside>
  </div>
</main>
<script>
const D={data};
const $=id=>document.getElementById(id);
const css=name=>getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const esc=value=>String(value).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');
const unique=values=>[...new Set(values)].sort((a,b)=>a.localeCompare(b));
const fmt=value=>value==null||!Number.isFinite(Number(value))?'—':Number(value).toPrecision(3);
const plotTheme=()=>({{font:{{color:css('--fg')}},paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)'}});
const group=$('group'),feature=$('feature'),view=$('view'),glmFeature=$('glmFeature');
const targetLabels=D.target?D.target.labels:[];
group.innerHTML=['All',...unique([...D.source.labels,...targetLabels])].map(value=>`<option value="${{esc(value)}}">${{esc(value)}}</option>`).join('');
const featureNames=Object.keys(D.source.features||{{}});
feature.innerHTML=(featureNames.length?featureNames:['unavailable']).map(value=>`<option value="${{esc(value)}}">${{esc(value)}}</option>`).join('');
if(featureNames.includes(D.defaultFeature))feature.value=D.defaultFeature;
const glmNames=Object.keys(D.glm||{{}});
glmFeature.innerHTML=(glmNames.length?glmNames:['unavailable']).map(value=>`<option value="${{esc(value)}}">${{esc(value)}}</option>`).join('');
if(!D.target)view.querySelector('option[value="overlay"]').disabled=view.querySelector('option[value="target"]').disabled=true;
$('run').textContent=D.runId;$('runStatus').textContent=D.runStatus;
$('stages').innerHTML=(D.stages.length?D.stages:[{{name:'untracked',status:'unknown'}}]).map(stage=>`<span class="stage ${{esc(stage.status)}}">${{esc(stage.name)}} · ${{esc(stage.status)}}${{stage.reusedFrom?' ↩ '+esc(stage.reusedFrom):''}}</span>`).join('');
function indices(dataset){{const selected=group.value;return dataset.labels.map((label,index)=>selected==='All'||label===selected?index:-1).filter(index=>index>=0);}}
function median(values){{const clean=values.filter(Number.isFinite).sort((a,b)=>a-b);if(!clean.length)return null;const middle=Math.floor(clean.length/2);return clean.length%2?clean[middle]:(clean[middle-1]+clean[middle])/2;}}
function pointTrace(dataset,idx,name,color,opacity=1){{return {{type:'scatter3d',mode:'markers',name,x:idx.map(i=>dataset.coords[i][0]),y:idx.map(i=>dataset.coords[i][1]),z:idx.map(i=>dataset.coords[i][2]),text:idx.map(i=>`${{dataset.obsNames[i]}}<br>${{D.groupby}}: ${{dataset.labels[i]}}`),hoverinfo:'text',marker:{{size:Number($('pointSize').value),color,opacity}}}};}}
function renderSpatial(){{
  const sourceIndex=indices(D.source),targetIndex=D.target?indices(D.target):[],key=feature.value,values=(D.source.features||{{}})[key]||new Array(D.source.labels.length).fill(0),mode=view.value,scale=Number($('vectorScale').value),traces=[];
  if(mode==='source'||mode==='flow'||mode==='overlay'){{const trace=pointTrace(D.source,sourceIndex,'Source',sourceIndex.map(i=>values[i]));trace.marker.colorscale='Viridis';trace.marker.showscale=true;trace.marker.colorbar={{title:key}};trace.text=sourceIndex.map(i=>`${{D.source.obsNames[i]}}<br>${{D.groupby}}: ${{D.source.labels[i]}}<br>${{key}}: ${{fmt(values[i])}}`);traces.push(trace);}}
  if((mode==='target'||mode==='overlay')&&D.target)traces.push(pointTrace(D.target,targetIndex,'Target',css('--target'),.75));
  if(mode==='flow'&&D.source.mapped){{traces.push({{type:'scatter3d',mode:'markers',name:'Mapped source',x:sourceIndex.map(i=>D.source.mapped[i][0]),y:sourceIndex.map(i=>D.source.mapped[i][1]),z:sourceIndex.map(i=>D.source.mapped[i][2]),hoverinfo:'skip',marker:{{size:Math.max(1,Number($('pointSize').value)-1),color:css('--mapped'),opacity:.42}}}});}}
  if(mode==='flow'&&D.source.vectors&&scale>0){{const step=Math.max(1,Math.ceil(sourceIndex.length/D.vectorLimit)),vectorIndex=sourceIndex.filter((_,j)=>j%step===0);traces.push({{type:'cone',name:'Morphogenesis vectors',x:vectorIndex.map(i=>D.source.coords[i][0]),y:vectorIndex.map(i=>D.source.coords[i][1]),z:vectorIndex.map(i=>D.source.coords[i][2]),u:vectorIndex.map(i=>D.source.vectors[i][0]*scale),v:vectorIndex.map(i=>D.source.vectors[i][1]*scale),w:vectorIndex.map(i=>D.source.vectors[i][2]*scale),sizemode:'absolute',sizeref:1,showscale:false,opacity:.5,colorscale:[[0,'#64748b'],[1,'#64748b']]}});}}
  Plotly.react('spatial',traces,{{...plotTheme(),margin:{{l:0,r:0,t:0,b:0}},scene:{{aspectmode:'data',xaxis:{{title:'x',gridcolor:css('--line')}},yaxis:{{title:'y',gridcolor:css('--line')}},zaxis:{{title:'z',gridcolor:css('--line')}}}},legend:{{orientation:'h'}}}},{{responsive:true,displaylogo:false}});
  const sourceShown=mode==='target'?0:sourceIndex.length,targetShown=(mode==='target'||mode==='overlay')?targetIndex.length:0;$('shown').textContent=`${{sourceShown+targetShown}} / ${{D.source.total+(D.target?D.target.total:0)}}`;
  const displacement=(D.source.features||{{}}).displacement||[];$('median').textContent=fmt(median(sourceIndex.map(i=>displacement[i]).filter(Number.isFinite)));
}}
function renderMapping(){{const selected=group.value,rows=D.mappingSummary,colors=rows.map(row=>selected==='All'||row.group===selected?css('--source'):css('--line'));Plotly.react('mapping',[{{type:'bar',orientation:'h',y:rows.map(row=>row.group),x:rows.map(row=>row.value),text:rows.map(row=>row.n?`n=${{row.n}}`:''),hovertemplate:'%{{y}}<br>median=%{{x:.3g}}<br>%{{text}}<extra></extra>',marker:{{color:colors}}}}],{{...plotTheme(),margin:{{l:110,r:8,t:8,b:38}},xaxis:{{title:'Mapped displacement',gridcolor:css('--line')}},yaxis:{{automargin:true}}}},{{responsive:true,displaylogo:false}});}}
function renderGlm(){{const rows=(D.glm||{{}})[glmFeature.value]||[];$('glmRows').innerHTML=rows.slice(0,25).map(row=>`<tr><td>${{esc(row.gene)}}</td><td>${{fmt(row.q)}}</td><td>${{fmt(row.score)}}</td></tr>`).join('')||'<tr><td colspan="3">No GLM table available</td></tr>';}}
[group,feature,view].forEach(control=>control.addEventListener('change',()=>{{renderSpatial();renderMapping();}}));glmFeature.addEventListener('change',renderGlm);
$('vectorScale').addEventListener('input',event=>{{$('vectorValue').textContent=event.target.value;renderSpatial();}});$('pointSize').addEventListener('input',event=>{{$('pointValue').textContent=event.target.value;renderSpatial();}});
$('ask').addEventListener('click',async()=>{{const prompt=`For tracked Spateo run ${{D.runId}}, inspect view=${{view.value}}, ${{D.groupby}}=${{group.value}}, feature=${{feature.value}}. Use the manifest and saved outputs to explain the pattern; distinguish display changes from parameters that require recomputation.`;try{{if(window.openai?.sendFollowUpMessage){{await window.openai.sendFollowUpMessage({{prompt,title:'Inspect Spateo 4D selection'}});$('status').textContent='Sent to ChatGPT.';}}else{{await navigator.clipboard.writeText(prompt);$('status').textContent='Prompt copied.';}}}}catch(error){{$('status').textContent='Could not send or copy the prompt.';}}}});
renderSpatial();renderMapping();renderGlm();
</script>
</body>
</html>"""


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    command.add_argument("--adata", type=Path, required=True)
    command.add_argument("--target-adata", type=Path)
    command.add_argument("--manifest", type=Path)
    command.add_argument("--mapping-summary", type=Path)
    command.add_argument("--metric-summary", type=Path)
    command.add_argument("--glm-dir", type=Path)
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--title")
    command.add_argument("--spatial-key", default="3d_align_spatial")
    command.add_argument("--target-spatial-key")
    command.add_argument("--groupby", default="anno")
    command.add_argument("--vector-key", default="V_cells_mapping")
    command.add_argument("--mapped-key", default="X_cells_mapping")
    command.add_argument("--features", nargs="+", default=list(FEATURE_CANDIDATES))
    command.add_argument("--default-feature", default="displacement")
    command.add_argument("--max-points", type=int, default=6000)
    command.add_argument("--max-target-points", type=int, default=4000)
    command.add_argument("--max-vectors", type=int, default=1200)
    command.add_argument("--seed", type=int, default=0)
    command.add_argument("--cdn", action="store_true")
    return command


def main() -> int:
    try:
        args = parser().parse_args()
        payload = build_payload(args)
        output = args.output.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(render_html(payload, args.cdn, args.title), encoding="utf-8")
        print(
            json.dumps(
                {
                    "output": str(output),
                    "run_id": payload["runId"],
                    "source_displayed": payload["source"]["displayed"],
                    "source_total": payload["source"]["total"],
                    "target_displayed": (
                        payload["target"]["displayed"] if payload["target"] else 0
                    ),
                    "features": sorted(payload["source"].get("features", {})),
                },
                indent=2,
            )
        )
        return 0
    except Exception as exc:  # noqa: BLE001 - CLI boundary
        print(
            json.dumps(
                {"error": type(exc).__name__, "message": str(exc)}, ensure_ascii=False
            ),
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
