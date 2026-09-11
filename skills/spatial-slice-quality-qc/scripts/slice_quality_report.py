"""Skill-local HTML reports for serial-slice quality control.

Scientific metrics and decisions come from spateo.preprocessing.slice_quality;
this module owns only presentation and collection-report orchestration.
"""

from __future__ import annotations

import html
import json
import os
import re
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Union

import numpy as np
import pandas as pd
from spateo.preprocessing.slice_quality import (
    _REPORT_STATISTICS,
    SliceSeriesResult,
    _jsonable,
    _natural_key,
    _report_measurement_context,
    _summarize_report_statistics,
    evaluate_slice_calls,
)


def render_slice_quality_report_from_outputs(
    input_dir: Union[str, Path],
    output_html: Optional[Union[str, Path]] = None,
    *,
    title: Optional[str] = None,
) -> str:
    """Render HTML locally from remote-produced CSV and display payload files."""
    input_path = Path(input_dir).expanduser().resolve()
    metrics_path = input_path / "slice_quality_metrics.csv"
    payload_path = input_path / "slice_quality_display_payload.json"
    if not metrics_path.exists() or not payload_path.exists():
        raise FileNotFoundError(
            "Rendering requires slice_quality_metrics.csv and slice_quality_display_payload.json"
        )
    metrics = pd.read_csv(metrics_path, dtype={"slice_id": str})
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    result = SliceSeriesResult(
        metrics=metrics,
        point_samples=payload.get("point_samples", {}),
        provenance=payload.get("provenance", {}),
    )
    destination = (
        Path(output_html).expanduser().resolve()
        if output_html is not None
        else input_path / "slice_quality_report.html"
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        _render_report(
            result,
            title=title or payload.get("title") or "Slice QC",
            ground_truth=payload.get("ground_truth") or None,
        ),
        encoding="utf-8",
    )
    from slice_quality_preregistration import attach_preregistration_html

    attach_preregistration_html(input_path, destination)
    return str(destination)


def render_binary_slice_quality_report_from_outputs(
    input_dir: Union[str, Path],
    output_html: Optional[Union[str, Path]] = None,
    *,
    title: Optional[str] = None,
) -> str:
    """Render every slice with a final operational keep/exclude action."""
    input_path = Path(input_dir).expanduser().resolve()
    metrics_path = input_path / "slice_quality_metrics.csv"
    binary_path = input_path / "slice_quality_binary_audit.csv"
    payload_path = input_path / "slice_quality_display_payload.json"
    if (
        not metrics_path.exists()
        or not binary_path.exists()
        or not payload_path.exists()
    ):
        raise FileNotFoundError(
            "Binary rendering requires slice_quality_metrics.csv, "
            "slice_quality_binary_audit.csv and slice_quality_display_payload.json"
        )
    all_metrics = pd.read_csv(metrics_path, dtype={"slice_id": str})
    binary = pd.read_csv(binary_path, dtype={"slice_id": str})
    final_values = binary["final_call"].astype("string")
    invalid = final_values.isna() | ~final_values.isin(["keep", "exclude"])
    if invalid.any():
        raise ValueError(
            "Complete binary rendering requires a keep/exclude final_call for every slice. "
            "Republish with unresolved_action='keep'."
        )
    final_rows = binary.copy()
    final_rows["recommendation"] = final_values.astype(str)
    review_resolved = (
        final_rows.get("decision_basis", pd.Series(dtype=str))
        .astype(str)
        .isin(
            [
                "review_resolved_keep",
                "adaptive_multidomain_resolution",
                "tiered_multidomain_resolution",
            ]
        )
    )
    final_rows.drop(
        columns=[
            "internal_recommendation",
            "threshold_band",
            "threshold_triage_call",
            "threshold_triage_reason",
            "review_resolution",
            "review_resolution_reason",
        ],
        errors="ignore",
        inplace=True,
    )
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    resolved_title = title or payload.get("title") or "Spateo referee"
    records = [_jsonable(record) for record in final_rows.to_dict("records")]
    order = all_metrics.sort_values("slice_index")["slice_id"].astype(str).tolist()
    page_payload = {
        "title": resolved_title,
        "records": records,
        "order": order,
        "points": _jsonable(payload.get("point_samples", {})),
        "provenance": _jsonable(payload.get("provenance", {})),
        "statisticDefinitions": _jsonable(_REPORT_STATISTICS),
        "statisticSummary": _jsonable(_summarize_report_statistics(final_rows)),
        "measurementContext": _jsonable(
            _report_measurement_context(payload.get("provenance", {}))
        ),
        "n_total": int(len(binary)),
        "n_certified": int(
            final_rows.get("decision_basis", pd.Series(dtype=str))
            .astype(str)
            .eq("independently_certified")
            .sum()
        ),
        "n_review_resolved": int(review_resolved.sum()),
        "n_operational_default": int(
            final_rows.get("decision_basis", pd.Series(dtype=str))
            .astype(str)
            .isin(["conservative_keep_default", "review_resolved_keep"])
            .sum()
        ),
    }
    payload_json = json.dumps(page_payload, ensure_ascii=False).replace("</", "<\\/")
    destination = (
        Path(output_html).expanduser().resolve()
        if output_html is not None
        else input_path / "slice_quality_report.html"
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    template = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title><style>
:root{--bg:#f4f7fb;--card:#fff;--ink:#122033;--muted:#607086;--line:#dce4ef;--keep:#1f9d73;--exclude:#d94c4c;--accent:#315efb}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}header{padding:24px 30px;background:linear-gradient(135deg,#13233b,#284b7a);color:#fff}header h1{margin:0 0 6px;font-size:25px}header p{margin:0;color:#cbd8e9}main{padding:20px;max-width:1700px;margin:auto}.cards{display:grid;grid-template-columns:repeat(5,minmax(140px,1fr));gap:12px;margin-bottom:14px}.card,.panel{background:var(--card);border:1px solid var(--line);border-radius:12px;box-shadow:0 2px 8px #17304b0d}.card{padding:14px}.card b{display:block;font-size:24px}.card span,.note{color:var(--muted)}.layout{display:grid;grid-template-columns:minmax(0,1.25fr) minmax(360px,.75fr);gap:14px}.panel{padding:15px;margin-bottom:14px}h2{font-size:16px;margin:0 0 12px}.toolbar{display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin-bottom:10px}select,input{font:inherit;border:1px solid #c9d4e2;border-radius:7px;padding:7px 9px;background:#fff}#timeline{width:100%;height:275px;border:1px solid var(--line);border-radius:8px;background:#fbfdff}.triad{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:8px}canvas{width:100%;height:230px;border:1px solid var(--line);border-radius:8px;background:#07101f}.caption{text-align:center;color:var(--muted);margin-top:4px;overflow-wrap:anywhere}.metrics{display:grid;grid-template-columns:repeat(2,1fr);gap:7px}.metric{padding:8px;background:#f7f9fc;border-radius:7px}.metric b{display:block}.badge{display:inline-block;padding:3px 8px;border-radius:999px;color:#fff;font-weight:650}.keep{background:var(--keep)}.exclude{background:var(--exclude)}.tablewrap{overflow:auto;max-height:520px}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:7px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}th{position:sticky;top:0;background:#edf3fa}tr.active{background:#edf3ff}.reason{white-space:normal;min-width:280px}.statgrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(245px,1fr));gap:10px}.statplot{border:1px solid var(--line);border-radius:9px;padding:10px;background:#fbfdff}.stathead{display:flex;justify-content:space-between;gap:8px;align-items:flex-start}.stathead b{font-size:13px}.statvalue{color:var(--accent);font-weight:700;white-space:nowrap}.statplot svg{width:100%;height:145px;display:block}.stathelp{color:var(--muted);font-size:11px;min-height:32px;margin:3px 0 0}@media(max-width:1000px){.layout{grid-template-columns:1fr}.cards{grid-template-columns:repeat(2,1fr)}}
</style></head><body><header><h1 id="title"></h1><p>Two-stage serial-section quality decisions · final output: keep / exclude</p></header><main><section class="cards" id="cards"></section><section class="panel"><h2>Descriptive statistics across all slices</h2><p class="note" id="statNote"></p><div class="statgrid" id="statGrid"></div></section><div class="layout"><div><section class="panel"><h2>Complete binary decision timeline</h2><svg id="timeline" viewBox="0 0 1100 275" preserveAspectRatio="none"></svg><p class="note">Stage 1 performs threshold triage. Only the internal uncertainty queue enters stage 2, where multi-domain, anatomical and 3/5/7-window evidence resolves the final action. The public result is always keep or exclude.</p></section><section class="panel"><h2>All slice decisions</h2><div class="toolbar"><select id="filter"><option value="all">All calls</option><option value="exclude">Exclude</option><option value="keep">Keep</option></select><input id="search" placeholder="Search slice"></div><div class="tablewrap"><table><thead><tr><th>#</th><th>slice</th><th>score</th><th>call</th><th>density</th><th>expression</th><th>damage</th><th>continuity</th><th>reason</th></tr></thead><tbody id="rows"></tbody></table></div></section></div><div><section class="panel"><h2>Local anatomical context</h2><div class="toolbar"><label>Slice <select id="sliceSelect"></select></label><select id="colorMode"><option value="depth">captured counts · shared scale</option><option value="plain">uniform</option></select></div><div class="triad"><div><canvas id="prev"></canvas><div class="caption" id="prevCap"></div></div><div><canvas id="mid"></canvas><div class="caption" id="midCap"></div></div><div><canvas id="next"></canvas><div class="caption" id="nextCap"></div></div></div></section><section class="panel"><p class="note">Absolute counts share the window color range when available. Legacy payloads without counts use within-slice relative depth and cannot compare capture levels.</p><h2>Selected evidence</h2><div id="selected"></div></section><section class="panel"><h2>Decision contract</h2><p>Initial <b>keep</b> and safe high-confidence <b>exclude</b> calls may pass directly. Only internally uncertain slices receive the fine screen. They become <b>exclude</b> only when all configured evidence checks pass; otherwise they become <b>keep</b>. No source H5AD is modified.</p><p class="note">The intermediate triage is retained only in the audit CSV. It never creates a third final user action.</p></section></div></div></main><script id="payload" type="application/json">__PAYLOAD__</script><script>
const D=JSON.parse(document.getElementById('payload').textContent),R=D.records,color={keep:'#1f9d73',exclude:'#d94c4c'};let selected=0;document.getElementById('title').textContent=D.title;const count=x=>R.filter(r=>r.recommendation===x).length;document.getElementById('cards').innerHTML=[['total slices',D.n_total],['classified',R.length],['keep',count('keep')],['exclude',count('exclude')],['fine-screen resolved',D.n_review_resolved]].map(([k,v])=>`<div class="card"><b style="color:${color[k]||''}">${v}</b><span>${k}</span></div>`).join('');
function statFmt(x){if(x===null||x===undefined||x==='')return 'NA';let v=Number(x);if(!Number.isFinite(v))return 'NA';let a=Math.abs(v);if(a>0&&(a<.01||a>=1e6))return v.toExponential(2);return new Intl.NumberFormat('en-US',{maximumFractionDigits:a<10?3:a<100?1:0}).format(v)}
function statistics(){let defs=(D.statisticDefinitions||[]).filter(d=>D.statisticSummary&&D.statisticSummary[d.key]),den=Math.max(D.n_total-1,1),w=330,h=145,pl=42,pr=8,pt=13,pb=24,iw=w-pl-pr,ih=h-pt-pb;document.getElementById('statNote').textContent=`${D.measurementContext.count_semantics}. ${D.measurementContext.cross_dataset_note} Observed spot spacing is not nominal platform resolution.`;document.getElementById('statGrid').innerHTML=defs.map((d,k)=>{let s=D.statisticSummary[d.key],vals=R.map((r,i)=>({r,i,v:r[d.key]===null||r[d.key]===undefined?NaN:Number(r[d.key])})).filter(o=>Number.isFinite(o.v)),lo=Math.min(...vals.map(o=>o.v)),hi=Math.max(...vals.map(o=>o.v));if(!(hi>lo)){lo-=Math.abs(lo||1)*.05;hi+=Math.abs(hi||1)*.05}let sy=v=>pt+(hi-v)*ih/(hi-lo),grid=[lo,(lo+hi)/2,hi].map(v=>`<line x1="${pl}" x2="${w-pr}" y1="${sy(v)}" y2="${sy(v)}" stroke="#e3e9f1"/><text x="2" y="${sy(v)+4}" font-size="9" fill="#607086">${statFmt(v)}</text>`).join(''),median=`<line x1="${pl}" x2="${w-pr}" y1="${sy(s.median)}" y2="${sy(s.median)}" stroke="#315efb" stroke-dasharray="4 4"/>`,pts=vals.map(o=>{let x=pl+(Number.isFinite(+o.r.slice_index)?+o.r.slice_index:o.i)*iw/den;return `<circle data-i="${o.i}" cx="${x}" cy="${sy(o.v)}" r="${o.i===selected?5:3.5}" fill="${color[o.r.recommendation]}" stroke="white" stroke-width="1.5"><title>${o.r.slice_id}: ${statFmt(o.v)} ${d.unit}</title></circle>`}).join('');return `<div class="statplot" title="${d.help}"><div class="stathead"><b>${d.label}</b><span class="statvalue">median ${statFmt(s.median)}</span></div><svg class="statSvg" data-k="${k}" viewBox="0 0 ${w} ${h}" preserveAspectRatio="none">${grid}${median}${pts}<text x="${pl}" y="${h-6}" font-size="9" fill="#607086">slice order →</text></svg><p class="stathelp">${d.help}</p></div>`}).join('');document.querySelectorAll('.statSvg circle').forEach(n=>n.onclick=()=>select(+n.dataset.i))}
function timeline(){let svg=document.getElementById('timeline'),w=1100,h=275,p=35,inner=w-2*p,den=Math.max(D.n_total-1,1),pts=R.map(r=>[p+(+r.slice_index)*inner/den,h-p-(+r.quality_anomaly_score)*(h-2*p)]);let grid=[0,.25,.5,.75,1].map(v=>`<line x1="${p}" x2="${w-p}" y1="${h-p-v*(h-2*p)}" y2="${h-p-v*(h-2*p)}" stroke="#dce4ef"/><text x="6" y="${h-p-v*(h-2*p)+4}" fill="#607086" font-size="11">${v.toFixed(2)}</text>`).join('');let dots=pts.map((q,i)=>`<circle data-i="${i}" cx="${q[0]}" cy="${q[1]}" r="${i===selected?7:5}" fill="${color[R[i].recommendation]}" stroke="white" stroke-width="2"><title>${R[i].slice_id}: ${Number(R[i].quality_anomaly_score).toFixed(3)} (${R[i].recommendation})</title></circle>`).join('');svg.innerHTML=grid+dots;svg.querySelectorAll('circle').forEach(n=>n.onclick=()=>select(+n.dataset.i))}
function table(){let f=document.getElementById('filter').value,q=document.getElementById('search').value.toLowerCase(),items=R.map((r,i)=>({r,i})).filter(o=>(f==='all'||o.r.recommendation===f)&&String(o.r.slice_id).toLowerCase().includes(q));let fmt=x=>x!==null&&x!==undefined&&Number.isFinite(+x)?(+x).toFixed(3):'NA';document.getElementById('rows').innerHTML=items.map(o=>`<tr data-i="${o.i}" class="${o.i===selected?'active':''}"><td>${o.r.slice_index}</td><td>${o.r.slice_id}</td><td>${fmt(o.r.quality_anomaly_score)}</td><td><span class="badge ${o.r.recommendation}">${o.r.recommendation}</span></td><td>${fmt(o.r.density_domain_score)}</td><td>${fmt(o.r.expression_domain_score)}</td><td>${fmt(o.r.damage_domain_score)}</td><td>${fmt(o.r.continuity_domain_score)}</td><td class="reason">${o.r.reason}</td></tr>`).join('');document.querySelectorAll('#rows tr').forEach(n=>n.onclick=()=>select(+n.dataset.i))}
let depthScale=[0,1];
function refreshDepthScale(i){let oi=D.order.indexOf(String(R[i].slice_id)),vals=D.order.slice(Math.max(0,oi-1),oi+2).flatMap(id=>D.points[id]?.counts||[]).filter(v=>typeof v==='number'&&Number.isFinite(v)).sort((a,b)=>a-b);depthScale=vals.length?[vals[Math.floor(.04*(vals.length-1))],vals[Math.floor(.96*(vals.length-1))]]:[0,1];}
function sharedDepth(v){return typeof v==='number'&&Number.isFinite(v)?Math.max(0,Math.min(1,(v-depthScale[0])/Math.max(depthScale[1]-depthScale[0],1e-9))):.5;}
function draw(canvas,id){let ctx=canvas.getContext('2d'),d=D.points[id];canvas.width=canvas.clientWidth*devicePixelRatio;canvas.height=canvas.clientHeight*devicePixelRatio;ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);let w=canvas.clientWidth,h=canvas.clientHeight;ctx.fillStyle='#07101f';ctx.fillRect(0,0,w,h);if(!d||!d.x.length)return;let minx=Math.min(...d.x),maxx=Math.max(...d.x),miny=Math.min(...d.y),maxy=Math.max(...d.y),s=Math.min((w-18)/Math.max(maxx-minx,1e-9),(h-18)/Math.max(maxy-miny,1e-9));for(let i=0;i<d.x.length;i++){let x=9+(d.x[i]-minx)*s,y=h-9-(d.y[i]-miny)*s,v=d.counts?sharedDepth(d.counts[i]):d.depth[i];ctx.fillStyle=document.getElementById('colorMode').value==='plain'?'#6db7ff':`rgb(${Math.round(45+210*v)},${Math.round(95+90*(1-v))},${Math.round(230-170*v)})`;ctx.globalAlpha=.78;ctx.beginPath();ctx.arc(x,y,1.6,0,Math.PI*2);ctx.fill()}ctx.globalAlpha=1}
function select(i){if(!R.length)return;selected=Math.max(0,Math.min(R.length-1,i));refreshDepthScale(selected);document.getElementById('sliceSelect').value=selected;let r=R[selected],oi=D.order.indexOf(String(r.slice_id)),ids=[D.order[oi-1],D.order[oi],D.order[oi+1]],can=['prev','mid','next'];can.forEach((c,j)=>{draw(document.getElementById(c),ids[j]);document.getElementById(c+'Cap').textContent=ids[j]??'—'});document.getElementById('selected').innerHTML=`<p><span class="badge ${r.recommendation}">${r.recommendation}</span> <b>${r.slice_id}</b> · score ${statFmt(r.quality_anomaly_score)}</p><div class="metrics"><div class="metric"><b>${statFmt(r.n_locations)}</b>locations / spots</div><div class="metric"><b>${statFmt(r.median_n_genes)}</b>median genes / spot</div><div class="metric"><b>${statFmt(r.median_total_counts)}</b>median captured counts / spot</div><div class="metric"><b>${statFmt(r.detected_genes)}</b>total detected genes</div><div class="metric"><b>${statFmt(r.cell_density)}</b>density</div><div class="metric"><b>${statFmt(r.median_nn_distance)}</b>observed spot spacing</div><div class="metric"><b>${statFmt(r.hole_fraction)}</b>hole fraction</div><div class="metric"><b>${statFmt(r.regional_low_depth_cluster_fraction)}</b>regional low depth</div></div><p>${r.binary_reason}</p>`;timeline();table();statistics()}
document.getElementById('sliceSelect').innerHTML=R.map((r,i)=>`<option value="${i}">${r.slice_id}</option>`).join('');document.getElementById('sliceSelect').onchange=e=>select(+e.target.value);document.getElementById('colorMode').onchange=()=>select(selected);document.getElementById('filter').onchange=table;document.getElementById('search').oninput=table;timeline();table();select(0);
</script></body></html>"""
    destination.write_text(
        template.replace("__TITLE__", html.escape(resolved_title)).replace(
            "__PAYLOAD__", payload_json
        ),
        encoding="utf-8",
    )
    from slice_quality_preregistration import attach_preregistration_html

    attach_preregistration_html(input_path, destination)
    return str(destination)


def _render_report(
    result: SliceSeriesResult,
    *,
    title: str,
    ground_truth: Optional[Mapping[str, bool]],
) -> str:
    records = [_jsonable(record) for record in result.metrics.to_dict("records")]
    payload = {
        "title": title,
        "records": records,
        "points": _jsonable(result.point_samples),
        "provenance": _jsonable(result.provenance),
        "statisticDefinitions": _jsonable(_REPORT_STATISTICS),
        "statisticSummary": _jsonable(_summarize_report_statistics(result.metrics)),
        "measurementContext": _jsonable(_report_measurement_context(result.provenance)),
        "groundTruth": _jsonable(ground_truth or {}),
        "evaluation": (
            _jsonable(evaluate_slice_calls(result.metrics, ground_truth))
            if ground_truth is not None
            else None
        ),
    }
    payload_json = json.dumps(payload, ensure_ascii=False).replace("</", "<\\/")
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<style>
:root{{--bg:#f4f7fb;--card:#fff;--ink:#122033;--muted:#607086;--line:#dce4ef;--keep:#1f9d73;--review:#e6a029;--exclude:#d94c4c;--accent:#315efb}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}}
header{{padding:24px 30px;background:linear-gradient(135deg,#13233b,#284b7a);color:white}}header h1{{margin:0 0 6px;font-size:25px}}header p{{margin:0;color:#cbd8e9}}
main{{padding:20px;max-width:1700px;margin:auto}}.cards{{display:grid;grid-template-columns:repeat(5,minmax(140px,1fr));gap:12px;margin-bottom:14px}}
.card,.panel{{background:var(--card);border:1px solid var(--line);border-radius:12px;box-shadow:0 2px 8px #17304b0d}}.card{{padding:14px}}.card b{{display:block;font-size:24px}}.card span{{color:var(--muted)}}
.layout{{display:grid;grid-template-columns:minmax(0,1.25fr) minmax(360px,.75fr);gap:14px}}.panel{{padding:15px;margin-bottom:14px}}h2{{font-size:16px;margin:0 0 12px}}
.toolbar{{display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin-bottom:10px}}select,input,button{{font:inherit;border:1px solid #c9d4e2;border-radius:7px;padding:7px 9px;background:white}}button{{cursor:pointer;color:#174b94}}label{{color:var(--muted)}}
#timeline{{width:100%;height:275px;border:1px solid var(--line);border-radius:8px;background:#fbfdff}}.triad{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:8px}}canvas{{width:100%;height:230px;border:1px solid var(--line);border-radius:8px;background:#07101f}}.caption{{text-align:center;color:var(--muted);margin-top:4px;overflow-wrap:anywhere}}
.metrics{{display:grid;grid-template-columns:repeat(2,1fr);gap:7px}}.metric{{padding:8px;background:#f7f9fc;border-radius:7px}}.metric b{{display:block}}.badge{{display:inline-block;padding:3px 8px;border-radius:999px;color:#fff;font-weight:650}}.keep{{background:var(--keep)}}.review{{background:var(--review)}}.exclude{{background:var(--exclude)}}
.statgrid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(245px,1fr));gap:10px}}.statplot{{border:1px solid var(--line);border-radius:9px;padding:10px;background:#fbfdff}}.stathead{{display:flex;justify-content:space-between;gap:8px;align-items:flex-start}}.stathead b{{font-size:13px}}.statvalue{{color:var(--accent);font-weight:700;white-space:nowrap}}.statplot svg{{width:100%;height:145px;display:block}}.stathelp{{color:var(--muted);font-size:11px;min-height:32px;margin:3px 0 0}}
.tablewrap{{overflow:auto;max-height:520px}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{padding:7px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}}th{{position:sticky;top:0;background:#edf3fa;cursor:pointer}}tr.active{{background:#edf3ff}}code{{background:#eef2f7;padding:2px 4px;border-radius:4px}}.note{{color:var(--muted);font-size:12px}}.reason{{white-space:normal;min-width:280px}}
@media(max-width:1000px){{.layout{{grid-template-columns:1fr}}.cards{{grid-template-columns:repeat(2,1fr)}}}}
</style></head><body>
<header><h1 id="title"></h1><p>Non-destructive triage before alignment · keep / review / exclude candidates · default centered local window</p></header>
<main><section class="cards" id="cards"></section>
<section class="panel"><h2>Descriptive statistics across slices</h2><p class="note" id="statNote"></p><div class="statgrid" id="statGrid"></div></section>
<div class="layout"><div>
<section class="panel"><h2>Slice anomaly timeline</h2><div class="toolbar"><label>Review ≥ <input id="reviewCut" type="range" min="0" max="1" step="0.01"></label><span id="reviewVal"></span><label>Exclude ≥ <input id="excludeCut" type="range" min="0" max="1" step="0.01"></label><span id="excludeVal"></span><button id="download">Export adjusted calls CSV</button></div><svg id="timeline" viewBox="0 0 1100 275" preserveAspectRatio="none"></svg><p class="note">Threshold sliders are for visual sensitivity analysis only. The saved metrics retain the reproducible detector calls.</p></section>
<section class="panel"><h2>Slice table</h2><div class="toolbar"><select id="filter"><option value="all">All calls</option><option value="exclude">Exclude</option><option value="review">Review</option><option value="keep">Keep</option></select><input id="search" placeholder="Search slice"></div><div class="tablewrap"><table><thead><tr><th data-key="slice_index">#</th><th data-key="slice_id">slice</th><th data-key="quality_anomaly_score">score</th><th data-key="recommendation">call</th><th data-key="density_domain_score">density</th><th data-key="expression_domain_score">expression</th><th data-key="damage_domain_score">damage</th><th data-key="continuity_domain_score">continuity</th><th>reason</th></tr></thead><tbody id="rows"></tbody></table></div></section>
</div><div>
<section class="panel"><h2>Triad / local-window review</h2><div class="toolbar"><label>Focal slice <select id="sliceSelect"></select></label><select id="colorMode"><option value="depth">captured counts · shared scale</option><option value="plain">uniform</option></select></div><div class="triad"><div><canvas id="prev"></canvas><div class="caption" id="prevCap"></div></div><div><canvas id="mid"></canvas><div class="caption" id="midCap"></div></div><div><canvas id="next"></canvas><div class="caption" id="nextCap"></div></div></div></section>
<section class="panel"><p class="note">Absolute counts share the window color range when available. Legacy payloads without counts use within-slice relative depth and cannot compare capture levels.</p><h2>Selected slice evidence</h2><div id="selected"></div></section>
<section class="panel"><h2>Interpretation contract</h2><p><b>Exclude</b> requires strong expression loss or corroborated evidence from multiple domains. Geometry-only anomalies, naturally tapering end sections, and one-sided endpoints are capped at manual review.</p><p class="note">Library size and detected genes are count-matrix proxies. Raw reads, duplication and saturation cannot be recovered from H5AD unless supplied separately. No source file is modified and no slice is automatically removed.</p><details><summary>Provenance</summary><pre id="prov" style="white-space:pre-wrap;font-size:11px"></pre></details></section>
</div></div></main>
<script id="payload" type="application/json">{payload_json}</script>
<script>
const D=JSON.parse(document.getElementById('payload').textContent), R=D.records; document.getElementById('title').textContent=D.title;
const reviewInput=document.getElementById('reviewCut'), excludeInput=document.getElementById('excludeCut');
reviewInput.value=D.provenance.config.review_threshold??0.38; excludeInput.value=D.provenance.config.exclude_threshold??0.64;
const initialReview=+reviewInput.value, initialExclude=+excludeInput.value;
const color={{keep:'#1f9d73',review:'#e6a029',exclude:'#d94c4c'}}; let selected=0, sortKey='slice_index', sortAsc=true;
function adjustedCall(r){{let a=+r.quality_anomaly_score, e=+excludeInput.value, v=+reviewInput.value;if(Math.abs(e-initialExclude)<1e-9&&Math.abs(v-initialReview)<1e-9)return r.recommendation;if(r.partial_structure_protection||r.window_context!=='two_sided')return a>=v?'review':'keep';return a>=e?'exclude':a>=v?'review':'keep'}}
function cards(){{let calls=R.map(adjustedCall), count=x=>calls.filter(y=>y===x).length;document.getElementById('cards').innerHTML=`<div class="card"><b>${{R.length}}</b><span>slices</span></div><div class="card"><b style="color:var(--keep)">${{count('keep')}}</b><span>keep</span></div><div class="card"><b style="color:var(--review)">${{count('review')}}</b><span>review</span></div><div class="card"><b style="color:var(--exclude)">${{count('exclude')}}</b><span>exclude candidates</span></div><div class="card"><b>${{D.provenance.selected_window}}</b><span>selected window</span></div>`;}}
function statFmt(x){{if(x===null||x===undefined||x==='')return 'NA';let v=Number(x);if(!Number.isFinite(v))return 'NA';let a=Math.abs(v);if(a>0&&(a<.01||a>=1e6))return v.toExponential(2);return new Intl.NumberFormat('en-US',{{maximumFractionDigits:a<10?3:a<100?1:0}}).format(v)}}
function statistics(){{let defs=(D.statisticDefinitions||[]).filter(d=>D.statisticSummary&&D.statisticSummary[d.key]),den=Math.max(R.length-1,1),w=330,h=145,pl=42,pr=8,pt=13,pb=24,iw=w-pl-pr,ih=h-pt-pb;document.getElementById('statNote').textContent=`${{D.measurementContext.count_semantics}}. ${{D.measurementContext.cross_dataset_note}} Observed spot spacing is not nominal platform resolution.`;document.getElementById('statGrid').innerHTML=defs.map((d,k)=>{{let s=D.statisticSummary[d.key],vals=R.map((r,i)=>({{r,i,v:r[d.key]===null||r[d.key]===undefined?NaN:Number(r[d.key])}})).filter(o=>Number.isFinite(o.v)),lo=Math.min(...vals.map(o=>o.v)),hi=Math.max(...vals.map(o=>o.v));if(!(hi>lo)){{lo-=Math.abs(lo||1)*.05;hi+=Math.abs(hi||1)*.05}}let sy=v=>pt+(hi-v)*ih/(hi-lo),grid=[lo,(lo+hi)/2,hi].map(v=>`<line x1="${{pl}}" x2="${{w-pr}}" y1="${{sy(v)}}" y2="${{sy(v)}}" stroke="#e3e9f1"/><text x="2" y="${{sy(v)+4}}" font-size="9" fill="#607086">${{statFmt(v)}}</text>`).join(''),median=`<line x1="${{pl}}" x2="${{w-pr}}" y1="${{sy(s.median)}}" y2="${{sy(s.median)}}" stroke="#315efb" stroke-dasharray="4 4"/>`,pts=vals.map(o=>{{let x=pl+o.i*iw/den;return `<circle data-i="${{o.i}}" cx="${{x}}" cy="${{sy(o.v)}}" r="${{o.i===selected?5:3.5}}" fill="${{color[adjustedCall(o.r)]}}" stroke="white" stroke-width="1.5"><title>${{o.r.slice_id}}: ${{statFmt(o.v)}} ${{d.unit}}</title></circle>`}}).join('');return `<div class="statplot" title="${{d.help}}"><div class="stathead"><b>${{d.label}}</b><span class="statvalue">median ${{statFmt(s.median)}}</span></div><svg class="statSvg" data-k="${{k}}" viewBox="0 0 ${{w}} ${{h}}" preserveAspectRatio="none">${{grid}}${{median}}${{pts}}<text x="${{pl}}" y="${{h-6}}" font-size="9" fill="#607086">slice order →</text></svg><p class="stathelp">${{d.help}}</p></div>`}}).join('');document.querySelectorAll('.statSvg circle').forEach(n=>n.onclick=()=>select(+n.dataset.i))}}
function timeline(){{let svg=document.getElementById('timeline'),w=1100,h=275,p=35,inner=w-2*p;let pts=R.map((r,i)=>[p+i*inner/Math.max(R.length-1,1),h-p-(+r.quality_anomaly_score)*(h-2*p)]);let grid=[0,.25,.5,.75,1].map(v=>`<line x1="${{p}}" x2="${{w-p}}" y1="${{h-p-v*(h-2*p)}}" y2="${{h-p-v*(h-2*p)}}" stroke="#dce4ef"/><text x="6" y="${{h-p-v*(h-2*p)+4}}" fill="#607086" font-size="11">${{v.toFixed(2)}}</text>`).join('');let line=`<polyline fill="none" stroke="#315efb" stroke-width="2" points="${{pts.map(p=>p.join(',')).join(' ')}}"/>`;let dots=pts.map((p,i)=>`<circle data-i="${{i}}" cx="${{p[0]}}" cy="${{p[1]}}" r="${{i===selected?7:5}}" fill="${{color[adjustedCall(R[i])]}}" stroke="white" stroke-width="2"><title>${{R[i].slice_id}}: ${{(+R[i].quality_anomaly_score).toFixed(3)}}</title></circle>`).join('');let cuts=[['reviewCut','#e6a029'],['excludeCut','#d94c4c']].map(([id,c])=>{{let v=+document.getElementById(id).value,y=h-p-v*(h-2*p);return `<line x1="${{p}}" x2="${{w-p}}" y1="${{y}}" y2="${{y}}" stroke="${{c}}" stroke-dasharray="6 5"/>`}}).join('');svg.innerHTML=grid+cuts+line+dots;svg.querySelectorAll('circle').forEach(n=>n.onclick=()=>select(+n.dataset.i));}}
function table(){{let f=document.getElementById('filter').value,q=document.getElementById('search').value.toLowerCase();let items=R.map((r,i)=>({{r,i}})).filter(o=>(f==='all'||adjustedCall(o.r)===f)&&String(o.r.slice_id).toLowerCase().includes(q));items.sort((a,b)=>{{let x=a.r[sortKey],y=b.r[sortKey];return (typeof x==='number'?x-y:String(x).localeCompare(String(y)))*(sortAsc?1:-1)}});document.getElementById('rows').innerHTML=items.map(o=>{{let r=o.r,c=adjustedCall(r),fmt=x=>x!==null&&x!==undefined&&Number.isFinite(+x)?(+x).toFixed(3):'NA';return `<tr data-i="${{o.i}}" class="${{o.i===selected?'active':''}}"><td>${{r.slice_index}}</td><td>${{r.slice_id}}</td><td>${{fmt(r.quality_anomaly_score)}}</td><td><span class="badge ${{c}}">${{c}}</span></td><td>${{fmt(r.density_domain_score)}}</td><td>${{fmt(r.expression_domain_score)}}</td><td>${{fmt(r.damage_domain_score)}}</td><td>${{fmt(r.continuity_domain_score)}}</td><td class="reason">${{r.reason}}</td></tr>`}}).join('');document.querySelectorAll('#rows tr').forEach(n=>n.onclick=()=>select(+n.dataset.i));}}
let depthScale=[0,1];
function refreshDepthScale(i){{let vals=R.slice(Math.max(0,i-1),i+2).flatMap(r=>D.points[r.slice_id]?.counts||[]).filter(v=>typeof v==='number'&&Number.isFinite(v)).sort((a,b)=>a-b);depthScale=vals.length?[vals[Math.floor(.04*(vals.length-1))],vals[Math.floor(.96*(vals.length-1))]]:[0,1];}}
function sharedDepth(v){{return typeof v==='number'&&Number.isFinite(v)?Math.max(0,Math.min(1,(v-depthScale[0])/Math.max(depthScale[1]-depthScale[0],1e-9))):.5;}}
function draw(canvas,id){{let ctx=canvas.getContext('2d'),d=D.points[id];canvas.width=canvas.clientWidth*devicePixelRatio;canvas.height=canvas.clientHeight*devicePixelRatio;ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);let w=canvas.clientWidth,h=canvas.clientHeight;ctx.fillStyle='#07101f';ctx.fillRect(0,0,w,h);if(!d||!d.x.length)return;let minx=Math.min(...d.x),maxx=Math.max(...d.x),miny=Math.min(...d.y),maxy=Math.max(...d.y),s=Math.min((w-18)/Math.max(maxx-minx,1e-9),(h-18)/Math.max(maxy-miny,1e-9));for(let i=0;i<d.x.length;i++){{let x=9+(d.x[i]-minx)*s,y=h-9-(d.y[i]-miny)*s,v=d.counts?sharedDepth(d.counts[i]):d.depth[i];ctx.fillStyle=document.getElementById('colorMode').value==='plain'?'#6db7ff':`rgb(${{Math.round(45+210*v)}},${{Math.round(95+90*(1-v))}},${{Math.round(230-170*v)}})`;ctx.globalAlpha=.78;ctx.beginPath();ctx.arc(x,y,1.6,0,Math.PI*2);ctx.fill();}}ctx.globalAlpha=1;}}
function select(i){{refreshDepthScale(i);selected=Math.max(0,Math.min(R.length-1,i));document.getElementById('sliceSelect').value=selected;let ids=[R[selected-1]?.slice_id,R[selected].slice_id,R[selected+1]?.slice_id],can=['prev','mid','next'];can.forEach((c,j)=>{{draw(document.getElementById(c),ids[j]);document.getElementById(c+'Cap').textContent=ids[j]??'—'}});let r=R[selected];document.getElementById('selected').innerHTML=`<p><span class="badge ${{adjustedCall(r)}}">${{adjustedCall(r)}}</span> <b>${{r.slice_id}}</b> · score ${{statFmt(r.quality_anomaly_score)}}</p><div class="metrics"><div class="metric"><b>${{statFmt(r.n_locations)}}</b>locations / spots</div><div class="metric"><b>${{statFmt(r.median_n_genes)}}</b>median genes / spot</div><div class="metric"><b>${{statFmt(r.median_total_counts)}}</b>median captured counts / spot</div><div class="metric"><b>${{statFmt(r.detected_genes)}}</b>total detected genes</div><div class="metric"><b>${{statFmt(r.cell_density)}}</b>density</div><div class="metric"><b>${{statFmt(r.median_nn_distance)}}</b>observed spot spacing</div><div class="metric"><b>${{statFmt(r.hole_fraction)}}</b>hole fraction</div><div class="metric"><b>${{statFmt(r.regional_low_depth_cluster_fraction)}}</b>regional low depth</div></div><p>${{r.reason}}</p>`;timeline();table();statistics();}}
function refresh(){{document.getElementById('reviewVal').textContent=(+reviewInput.value).toFixed(2);document.getElementById('excludeVal').textContent=(+excludeInput.value).toFixed(2);if(+reviewInput.value>+excludeInput.value)reviewInput.value=excludeInput.value;cards();timeline();table();select(selected)}}
document.getElementById('sliceSelect').innerHTML=R.map((r,i)=>`<option value="${{i}}">${{r.slice_id}}</option>`).join('');document.getElementById('sliceSelect').onchange=e=>select(+e.target.value);document.getElementById('colorMode').onchange=()=>select(selected);reviewInput.oninput=refresh;excludeInput.oninput=refresh;document.getElementById('filter').onchange=table;document.getElementById('search').oninput=table;
document.querySelectorAll('th[data-key]').forEach(n=>n.onclick=()=>{{if(sortKey===n.dataset.key)sortAsc=!sortAsc;else{{sortKey=n.dataset.key;sortAsc=true}}table()}});
document.getElementById('download').onclick=()=>{{let rows=['slice_id,score,adjusted_call'].concat(R.map(r=>`${{JSON.stringify(r.slice_id)}},${{r.quality_anomaly_score}},${{adjustedCall(r)}}`));let a=document.createElement('a');a.href=URL.createObjectURL(new Blob([rows.join('\\n')],{{type:'text/csv'}}));a.download='adjusted_slice_calls.csv';a.click();URL.revokeObjectURL(a.href)}};document.getElementById('prov').textContent=JSON.stringify(D.provenance,null,2);refresh();
</script></body></html>"""


def write_slice_quality_collection_report(
    run_dirs: Sequence[Union[str, Path]],
    output_html: Union[str, Path],
    *,
    title: str = "Spateo referee",
    binary_only: bool = False,
    catalog_records: Optional[Sequence[Mapping[str, Any]]] = None,
) -> dict[str, Any]:
    """Bundle multiple completed slice-QC runs into one interactive HTML index.

    Each run directory must contain ``slice_quality_metrics.csv`` and
    ``slice_quality_manifest.json``.  Synthetic benchmark directories may also
    contain ``synthetic_benchmark_evaluation.json`` and
    ``paired_detection_comparison.csv``; their paired expression-detection
    evidence is surfaced separately from real-data calls.
    """
    output_path = Path(output_html).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    datasets: list[dict[str, Any]] = []
    candidates: list[dict[str, Any]] = []
    benchmarks: list[dict[str, Any]] = []
    binary_policies: dict[str, dict[str, Any]] = {}
    for value in run_dirs:
        run_dir = Path(value).expanduser().resolve()
        metrics_path = run_dir / "slice_quality_metrics.csv"
        manifest_path = run_dir / "slice_quality_manifest.json"
        report_path = run_dir / "slice_quality_report.html"
        if (
            not metrics_path.exists()
            or not manifest_path.exists()
            or not report_path.exists()
        ):
            raise FileNotFoundError(f"Incomplete slice-QC run directory: {run_dir}")
        metrics = pd.read_csv(metrics_path, dtype={"slice_id": str})
        descriptive_metrics = metrics.copy()
        published_count = int(len(metrics))
        withheld_count = 0
        if binary_only:
            binary_path = run_dir / "slice_quality_binary_audit.csv"
            if not binary_path.exists():
                raise FileNotFoundError(
                    f"Binary collection mode requires a published audit table: {binary_path}"
                )
            binary = pd.read_csv(binary_path, dtype={"slice_id": str})
            final_values = binary["final_call"].astype("string")
            invalid = final_values.isna() | ~final_values.isin(["keep", "exclude"])
            if invalid.any():
                raise ValueError(
                    f"Complete binary collection requires a keep/exclude final_call for every slice: {run_dir}"
                )
            published_count = int(len(binary))
            withheld_count = 0
            metrics = binary.copy()
            metrics["recommendation"] = final_values.astype(str)
            application_path = run_dir / "binary_policy_application.json"
            if application_path.exists():
                application = json.loads(application_path.read_text(encoding="utf-8"))
                policy_record = dict(application.get("policy", {}))
                policy_id = str(policy_record.get("calibration_id", "unknown"))
                binary_policies[policy_id] = policy_record
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        name = run_dir.name
        config = manifest.get("config", {})
        source_records = manifest.get("sources", [])
        source_path = str(source_records[0].get("path", "")) if source_records else ""
        raw_group = str(
            manifest.get("collection_group")
            or config.get("collection_group")
            or run_dir.parent.name
            or "Other"
        )
        group_aliases = {
            "cs8": "Human embryo",
            "human_embryo": "Human embryo",
            "drosophila": "Drosophila",
            "planarian": "Planarian",
            "mouse": "Mouse organs",
            "mouse_brain": "Mouse organs",
            "zebrafish": "Zebrafish",
            "artista": "Axolotl / Artista",
            "cerebellum": "Cross-species cerebellum",
            "expression_benchmarks": "Validation benchmarks",
            "artificial_unaligned_mild_20260902": "Validation benchmarks",
        }
        normalized_group = re.sub(r"[^a-z0-9]+", "_", raw_group.lower()).strip("_")
        group = group_aliases.get(
            normalized_group, raw_group.replace("_", " ").strip().title()
        )
        if group == "Cross-species cerebellum":
            match = re.match(r"(Marmoset|Macaque|Mouse)", name, flags=re.IGNORECASE)
            subgroup = (
                f"{match.group(1).title()} sections" if match else "Other sections"
            )
        elif group == "Drosophila":
            subgroup = "Embryo" if re.match(r"(?:\d|E\d)", name) else "Larva / pupa"
        else:
            subgroup = str(
                manifest.get("collection_subgroup")
                or config.get("collection_subgroup")
                or ""
            )
        try:
            report_href = report_path.relative_to(output_path.parent).as_posix()
        except ValueError:
            report_href = Path(
                os.path.relpath(report_path, output_path.parent)
            ).as_posix()
        counts = metrics["recommendation"].astype(str).value_counts().to_dict()
        selected_window = manifest.get("selected_window")
        if selected_window is None:
            window_values = descriptive_metrics.get(
                "window_size", pd.Series(dtype=float)
            )
            selected_window = window_values.iloc[0] if len(window_values) else 3
        score_values = pd.to_numeric(
            metrics["quality_anomaly_score"], errors="coerce"
        ).dropna()
        record = {
            "name": name,
            "group": group,
            "subgroup": subgroup,
            "title": str(config.get("report_title", name)),
            "source": source_path,
            "n_slices": int(published_count + withheld_count),
            "published": published_count,
            "withheld": withheld_count,
            "coverage": float(
                published_count / max(published_count + withheld_count, 1)
            ),
            "selected_window": int(selected_window),
            "keep": int(counts.get("keep", 0)),
            "review": int(counts.get("review", 0)),
            "exclude": int(counts.get("exclude", 0)),
            "max_score": float(score_values.max()) if len(score_values) else None,
            "report_href": report_href,
            "statistics": _jsonable(_summarize_report_statistics(descriptive_metrics)),
            "measurement_context": _jsonable(_report_measurement_context(manifest)),
            "timeline": _jsonable(
                metrics[
                    (["slice_index"] if "slice_index" in metrics else [])
                    + ["slice_id", "quality_anomaly_score", "recommendation"]
                ].to_dict("records")
            ),
        }
        datasets.append(record)
        for row in metrics.loc[
            metrics["recommendation"].astype(str).ne("keep")
        ].to_dict("records"):
            candidates.append(
                {
                    "dataset": name,
                    "group": group,
                    "subgroup": subgroup,
                    "slice_id": str(row["slice_id"]),
                    "score": _jsonable(row.get("quality_anomaly_score")),
                    "recommendation": str(row.get("recommendation", "review")),
                    "expression": _jsonable(row.get("expression_domain_score")),
                    "density": _jsonable(row.get("density_domain_score")),
                    "partial_protection": bool(
                        row.get("partial_structure_protection", False)
                    ),
                    "reason": str(row.get("reason", "")),
                    "report_href": report_href,
                }
            )

        evaluation_path = run_dir / "synthetic_benchmark_evaluation.json"
        paired_path = run_dir / "paired_detection_comparison.csv"
        if evaluation_path.exists() and paired_path.exists():
            evaluation = json.loads(evaluation_path.read_text(encoding="utf-8"))
            paired_summary = evaluation.get("paired_evaluation", {})
            paired = pd.read_csv(paired_path)
            injected_rows = paired.loc[paired["injected"].astype(bool)].to_dict(
                "records"
            )
            benchmarks.append(
                {
                    "dataset": name,
                    "group": group,
                    "report_href": report_href,
                    "selected_window": evaluation.get(
                        "selected_window_locked_from_baseline"
                    ),
                    "summary": _jsonable(paired_summary),
                    "injected": _jsonable(injected_rows),
                }
            )

    if not datasets:
        raise ValueError("At least one completed run directory is required")
    datasets.sort(
        key=lambda item: (
            _natural_key(item["group"]),
            _natural_key(item["subgroup"]),
            _natural_key(item["name"]),
        )
    )
    candidates.sort(
        key=lambda item: (
            _natural_key(item["group"]),
            _natural_key(item["dataset"]),
            -float(item["score"] or 0),
        )
    )
    catalog = [_jsonable(dict(item)) for item in (catalog_records or [])]
    catalog.sort(
        key=lambda item: (
            _natural_key(str(item.get("group", "Other"))),
            _natural_key(str(item.get("cohort", ""))),
            _natural_key(str(item.get("filename", ""))),
        )
    )
    payload = {
        "title": title,
        "binary_only": bool(binary_only),
        "statistic_definitions": _jsonable(_REPORT_STATISTICS),
        "datasets": datasets,
        "candidates": candidates,
        "benchmarks": benchmarks,
        "binary_policies": list(binary_policies.values()),
        "catalog": catalog,
    }
    payload_json = json.dumps(_jsonable(payload), ensure_ascii=False).replace(
        "</", "<\\/"
    )
    template = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title><style>
:root{--bg:#f3f6fa;--card:#fff;--ink:#132238;--muted:#627287;--line:#dce4ee;--keep:#178963;--review:#df9820;--exclude:#d54848;--accent:#315efb}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
header{padding:27px 32px;background:linear-gradient(125deg,#12243d,#285786);color:#fff}header h1{margin:0 0 6px;font-size:27px}header p{margin:0;color:#cfdbeb}
main{max-width:1680px;margin:auto;padding:20px}.cards{display:grid;grid-template-columns:repeat(5,minmax(130px,1fr));gap:12px}.card,.panel,.dataset{background:var(--card);border:1px solid var(--line);border-radius:12px;box-shadow:0 2px 8px #17304b0d}.card{padding:14px}.card b{font-size:24px;display:block}.card span,.muted{color:var(--muted)}
.panel{padding:16px;margin-top:14px}h2{font-size:17px;margin:0 0 12px}.datasetgroup{margin:16px 0 6px}.datasetgroup:first-child{margin-top:0}.datasetgroup h3{font-size:16px;margin:0 0 9px;padding-bottom:7px;border-bottom:2px solid #dce6f3}.subgroup{color:var(--muted);font-size:12px;margin:8px 0 6px}.datasetgrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(265px,1fr));gap:10px}.dataset{padding:13px;cursor:pointer}.dataset:hover{border-color:#9bb4dc}.dataset h4{font-size:14px;margin:0 0 7px}.counts{display:flex;gap:7px;flex-wrap:wrap}.pill,.badge,.status{border-radius:999px;padding:3px 8px;font-weight:650}.pill{background:#edf2f8}.badge{color:#fff}.keep{background:var(--keep)}.review{background:var(--review)}.exclude{background:var(--exclude)}.cataloggroup{margin-top:13px}.cataloggroup h3{margin:0 0 7px}.cataloggroup details{margin:7px 0;border:1px solid var(--line);border-radius:8px;padding:8px 10px;background:#fbfdff}.catalogfiles{columns:3;column-gap:24px;margin:8px 0 0;padding-left:20px;font-size:12px}.catalogfiles li{break-inside:avoid;margin:3px 0}.serial_qc,.serial_qc_grouped{background:#dff5ec;color:#116b4e}.not_applicable_single_slice{background:#eef2f7;color:#536273}.scan_failed,.missing_or_incomplete{background:#fde5e5;color:#a62f2f}.downloaded_unclassified{background:#fff0cf;color:#895d00}
.toolbar{display:flex;gap:9px;flex-wrap:wrap;margin-bottom:10px}select,input{font:inherit;padding:7px 9px;border:1px solid #c9d4e2;border-radius:7px;background:#fff}svg{width:100%;height:260px;border:1px solid var(--line);border-radius:8px;background:#fbfdff}.tablewrap{overflow:auto;max-height:570px}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:7px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}th{position:sticky;top:0;background:#edf3fa}.reason{white-space:normal;min-width:300px}a{color:#245db8;text-decoration:none}a:hover{text-decoration:underline}.bench{border-left:4px solid #6d5ee7;padding:12px;background:#f8f7ff;border-radius:8px;margin:9px 0}.delta{font-variant-numeric:tabular-nums}.ok{color:var(--keep);font-weight:700}.comparegrid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}.comparepanel{border:1px solid var(--line);border-radius:10px;padding:12px;background:#fbfdff;min-width:0}.comparehead{display:flex;justify-content:space-between;gap:8px;margin-bottom:8px}.comparehead b{font-size:13px}.comparebody{max-height:390px;overflow:auto;padding-right:4px}.comparerow{display:grid;grid-template-columns:minmax(110px,1.15fr) minmax(110px,2fr) 86px;gap:8px;align-items:center;min-height:27px;border-top:1px solid #edf1f6;cursor:pointer}.comparerow:hover{background:#f1f6ff}.comparename{overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-size:11px}.comparetrack{height:12px;background:#e7edf5;border-radius:999px;position:relative}.compareiqr{height:6px;top:3px;position:absolute;background:#94b1df;border-radius:999px}.comparedot{width:9px;height:9px;top:1.5px;position:absolute;background:var(--accent);border:2px solid #fff;border-radius:50%;box-shadow:0 0 0 1px #315efb}.comparevalue{text-align:right;font-variant-numeric:tabular-nums;font-size:11px}.scale{font-size:10px;color:var(--muted);font-weight:500}
@media(max-width:800px){.cards{grid-template-columns:repeat(2,1fr)}.comparegrid{grid-template-columns:1fr}.comparerow{grid-template-columns:minmax(100px,1fr) minmax(90px,1.4fr) 78px}}
</style></head><body><header><h1 id="title"></h1><p>Cross-species serial-section quality benchmark · datasets grouped by biological system and data type</p></header><main>
<section class="cards" id="cards"></section>
<section class="panel" id="policyPanel"><h2>Complete binary decision policy</h2><div id="policySummary"></div></section>
<section class="panel"><h2>Datasets by type</h2><div id="datasetGroups"></div></section>
<section class="panel" id="statisticsPanel"><h2>Dataset-level descriptive statistics</h2><p class="muted">Each row shows the median across slices; the horizontal band is the slice interquartile range. Panels use independent scales. Captured counts are an H5AD count/UMI proxy rather than raw read depth; coordinate-dependent metrics should be compared only when platform and coordinate units match.</p><div class="toolbar"><label>Biological type <select id="statisticsGroup"><option value="all">All data types</option></select></label></div><div class="comparegrid" id="statisticsComparison"></div></section>
<section class="panel" id="catalogPanel"><h2>Official file catalogue</h2><p class="muted">Every requested official file is grouped by biological type. Serial QC is performed only when a valid ordered slice series is present.</p><div id="catalog"></div></section>
<section class="panel"><h2>Selected dataset timeline</h2><div class="toolbar"><select id="datasetSelect"></select><a id="openReport" target="_self">Open full interactive report →</a></div><svg id="timeline" viewBox="0 0 1200 260" preserveAspectRatio="none"></svg></section>
<section class="panel" id="benchmarkPanel"><h2>Paired low-expression detection</h2><p class="muted">The same slices are compared before and after count thinning; native baseline candidates are not relabelled as synthetic false positives.</p><div id="benchmarks"></div></section>
<section class="panel"><h2 id="queueTitle">Review queue</h2><div class="toolbar"><select id="callFilter"><option value="all">All non-keep</option><option value="exclude">Exclude candidates</option><option value="review">Review</option></select><select id="groupFilter"><option value="all">All types</option></select><select id="datasetFilter"><option value="all">All datasets</option></select><input id="search" placeholder="Search slice or reason"></div><div class="tablewrap"><table><thead><tr><th>type</th><th>dataset</th><th>slice</th><th>score</th><th>call</th><th>expression</th><th>density</th><th>partial protected</th><th>reason</th></tr></thead><tbody id="rows"></tbody></table></div></section>
</main><script id="payload" type="application/json">__PAYLOAD__</script><script>
const D=JSON.parse(document.getElementById('payload').textContent), colors={keep:'#178963',review:'#df9820',exclude:'#d54848'};
document.getElementById('title').textContent=D.title;
const sums=k=>D.datasets.reduce((a,d)=>a+(d[k]||0),0), cards=D.binary_only?[['datasets',D.datasets.length],['slices',sums('n_slices')],['classified',sums('published')],['keep',sums('keep')],['exclude',sums('exclude')]]:[['datasets',D.datasets.length],['slices',sums('n_slices')],['keep',sums('keep')],['review',sums('review')],['exclude',sums('exclude')]];
document.getElementById('cards').innerHTML=cards.map(([k,v])=>`<div class="card"><b style="color:${colors[k]||''}">${v}</b><span>${k}</span></div>`).join('');
const pp=document.getElementById('policyPanel');if(!D.binary_only||!D.binary_policies.length)pp.style.display='none';else document.getElementById('policySummary').innerHTML=D.binary_policies.map(p=>{const b=p.benchmark_summary||{},c=b.calibration||{},ext=b.adaptive_extension||{},tiers=p.review_exclusion_tiers||[],legacy=Number.isFinite(Number(p.adaptive_exclude_min_score)),resolver=tiers.length>0||legacy,basePassed=(b.base_policy_validation_passed??b.validation_passed)===true,resolverPassed=b.adaptive_extension_validation_passed===true,windows=(ext.windows||[3,5,7]).join('/');let decision=`stage 1: keep ≤ ${Number(p.keep_max_score).toFixed(3)} · direct exclude ≥ ${Number(p.exclude_min_score).toFixed(3)} only with detector and anatomical guardrail agreement`;if(tiers.length){let bands=tiers.map(t=>{const range=`${t.name} ${Number(t.min_score).toFixed(3)}–${t.max_score==null?'up':Number(t.max_score).toFixed(3)}`;return t.enable_exclude===false?`${range}: calibrated keep-only band (no safe incremental exclusion rule)`:`${range}: ≥${t.min_corroborating_domains} domains, severity ≥${Number(t.severe_domain_threshold).toFixed(2)}`}).join(' · ');decision+=` · stage 2 assesses the full review interval with ${windows}-slice stability (${bands})`}else if(legacy)decision+=` · legacy review resolver ≥ ${Number(p.adaptive_exclude_min_score).toFixed(3)} uses two-sided, unprotected, severe multi-domain evidence stable at ${windows} slices`;decision+=' · all slices failing an enabled exclusion route are kept';let validation=resolver?`base policy validation ${basePassed?'passed':'not passed'} · review resolver ${resolverPassed?'independent validation passed':'not independently validated'}`:`validation ${basePassed?'passed':'not passed'}`,cal=ext.calibration_evaluation||{},val=ext.validation_evaluation||{};return `<div class="bench"><b>${p.calibration_id}</b><p>${decision}</p><p class="muted">${validation}. Every slice receives one final action. Base calibration datasets ${(b.calibration_dataset_ids||[]).length} · base unseen validation datasets ${(b.validation_dataset_ids||[]).length} · review calibration trials ${(cal.positive_trials||0)+(cal.negative_trials||0)||'—'} · review unseen-validation trials ${(val.positive_trials||0)+(val.negative_trials||0)||'—'}</p></div>`}).join('');
const ds=document.getElementById('datasetSelect'), groups=[...new Set(D.datasets.map(d=>d.group))];ds.innerHTML=groups.map(g=>`<optgroup label="${g}">${D.datasets.map((d,i)=>d.group===g?`<option value="${i}">${d.name}</option>`:'').join('')}</optgroup>`).join('');
document.getElementById('groupFilter').innerHTML+=groups.map(g=>`<option value="${g}">${g}</option>`).join('');
document.getElementById('datasetFilter').innerHTML+=[...D.datasets].map(d=>`<option value="${d.name}">${d.name}</option>`).join('');
document.getElementById('statisticsGroup').innerHTML+=groups.map(g=>`<option value="${g}">${g}</option>`).join('');
function statFmt(x){if(x===null||x===undefined||x==='')return 'NA';let v=Number(x);if(!Number.isFinite(v))return 'NA';let a=Math.abs(v);if(a>0&&(a<.01||a>=1e6))return v.toExponential(2);return new Intl.NumberFormat('en-US',{maximumFractionDigits:a<10?3:a<100?1:0}).format(v)}
function renderStatistics(){let group=document.getElementById('statisticsGroup').value,sets=D.datasets.filter(d=>group==='all'||d.group===group),panels=(D.statistic_definitions||[]).map(def=>{let rows=sets.map(d=>({d,s:d.statistics&&d.statistics[def.key]})).filter(x=>x.s&&Number.isFinite(+x.s.median));if(!rows.length)return '';rows.sort((a,b)=>b.s.median-a.s.median);let raw=rows.flatMap(x=>[+x.s.q25,+x.s.median,+x.s.q75]).filter(Number.isFinite),positive=raw.filter(v=>v>0),lo=Math.min(...raw),hi=Math.max(...raw),logScale=positive.length===raw.length&&Math.max(...positive)/Math.max(Math.min(...positive),1e-300)>100;if(logScale){lo=Math.log10(Math.min(...positive));hi=Math.log10(Math.max(...positive))}else{lo=Math.min(0,lo)}if(!(hi>lo))hi=lo+1;let pos=v=>Math.max(0,Math.min(100,100*((logScale?Math.log10(Math.max(+v,1e-300)):+v)-lo)/(hi-lo)));let body=rows.map(x=>{let a=pos(x.s.q25),b=pos(x.s.q75),m=pos(x.s.median),context=x.d.measurement_context?.count_semantics||'';return `<div class="comparerow" data-href="${x.d.report_href}" title="${x.d.group}${x.d.subgroup?' · '+x.d.subgroup:''} · ${context}"><span class="comparename">${x.d.name}</span><span class="comparetrack"><span class="compareiqr" style="left:${a}%;width:${Math.max(2,b-a)}%"></span><span class="comparedot" style="left:calc(${m}% - 4px)"></span></span><span class="comparevalue">${statFmt(x.s.median)}</span></div>`}).join('');return `<div class="comparepanel"><div class="comparehead"><b title="${def.help}">${def.label}</b><span class="scale">${rows.length} datasets · ${logScale?'log':'linear'}</span></div><div class="comparebody">${body}</div></div>`}).filter(Boolean);document.getElementById('statisticsComparison').innerHTML=panels.join('')||'<p class="muted">No compatible statistics are available for this selection.</p>';document.querySelectorAll('.comparerow').forEach(n=>n.onclick=()=>{window.location.href=n.dataset.href})}
document.getElementById('datasetGroups').innerHTML=groups.map(g=>{let subsets=[...new Set(D.datasets.filter(d=>d.group===g).map(d=>d.subgroup||''))];return `<div class="datasetgroup"><h3>${g} <span class="muted">(${D.datasets.filter(d=>d.group===g).length})</span></h3>`+subsets.map(s=>`${s?`<div class="subgroup">${s}</div>`:''}<div class="datasetgrid">${D.datasets.map((d,i)=>d.group===g&&(d.subgroup||'')===s?`<div class="dataset" data-i="${i}"><h4>${d.name}</h4><div class="counts"><span class="pill">${d.n_slices} slices</span><span class="pill">window ${d.selected_window}</span><span class="pill">max ${Number(d.max_score).toFixed(3)}</span>${D.binary_only?`<span class="pill">${Math.round(100*d.coverage)}% classified</span><span class="badge keep">${d.keep}</span><span class="badge exclude">${d.exclude}</span>`:`<span class="badge keep">${d.keep}</span><span class="badge review">${d.review}</span><span class="badge exclude">${d.exclude}</span>`}</div></div>`:'').join('')}</div>`).join('')+`</div>`}).join('');
const cp=document.getElementById('catalogPanel'),statusLabel={serial_qc:'serial QC',serial_qc_grouped:'serial QC (grouped)',not_applicable_single_slice:'single slice · N/A',scan_failed:'scan failed',missing_or_incomplete:'missing/incomplete',downloaded_unclassified:'downloaded · unclassified'};if(!D.catalog.length)cp.style.display='none';else{let cg=[...new Set(D.catalog.map(x=>x.group))];document.getElementById('catalog').innerHTML=cg.map(g=>{let items=D.catalog.filter(x=>x.group===g),coh=[...new Set(items.map(x=>x.cohort||''))];return `<div class="cataloggroup"><h3>${g} <span class="muted">(${items.length} files)</span></h3>`+coh.map(c=>{let f=items.filter(x=>(x.cohort||'')===c),counts=[...new Set(f.map(x=>x.qc_status))].map(s=>`<span class="status ${s}">${statusLabel[s]||s}: ${f.filter(x=>x.qc_status===s).length}</span>`).join(' ');return `<details><summary>${c||'Other'} · ${f.length} files · ${counts}</summary><ul class="catalogfiles">${f.map(x=>`<li title="${x.local_path}">${x.filename} <span class="status ${x.qc_status}">${statusLabel[x.qc_status]||x.qc_status}</span></li>`).join('')}</ul></details>`}).join('')+`</div>`}).join('')}
function timeline(i){const d=D.datasets[i],R=d.timeline,svg=document.getElementById('timeline'),w=1200,h=260,p=36,iw=w-2*p,ih=h-2*p,pts=R.map((r,j)=>[p+(Number.isFinite(+r.slice_index)?+r.slice_index:j)*iw/Math.max(1,d.n_slices-1),h-p-Number(r.quality_anomaly_score)*ih]);let grid=[0,.25,.5,.75,1].map(v=>`<line x1="${p}" x2="${w-p}" y1="${h-p-v*ih}" y2="${h-p-v*ih}" stroke="#dce4ee"/><text x="5" y="${h-p-v*ih+4}" font-size="11" fill="#627287">${v.toFixed(2)}</text>`).join('');let line=`<polyline points="${pts.map(x=>x.join(',')).join(' ')}" fill="none" stroke="#315efb" stroke-width="2"/>`;let dots=pts.map((x,j)=>`<circle cx="${x[0]}" cy="${x[1]}" r="5" fill="${colors[R[j].recommendation]}"><title>${R[j].slice_id}: ${Number(R[j].quality_anomaly_score).toFixed(3)} (${R[j].recommendation})</title></circle>`).join('');svg.innerHTML=grid+line+dots;document.getElementById('openReport').href=d.report_href;ds.value=i}
ds.onchange=e=>timeline(+e.target.value);document.querySelectorAll('.dataset').forEach(n=>n.onclick=()=>{timeline(+n.dataset.i);document.getElementById('timeline').scrollIntoView({behavior:'smooth',block:'center'})});
function rows(){const cf=document.getElementById('callFilter').value,gf=document.getElementById('groupFilter').value,df=document.getElementById('datasetFilter').value,q=document.getElementById('search').value.toLowerCase();let r=D.candidates.filter(x=>(cf==='all'||x.recommendation===cf)&&(gf==='all'||x.group===gf)&&(df==='all'||x.dataset===df)&&(`${x.slice_id} ${x.reason}`).toLowerCase().includes(q));r.sort((a,b)=>b.score-a.score);document.getElementById('rows').innerHTML=r.map(x=>`<tr><td>${x.group}</td><td><a href="${x.report_href}">${x.dataset}</a></td><td>${x.slice_id}</td><td>${Number(x.score).toFixed(3)}</td><td><span class="badge ${x.recommendation}">${x.recommendation}</span></td><td>${Number(x.expression).toFixed(3)}</td><td>${Number(x.density).toFixed(3)}</td><td>${x.partial_protection?'yes':'no'}</td><td class="reason">${x.reason}</td></tr>`).join('')}
['callFilter','groupFilter','datasetFilter','search'].forEach(id=>document.getElementById(id).addEventListener(id==='search'?'input':'change',rows));
document.getElementById('statisticsGroup').onchange=renderStatistics;
if(D.binary_only){document.getElementById('queueTitle').textContent='Automatic exclusions';document.getElementById('callFilter').innerHTML='<option value="exclude">Exclude</option>'}
const bp=document.getElementById('benchmarkPanel');if(!D.benchmarks.length)bp.style.display='none';else document.getElementById('benchmarks').innerHTML=D.benchmarks.map(b=>{const s=b.summary,items=b.injected.map(x=>`<tr><td>${x.slice_id}</td><td>${x.artifact}</td><td>${x.baseline_recommendation}</td><td>${x.simulated_recommendation}</td><td class="delta">${Number(x.baseline_score).toFixed(3)} → ${Number(x.simulated_score).toFixed(3)} (+${Number(x.score_delta).toFixed(3)})</td></tr>`).join('');return `<div class="bench"><b>${b.dataset}</b> · window ${b.selected_window} · <span class="ok">${s.injected_detected_after}/${s.injected_total} detected</span> · median score Δ ${Number(s.injected_median_score_delta).toFixed(3)} · non-injected newly flagged ${s.noninjected_newly_flagged}<table><thead><tr><th>slice</th><th>artifact</th><th>before</th><th>after</th><th>score change</th></tr></thead><tbody>${items}</tbody></table><a href="${b.report_href}">Open synthetic report →</a></div>`}).join('');
timeline(0);rows();renderStatistics();
</script></body></html>"""
    rendered = template.replace("__TITLE__", html.escape(title)).replace(
        "__PAYLOAD__", payload_json
    )
    output_path.write_text(rendered, encoding="utf-8")
    dataset_csv = output_path.with_name(f"{output_path.stem}_dataset_summary.csv")
    candidate_csv = output_path.with_name(f"{output_path.stem}_candidate_slices.csv")
    benchmark_csv = output_path.with_name(
        f"{output_path.stem}_expression_benchmarks.csv"
    )
    dataset_rows: list[dict[str, Any]] = []
    for item in datasets:
        row = {
            key: value
            for key, value in item.items()
            if key not in {"timeline", "statistics", "measurement_context"}
        }
        row["count_semantics"] = item.get("measurement_context", {}).get(
            "count_semantics"
        )
        for metric_key, summary in item.get("statistics", {}).items():
            for summary_key in ("median", "q25", "q75", "n_slices"):
                row[f"{metric_key}_{summary_key}"] = summary.get(summary_key)
        dataset_rows.append(row)
    pd.DataFrame(dataset_rows).to_csv(dataset_csv, index=False)
    pd.DataFrame(candidates).to_csv(candidate_csv, index=False)
    benchmark_rows: list[dict[str, Any]] = []
    for benchmark in benchmarks:
        summary = benchmark["summary"]
        for injected in benchmark["injected"]:
            benchmark_rows.append(
                {
                    "dataset": benchmark["dataset"],
                    "selected_window": benchmark["selected_window"],
                    "slice_id": injected.get("slice_id"),
                    "artifact": injected.get("artifact"),
                    "baseline_recommendation": injected.get("baseline_recommendation"),
                    "simulated_recommendation": injected.get(
                        "simulated_recommendation"
                    ),
                    "baseline_score": injected.get("baseline_score"),
                    "simulated_score": injected.get("simulated_score"),
                    "score_delta": injected.get("score_delta"),
                    "detected_after": injected.get("detected_after"),
                    "noninjected_newly_flagged_for_run": summary.get(
                        "noninjected_newly_flagged"
                    ),
                }
            )
    pd.DataFrame(benchmark_rows).to_csv(benchmark_csv, index=False)
    return {
        "report": str(output_path),
        "dataset_summary_csv": str(dataset_csv),
        "candidate_slices_csv": str(candidate_csv),
        "expression_benchmarks_csv": str(benchmark_csv),
        "datasets": int(len(datasets)),
        "catalog_files": int(len(catalog)),
        "slices": int(sum(item["n_slices"] for item in datasets)),
        "keep": int(sum(item["keep"] for item in datasets)),
        "review": int(sum(item["review"] for item in datasets)),
        "exclude": int(sum(item["exclude"] for item in datasets)),
        "published": int(sum(item["published"] for item in datasets)),
        "withheld": int(sum(item["withheld"] for item in datasets)),
        "benchmarks": int(len(benchmarks)),
    }


__all__ = [
    "SliceQCConfig",
    "SliceSeriesResult",
    "HighConfidencePolicy",
    "add_multiscale_exclusion_evidence",
    "apply_high_confidence_policy",
    "calculate_slice_quality",
    "scan_h5ad_series",
    "simulate_slice_quality_artifacts",
    "evaluate_paired_simulation",
    "evaluate_slice_calls",
    "write_slice_quality_outputs",
    "write_high_confidence_outputs",
    "render_slice_quality_report_from_outputs",
    "write_slice_quality_collection_report",
]


__all__ = [
    "render_slice_quality_report_from_outputs",
    "render_binary_slice_quality_report_from_outputs",
    "write_slice_quality_collection_report",
]
