/* Extension of the existing triad, table, timeline and export report. */
(()=>{
const P=JSON.parse(document.getElementById('preregistrationPayload').textContent);
// Reuse the report's point arrays; accept older self-contained extension data.
if(!P.point_samples&&typeof D!=='undefined')P.point_samples=D.points;
const M=P.preregistration;
if(!M)return;
const toolbar=document.getElementById('sliceSelect').closest('.toolbar');
const mode=document.createElement('select');mode.id='coordinateMode';
mode.innerHTML='<option value="display">Display preregistration · shared frame</option><option value="raw">Original coordinates · shared raw axes</option>';
toolbar.appendChild(mode);
const info=document.createElement('div');info.id='roiInfo';info.className='note';
const panel=toolbar.closest('.panel');panel.appendChild(info);
const fitNote=document.createElement('p');fitNote.className='note';
fitNote.textContent=`Display fit: rotation + translation; up to ${M.parameters.sample_size} fitting points per slice; seed ${M.parameters.seed}. ${M.parameters.annotation_weight?`Cell-type candidate weight ${M.parameters.annotation_weight}; annotation: ${[...new Set(Object.values(M.annotation_sources||{}))].join(', ')||'unavailable'}.`:'Geometry-only candidate selection.'}`;
panel.appendChild(fitNote);
const controls=document.createElement('div');controls.className='toolbar';
controls.innerHTML='<button id="clearRoi">Clear ROI</button><button id="exportRoi">Export ROI + sampled point IDs</button>';
panel.appendChild(controls);
const caution=document.createElement('p');caution.style.color='#a34f00';caution.textContent='Display only. Same ROI does not guarantee the same anatomy. No scaling, reflection or local deformation. Selection counts below refer to displayed samples; replay the exported ROI on the server for all points.';panel.appendChild(caution);
let bounds=null, drag=null, maps={}, active=[];
const sliceOrder=(typeof D!=='undefined'&&Array.isArray(D.order)?D.order:R.map(r=>r.slice_id)).map(String);
const points=id=>P.point_samples[id];
const coords=d=>mode.value==='display'?[d.display_x,d.display_y]:[d.x,d.y];
function windowBounds(){
 const focal=sliceOrder.indexOf(String(R[selected]?.slice_id));
 active=sliceOrder.slice(Math.max(0,focal-1),focal+2);
 let xmin=Infinity,xmax=-Infinity,ymin=Infinity,ymax=-Infinity;
 active.forEach(id=>{const d=points(id);if(!d)return;const [x,y]=coords(d);for(let i=0;i<x.length;i++){if(!Number.isFinite(x[i])||!Number.isFinite(y[i]))continue;xmin=Math.min(xmin,x[i]);xmax=Math.max(xmax,x[i]);ymin=Math.min(ymin,y[i]);ymax=Math.max(ymax,y[i]);}});
 return [xmin,xmax,ymin,ymax];
}
function selection(d){if(!bounds||mode.value!=='display')return [];return d.display_x.map((x,i)=>i).filter(i=>d.display_x[i]>=bounds[0]&&d.display_x[i]<=bounds[1]&&d.display_y[i]>=bounds[2]&&d.display_y[i]<=bounds[3]);}
function updateInfo(){
 const details=active.map(id=>{const d=points(id),s=M.slices[id];return `${id}: ${s.status}; ${s.warnings.join(', ')||'anatomical accuracy unvalidated'}; selected ${selection(d).length}/${d.n_displayed} shown (${d.n_total} total); raw ${d.coordinate_source}`;});
 info.textContent=`Frame: ${mode.value==='display'?M.frame_id:'raw coordinates (not comparable anatomy)'} | reference: ${M.reference_slice}. ${details.join(' | ')}`;
}
draw=function(canvas,id){
 const ctx=canvas.getContext('2d'),d=points(id),w=canvas.clientWidth,h=canvas.clientHeight,ratio=devicePixelRatio;
 canvas.width=w*ratio;canvas.height=h*ratio;ctx.setTransform(ratio,0,0,ratio,0,0);ctx.fillStyle='#07101f';ctx.fillRect(0,0,w,h);
 if(!d)return;const [x,y]=coords(d),b=windowBounds();
 const scale=Math.min((w-20)/Math.max(b[1]-b[0],1e-9),(h-20)/Math.max(b[3]-b[2],1e-9));
 const ox=(w-scale*(b[1]-b[0]))/2,oy=(h-scale*(b[3]-b[2]))/2;
 const px=x=>ox+(x-b[0])*scale,py=y=>h-oy-(y-b[2])*scale;
 maps[canvas.id]={toData:(cx,cy)=>[b[0]+(cx-ox)/scale,b[2]+(h-oy-cy)/scale]};
 const values=active.flatMap(s=>points(s)?.counts||[]).filter(Number.isFinite).sort((a,b)=>a-b);
 const lo=values[Math.floor(.04*(values.length-1))]||0,hi=values[Math.floor(.96*(values.length-1))]||1;
 const hits=new Set(selection(d));
 for(let i=0;i<x.length;i++){
  const v=Number.isFinite(d.counts[i])?Math.max(0,Math.min(1,(d.counts[i]-lo)/Math.max(hi-lo,1e-9))):.5;
  ctx.fillStyle=hits.has(i)?'#ffe178':document.getElementById('colorMode').value==='plain'?'#6db7ff':`rgb(${Math.round(45+210*v)},${Math.round(95+90*(1-v))},${Math.round(230-170*v)})`;
  ctx.globalAlpha=bounds&&!hits.has(i)&&mode.value==='display'?.28:.82;
  ctx.beginPath();ctx.arc(px(x[i]),py(y[i]),hits.has(i)?2.1:1.4,0,2*Math.PI);ctx.fill();
 }
 ctx.globalAlpha=1;
 if(bounds&&mode.value==='display'){ctx.strokeStyle='#ffe178';ctx.lineWidth=2;ctx.strokeRect(px(bounds[0]),py(bounds[3]),(bounds[1]-bounds[0])*scale,(bounds[3]-bounds[2])*scale);}
 ctx.fillStyle='#ffffff';ctx.font='10px sans-serif';ctx.fillText(values.length?`counts: ${lo.toPrecision(3)}–${hi.toPrecision(3)} (shared)`:'Expression capture unavailable (annotation matrix)',6,h-5);
 updateInfo();
};
const oldSelect=select;
select=function(i){oldSelect(i);updateInfo();};
const mid=document.getElementById('mid');mid.style.touchAction='none';
mid.addEventListener('pointerdown',e=>{if(mode.value!=='display')return;const r=mid.getBoundingClientRect();drag=maps.mid.toData(e.clientX-r.left,e.clientY-r.top);mid.setPointerCapture(e.pointerId);});
mid.addEventListener('pointermove',e=>{if(!drag)return;const r=mid.getBoundingClientRect(),q=maps.mid.toData(e.clientX-r.left,e.clientY-r.top);bounds=[Math.min(drag[0],q[0]),Math.max(drag[0],q[0]),Math.min(drag[1],q[1]),Math.max(drag[1],q[1])];select(selected);});
mid.addEventListener('pointerup',()=>{drag=null;});
mode.onchange=()=>select(selected);
document.getElementById('clearRoi').onclick=()=>{bounds=null;select(selected);};
document.getElementById('exportRoi').onclick=()=>{
 if(!bounds||mode.value!=='display'){info.textContent='Draw a rectangle on the middle slice in display mode first.';return;}
 const exportData={frame_id:M.frame_id,coordinate_mode:'display',bounds,slice_ids:active,reference_slice:M.reference_slice,
  selection_scope:'displayed sample only; use export-roi for exact full-data replay',
  selected_points:Object.fromEntries(active.map(id=>{const d=points(id),s=M.slices[id];return [id,{status:s.status,warnings:s.warnings,source_file:d.source_file,
    raw_roi_polygon:[[bounds[0],bounds[2]],[bounds[1],bounds[2]],[bounds[1],bounds[3]],[bounds[0],bounds[3]]].map(([x,y])=>{const t=s.inverse_matrix;return [t[0][0]*x+t[0][1]*y+t[0][2],t[1][0]*x+t[1][1]*y+t[1][2]];}),
    points:selection(d).map(i=>({obs_name:d.obs_names[i],source_row:d.source_row[i],raw_xy:[d.x[i],d.y[i]],display_xy:[d.display_x[i],d.display_y[i]]}))}];}))};
 const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([JSON.stringify(exportData,null,2)],{type:'application/json'}));a.download='display_roi.json';a.click();URL.revokeObjectURL(a.href);
};
window.refereeROI={setBounds:b=>{bounds=b;select(selected);},getBounds:()=>bounds,selectedCounts:()=>Object.fromEntries(active.map(id=>[id,selection(points(id)).length])),frameId:M.frame_id};
select(selected);window.addEventListener('resize',()=>select(selected));
})();
