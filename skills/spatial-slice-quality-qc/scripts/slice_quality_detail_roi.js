/* Adapter for the existing multi-panel binary appendix; no QC recomputation. */
(()=>{
const M=D.preregistration;if(!M)return;
const toolbar=document.getElementById('sliceSelect').closest('.toolbar');
const mode=document.createElement('select');mode.id='coordinateMode';
mode.innerHTML='<option value="display">Display preregistration</option><option value="raw">Original coordinates</option>';toolbar.appendChild(mode);
for(const [id,label] of [['clearDetailRoi','Clear ROI'],['exportDetailRoi','Export ROI + point IDs']]){const b=document.createElement('button');b.id=id;b.textContent=label;toolbar.appendChild(b);}
const info=document.createElement('p');info.id='detailRoiInfo';info.className='note';toolbar.parentElement.appendChild(info);
const note=document.createElement('p');note.className='note';
note.textContent=`Frame: ${M.frame_id}; reference: ${M.reference_slice}. Drag on the focal slice in any evidence panel. Same display bounds do not guarantee the same anatomy. Display fitting never changes QC. ${D.evidence_source.note||''}`;
toolbar.parentElement.appendChild(note);
const points=D.points;for(const p of Object.values(points)){p.raw_x=p.x;p.raw_y=p.y;}
let bounds=null,drag=null,maps={},active=[];
const oldDraw=drawCanvas,oldEvidence=drawEvidence;
function selected(id){const p=points[id];return bounds&&mode.value==='display'?p.display_x.map((x,i)=>i).filter(i=>p.display_x[i]>=bounds[0]&&p.display_x[i]<=bounds[1]&&p.display_y[i]>=bounds[2]&&p.display_y[i]<=bounds[3]):[];}
function updateInfo(){active=currentIds();info.textContent=active.map(id=>{const p=points[id],s=M.slices[id];return `${id} · ${byId[id].final_call}; ${s.status}; ${s.warnings.join(', ')||'no registration warning'}; selected ${selected(id).length}/${p.displayed} displayed of ${p.available}; raw ${p.coordinate_source}`;}).join(' | ');}
function toData(canvas,event){const m=maps[canvas.id],r=canvas.getBoundingClientRect();return [(event.clientX-r.left-m.ox)/m.s+m.e[0],(m.h-m.oy-(event.clientY-r.top))/m.s+m.e[2]];}
function redraw(){drawEvidence();}
drawCanvas=function(canvas,id,extent,metric,lo,hi,focus,region){
 oldDraw(canvas,id,extent,metric,lo,hi,focus,mode.value==='display'?region:null);
 const p=points[id];if(!p?.x.length)return;
 const w=canvas.clientWidth,h=canvas.clientHeight,e=extent||[Math.min(...p.x),Math.max(...p.x),Math.min(...p.y),Math.max(...p.y)],s=Math.min((w-18)/Math.max(e[1]-e[0],1e-9),(h-18)/Math.max(e[3]-e[2],1e-9)),ox=(w-(e[1]-e[0])*s)/2,oy=(h-(e[3]-e[2])*s)/2;
 maps[canvas.id]={e,s,ox,oy,h,canvas,args:[id,extent,metric,lo,hi,focus,region]};const ctx=canvas.getContext('2d'),px=x=>ox+(x-e[0])*s,py=y=>h-oy-(y-e[2])*s;
 if(bounds&&mode.value==='display'){ctx.strokeStyle='#ffe178';ctx.lineWidth=1.5;ctx.setLineDash([]);ctx.strokeRect(px(bounds[0]),py(bounds[3]),(bounds[1]-bounds[0])*s,(bounds[3]-bounds[2])*s);for(const i of selected(id)){ctx.beginPath();ctx.arc(px(p.x[i]),py(p.y[i]),3.3,0,2*Math.PI);ctx.stroke();}}
 if(!focus)return;canvas.style.touchAction='none';
 canvas.onpointerdown=ev=>{if(mode.value!=='display')return;drag={canvasId:canvas.id,start:toData(canvas,ev)};canvas.setPointerCapture(ev.pointerId);};
 // Events continue on the captured canvas even when the report redraw replaces it.
 canvas.onpointermove=ev=>{if(!drag)return;const q=toData(canvas,ev),a=drag.start;bounds=[Math.min(a[0],q[0]),Math.max(a[0],q[0]),Math.min(a[1],q[1]),Math.max(a[1],q[1])];for(const m of Object.values(maps))drawCanvas(m.canvas,...m.args);updateInfo();};
 canvas.onpointerup=()=>{drag=null;redraw();};canvas.onpointercancel=()=>{drag=null;};
};
drawEvidence=function(){maps={};for(const p of Object.values(points)){p.x=mode.value==='display'?p.display_x:p.raw_x;p.y=mode.value==='display'?p.display_y:p.raw_y;}if(mode.value==='display')document.getElementById('extentSelect').value='shared';oldEvidence();if(mode.value==='raw')document.getElementById('regionNote').textContent='Automatic cross-slice deficit boxes are hidden in the raw coordinate frame.';updateInfo();};
mode.onchange=()=>{bounds=null;redraw();};
// Previous handlers captured the original function; explicitly reconnect them.
for(const id of ['sliceSelect','windowSelect','expressionSelect','extentSelect'])document.getElementById(id).onchange=redraw;
document.getElementById('clearDetailRoi').onclick=()=>{bounds=null;redraw();};
document.getElementById('exportDetailRoi').onclick=()=>{
 if(!bounds||mode.value!=='display'){info.textContent='Draw an ROI in display mode first.';return;}
 const data={frame_id:M.frame_id,coordinate_mode:'display',bounds,slice_ids:active,reference_slice:M.reference_slice,selection_scope:'displayed sample; replay export-roi on full caches',selected_points:{}};
 for(const id of active){const p=points[id],s=M.slices[id],t=s.inverse_matrix;data.selected_points[id]={source_file:p.source_file,status:s.status,warnings:s.warnings,raw_roi_polygon:[[bounds[0],bounds[2]],[bounds[1],bounds[2]],[bounds[1],bounds[3]],[bounds[0],bounds[3]]].map(([x,y])=>[t[0][0]*x+t[0][1]*y+t[0][2],t[1][0]*x+t[1][1]*y+t[1][2]]),points:selected(id).map(i=>({source_row:p.source_row[i],obs_name:p.obs_names[i],raw_xy:[p.raw_x[i],p.raw_y[i]],display_xy:[p.display_x[i],p.display_y[i]]}))};}
 const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([JSON.stringify(data,null,2)],{type:'application/json'}));a.download='display_roi.json';a.click();URL.revokeObjectURL(a.href);
};
window.refereeDetailROI={frameId:M.frame_id,getBounds:()=>bounds,setBounds:b=>{bounds=b;redraw();},selectedCounts:()=>Object.fromEntries(currentIds().map(id=>[id,selected(id).length]))};
redraw();
})();
