// Actual detailed report scripts; DOM/canvas substitutes, not browser rendering.
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const html=fs.readFileSync(process.argv[2],'utf8'),ids={},downloads=[];
function node(id=''){
 const e={style:{},value:'',textContent:'',clientWidth:300,clientHeight:200,events:{},dataset:{},
  classList:{contains:()=>true,toggle(){}},appendChild(){},querySelectorAll:()=>[],setAttribute(){},
  closest:s=>s==='.toolbar'?toolbar:panel,addEventListener(n,f){this.events[n]=f;},
  getBoundingClientRect:()=>({left:0,top:0}),setPointerCapture(){},click(){if(this.onclick)this.onclick();},
  getContext:()=>new Proxy({},{get:()=>()=>{},set:()=>true})};
 Object.defineProperty(e,'id',{get(){return this._id;},set(v){this._id=v;ids[v]=this;}});
 Object.defineProperty(e,'innerHTML',{get(){return this._html||'';},set(v){this._html=v;if(v.includes('value="display"'))this.value='display';}});
 e.id=id;return e;
}
const toolbar=node('toolbar'),panel=node('panel');toolbar.parentElement=panel;
const document={getElementById:id=>ids[id]||node(id),querySelectorAll:()=>[],querySelector:()=>node(),createElement:()=>node()};
for(const [id,v] of [['extentSelect','shared'],['expressionSelect','log_captured_counts'],['callFilter','all']])document.getElementById(id).value=v;
const context={document,window:{devicePixelRatio:1,addEventListener(){}},requestAnimationFrame:f=>f(),Blob,URL:{createObjectURL:b=>{downloads.push(b);return 'test';},revokeObjectURL(){}}};vm.createContext(context);
for(const [,attrs,code] of html.matchAll(/<script([^>]*)>([\s\S]*?)<\/script>/g)){
 if(attrs.includes('application/json'))document.getElementById(attrs.match(/id="([^"]+)"/)[1]).textContent=code;
 else vm.runInContext(code,context,{filename:process.argv[2]});
}
assert.strictEqual(vm.runInContext('fmt(null)',context),'NA');
const D=JSON.parse(document.getElementById('payload').textContent);assert(D.records.every(r=>['keep','exclude'].includes(r.final_call)));
assert.strictEqual(new Set(D.records.map(r=>r.slice_id)).size,D.order.length);
for(const windowSize of [3,5]){
 document.getElementById('windowSelect').value=String(windowSize);document.getElementById('windowSelect').onchange();
 const active=JSON.parse(vm.runInContext('JSON.stringify(currentIds())',context));assert.strictEqual(active.length,Math.min(windowSize,D.order.length));
 for(const id of ['kdeWindow','expressionWindow','componentWindow'])assert.strictEqual((document.getElementById(id).innerHTML.match(/<canvas /g)||[]).length,active.length);
}
const finite=D.statistics.reduce((n,s)=>n+D.records.filter(r=>r[s.key]!=null&&Number.isFinite(Number(r[s.key]))).length,0);
assert.strictEqual((document.getElementById('plots').innerHTML.match(/data-id=/g)||[]).length,finite);
async function finish(){
 let selectedCount=null;
 if(D.preregistration){
  const mid=D.order[Math.floor(D.order.length/2)];document.getElementById('sliceSelect').value=mid;document.getElementById('sliceSelect').onchange();
  const active=JSON.parse(vm.runInContext('JSON.stringify(currentIds())',context)),idx=active.indexOf(mid),canvas=document.getElementById('kdeWindow'+idx);
  canvas.onpointerdown({clientX:80,clientY:40,pointerId:1});canvas.onpointermove({clientX:210,clientY:160});canvas.onpointerup();
  const bounds=Array.from(context.window.refereeDetailROI.getBounds());assert(bounds[0]<bounds[1]&&bounds[2]<bounds[3]);
  const expected=Object.fromEntries(active.map(id=>[id,D.points[id].display_x.filter((x,i)=>x>=bounds[0]&&x<=bounds[1]&&D.points[id].display_y[i]>=bounds[2]&&D.points[id].display_y[i]<=bounds[3]).length]));
  assert.deepStrictEqual(JSON.parse(JSON.stringify(context.window.refereeDetailROI.selectedCounts())),expected);selectedCount=expected;
  document.getElementById('exportDetailRoi').onclick();assert.strictEqual(downloads.length,1);
  const data=JSON.parse(await downloads[0].text());if(process.argv[3])fs.writeFileSync(process.argv[3],JSON.stringify(data,null,2));assert.strictEqual(data.frame_id,D.preregistration.frame_id);
  for(const id of active){const p=D.points[id],entry=data.selected_points[id],m=D.preregistration.slices[id].matrix;
   assert.strictEqual(entry.points.length,expected[id]);
   entry.points.forEach(point=>{const i=p.source_row.indexOf(point.source_row);assert.strictEqual(point.obs_name,p.obs_names[i]);assert.deepStrictEqual(point.raw_xy,[p.x[i],p.y[i]]);});
   entry.raw_roi_polygon.forEach(([x,y],i)=>{const corner=[[bounds[0],bounds[2]],[bounds[1],bounds[2]],[bounds[1],bounds[3]],[bounds[0],bounds[3]]][i];assert(Math.abs(m[0][0]*x+m[0][1]*y+m[0][2]-corner[0])<1e-7);assert(Math.abs(m[1][0]*x+m[1][1]*y+m[1][2]-corner[1])<1e-7);});
  }
  for(const sid of D.order){document.getElementById('sliceSelect').value=sid;document.getElementById('sliceSelect').onchange();for(const w of D.preregistration.slices[sid].warnings)assert(document.getElementById('detailRoiInfo').textContent.includes(w));}
  document.getElementById('coordinateMode').value='raw';document.getElementById('coordinateMode').onchange();assert.strictEqual(context.window.refereeDetailROI.getBounds(),null);
  assert(vm.runInContext('Object.values(D.points).every(p=>p.x===p.raw_x&&p.y===p.raw_y)',context));
 }
 console.log(JSON.stringify({html:process.argv[2],records:D.records.length,binary_only:true,three_evidence_panels:true,statistics:true,windows:[3,5],selected_sample_counts:selectedCount,browser_rendering_tested:false}));
}
finish().catch(e=>{console.error(e);process.exitCode=1;});
