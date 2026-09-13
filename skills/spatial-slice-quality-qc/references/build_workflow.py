"""Build an Illustrator-editable branched SVG from the current Referee policy.
No raster images, foreignObject, CSS layout, external fonts, or SVG markers.
Each panel and object is a named group; text remains live text.
"""
from pathlib import Path
import argparse, json, math
from xml.sax.saxutils import escape

p=argparse.ArgumentParser();p.add_argument('--policy',type=Path,default=Path(__file__).resolve().parents[1]/'policies/joint_review_v2.json');p.add_argument('--output',type=Path,default=Path(__file__).with_name('workflow.svg' if Path(__file__).parent.name=='references' else 'spateo_referee_complete_editable.svg'));a=p.parse_args()
P=json.loads(a.policy.read_text())['policy'];L,H=P['review_exclusion_tiers'];K=P['keep_max_score'];E=P['exclude_min_score'];B=L['max_score']
W,HEIGHT=3900,3200
C={'ink':'#183047','muted':'#56697D','line':'#657B91','blue':'#EAF2FB','teal':'#E3F3F0','purple':'#EEE9F7','amber':'#FFF3D9','red':'#FBE6E8','green':'#E5F3E8','white':'#FFFFFF'}
parts=[f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" width="{W}pt" height="{HEIGHT}pt" viewBox="0 0 {W} {HEIGHT}" version="1.1">
<title>Spateo Referee — signal-dependent QC, window handling and binary decisions</title>
<desc>Editable vector workflow. Joint review v2. Panels A–F contain current implementation branches; code function names in object descriptions provide traceability. All labels are live Arial text. No raster or embedded HTML.</desc>
<rect x="0" y="0" width="3900" height="3200" fill="#FFFFFF"/>
''']
objects=[];panel=None;ox=oy=0;serial=0

def txt(x,y,s,size=24,bold=False,anchor='middle',color='ink'):
 parts.append(f'<text x="{x}" y="{y}" font-family="Arial" font-size="{size}" font-weight="{700 if bold else 400}" text-anchor="{anchor}" fill="{C.get(color,color)}">{escape(str(s))}</text>')

def group(id,source=''):
 parts.append(f'<g id="{id}"><title>{escape(id.replace("_"," "))}</title><desc>{escape(source)}</desc>')

def box(id,x,y,w,h,title,lines=(),fill='blue',size=24,source=''):
 group(id,source);parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="{C[fill]}" stroke="{C["line"]}" stroke-width="1.8"/>')
 txt(x+w/2,y+34,title,26,True)
 for i,s in enumerate(lines):txt(x+w/2,y+67+29*i,s,size)
 parts.append('</g>');objects.append(dict(panel=panel,id=id,type='box',x=x,y=y,w=w,h=h,source=source))

def diamond(id,cx,cy,w,h,lines,source=''):
 group(id,source);parts.append(f'<path d="M {cx} {cy-h/2} L {cx+w/2} {cy} L {cx} {cy+h/2} L {cx-w/2} {cy} Z" fill="{C["purple"]}" stroke="{C["line"]}" stroke-width="2"/>')
 for i,s in enumerate(lines):txt(cx,cy+(i-(len(lines)-1)/2)*27+8,s,23,True)
 parts.append('</g>');objects.append(dict(panel=panel,id=id,type='decision',x=cx-w/2,y=cy-h/2,w=w,h=h,source=source))

def arrow(pts,label='',at=None,color='line',dash=False):
 global serial
 serial+=1;col=C.get(color,color);group(f'connector_{serial}')
 points=' '.join(f'{x},{y}' for x,y in pts)
 parts.append(f'<polyline points="{points}" fill="none" stroke="{col}" stroke-width="2.5" stroke-linejoin="round"'+(' stroke-dasharray="9 7"' if dash else '')+'/>')
 (x0,y0),(x1,y1)=pts[-2:];theta=math.atan2(y1-y0,x1-x0);length=13;half=5
 bx=x1-length*math.cos(theta);by=y1-length*math.sin(theta)
 parts.append(f'<path d="M{x1},{y1} L{bx-half*math.sin(theta)},{by+half*math.cos(theta)} L{bx+half*math.sin(theta)},{by-half*math.cos(theta)} Z" fill="{col}"/>')
 if label:
  x,y=at;tw=len(label)*12+16;parts.append(f'<rect x="{x-tw/2}" y="{y-23}" width="{tw}" height="30" fill="white"/>');txt(x,y,label,22,True,color='muted')
 parts.append('</g>')

def badge(x,y,label,fill='teal'):
 parts.append(f'<circle cx="{x}" cy="{y}" r="27" fill="{C[fill]}" stroke="{C["line"]}" stroke-width="2"/>');txt(x,y+8,label,24,True)

def begin(letter,x,y,title,subtitle):
 global panel,ox,oy
 panel=letter;ox=x;oy=y
 parts.append(f'<g id="panel_{letter}" inkscape:groupmode="layer" inkscape:label="{letter} - {escape(title)}" transform="translate({x},{y})">')
 parts.append('<rect x="0" y="0" width="1200" height="1410" rx="20" fill="#FFFFFF" stroke="#C3D0DD" stroke-width="2"/>')
 parts.append('<path d="M20 0 H1180 Q1200 0 1200 20 V112 H0 V20 Q0 0 20 0" fill="#F1F5F9"/>')
 badge(52,47,letter,'blue');txt(101,45,title,32,True,'start');txt(101,83,subtitle,22,False,'start','muted')

def end():parts.append('</g>')

def note(y,lines):
 for i,s in enumerate(lines):txt(600,y+27*i,s,21,color='muted')

# Overall title and non-linear panel routing.
group('title_and_reading_key');txt(60,66,'Spateo Referee',49,True,'start');txt(60,110,'Signal-dependent slice QC • sliding windows • evidence aggregation • auditable keep / exclude',28,False,'start','muted')
txt(3840,54,'Joint review v2  |  13 September 2026',24,True,'end');txt(3840,94,'Metric-stress tested; independent biological validation incomplete',23,False,'end','muted');parts.append('</g>')

# A: three input branches plus distinct signal routes.
begin('A',60,160,'Input contracts & signal routing','One biological series at a time; preserve raw coordinates and identities')
box('a_input',360,138,480,92,'Input format',['H5AD or ordered slice files'],source='scan_h5ad_series; scan_h5ad_collection')
for x,title,lines in [(30,'Slice labels',['Use selected obs field']),(430,'Discrete z',['Derive z-based slice IDs']),(830,'One file / slice',['Natural sort or manifest'])]:
 box('a_'+title.replace(' ','_').replace('/',''),x,290,340,106,title,lines,size=22,source='_resolve_slice_labels; _resolve_order')
 arrow([(600,230),(600,258),(x+170,258),(x+170,290)])
for x in [200,600,1000]:arrow([(x,396),(x,427),(600,427),(600,454)])
box('a_contract',250,454,700,107,'Resolve slice order + original x/y',['Record units, selected layer and source hashes'],source='_resolve_xy; _resolve_order; _extract_one_file')
arrow([(950,505),(1145,505),(1145,584)],'invalid',(1048,485))
box('a_stop',875,584,290,100,'Stop affected run',['Do not invent fields'],fill='red',size=21)
arrow([(600,561),(600,605)]);diamond('a_signals',600,690,470,170,['Which signals','are available?'],'_extract_one_file')
# These are concurrent availability routes, not mutually exclusive.
arrow([(365,690),(165,690),(165,813)],'coordinates',(248,673))
arrow([(600,775),(600,813)],'matrix',(660,793))
arrow([(835,690),(1050,690),(1050,813)],'annotations',(1050,758))
box('a_geometry',25,813,285,164,'Spatial geometry',['Finite x/y points','Hull, NN, holes, graph','Degenerate → limited metrics'],fill='teal',size=20,source='_geometry_metrics: n<2 returns NaN geometry and conservative defaults; invalid x/y excluded from geometry')
box('a_labels',890,813,285,164,'Cell-type labels',['Present → type proportions','Absent → no type evidence','Not an expression matrix'],fill='teal',size=20,source='_celltype_js_profiles; _extract_one_file')
diamond('a_matrix_type',600,890,485,140,['Selected matrix','semantics?'],'_choose_layer; _is_count_like; blind_input_contract.expression_representation')
for x,title,lines,fill in [
 (20,'Measured counts',['Row sums / nonzero genes','Override stale obs summaries','Capture + expression profile'],'blue'),
 (425,'Non-count matrix',['Use obs counts / genes if present','Else matrix-derived proxy','Record limitations; profile retained'],'amber'),
 (830,'Declared one-hot X',['Disable capture + mito proxies','Skip expression pseudobulk','Retain geometry + actual labels'],'amber')]:
 box('a_semantic_'+str(x),x,1055,350,163,title,lines,fill,size=20,source='_extract_one_file; _align_profiles')
arrow([(357,890),(340,890),(340,1005),(195,1005),(195,1055)],'counts',(240,1038))
arrow([(600,960),(600,1055)],'non-count',(689,1014))
arrow([(843,890),(861,890),(861,1005),(1005,1005),(1005,1055)],'one-hot',(1019,1038))
box('a_output',240,1280,720,84,'Per-slice metrics + availability + provenance',[],fill='teal',source='_extract_one_file')
for x in [195,600,1005]:arrow([(x,1218),(x,1250),(600,1250),(600,1280)])
arrow([(25,960),(10,960),(10,1322),(240,1322)]);arrow([(1175,960),(1190,960),(1190,1322),(960,1322)])
end()

# B: actual sliding-window illustration and endpoint decisions.
begin('B',1350,160,'Sliding windows & boundaries','Focal slice is excluded from its own local reference')
diamond('b_window_mode',600,199,400,125,['Window mode?'],'SliceQCConfig; _choose_window')
box('b_fixed',25,319,455,128,'Fixed window',['Default W = 3; odd width ≥3','Keep selected width in provenance'],size=22)
box('b_auto',620,319,555,156,'Automatic choice',['Candidates 3 / 5 / 7, limited by series length','Test ≤8 deterministic target positions','Maximize recorded perturbation utility'],size=22,source='_choose_window: utility=recall+.35*delta-.35*baseline_flags-.012*(W-3); no candidate -> W=3')
arrow([(400,199),(253,199),(253,319)],'fixed',(300,270));arrow([(800,199),(897,199),(897,319)],'auto',(927,270))
arrow([(253,447),(253,511),(600,511)]);arrow([(897,475),(897,511),(600,511)])
group('b_sliding_window_example','_local_expected uses positional slice order, not physical z spacing; focal excluded; no wrap or padding')
txt(600,562,'Move the focal index i through the ordered series',24,True)
for j in range(9):
 x=141+j*106;fill='purple' if j==4 else ('blue' if 2<=j<=6 else 'white')
 parts.append(f'<rect x="{x}" y="590" width="78" height="72" rx="6" fill="{C[fill]}" stroke="{C["line"]}" stroke-width="2"/>');txt(x+39,634,['…','i−3','i−2','i−1','i','i+1','i+2','i+3','…'][j],23,j==4)
parts.append('<rect x="340" y="579" width="522" height="95" rx="10" fill="none" stroke="#3E8B94" stroke-width="3" stroke-dasharray="10 7"/>')
arrow([(375,704),(825,704)],'example W = 5; advance i →',(600,735));parts.append('</g>')
diamond('b_context',600,832,560,130,['Finite neighbors for this metric?'],'_local_expected: finite neighbors within truncated window')
arrow([(600,748),(600,767)])
for x,title,lines,fill in [
 (25,'Both sides',['Exclude focal slice','Median pairwise-slope trend','Estimate reference at i'],'blue'),
 (425,'One side only',['Truncate at series boundary','Neighbor median reference','No padding or wraparound'],'amber'),
 (825,'No valid neighbor',['Reference = unavailable','Metric anomaly = 0 internally','0 means no evidence, not health'],'amber')]:
 box('b_context_'+str(x),x,993,350,163,title,lines,fill,size=21,source='_local_expected; _metric_anomaly')
arrow([(320,832),(200,832),(200,993)],'two-sided',(235,944));arrow([(600,897),(600,993)],'one-sided',(680,956));arrow([(880,832),(1000,832),(1000,993)],'none',(1014,944))
for x in [200,600,1000]:arrow([(x,1156),(x,1190),(600,1190),(600,1225)])
box('b_support',85,1225,1030,121,'Primary + robust support window',['Score both; take larger metric anomaly. Support width ≥W, usually up to 9.','Keep primary context; score_confidence = 1 if supported on both sides, else 0.72.'],fill='teal',size=21,source='_score_metrics: support=max(W, max(3, min(9, largest odd <=N))); internal=np.any(two_sided_matrix)')
note(1380,['Different specimens never share neighbors; missing slice IDs do not imply poor tissue.'])
end()

# C: directional normalization split then recombination.
begin('C',2640,160,'Metric standardization','Directional effect size + robust residual; no probability interpretation')
box('c_input',270,144,660,90,'Observed value v + local expectation e',[],source='_metric_anomaly')
diamond('c_direction',600,337,430,132,['Adverse direction?'],'_LOW_BAD_RULES; _HIGH_BAD_RULES')
arrow([(600,234),(600,271)])
box('c_low',25,480,540,193,'Lower is worse',['Counts, genes, density, tissue amount,','complexity, largest component fraction','m = max((e − v) / |e|, 0)'],size=22,source='_LOW_BAD_RULES; EPS prevents zero division')
box('c_high',635,480,540,193,'Higher is worse',['Holes, fragmentation, zero fraction,','NN tail, low-capture regions, mito','m = max(v − e, 0)'],size=22,source='_HIGH_BAD_RULES')
arrow([(385,337),(295,337),(295,480)],'decrease',(295,431));arrow([(815,337),(905,337),(905,480)],'increase',(905,431))
arrow([(295,673),(295,704),(600,704)]);arrow([(905,673),(905,704),(600,704)])
box('c_effect',25,755,540,162,'Effect-size score R',['R = clip((m − mild)/(severe − mild), 0, 1)','Metric-specific mild / severe limits','See parameter table in methods'],size=21,source='_metric_anomaly; _LOW_BAD_RULES; _HIGH_BAD_RULES')
box('c_z',635,755,540,179,'Robust residual score Z',['Residual: e−v (low) or v−e (high)','z = positive centered residual / scale','Z = clip((z − 1.5)/3, 0, 1)','Scale: 1.4826×MAD; fallback SD, then 1'],size=21,source='_robust_scale:1.4826*MAD; fallback std then1; _metric_anomaly')
arrow([(600,704),(295,704),(295,755)]);arrow([(600,704),(905,704),(905,755)])
box('c_small_high',635,951,540,128,'High-direction safeguard',['Z ← Z × clip(m / mild, 0, 1)','Tiny absolute effects cannot dominate'],fill='amber',size=22,source='_metric_anomaly: only direction=high')
arrow([(905,934),(905,951)],'high only',(1030,943))
arrow([(1175,845),(1187,845),(1187,1098),(1070,1098),(1070,1119)],'low: no rescaling',(1020,1089))
box('c_merge',35,1119,1130,174,'Combine, downweight and retain provenance',['A = max(R, Z) × context factor × agreement penalty','Context factor: 1 (two-sided) or 0.72 (one-sided)','Agreement penalty = clip(exp(−left/right disagreement / (2.5 × global scale)), .35, 1)'],fill='teal',size=21,source='_metric_anomaly; invalid observation/reference gives zero anomaly evidence')
arrow([(295,917),(295,1119)]);arrow([(905,1079),(905,1119)])
note(1334,['For counts / genes, retain the larger direct relative-drop score (same context penalty).','Primary and support scores are combined by max; missing evidence is not a normality claim.','Profile / type continuity uses its own divergence mapping before domain aggregation.'])
end()

# D: parallel domains then protection with different outcomes.
begin('D',60,1640,'Four domains → aggregate score','Signals may contribute to more than one domain; domains are not independent')
for x,y,title,lines,source in [
 (25,150,'Density / amount D',['Density .28; amount .10; holes .27','NN tail .13; fragments .12; largest .10','Max(weighted mean, top-two evidence)'],'_score_metrics: density_domain_score'),
 (625,150,'Expression capture X',['Counts .27; genes .24; zeros .14','Complexity .11; gene fraction .08','Local / regional low capture .08 / .08','Max(weighted mean, top-two evidence)'],'_score_metrics: expression_domain_score'),
 (25,402,'Damage proxies M',['Mito .65 + fragmentation .35','Weighted available signals'],'_score_metrics: damage_domain_score'),
 (625,402,'Continuity C',['Expression profile .58; cell types .22','Holes .20; weighted available signals','Unavailable profile / labels add no evidence'],'_profile_continuity; _celltype_js_profiles; _score_metrics')]:
 box('d_domain_'+title.split()[0],x,y,550,205 if y==150 else 155,title,lines,fill='blue' if y==150 else 'teal',size=22,source=source)
box('d_available',25,602,1150,92,'Within-domain rules',['Weighted mean renormalizes finite signals; top-two = .65 × strongest + .35 × second'],fill='white',size=22,source='_weighted_available; _top_two_evidence')
for x,y in [(300,355),(900,355)]:arrow([(x,y),(600,y),(600,582),(x,582),(x,602)])
arrow([(300,557),(300,602)]);arrow([(900,557),(900,602)])
box('d_rawscore',25,733,1150,144,'Aggregate anomaly score S₀',['S₀ = max(.31D + .39X + .12M + .18C, clip(.48 + .78(max(D,X,M,C) − .68), 0, 1))','Supporting-domain count = number of domains ≥.45'],fill='purple',size=21,source='_score_metrics: raw_score; corroborating_domains')
arrow([(600,694),(600,733)])
diamond('d_protection',600,969,495,125,['Anatomical protection?'],'_score_metrics: density_only OR taper_like')
arrow([(600,877),(600,906.5)])
box('d_guard_rule',25,1068,680,166,'Protection conditions (either suffices)',['D≥.48, X<.32, M<.38, C<.45; OR','amount anomaly≥.45, density anomaly<.30,','hole anomaly<.30, X<.35'],fill='amber',size=22,source='_score_metrics: density_only; taper_like')
box('d_guard_yes',785,1070,390,130,'Protected',['S ← min(S₀, .62)','Block automatic exclusion'],fill='amber',size=22)
arrow([(847.5,969),(980,969),(980,1070)],'yes',(1020,1022));arrow([(352.5,969),(12,969),(12,1286),(60,1286)],'no: S = S₀',(141,951))
box('d_endpoint',60,1260,1080,110,'Final score and boundary status',['If one-sided: S ← .82 × S. Otherwise keep S. Clip final score to [0,1].','Pass S, domain scores, context and protection flag to E.'],fill='teal',size=22,source='_score_metrics: final_score; window_context; score_confidence')
arrow([(980,1200),(980,1260)]);arrow([(365,1068),(365,1045),(460,1045),(460,1004)],dash=True)
end()

# E: multiple branches, reviewed tiers, AND gate and binary exits.
begin('E',1350,1640,'Two-stage filtering','Only keep / exclude are public; all review routes remain in the audit')
box('e_detector',25,145,1150,185,'Internal detector (distinct from final thresholds)',['With two sides + no protection: exclude if (S≥.64 AND ≥2 strong domains),','OR X≥.78, OR (D≥.82 AND holes≥.55 AND C≥.42).','Else review if S≥.38, any domain≥.45, any severe domain≥.78,','any metric≥.75 or protection; otherwise detector keep.'],size=21,source='_score_metrics: recommendations')
diamond('e_stage1',600,424,520,124,[f'Stage 1: compare S',f'K = {K:.3f}     E = {E:.3f}'],'apply_high_confidence_policy')
arrow([(600,330),(600,362)])
box('e_keep_candidate',25,555,300,105,'Keep candidate',[f'S ≤ {K:.3f}'],fill='green',size=23)
box('e_middle',405,555,350,105,'Internal review',[f'{K:.3f} < S < {E:.3f}'],fill='purple',size=23)
box('e_high',835,540,340,143,'Direct-exclude gate',['Detector exclude?','Two-sided? No protection?'],fill='purple',size=22)
arrow([(340,424),(175,424),(175,555)],'low',(219,498));arrow([(600,486),(600,555)],'middle',(670,523));arrow([(860,424),(1005,424),(1005,540)],'high',(1041,498))
arrow([(835,610),(790,610),(790,696),(600,696)],'any fails',(709,685))
diamond('e_review_split',600,789,500,125,[f'Review score < {B:.3f}?'],'_validate_review_tiers; _evaluate_review_tier_row')
arrow([(580,660),(580,696),(600,696),(600,726.5)])
box('e_low_tier',55,922,495,132,'Low-score review',[f'{K:.3f} < S < {B:.3f}',f'≥{L["min_corroborating_domains"]} domains; strongest ≥{L["severe_domain_threshold"]:.2f}'],fill='purple',size=23)
box('e_high_tier',650,922,495,132,'Higher-score review',[f'S ≥ {B:.3f}; includes downgraded highs',f'≥{H["min_corroborating_domains"]} domains; strongest ≥{H["severe_domain_threshold"]:.2f}'],fill='purple',size=22)
arrow([(350,789),(302,789),(302,922)],'yes',(310,885));arrow([(850,789),(897,789),(897,922)],'no',(910,885))
box('e_all_gates',140,1098,920,150,'AND: every required evidence condition passes',['Detector review/exclude + two-sided + no protection + confidence≥.90','Every available 3/5/7 window: same guards, domain count / severity,','and selected tier floor (.129 or .540); no window details → fail.'],fill='amber',size=21,source='_evaluate_review_tier_row: primary confidence only; all-window support=1.0; upper tier bound applies to primary only')
arrow([(302,1054),(302,1098)]);arrow([(897,1054),(897,1098)])
box('e_keep',125,1310,380,68,'KEEP',fill='green');box('e_exclude',695,1310,380,68,'EXCLUDE',fill='red')
arrow([(380,1248),(380,1310)],'any fails',(287,1289));arrow([(820,1248),(820,1310)],'all pass',(913,1289))
arrow([(175,660),(13,660),(13,1344),(125,1344)])
arrow([(1175,610),(1188,610),(1188,1344),(1075,1344)],'all pass',(1090,873))
note(1400,['Non-finite score / insufficient evidence → conservative keep + recorded diagnostics.'])
end()

# F: separate display lane, decision audit and ROI output.
begin('F',2640,1640,'Audit, display & exact ROI','Display transforms never feed back into QC scores or decisions')
box('f_audit',25,148,1150,140,'Save one audit row per original slice',['Raw score, detector call, stage-1 route, stage-2 tier, window checks, failed gates, final_call','Freeze policy + hashes; experimental_policy clears certification. Keep ≠ proven healthy.'],fill='teal',size=21,source='write_high_confidence_outputs; run_publish')
arrow([(600,288),(600,348)])
diamond('f_prereg_request',600,421,530,146,['Shared-coordinate','comparison requested?'],'run_preregister')
box('f_raw_only',30,575,355,129,'Raw-only view',['Keep original coordinates','No shared-anatomy claim'],fill='amber',size=22)
box('f_rigid',465,550,710,215,'Display-only rigid preregistration',['Separate middle-slice reference for each dataset','Seed13; ≤1000 points; 24 angles; 50 iterations; trim80%','Translation + rotation; no scale / mirror / deformation','Compose T(child→ref) = T(parent→ref) × H(child→parent)'],fill='blue',size=21,source='slice_preregistration.py; slice_quality_preregistration.py')
arrow([(335,421),(207,421),(207,575)],'no',(248,524));arrow([(865,421),(980,421),(980,550)],'yes',(1015,516))
diamond('f_fit_status',825,866,525,134,['Fit / reference chain','reliable enough to compare?'],'slice_preregistration.py: failure / ambiguity warnings propagate along composed chain')
arrow([(825,765),(825,799)])
box('f_warn',465,1002,335,128,'Uncertain / failed',['Save status + visible warning','Do not turn fit error into QC'],fill='amber',size=20)
box('f_fit_ok',840,1002,335,128,'Fitted coordinates',['Save matrices + inverses','Not anatomical ground truth'],fill='teal',size=20)
arrow([(562.5,866),(510,866),(510,970),(632,970),(632,1002)],'no',(616,957));arrow([(1087.5,866),(1147,866),(1147,970),(1007,970),(1007,1002)],'yes',(1020,957))
box('f_viewer',25,1180,1150,85,'Current viewer: KDE • capture • components • statistics • binary calls',[],fill='teal',source='slice_quality_visualization.py: render_binary_slice_quality_appendix; raw/display toggle; linked 3/5 slice windows')
arrow([(207,704),(207,1180)]);arrow([(632,1130),(632,1180)]);arrow([(1007,1130),(1007,1180)])
box('f_roi',25,1304,1150,82,'ROI: display rectangle → exact selected points → inverse raw polygon',[],fill='blue',source='slice_quality_detail_roi.js; export-roi: frame_id, source_row, obs_name, full-cache hashes; never use raw axis-aligned bounding box as exact inverse')
arrow([(600,1265),(600,1304)]);end()

# Reading links between panels; short connections avoid a single vertical trunk.
group('panel_routing')
for x,y in [(1305,865),(2595,865),(1305,2345),(2595,2345)]:
 arrow([(x-35,y),(x+35,y)],color='#3E8B94')
# C exits to a labeled continuation rather than crossing the entire figure.
badge(3760,1600,'D','purple');arrow([(3760,1570),(3760,1573)])
badge(125,1600,'C','purple');arrow([(125,1627),(125,1640)])
txt(3420,1608,'Continue in D',22,True,color='muted');txt(350,1608,'From C: normalized evidence',22,True,color='muted')
parts.append('</g>')
group('figure_notes')
txt(60,3103,'READING KEY',22,True,'start');txt(250,3103,'A → B → C → D → E → F; availability branches in A may operate in parallel.',23,False,'start')
txt(60,3142,'Frozen policy:',22,True,'start');txt(250,3142,f'K={K:.3f}; split={B:.3f}; E={E:.3f}; strongest domain low/high={L["severe_domain_threshold"]:.2f}/{H["severe_domain_threshold"]:.2f}. Local references can adapt; these policy values do not.',23,False,'start')
txt(60,3179,'Traceability:',22,True,'start');txt(250,3179,'Named SVG groups contain source-function descriptions. Full metric limits, score mappings and validation scope: companion methods document.',23,False,'start')
parts.append('</g></svg>')
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text('\n'.join(parts))
a.output.with_suffix('.objects.json').write_text(json.dumps(objects,indent=2))
print(a.output)
