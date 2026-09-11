"""Maintain one end-to-end SVG; thresholds are read from the packaged policy."""
from pathlib import Path
from xml.sax.saxutils import escape
import json

ROOT=Path(__file__).resolve().parents[1]
p=json.loads((ROOT/'policies/joint_review_v2.json').read_text())['policy']
low,high=p['review_exclusion_tiers'];K=p['keep_max_score'];E=p['exclude_min_score'];B=low['max_score']
W,H=2400,6050
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
'<defs><marker id="arrow" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto"><path d="M0,0 L10,4 L0,8z" fill="#415165"/></marker></defs>',
'<rect width="100%" height="100%" fill="white"/>']
colors={'data':'#e7eef6','metric':'#e8edf9','decision':'#eae3f5','keep':'#dff1e7','exclude':'#f7e0e2','display':'#e1f1f2','warning':'#fff0cd'}

def text(x,y,value,size=26,weight=400,anchor='middle',color='#193047'):
 svg.append(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="Arial,sans-serif" font-size="{size}" font-weight="{weight}" fill="{color}">{escape(value)}</text>')

def box(x,y,w,h,title,lines=(),kind='data'):
 svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="14" fill="{colors[kind]}" stroke="#556579" stroke-width="2"/>')
 text(x+w/2,y+40,title,29,700)
 for i,line in enumerate(lines):text(x+w/2,y+80+32*i,line,24)

def diamond(cx,cy,w,h,lines):
 svg.append(f'<path d="M{cx} {cy-h/2} L{cx+w/2} {cy} L{cx} {cy+h/2} L{cx-w/2} {cy} Z" fill="{colors["decision"]}" stroke="#556579" stroke-width="2"/>')
 for i,line in enumerate(lines):text(cx,cy+(i-(len(lines)-1)/2)*29+8,line,24,700)

def edge(points,label='',label_xy=None):
 svg.append('<polyline points="'+' '.join(f'{x},{y}' for x,y in points)+'" fill="none" stroke="#415165" stroke-width="3" marker-end="url(#arrow)"/>')
 if label:
  x,y=label_xy or points[len(points)//2];text(x,y,label,24,700,color='#76531a')

def heading(y,title):
 text(50,y,title,32,700,'start');svg.append(f'<line x1="50" y1="{y+16}" x2="2350" y2="{y+16}" stroke="#c7d0db"/>')

text(1200,65,'Spateo Referee · complete QC and viewer workflow',44,700)
text(1200,112,'Joint review v2 · metric-stress tested · independent biological validation incomplete',27,color='#95681e')
heading(160,'1  Input contracts and observable evidence')
box(920,205,560,150,'Ordered biological series',['H5AD / one file per slice / dataset manifest','Never join specimen boundaries'])
box(920,415,560,150,'Resolve and verify input',['Slice ID + order + raw x/y + units','Expression layer + source / identity hashes'])
box(1820,415,500,150,'Invalid input? Stop affected run',['Missing identity / coordinate contracts','Do not invent slices, keys or values'],'warning')
box(70,405,730,180,'Expression semantics branch',['Count-like matrix: compute captured counts / genes','Normalized or one-hot: mark capture unavailable','Missing fields do not become zero evidence'],'warning')
edge([(1200,355),(1200,415)]);edge([(1480,490),(1820,490)],'invalid',(1650,476));edge([(920,490),(800,490)])
metrics=[('Density / amount',['Locations, area, convex-hull density','NN-tail ratio, holes, components','Observed tissue amount ≠ quality alone']),('Expression capture',['Counts, genes, zero fraction, complexity','Local low-capture fraction and clustering','Selected matrix overrides stale obs summaries']),('Damage proxies',['Fragmentation + mitochondrial fraction','Largest connected-tissue fraction','Not a specific experimental diagnosis']),('Cross-slice continuity',['Expression profiles, type composition, holes','Requires the corresponding measurements','Coordinate display fitting is not input evidence'])]
for x,(title,lines) in zip([40,635,1230,1825],metrics):
 box(x,700,535,225,title,lines,'metric');edge([(1200,565),(1200,635),(x+267,635),(x+267,700)])
for x in [307,902,1497,2092]:edge([(x,925),(x,985),(1200,985),(1200,1100)])
heading(1055,'2  Local expectations, anomaly scores and protection')
box(540,1100,1320,210,'Local expectation and metric anomalies',[
 'Primary window: default 3; optional auto compares 3 / 5 / 7; support trend up to 9',
 'Leave focal slice out; two-sided local trend, otherwise one-sided reference',
 'Directional effect size + robust MAD residual; metric-specific limits: methods table'], 'metric')
box(540,1370,1320,245,'Aggregate observed evidence',[
 'Four domains: weighted available metrics and top-two corroborating signals',
 'Weighted sum: density 0.31 + expression 0.39 + damage 0.12 + continuity 0.18',
 'Score = max(weighted sum, clip(0.48 + 0.78 × (maximum domain − 0.68), 0, 1))',
 'Apply anatomical protection cap; one-sided endpoint score × 0.82; retain warnings'], 'metric')
edge([(1200,1310),(1200,1370)])
box(540,1670,1320,300,'Internal detector branches (not final policy calls)',[
 'Exclude requires two-sided + no protection, and at least one of:',
 '(S ≥0.64 AND ≥2 domains at 0.45) OR expression ≥0.78',
 'OR (density ≥0.82 AND holes ≥0.55 AND continuity ≥0.42)',
 'Otherwise review if S ≥0.38 OR any domain ≥0.45 OR any metric ≥0.75',
 'OR a severe domain / protection flag; otherwise detector keep'], 'decision')
edge([(1200,1615),(1200,1670)])
svg.append('<g transform="translate(0,370)">')
heading(1670,'3  Frozen policy scope and stage-1 routing')
box(540,1720,1320,155,'Check policy and scope',[
 'Bundled joint policy: experimental_policy; require passed frozen metric-stress evaluation',
 'No raw-matrix / biological certification is inferred; retain policy ID and hash'], 'warning')
edge([(1200,1600),(1200,1720)])
diamond(1200,2015,720,200,[f'Compare finite score S',f'K = {K:.3f}; E = {E:.3f}'])
edge([(1200,1875),(1200,1915)])
box(70,2180,535,175,'Stage-1 keep candidate',[
 'Operational keep under complete binary mode',
 'Retain detector / guard agreement in audit',
 'Keep is not proof of defect-free tissue'],'keep')
box(1780,2180,550,190,'Direct exclusion gate',[
 'Detector = exclude AND two-sided context',
 'AND no anatomical protection',
 'Every condition must pass'],'decision')
edge([(840,2015),(337,2015),(337,2180)],f'S ≤ {K:.3f}',(580,1995));edge([(1560,2015),(2055,2015),(2055,2180)],f'S ≥ {E:.3f}',(1820,1995))
box(870,2450,660,145,'Internal review queue',[
 'Middle scores + unsafe high-score candidates',
 'Review is not a third final viewer state'],'decision')
edge([(1200,2115),(1200,2450)],f'{K:.3f} < S < {E:.3f}',(1420,2280));edge([(1780,2300),(1600,2300),(1600,2520),(1530,2520)],'gate fails',(1640,2270))
diamond(1200,2725,640,170,[f'S < {B:.3f}?','Select one evidence band'])
edge([(1200,2595),(1200,2640)])
box(590,2910,560,175,'Low-score evidence band',[
 f'{K:.3f} < S < {B:.3f}',
 f'Domains ≥ {low["min_corroborating_domains"]}; strongest domain ≥ {low["severe_domain_threshold"]:.2f}',
 'Higher exclusion-evidence threshold'],'decision')
box(1260,2910,560,175,'Higher-score evidence band',[
 f'S ≥ {B:.3f}, only if still in review',
 f'Domains ≥ {high["min_corroborating_domains"]}; strongest domain ≥ {high["severe_domain_threshold"]:.2f}',
 'Includes high-score gate downgrades'],'decision')
edge([(880,2725),(870,2725),(870,2910)],'yes',(810,2835));edge([(1520,2725),(1540,2725),(1540,2910)],'no',(1600,2835))
box(540,3170,1320,210,'Check every evidence condition',[
 'Primary: detector review/exclude; two-sided; no protection; confidence field ≥0.90',
 'Supporting domains have score ≥0.45; apply the selected maximum-domain requirement',
 'Every available 3/5/7 window: same detector, guards, domains, severity and band floor',
 'No window evidence ⇒ fail. Required support =100%. No keep-only shortcut.'],'decision')
edge([(870,3085),(870,3120),(1200,3120),(1200,3170)]);edge([(1540,3085),(1540,3120),(1200,3120)])
diamond(1200,3480,620,150,['All selected-band','conditions pass?'])
edge([(1200,3380),(1200,3405)])
box(70,3550,535,135,'KEEP',['Save failed conditions / conservative retention'],'keep')
box(1795,3550,535,135,'EXCLUDE',['Save direct or review evidence route'],'exclude')
edge([(890,3480),(700,3480),(700,3615),(605,3615)],'no',(780,3465));edge([(1510,3480),(1660,3480),(1660,3615),(1795,3615)],'yes',(1600,3465))
edge([(337,2355),(337,3550)]);edge([(2330,2275),(2360,2275),(2360,3520),(2060,3520),(2060,3550)],'direct gate passes',(2140,3410))
box(540,3790,1320,170,'Write complete decision audit',[
 'One row per original slice: score, stage-1 triage, tier, all window checks and final_call',
 'Policy / input hashes + experimental scope; public CSV and HTML contain keep/exclude',
 'Non-finite score or insufficient evidence: conservative keep + explicit diagnostics'],'data')
edge([(337,3685),(337,3875),(540,3875)]);edge([(2060,3685),(2060,3875),(1860,3875)])
heading(4030,'4  Optional display preregistration and exact ROI; no feedback into QC')
diamond(1200,4190,680,170,['Shared-coordinate ROI','comparison requested?'])
edge([(1200,3960),(1200,4105)])
box(60,4360,1110,245,'Fit display-only rigid transforms',[
 'Original x/y kept; middle slice anchors each separate dataset',
 'Seed 13; ≤1000 fit points; 24 angles; ≤50 iterations; closest 80% matches',
 'Optional observed-label ranking (weight 0.5); no reference coordinates for fitting',
 'No scale, reflection, nonrigid deformation or smoothing',
 'Compose T(child→reference) = T(parent→reference) × H(child→parent)'],'display')
box(1550,4360,780,180,'Raw-coordinate display only',[
 'Preserve available evidence views',
 'Do not imply a common anatomical frame',
 'Shared-frame ROI requires a saved transform'],'warning')
edge([(860,4190),(615,4190),(615,4360)],'yes',(720,4170));edge([(1540,4190),(1940,4190),(1940,4360)],'no',(1740,4170))
box(60,4660,1110,180,'Persist transforms and point provenance',[
 'Reference, pair edges, matrices, inverses, seeds, parameters and source hashes',
 'Full raw/display caches retain source row and obs identity',
 'Failed / ambiguous edges and blocked chains remain visible'],'display')
edge([(615,4605),(615,4660)])
box(1410,4660,920,180,'Fit failure / ambiguity / drift branch',[
 'Store failure status; warn on affected slices and downstream chains',
 'Do not interpret misalignment as a low-quality slice',
 'Equal display range does not ensure equal anatomy'],'warning')
edge([(1170,4750),(1410,4750)],'warning',(1290,4730));edge([(1940,4540),(2375,4540),(2375,4910),(1600,4910),(1600,4950)])
box(540,4950,1320,220,'Render the current detailed viewer + collection',[
 'Overview · Statistics · Slice evidence · All slices · Methods',
 'KDE / expression capture / connected components in linked 3- or 5-slice windows',
 'Raw/display switch; common extents; QC call, band requirements and fit status',
 'Display samples and unavailable measures disclosed; no point colors alter decisions'],'display')
edge([(615,4840),(615,4890),(1200,4890),(1200,4950)]);edge([(1870,4840),(1870,4890),(1200,4890)])
box(540,5240,1320,215,'Select and export a traceable ROI',[
 'Drag on focal slice; test transformed points against one display-frame rectangle',
 'Show corresponding neighbors; export frame_id, slice IDs and selected source identities',
 'Inverse-transform rectangle corners to a raw polygon, not an axis-aligned bounding box',
 'Full-point replay verifies frame/cache hashes; original H5AD is never overwritten'],'display')
edge([(1200,5170),(1200,5240)])
text(1200,5530,'Implementation and every metric limit: methods.md · policy values: ../policies/joint_review_v2.json',25)
text(1200,5570,'Editable source: build_workflow.py · Generated from actual packaged policy · No biological accuracy claim',24)
svg.append('</g></svg>')
(ROOT/'references/workflow.svg').write_text('\n'.join(svg))
