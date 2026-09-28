"""Verify exported paths, IDs and comparison provenance independently of the UI."""
import hashlib,json,sys,re
from pathlib import Path
import numpy as np
html=Path(sys.argv[1]);output=Path(sys.argv[2]);text=html.read_text()
payload=json.loads(text.rsplit('const D=',1)[1].split(',$=id=>',1)[0])
assert payload['capabilities']==['alignment','mapping','field','trajectory','features']
assert len(payload['registrationComparisons'])==2
for k,c in enumerate(payload['alignment']):
    assert len(c['ids'])==len(set(c['ids']))==c['total']
    for v in payload['registrationComparisons']:
        assert len(v['coordinates'][k])==len(c['ids'])
        assert v['qc']['xyz']['dorsoventral_z']['opposite_signs']
assert not payload['registrationReview']['xyz']['dorsoventral_z']['opposite_signs']
a=payload['analysis'][0];ids={cell:i for i,cell in enumerate(a['ids'])}
for cell,path in zip(payload['trajectories']['ids'],payload['trajectories']['paths']):
    np.testing.assert_allclose(path[0],a['xyz'][ids[cell]],rtol=0,atol=0)
assert all(len(v)==len(a['ids']) for v in a['features'].values())
assert not re.search(r'<script[^>]+src=[\"\']https?://',text,re.I)
result={'passed':True,'html_sha256':hashlib.sha256(html.read_bytes()).hexdigest(),
        'full_cells':[c['total'] for c in payload['alignment']],
        'displayed_paths':len(payload['trajectories']['ids']),
        'comparison_count':len(payload['registrationComparisons']),
        'checks':['full unique cell identities','historical orientation flags','corrected orientation flag',
                  'every displayed trajectory starts at its named seed','feature lengths match source IDs',
                  'no Plotly CDN dependency']}
output.write_text(json.dumps(result,indent=2));print(result)
