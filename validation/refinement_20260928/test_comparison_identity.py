"""Registration comparison must reject a same-size, wrong-coordinate specimen."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'skills/spateo-4d-pipeline/scripts'))
from viewer_payload import attach_comparison

def test_comparison_identity_and_coordinate_guard(tmp_path):
    xyz=np.array([[0.,0.,0.],[1.,2.,3.],[2.,4.,6.]])
    payload={'alignment':[{'ids':['a','b','c'],'before':xyz.tolist(),'groups':['CNS']*3} for _ in range(2)]}
    arrays={}
    for k in (1,2):
        arrays.update({f'stage{k}_ids':np.array(['c','a','b']),f'stage{k}_before':xyz[[2,0,1]],
          f'stage{k}_benchmark_rigid':xyz[[2,0,1]]+k,f'stage{k}_labels':np.array(['CNS']*3)})
    np.savez(tmp_path/'coordinates.npz',**arrays)
    record={'status':'completed','coordinates_sha256':hashlib.sha256((tmp_path/'coordinates.npz').read_bytes()).hexdigest(),
            'inputs':[],'parameters':{},'spateo_root':'test'}
    (tmp_path/'comparison.json').write_text(json.dumps(record))
    attach_comparison(payload,'permuted rows',tmp_path)
    assert payload['registrationComparisons'][0]['coordinates'][0]==(xyz+1).tolist()
    payload['alignment'][0]['before'][0][0]=20
    with pytest.raises(AssertionError):attach_comparison(payload,'wrong specimen',tmp_path)
