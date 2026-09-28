"""A real failed-CI case must export gaps without fabricating finite bounds."""
import importlib.util,json,sys
from pathlib import Path
scripts=Path(__file__).resolve().parents[2]/'skills/spateo-4d-pipeline/scripts'
sys.path.insert(0,str(scripts))
from features_stage import finite_values

def test_nonfinite_fit_bounds_export_as_gaps():
    values=finite_values([1.,float('nan'),float('inf'),-float('inf'),0.])
    assert json.loads(json.dumps(values,allow_nan=False))==[1.,None,None,None,0.]
