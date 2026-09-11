"""Execute CLI against this checkout while using the installed Spateo runtime.
No installed files are changed. Usage: python activate_candidate.py scan ...
"""
import importlib.util
import os
from pathlib import Path
import runpy
import sys


def activate():
    import spateo
    import spateo.preprocessing as pp
    # Older server installations lack the current public preprocessing alias.
    if not hasattr(spateo, 'pp'):
        spateo.pp = pp
    bundled = Path(__file__).resolve().parents[1] / 'runtime'
    default_source = bundled if (bundled/'spateo/preprocessing/slice_quality.py').is_file() else Path(__file__).resolve().parents[3]
    root=Path(os.environ.get('SPATEO_REFEREE_SOURCE', default_source)).expanduser().resolve()
    if not (root/'spateo/preprocessing/slice_preregistration.py').is_file():
        raise FileNotFoundError('Set SPATEO_REFEREE_SOURCE to the updated Referee checkout, or run this wrapper from that checkout.')
    for short in ['slice_quality','slice_preregistration']:
        name='spateo.preprocessing.'+short
        spec=importlib.util.spec_from_file_location(name,root/'spateo/preprocessing'/f'{short}.py')
        mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod)
        setattr(pp,short,mod)
        for key in dir(mod):
            if getattr(getattr(mod,key), "__module__", "") == name:
                setattr(pp,key,getattr(mod,key))
                if hasattr(spateo,"pp"): setattr(spateo.pp,key,getattr(mod,key))
    return root

if __name__=='__main__':
    activate()
    runpy.run_path(str(Path(__file__).with_name('run_slice_quality_qc.py')),run_name='__main__')
