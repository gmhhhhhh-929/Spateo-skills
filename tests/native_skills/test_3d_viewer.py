import importlib.util
import json
from pathlib import Path

import numpy as np
import pyvista as pv
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'skills/spateo-3d-pipeline/subskills/spateo-render-3d-viewer/scripts/build_viewer.py'
spec = importlib.util.spec_from_file_location('mesh_review_viewer', SCRIPT)
viewer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(viewer)


@pytest.fixture
def fixture_config(tmp_path):
    sphere = pv.Sphere(radius=1, theta_resolution=16, phi_resolution=16).triangulate()
    sphere.save(tmp_path / 'tissue.vtk')
    pc = pv.PolyData(np.array([[0.,0,0],[.1,0,0],[3,0,0],[0,.2,0]]))
    pc['obs_index'] = np.array(['a','b','c','d'])
    pc['tissue'] = np.array(['CNS','CNS','CNS','Gut'])
    pc.save(tmp_path / 'pc.vtk')
    config = {'datasets': [{'name': 'Stage 1', 'point_cloud': 'pc.vtk', 'label_key': 'tissue',
                           'point_budget': 2, 'layers': [{'name': 'CNS', 'mesh': 'tissue.vtk', 'values': ['CNS']}]}]}
    return tmp_path, config


def test_source_geometry_ids_and_unsampled_morphology(fixture_config):
    base, config = fixture_config
    before = {p: viewer.digest(p) for p in base.glob('*.vtk')}
    payload, audits = viewer.build_payload(config, base, {})
    layer = payload['datasets'][0]['layers'][0]
    original = pv.read(base / 'tissue.vtk')
    assert np.array_equal(layer['mesh']['points'], original.points)
    assert np.array_equal(layer['mesh']['faces'], original.faces.reshape(-1,4)[:,1:])
    assert layer['cells']['ids'] == ['a','c']  # Includes the outside source cell, not clipped.
    assert layer['features']['Cells_inside'] == 2
    assert layer['features']['cell_density'] == pytest.approx(2 / original.volume)
    assert audits[0]['layers'][0]['source_cells'] == 3
    assert all(viewer.digest(p) == h for p,h in before.items())


def test_native_manifest_relative_paths_and_body_selection(fixture_config):
    base, config = fixture_config
    manifest = {'meshes': [{'name':'Body','role':'body','vtk':'tissue.vtk','color':'#BFC7CD','selected_values':[]}]}
    (base / 'mesh.json').write_text(json.dumps(manifest))
    dataset = config['datasets'][0];dataset.pop('layers');dataset['mesh_manifest']='mesh.json'
    payload, _ = viewer.build_payload(config, base, {})
    body = payload['datasets'][0]['layers'][0]
    assert body['opacity'] == .15
    assert body['cells']['source_count'] == 4


@pytest.mark.parametrize('change', ['missing_label','missing_value','duplicate_ids','nonfinite','budget','duplicate_name'])
def test_invalid_inputs_fail(fixture_config, change):
    base, config = fixture_config
    dataset = config['datasets'][0]
    if change == 'missing_label': dataset['label_key']='absent'
    elif change == 'missing_value': dataset['layers'][0]['values']=['absent']
    elif change == 'budget': dataset['point_budget']=0
    elif change == 'duplicate_name': dataset['layers'] *= 2
    else:
        pc = pv.read(base/'pc.vtk')
        if change == 'duplicate_ids': pc['obs_index']=np.array(['a']*pc.n_points)
        else: pc.points[0,0]=np.nan
        pc.save(base/'pc.vtk')
    with pytest.raises(ValueError):viewer.build_payload(config, base, {})


def test_open_surface_suppresses_volume(fixture_config):
    base, config = fixture_config
    pv.Plane().triangulate().save(base/'tissue.vtk')
    payload, audits = viewer.build_payload(config, base, {})
    layer=payload['datasets'][0]['layers'][0]
    assert layer['features']['Volume'] is None
    assert layer['features']['V/SA_ratio'] is None
    assert 'cell_density' not in layer['features']
    assert audits[0]['layers'][0]['warnings']


def test_html_offline_escape_and_no_overwrite(fixture_config):
    base, config = fixture_config
    malicious = '</script><script>alert(1)</script>'
    config['datasets'][0]['layers'][0]['name']=malicious
    path=base/'config.json';path.write_text(json.dumps(config))
    output=base/'review';report=viewer.build(path,output)
    html=(output/'index.html').read_text()
    embedded=html.split('<script type="application/json" id="model-data">')[1].split('</script>')[0]
    assert malicious not in embedded
    assert json.loads(embedded)['datasets'][0]['layers'][0]['name']==malicious
    assert '<script src=' not in html
    assert '/*__PLOTLY__*/' not in html
    assert report['source_files_unchanged']
    assert report['html_sha256']==viewer.digest(output/'index.html')
    assert str(base) not in embedded
    with pytest.raises(FileExistsError):viewer.build(path,output)


def test_point_only_and_mesh_only(fixture_config):
    base, config=fixture_config
    config['datasets'][0]['layers']=[{'name':'Cells','cells':'pc.vtk'},{'name':'Surface','mesh':'tissue.vtk'}]
    payload,_=viewer.build_payload(config,base,{})
    points,mesh=payload['datasets'][0]['layers']
    assert points['mesh'] is None and points['cells'] is not None
    assert mesh['cells'] is None and 'cell_density' not in mesh['features']
