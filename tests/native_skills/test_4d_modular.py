"""Behavioral tests for independent entry, VTK IDs, full branches and honest viewer routes."""
import copy
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pytest
import pipeline_runtime as rt
from alignment_stage import attach_pointcloud
from viewer_payload import build_payload
from test_4d import inputs, completed


def test_alignment_only_has_only_registration_capability(inputs):
    root,path,c=inputs
    r=rt.execute(path,root/'routes',run_id='alignment_only',until='alignment')
    m=rt.read_json(r['manifest']);payload=rt.read_json(m['stages']['dashboard']['outputs']['viewer_payload']['path'])
    assert payload['capabilities']==['alignment']
    assert payload['analysis'] is None and payload['trajectories'] is None
    assert all(m['stages'][s]['status']=='skipped' for s in ['mapping','morphofield','trajectory','metrics','gp'])
    assert payload['alignment'][0]['total']==40


def test_vtk_import_matches_ids_not_row_order(inputs,completed,tmp_path):
    import pyvista as pv
    root,_,c=inputs
    r,m=completed
    a=ad.read_h5ad(root/'stage1.h5ad')
    pc=pv.PolyData(a.obsm['spatial_3d'][::-1]+[1,2,3])
    pc.point_data['obs_index']=a.obs_names[::-1].to_numpy(dtype=str)
    p=tmp_path/'permuted.vtk';pc.save(p)
    attach_pointcloud(a,str(p),'spatial_3d')
    np.testing.assert_allclose(a.obsm['spatial_3d'],a.obsm['spatial_3d_before_vtk']+[1,2,3])
    bad_ids=pc.point_data['obs_index'].copy().astype(object);bad_ids[0]='wrong-id'
    pc.point_data['obs_index']=np.asarray(bad_ids,dtype=str);pc.save(p)
    with pytest.raises(ValueError,match='exact match'):attach_pointcloud(a,str(p),'spatial_3d')


def test_import_aligned_pair_does_not_register(inputs,completed):
    root,_,c=inputs
    _,m=completed;c=copy.deepcopy(c)
    pair=m['stages']['alignment']['outputs']
    c['inputs'].update(stage1=pair['stage1_h5ad']['path'],stage2=pair['stage2_h5ad']['path'])
    c['workflow'].update(entry='aligned',frame_id='base',until='trajectory')
    p=root/'from_aligned.json';rt.write_json(p,c)
    r=rt.execute(p,root/'routes',run_id='from_aligned');n=rt.read_json(r['manifest'])
    imported,_=rt.load_pair(n['stages']['alignment']['outputs']);original,_=rt.load_pair(pair)
    np.testing.assert_array_equal(imported.obsm['aligned'],original.obsm['aligned'])
    payload=rt.read_json(n['stages']['dashboard']['outputs']['viewer_payload']['path'])
    assert payload['capabilities']==['alignment','mapping','field','trajectory']
    assert payload['trajectories']['total']==40
    assert np.asarray(payload['trajectories']['paths']).shape==(40,6,3)


def test_import_native_field_features_and_jacobian(inputs,completed):
    root,_,c=inputs
    _,m=completed;c=copy.deepcopy(c)
    pair=m['stages']['morphofield']['outputs']
    c['inputs'].update(stage1=pair['stage1_h5ad']['path'],stage2=pair['stage2_h5ad']['path'])
    c['workflow'].update(entry='field',frame_id='base')
    c['metrics']['selected'] += ['jacobian_xy','velocity_x']
    c['metrics']['glm_metrics'] += ['jacobian_frobenius']
    c['metrics']['glm_genes']='*'
    c['gp']['enabled']=False
    p=root/'from_field.json';rt.write_json(p,c)
    r=rt.execute(p,root/'routes',run_id='from_field');n=rt.read_json(r['manifest'])
    assert n['stages']['mapping']['status']=='imported'
    assert n['stages']['morphofield']['status']=='imported'
    a,_=rt.load_pair(n['stages']['metrics']['outputs'])
    J=a.uns['jacobian'];assert J.shape==(3,3,40)
    np.testing.assert_allclose(a.obs['jacobian_xy'],J[0,1])
    summary=rt.read_json(n['stages']['metrics']['outputs']['glm_summary']['path'])
    assert summary['jacobian_frobenius']['tested']==12
    payload=rt.read_json(n['stages']['dashboard']['outputs']['viewer_payload']['path'])
    assert 'features' in payload['capabilities']
    assert 'jacobian_frobenius' in payload['analysis'][0]['features']
    assert payload['trajectories']['ids']==a.obs_names.tolist()
    for path,xyz in zip(payload['trajectories']['paths'],a.obsm['aligned']):np.testing.assert_allclose(path[0],xyz)


def test_new_config_rejects_ignored_backend_keyword(inputs,tmp_path):
    _,_,c=inputs;c=copy.deepcopy(c);c['morphofield']['MaxIter']=100
    p=tmp_path/'bad.json';rt.write_json(p,c)
    with pytest.raises(ValueError,match='Unknown morphofield'):rt.canonical_config(p)


def test_missing_frame_contract_rejected(inputs,tmp_path):
    _,_,c=inputs;c=copy.deepcopy(c);c['workflow']['entry']='aligned'
    p=tmp_path/'bad.json';rt.write_json(p,c)
    with pytest.raises(ValueError,match='frame_id'):rt.canonical_config(p)


def test_imported_frame_mismatch_fails(inputs,completed):
    root,_,c=inputs
    _,m=completed;c=copy.deepcopy(c)
    pair=m['stages']['alignment']['outputs']
    c['inputs'].update(stage1=pair['stage1_h5ad']['path'],stage2=pair['stage2_h5ad']['path'])
    c['workflow'].update(entry='aligned',frame_id='different-recorded-frame')
    path=root/'wrong_frame.json';rt.write_json(path,c)
    with pytest.raises(ValueError,match='frame metadata'):
        rt.execute(path,root/'routes',run_id='wrong_frame')


def test_viewer_template_change_rebuilds_only_display(inputs,completed):
    _,_,c=inputs
    _,m=completed
    runtime=copy.deepcopy(m['implementation']);runtime['viewer_asset_sha256']='new-template'
    plan=rt.plan(c,m,input_hashes=m['inputs'],runtime=runtime)
    assert plan['recompute_stages']==['dashboard']
