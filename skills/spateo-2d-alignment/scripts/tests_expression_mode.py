#!/usr/bin/env python3
"""Synthetic no-annotation expression/profile contracts with a fake Spateo backend."""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import types

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SKILL=Path(__file__).resolve().parents[1]
PIPELINE=SKILL/'pipelines/continuity-guided'
sys.path.insert(0,str(PIPELINE))
import engine
import _core as core
import validate_inputs as validator

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def parse(directory,output,*extra,profile='generalized'):
    previous=sys.argv;sys.argv=['engine.py','--stage','synthetic','--slice-dir',str(directory),
        '--output-dir',str(output),'--representation','expression-pca','--device','cpu',
        *(['--profile',profile] if profile is not None else []),*extra]
    try:return engine.parse_args()
    finally:sys.argv=previous

def sync_hashes(parent):
    path=parent/'expression_pca_manifest.json';m=json.loads(path.read_text())
    for entry in m['slices']:
        p=parent/entry['path'];obj=ad.read_h5ad(p)
        entry.update(output_h5ad_sha256=validator.sha(p),features_sha256=validator.feature_sha(obj.obsm['X_pca']),
                     cell_ids_sha256=validator.cell_ids_sha(obj.obs_names),n_cells=obj.n_obs)
    path.write_text(json.dumps(m,indent=2))

def expect_error(call,text):
    try:call()
    except (ValueError,RuntimeError) as error:
        assert text in str(error),(text,str(error));return
    raise AssertionError('Expected rejection: '+text)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report',type=Path);options=p.parse_args()
    checks=[];with_calls=[]
    with tempfile.TemporaryDirectory(prefix='expression-no-annotation-') as temp:
        root=Path(temp);inputs=root/'blind';inputs.mkdir();rng=np.random.default_rng(938)
        counts=rng.poisson(3.,size=(160,40)).astype(np.float32)
        all_ids=np.asarray([f'cell_{i}' for i in range(160)])
        # Expression source includes forbidden reference keys; the preparation
        # script must access only allowed expression/identity datasets.
        source=ad.AnnData(X=sparse.csr_matrix(counts),obs=pd.DataFrame(index=all_ids),
                         var=pd.DataFrame(index=[f'g{i}' for i in range(40)]))
        source.obsm['spatial_3d']=rng.normal(size=(160,3))
        source_path=root/'expression_source.h5ad';source.write_h5ad(source_path)
        offset=0
        for i,n in enumerate([40,40,80]):
            obj=ad.AnnData(X=sparse.csr_matrix((n,1),dtype=np.float32),obs=pd.DataFrame({'physical_z':i*17.},index=all_ids[offset:offset+n]))
            obj.obsm['spatial']=rng.normal(size=(n,2)).astype(np.float32)+[i*.2,0]
            obj.write_h5ad(inputs/f'CS01_SL{i+1:03d}_YSYNTH.Spatial.h5ad');offset+=n
        prep=module('expression_preparation',SKILL/'scripts/prepare_expression_pca.py')
        prepared=root/'prepared'
        with contextlib.redirect_stdout(io.StringIO()):
            prep.prepare(argparse.Namespace(slice_dir=inputs,expression_source=source_path,expression_layer='X',
                matrix_state='counts',output_dir=prepared,n_hvg=30,n_pcs=5,seed=938,annotation_key=None,physical_z_key='physical_z'))
        slices_dir=prepared/'slice_h5ad';summary,_=validator.validate_inputs(slices_dir,'expression-pca',filename_style='continuity')
        assert summary['annotation_qc']=='off' and summary['shared_expression_basis_provenance']['status']=='verified'
        checks.append('real joint-PCA preparer and validator agree on source, shared basis and per-file hashes')
        for path in slices_dir.glob('*.h5ad'):assert 'anno' not in ad.read_h5ad(path).obs
        original_hashes={str(path):validator.sha(path) for path in slices_dir.glob('*.h5ad')}
        expected={}
        for path in slices_dir.glob('*.h5ad'):
            obj=ad.read_h5ad(path)
            expected.update({str(cid):row.copy() for cid,row in zip(obj.obs_names,obj.obsm['X_pca'])})
        args=parse(slices_dir,root/'run','--pca-provenance',str(prepared/'expression_pca_manifest.json'))
        slices=core.load_slices(args);models=engine.load_models(slices,args)
        assert all(set(s.annotations)=={core.QC_UNLABELED} for s in slices)
        assert all(set(m.obs['anno'])=={core.QC_UNLABELED} for m in models)
        assert len(core.QC_UNLABELED)>1 and all(s.annotations.dtype.itemsize>=4*len(core.QC_UNLABELED) for s in slices)
        checks.append('no annotation is required; full non-biological sentinel exists only in memory')
        expect_error(lambda:validator.validate_inputs(slices_dir,'expression-pca',filename_style='continuity',annotation_qc='provided'),'requires obs')
        for flag in ['--repair-policy','--priority-annotation']:
            value='learned_probe' if flag=='--repair-policy' else 'ROD'
            with contextlib.redirect_stderr(io.StringIO()):
                try:parse(slices_dir,root/'invalid',flag,value)
                except SystemExit as error:assert error.code==2
                else:raise AssertionError('Generalized accepted legacy bypass/priority')
        checks.append('provided QC missing annotations fails; generalized disallows learned_probe and fixed biological priority')
        fake_torch=types.ModuleType('torch');fake_torch.manual_seed=lambda _:None;fake_torch.cuda=types.SimpleNamespace(is_available=lambda:False)
        fake_spateo=types.ModuleType('spateo')
        def fake_align(*,models,**kwargs):
            for m in models:
                for cid,row in zip(m.obs_names,m.obsm['X_pca']):np.testing.assert_array_equal(row,expected[str(cid)])
            with_calls.append({'nn_init':kwargs['nn_init'],'init_transform':kwargs['init_transform'],'rep_field':kwargs['rep_field']})
            output=[m.copy() for m in models]
            for m in output:m.obsm['align_spatial']=np.asarray(m.obsm['spatial']).copy()
            return output,[np.zeros((left.n_obs,right.n_obs)) for left,right in zip(models[:-1],models[1:])]
        fake_spateo.align=types.SimpleNamespace(morpho_align=fake_align)
        previous={name:sys.modules.get(name) for name in ['spateo','torch']};sys.modules.update(spateo=fake_spateo,torch=fake_torch)
        try:
            with contextlib.redirect_stdout(io.StringIO()):engine.main(args)
            qc=json.loads((args.output_dir/'blind_qc_summary.json').read_text())
            frozen=json.loads((args.output_dir/'frozen_manifest.json').read_text())
            assert with_calls[0]['nn_init'] is False and with_calls[0]['rep_field']=='obsm'
            assert 'no biological annotation' in qc['blind_contract']['annotation_source']
            assert any('annotation_qc_disabled' in r.get('reason','') for r in qc['operations'])
            assert frozen['initial_nn_init'] is False and frozen['anonymous_policy_sha256']==engine.ANONYMOUS_POLICY_SHA256
            assert frozen['output']['n_cells']==160
            assert frozen['parameters']['pca_provenance']==str(prepared/'expression_pca_manifest.json')
            coordinates=pd.read_csv(args.output_dir/'aligned_coordinates.csv.gz',keep_default_na=False)
            assert set(coordinates['anno'])=={''} and set(coordinates['cell_id'])==set(all_ids)
            checks.append('no-annotation expression end-to-end reaches initial and geometry Spateo calls unchanged; terminal label repair skips with reason')
            before=len(with_calls)
            legacy=parse(slices_dir,root/'legacy','--postprocess','none',profile=None)
            with contextlib.redirect_stdout(io.StringIO()):engine.main(legacy)
            assert with_calls[before]['nn_init'] is True and legacy.priority_annotation=='ROD'
            checks.append('omitted profile preserves legacy default, nn_init=True and legacy priority')
            for representation in ['spatial-only','annotation-onehot']:
                derived=root/representation;shutil.copytree(prepared,derived)
                for path in (derived/'slice_h5ad').glob('*.h5ad'):
                    obj=ad.read_h5ad(path)
                    if representation=='spatial-only':
                        obj.obsm['X_pca']=np.ones((obj.n_obs,30),dtype=np.float32)
                    else:
                        obj.obs['anno']=np.where(np.arange(obj.n_obs)%2,'class_A','class_B')
                        obj.obsm['X_pca']=np.eye(2,dtype=np.float32)[np.arange(obj.n_obs)%2]
                    obj.write_h5ad(path)
                    expected.update({str(cid):row.copy() for cid,row in zip(obj.obs_names,obj.obsm['X_pca'])})
                mode_args=parse(derived/'slice_h5ad',root/(representation+'_run'),'--representation',representation,'--postprocess','none')
                before=len(with_calls)
                with contextlib.redirect_stdout(io.StringIO()):engine.main(mode_args)
                mode_qc=json.loads((mode_args.output_dir/'blind_qc_summary.json').read_text())
                assert with_calls[before]['nn_init'] is False
                assert mode_qc['parameters']['annotation_qc']==('provided' if representation=='annotation-onehot' else 'off')
                assert mode_qc['n_cells']==160
            checks.append('spatial-only without annotation and existing annotation-onehot mode both complete with unchanged features and correct QC defaults')
        finally:
            for name,value in previous.items():
                if value is None:sys.modules.pop(name,None)
                else:sys.modules[name]=value
        assert all(validator.sha(path)==digest for path,digest in original_hashes.items())
        checks.append('source prepared files unchanged after both end-to-end profiles')
        annotated=root/'annotated';shutil.copytree(prepared,annotated)
        for path in (annotated/'slice_h5ad').glob('*.h5ad'):
            obj=ad.read_h5ad(path);obj.obs['anno']=np.where(np.arange(obj.n_obs)%2,'low_quality_A','low_quality_B');obj.write_h5ad(path)
        sync_hashes(annotated)
        off=parse(annotated/'slice_h5ad',root/'unused');assert all(set(s.annotations)=={core.QC_UNLABELED} for s in core.load_slices(off))
        provided=parse(annotated/'slice_h5ad',root/'unused2','--annotation-qc','provided')
        validator.validate_inputs(provided.slice_dir,'expression-pca',filename_style='continuity',annotation_qc='provided')
        actual=core.load_slices(provided);engine.load_models(actual,provided)
        assert set(actual[0].annotations)=={'low_quality_A','low_quality_B'}
        checks.append('existing annotations ignored by default; explicit provided QC uses them without requiring label-constant PCA')
        for mode in ['onehot','constant','zero_norm','basis_mismatch','feature_hash','missing_manifest']:
            dst=root/mode;shutil.copytree(prepared,dst);manifest_path=dst/'expression_pca_manifest.json';m=json.loads(manifest_path.read_text())
            first=next((dst/'slice_h5ad').glob('*.h5ad'))
            if mode in {'onehot','constant','zero_norm'}:
                obj=ad.read_h5ad(first);v=obj.obsm['X_pca'].copy()
                if mode=='onehot':v[:]=0;v[:,0]=1
                elif mode=='constant':v[:]=v[0]
                else:v[0]=0
                obj.obsm['X_pca']=v;obj.write_h5ad(first)
            elif mode=='basis_mismatch':m['slices'][0]['basis_sha256']='0'*64;manifest_path.write_text(json.dumps(m))
            elif mode=='feature_hash':m['slices'][0]['features_sha256']='0'*64;manifest_path.write_text(json.dumps(m))
            else:manifest_path.unlink()
            text={'onehot':'one-hot','constant':'constant','zero_norm':'zero-norm','basis_mismatch':'bases differ','feature_hash':'feature hash','missing_manifest':'requires expression_pca_manifest'}[mode]
            expect_error(lambda:validator.validate_inputs(dst/'slice_h5ad','expression-pca',filename_style='continuity'),text)
        checks.append('one-hot disguised PCA, constant, zero-norm, distinct per-slice basis, bad feature hash and missing provenance all rejected')
        pair=module('pairwise_validator_expression',SKILL/'pipelines/pairwise-rigid/validate_inputs.py')
        pair.validate_inputs(slices_dir,'expression-pca')
        checks.append('same generated PCA output also passes pairwise filename/provenance contract')
    report={'status':'passed','checks':checks,'real_alignment_executed':False,'scope':'synthetic counts/PCA and fake Spateo calls; no production data',
            'anonymous_policy_sha256':validator.sha(PIPELINE/'_universal_policy.py')}
    if options.report:options.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
