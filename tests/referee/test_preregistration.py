import numpy as np
import pytest
from scipy import sparse
from spateo.preprocessing.slice_preregistration import (PreregistrationConfig, register_pair,
 register_series, transform_points, rigid_fit, roi_mask, inverse_roi_polygon, validate_rigid)
from spateo.preprocessing.slice_quality import _matrix_row_nnz


def shape():
 rng=np.random.default_rng(123)
 return np.vstack([rng.normal([0,0],[2,.3],(220,2)),rng.normal([1,1],[.3,.8],(100,2))])


def test_rigid_recovery_and_reproducibility():
 a=shape(); angle=1.2;t=np.array([[np.cos(angle),-np.sin(angle),9],[np.sin(angle),np.cos(angle),-4],[0,0,1]])
 b=transform_points(a,t); cfg=PreregistrationConfig(angles=24,max_iter=70)
 fitted,info=register_pair(a,b,cfg); again,_=register_pair(a,b,cfg)
 assert np.allclose(transform_points(a,fitted),b,atol=1e-4)
 assert np.array_equal(fitted,again)
 assert info['normalized_median_residual']<1e-5


def test_annotations_disambiguate_symmetric_rigid_candidates():
 half=shape()+[4,1];a=np.vstack([half,-half]);labels=np.repeat(['head','tail'],len(half))
 t=np.array([[0,-1,9],[1,0,-4],[0,0,1]]);b=transform_points(a,t)
 fitted,info=register_pair(a,b,PreregistrationConfig(angles=8,max_iter=30,annotation_weight=.5),labels,labels)
 assert info['annotation_used'] and info['annotation_coverage']==1
 assert np.allclose(transform_points(a,fitted),b,atol=1e-5)
 _,missing=register_pair(a,b,PreregistrationConfig(angles=4,max_iter=5,annotation_weight=.5))
 assert 'annotation_support_unavailable' in missing['warnings']


def test_composition_inverse_and_exact_rotated_roi():
 a=shape();t=rigid_fit(a,transform_points(a,np.array([[0,-1,3],[1,0,-2],[0,0,1]])))
 b=np.array([[np.cos(.6),-np.sin(.6),1],[np.sin(.6),np.cos(.6),2],[0,0,1]])
 assert np.allclose(transform_points(transform_points(a,t),b),transform_points(a,b@t))
 assert np.allclose(transform_points(transform_points(a,t),np.linalg.inv(t)),a)
 bounds=[-.5,1,-.25,.8];polygon=inverse_roi_polygon(b,bounds)
 assert np.allclose(transform_points(polygon,b),[[-.5,-.25],[1,-.25],[1,.8],[-.5,.8]])
 mask=roi_mask(a,b,bounds); transformed=transform_points(a,b)
 assert np.array_equal(mask,((transformed>=[-.5,-.25])&(transformed<=[1,.8])).all(1))
 bbox=((a>=polygon.min(0))&(a<=polygon.max(0))).all(1)
 assert np.any(bbox&~mask)  # bbox must over-select some points


def test_failure_chain_and_dataset_isolation():
 a=shape();s=[{'slice_id':str(i),'raw_xy':a.copy()} for i in range(5)]
 s[1]['raw_xy']=np.zeros((3,2))
 one=register_series(s,'A',PreregistrationConfig(angles=2,max_iter=2))
 two=register_series(s,'B',PreregistrationConfig(angles=2,max_iter=2))
 assert one['reference_slice']=='2';assert one['slices']['1']['status']=='failed'
 assert one['slices']['0']['status']=='blocked';assert one['frame_id']!=two['frame_id']
 assert np.array_equal(s[2]['raw_xy'],a)
 s[2]['raw_xy']=np.zeros((3,2))
 invalid=register_series(s,'invalid_reference',PreregistrationConfig(angles=2,max_iter=2))
 assert invalid['slices']['2']['status']=='failed'
 assert invalid['slices']['3']['status']=='blocked'


def test_no_scale_no_reflection_and_sparse_zeros():
 for t in [np.diag([2,1,1]),np.diag([-1,1,1])]:
  with pytest.raises(ValueError):validate_rigid(t)
 x=sparse.csr_matrix(([0.,2.,0.,3.],([0,0,1,1],[0,1,0,1])),shape=(2,2));copy=x.copy()
 assert np.array_equal(_matrix_row_nnz(x),[1,1]);assert np.array_equal(copy.data,x.data)


def test_selected_counts_override_stale_obs_and_onehot_is_missing(tmp_path):
 from anndata import AnnData
 import pandas as pd
 from spateo.preprocessing.slice_quality import scan_h5ad_series,SliceQCConfig
 rng=np.random.default_rng(3);x=sparse.csr_matrix(rng.poisson(2,(180,12)))
 ad=AnnData(x,obs=pd.DataFrame({'slice_id':np.repeat(['S1','S2','S3'],60)},index=[str(i) for i in range(180)]))
 ad.obsm['spatial']=rng.normal(size=(180,2));ad.layers['counts']=x
 ad.obs['total_counts']=9999;ad.obs['n_genes_by_counts']=9999
 p=tmp_path/'counts.h5ad';ad.write_h5ad(p)
 result=scan_h5ad_series([p],spatial_key='spatial',layer='counts')
 assert result.metrics.median_total_counts.max()<100
 assert result.provenance['sources'][0]['obs_totals_match_selected_matrix'] is False
 ad.uns['blind_input_contract']={'expression_representation':"global obs['anno'] one-hot"}
 p=tmp_path/'onehot.h5ad';ad.write_h5ad(p)
 result=scan_h5ad_series([p],spatial_key='spatial',layer='X')
 assert result.metrics.median_total_counts.isna().all()
 assert result.metrics.detected_genes.isna().all()
 assert result.profiles is None
 assert not result.metrics.expression_capture_available.any()


def test_saved_roi_roundtrip_and_report_modes(tmp_path):
 import json,sys
 from pathlib import Path
 from anndata import AnnData
 import pandas as pd
 from spateo.preprocessing.slice_quality import (scan_h5ad_series,write_slice_quality_outputs,
     HighConfidencePolicy,write_high_confidence_outputs)
 sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'skills/spatial-slice-quality-qc/scripts'))
 from slice_quality_preregistration import preregister_outputs,export_roi
 from slice_quality_report import render_slice_quality_report_from_outputs,render_binary_slice_quality_report_from_outputs,write_slice_quality_collection_report
 paths=[];a=shape()
 for i in range(3):
  ad=AnnData(sparse.csr_matrix(np.ones((len(a),10))*(i+1)));ad.obsm['spatial']=a+[i,2*i]
  ad.obs_names=[f'p{j}' for j in range(len(a))];p=tmp_path/f'slice{i}.h5ad';ad.write_h5ad(p);paths.append(p)
 result=scan_h5ad_series(paths,spatial_key='spatial');run=tmp_path/'run';write_slice_quality_outputs(result,run,write_display_payload=True)
 before=(run/'slice_quality_metrics.csv').read_bytes()
 reg=preregister_outputs(run,'test',PreregistrationConfig(angles=4,max_iter=10),120)
 assert before==(run/'slice_quality_metrics.csv').read_bytes()
 roi={'frame_id':reg['frame_id'],'bounds':[-1,4,0,4]};p=tmp_path/'roi.json';p.write_text(json.dumps(roi))
 summary=export_roi(run,p,tmp_path/'roi.csv');table=pd.read_csv(tmp_path/'roi.csv')
 assert len(table)==sum(summary['counts'].values())>0
 assert table.source_row.max()<len(a)
 render_slice_quality_report_from_outputs(run)
 write_high_confidence_outputs(result.metrics,HighConfidencePolicy(.1,.7,unresolved_action='keep'),run)
 render_binary_slice_quality_report_from_outputs(run,run/'binary.html')
 for page in [run/'slice_quality_report.html',run/'binary.html']:
  text=page.read_text();assert 'coordinateMode' in text and 'source_row' in text and 'Same ROI does not guarantee' in text
 report=write_slice_quality_collection_report([run],tmp_path/'index.html');assert report['slices']==3
 roi['frame_id']='wrong';p.write_text(json.dumps(roi))
 with pytest.raises(ValueError,match='frame mismatch'):export_roi(run,p,tmp_path/'wrong.csv')


def test_saved_viewer_preserves_numeric_slice_identity(tmp_path):
 import json,sys,re
 from pathlib import Path
 import pandas as pd
 from anndata import AnnData
 from spateo.preprocessing.slice_quality import scan_h5ad_series,write_slice_quality_outputs
 sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'skills/spatial-slice-quality-qc/scripts'))
 from slice_quality_preregistration import preregister_outputs
 from slice_quality_report import render_slice_quality_report_from_outputs
 ad=AnnData(sparse.csr_matrix(np.ones((120,5))),obs=pd.DataFrame({'slice_id':np.repeat(['001','002','003'],40)},index=[f'p{i}' for i in range(120)]))
 ad.obsm['spatial']=np.random.default_rng(7).normal(size=(120,2));source=tmp_path/'numeric.h5ad';ad.write_h5ad(source)
 result=scan_h5ad_series([source],slice_key='slice_id',spatial_key='spatial');run=tmp_path/'numeric_qc';write_slice_quality_outputs(result,run,write_display_payload=True)
 preregister_outputs(run,'numeric',PreregistrationConfig(angles=2,max_iter=2),40)
 page=Path(render_slice_quality_report_from_outputs(run)).read_text()
 payload=json.loads(re.search(r'<script id="payload" type="application/json">(.*?)</script>',page,re.S)[1])
 assert [r['slice_id'] for r in payload['records']]==['001','002','003']
 assert set(payload['points'])=={'001','002','003'}


def test_new_input_scope_preserves_frozen_decisions_without_certification(tmp_path):
 import json
 import pandas as pd
 from spateo.preprocessing.slice_quality import HighConfidencePolicy,write_high_confidence_outputs,scan_h5ad_series
 from anndata import AnnData
 a=AnnData(sparse.csr_matrix(np.ones((120,8))))
 a.obs['slices']=np.repeat(['001','002','003'],40)
 a.obsm['spatial']=np.random.default_rng(24).normal(size=(120,2))
 p=tmp_path/'input.h5ad';a.write_h5ad(p)
 result=scan_h5ad_series([p],slice_key='slices',spatial_key='spatial')
 policy=HighConfidencePolicy(.1,.7,unresolved_action='keep')
 write_high_confidence_outputs(result.metrics,policy,tmp_path/'old')
 write_high_confidence_outputs(result.metrics,policy,tmp_path/'new',application_scope='new_input_unvalidated')
 old=pd.read_csv(tmp_path/'old/slice_quality_binary_audit.csv');new=pd.read_csv(tmp_path/'new/slice_quality_binary_audit.csv')
 assert new.final_call.tolist()==old.final_call.tolist()
 assert new.certified_call.isna().all()
 assert 'independently_certified' not in new.decision_basis.tolist()
 assert new.policy_gate_basis.tolist()==old.decision_basis.tolist()
 assert json.loads((tmp_path/'new/binary_policy_application.json').read_text())['certified']==0


def test_detailed_preregistered_report_reuses_raw_cache_and_identity(tmp_path):
 import json,sys,re
 from pathlib import Path
 from anndata import AnnData
 from spateo.preprocessing.slice_quality import scan_h5ad_series,write_slice_quality_outputs,write_high_confidence_outputs,HighConfidencePolicy
 sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'skills/spatial-slice-quality-qc/scripts'))
 from slice_quality_preregistration import preregister_outputs
 from slice_quality_visualization import render_binary_slice_quality_appendix,_component_ranks
 a=AnnData(sparse.csr_matrix(np.ones((120,5))))
 a.obs['slices']=np.repeat(['001','002','003'],40);a.obsm['spatial']=np.random.default_rng(32).normal(size=(120,2))
 p=tmp_path/'input.h5ad';a.write_h5ad(p)
 result=scan_h5ad_series([p],slice_key='slices',spatial_key='spatial');root=tmp_path/'run'
 write_slice_quality_outputs(result,root,write_display_payload=True)
 m=preregister_outputs(root,'fixture',PreregistrationConfig(angles=4,max_iter=5),30)
 write_high_confidence_outputs(result.metrics,HighConfidencePolicy(.1,.7,unresolved_action='keep'),root,application_scope='new_input_unvalidated')
 output=tmp_path/'report';summary=render_binary_slice_quality_appendix(root,output,language='en',max_points_per_slice=20)
 page=(output/'index.html').read_text();d=json.loads(re.search(r'<script id="payload" type="application/json">(.*?)</script>',page,re.S).group(1))
 assert d['order']==['001','002','003']
 assert set(r['final_call'] for r in d['records'])<={'keep','exclude'}
 assert summary['connected_component_maps'] and summary['absolute_expression_maps']
 assert 'not independently validated on this input' in page
 for sid,points in d['points'].items():
  raw=np.column_stack([points['x'],points['y']]);display=np.column_stack([points['display_x'],points['display_y']])
  assert np.allclose(transform_points(raw,m['slices'][sid]['matrix']),display)
  assert points['available']==40 and points['displayed']==20
  assert len(points['component_rank'])==20 and points['captured_counts']==[5.]*20
  assert np.allclose(a.obsm['spatial'][points['source_row']],raw)
  assert a.obs_names[points['source_row']].tolist()==points['obs_names']
 # Cached evidence is immutable and checksum checked even when the output is a new viewer.
 cache=root/'display_preregistration'/m['slices']['001']['cache_file'];cache.write_bytes(cache.read_bytes()+b'changed')
 with pytest.raises(ValueError,match='checksum'):render_binary_slice_quality_appendix(root,tmp_path/'invalid')


def test_keep_only_band_does_not_execute_strict_exclusion_checks():
 import json,pandas as pd
 from spateo.preprocessing.slice_quality import HighConfidencePolicy,ReviewEvidenceTier,apply_high_confidence_policy
 from slice_quality_visualization import build_directional_slice_explanations
 scores=[.129,.129001,.539999,.54]
 rows=[]
 for i,s in enumerate(scores):
  rows.append(dict(slice_id=str(i),quality_anomaly_score=s,recommendation='review',window_context='two_sided',partial_structure_protection=False,score_confidence=1.,corroborating_domains=4,density_domain_score=1.,expression_domain_score=1.,damage_domain_score=1.,continuity_domain_score=1.,adaptive_window_details=json.dumps([dict(window=w,score=.6,detector_call='review',window_context='two_sided',partial_structure_protection=False,corroborating_domains=4,maximum_domain_score=1.) for w in [3,5,7]])))
 policy=HighConfidencePolicy(.129,.7,unresolved_action='keep',review_exclusion_tiers=(ReviewEvidenceTier(name='mid_low_keep_only',min_score=.129,max_score=.54,enable_exclude=False),ReviewEvidenceTier(name='high',min_score=.54,max_score=None,severe_domain_threshold=.8)))
 out=apply_high_confidence_policy(pd.DataFrame(rows),policy)
 assert out.final_call.tolist()==['keep','keep','keep','exclude']
 assert out.review_resolution_reason.iloc[1:3].str.contains('calibrated keep-only score band').all()
 reasons=build_directional_slice_explanations(out,language='en')
 assert 'not executed' in reasons.user_reason.iloc[1]


def test_every_review_including_downgraded_high_scores_has_stage2_record():
 import json,pandas as pd
 from spateo.preprocessing.slice_quality import HighConfidencePolicy,ReviewEvidenceTier,apply_high_confidence_policy
 rows=[]
 for sid,score,context,protected,detector in [('low',.3,'two_sided',False,'review'),('enabled',.6,'two_sided',False,'review'),('terminal',.8,'one_sided',False,'exclude'),('protected',.8,'two_sided',True,'exclude'),('unsupported',.8,'two_sided',False,'keep')]:
  rows.append(dict(slice_id=sid,quality_anomaly_score=score,recommendation=detector,window_context=context,partial_structure_protection=protected,score_confidence=1.,corroborating_domains=4,density_domain_score=1.,expression_domain_score=1.,damage_domain_score=1.,continuity_domain_score=1.,adaptive_window_details=json.dumps([dict(window=w,score=score,detector_call=detector,window_context=context,partial_structure_protection=protected,corroborating_domains=4,maximum_domain_score=1.) for w in [3,5,7]])))
 policy=HighConfidencePolicy(.129,.7,unresolved_action='keep',review_exclusion_tiers=(ReviewEvidenceTier(name='keep_only',min_score=.129,max_score=.54,enable_exclude=False),ReviewEvidenceTier(name='high',min_score=.54,max_score=None,severe_domain_threshold=.8)))
 out=apply_high_confidence_policy(pd.DataFrame(rows),policy)
 assert out.threshold_triage_call.eq('review').all()
 assert out.review_resolution_tier.tolist()==['keep_only','high','high','high','high']
 assert out.final_call.tolist()==['keep','exclude','keep','keep','keep']
 assert out.review_resolution_reason.str.len().gt(10).all()
 assert out.review_tier_diagnostics.map(lambda s:len(json.loads(s))>0).all()
