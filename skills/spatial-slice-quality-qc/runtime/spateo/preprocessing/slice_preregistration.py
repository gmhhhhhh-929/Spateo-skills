"""Display-only SE(2) registration. Never feeds registered geometry into QC.

Homogeneous column convention: display = (T @ [raw_x, raw_y, 1])[:2].
The middle slice is the dataset reference; pair transforms map child to parent.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class PreregistrationConfig:
    sample_size: int = 1000
    seed: int = 13
    angles: int = 24
    max_iter: int = 50
    trim_fraction: float = 0.8
    ambiguity_ratio: float = 1.05
    residual_warning: float = 0.08
    min_points: int = 20
    annotation_weight: float = 0.0


def transform_points(xy, matrix):
    xy = np.asarray(xy, dtype=float)
    t = np.asarray(matrix, dtype=float)
    return xy @ t[:2, :2].T + t[:2, 2]


def validate_rigid(matrix):
    t = np.asarray(matrix, dtype=float)
    if t.shape != (3, 3) or not np.isfinite(t).all():
        raise ValueError('Expected finite 3 x 3 rigid matrix')
    if not np.allclose(t[2], [0, 0, 1], atol=1e-9):
        raise ValueError('Invalid homogeneous transform')
    r = t[:2, :2]
    if not np.allclose(r.T @ r, np.eye(2), atol=1e-7) or not np.isclose(np.linalg.det(r), 1, atol=1e-7):
        raise ValueError('Scaling, reflection or nonrigid transform is forbidden')
    return t


def rigid_fit(source, target):
    a, b = np.asarray(source), np.asarray(target)
    ac, bc = a.mean(0), b.mean(0)
    u, _, vt = np.linalg.svd((a-ac).T @ (b-bc))
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1] *= -1
        r = vt.T @ u.T
    t = np.eye(3)
    t[:2, :2], t[:2, 2] = r, bc-r @ ac
    return validate_rigid(t)


def sample_indices(n, maximum, seed, identity):
    salt = int(hashlib.sha256(str(identity).encode()).hexdigest()[:8], 16)
    return np.sort(np.random.default_rng(seed+salt).choice(n, min(n, maximum), replace=False))


def roi_mask(raw_xy, matrix, bounds):
    """Exact rectangle test in display frame; no raw-axis bounding-box shortcut."""
    xmin, xmax, ymin, ymax = map(float, bounds)
    if not np.isfinite([xmin,xmax,ymin,ymax]).all() or xmin > xmax or ymin > ymax:
        raise ValueError('Invalid ROI bounds')
    xy = transform_points(raw_xy, validate_rigid(matrix))
    return (np.isfinite(xy).all(1) & (xy[:,0]>=xmin) & (xy[:,0]<=xmax)
            & (xy[:,1]>=ymin) & (xy[:,1]<=ymax))


def inverse_roi_polygon(matrix, bounds):
    x0,x1,y0,y1 = bounds
    return transform_points([[x0,y0],[x1,y0],[x1,y1],[x0,y1]], np.linalg.inv(validate_rigid(matrix)))


def _residual(a, b):
    # Symmetric, untrimmed median distance retains a check independent of fitting trim.
    da = cKDTree(b).query(a)[0]
    db = cKDTree(a).query(b)[0]
    return float((np.median(da)+np.median(db))/2)


def register_pair(source, target, config=None, source_labels=None, target_labels=None):
    cfg = config or PreregistrationConfig()
    a,b = np.asarray(source,dtype=float),np.asarray(target,dtype=float)
    if min(len(a),len(b)) < cfg.min_points or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Insufficient finite points for rigid fitting')
    if min(np.linalg.matrix_rank(a-a.mean(0)),np.linalg.matrix_rank(b-b.mean(0))) < 2:
        raise ValueError('Degenerate/collinear geometry')
    if cfg.angles < 1 or cfg.max_iter < 1 or not 0 < cfg.trim_fraction <= 1:
        raise ValueError('Invalid rigid fitting configuration')
    scale = float(np.linalg.norm(np.quantile(b,.95,axis=0)-np.quantile(b,.05,axis=0)))
    if scale <= 0: raise ValueError('Zero tissue extent')
    if not 0 <= cfg.annotation_weight <= 1: raise ValueError('annotation_weight must lie in [0, 1]')
    typed=[];coverage=0.0
    if cfg.annotation_weight and source_labels is not None and target_labels is not None:
        la,lb=np.asarray(source_labels,dtype=str),np.asarray(target_labels,dtype=str)
        if len(la)!=len(a) or len(lb)!=len(b):raise ValueError('Annotation/coordinate length mismatch')
        common=[k for k in np.intersect1d(la,lb) if k.lower() not in ('nan','none','na','unknown','unassigned','') and min(np.sum(la==k),np.sum(lb==k))>=5]
        coverage=min(float(np.isin(la,common).mean()),float(np.isin(lb,common).mean()))
        if len(common)>=2 and coverage>=.5:typed=[(la==k,lb==k) for k in common]
    tree=cKDTree(b); candidates=[]
    for angle in np.arange(cfg.angles)*2*np.pi/cfg.angles:
        r=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        t=np.eye(3); t[:2,:2]=r; t[:2,2]=np.median(b,0)-r@np.median(a,0)
        for iteration in range(cfg.max_iter):
            moved=transform_points(a,t)
            distances,idx=tree.query(moved)
            keep=distances<=np.quantile(distances,cfg.trim_fraction)
            if keep.sum()<3: break
            update=rigid_fit(moved[keep],b[idx[keep]])
            t=update@t
            if np.max(np.linalg.norm(transform_points(moved,update)-moved,axis=1))<scale*1e-6: break
        moved=transform_points(a,t);geometry=_residual(moved,b)/scale
        annotation=float(np.mean([_residual(moved[ia],b[ib])/scale for ia,ib in typed])) if typed else None
        score=(1-cfg.annotation_weight)*geometry+cfg.annotation_weight*annotation if typed else geometry
        candidates.append((score,t,iteration+1,geometry,annotation))
    candidates.sort(key=lambda x:x[0]); score,t,it,geometry,annotation=candidates[0]
    theta=np.arctan2(t[1,0],t[0,0])
    alternatives=[v for v in candidates[1:] if abs(np.arctan2(np.sin(np.arctan2(v[1][1,0],v[1][0,0])-theta),np.cos(np.arctan2(v[1][1,0],v[1][0,0])-theta)))>np.pi/6]
    ratio=alternatives[0][0]/max(score,1e-12) if alternatives else None
    warnings=[]
    if ratio is not None and ratio<cfg.ambiguity_ratio: warnings.append('orientation_ambiguous')
    if geometry>cfg.residual_warning: warnings.append('low_geometric_overlap')
    if cfg.annotation_weight and not typed:warnings.append('annotation_support_unavailable')
    if annotation is not None and annotation>cfg.residual_warning:warnings.append('annotation_pattern_mismatch')
    return validate_rigid(t), {'normalized_median_residual':geometry,
        'normalized_selection_cost':score,'annotation_residual':annotation,
        'annotation_used':bool(typed),'annotation_common_labels':common if typed else [],'annotation_coverage':coverage,
        'raw_normalized_median_residual':_residual(a,b)/scale,
        'alternative_orientation_ratio':ratio,'iterations':it,'warnings':warnings,
        'candidate_count':cfg.angles,'status':'warning' if warnings else 'estimated',
        'accuracy_claim':'Geometric fit only; anatomical correspondence unvalidated'}


def register_series(slices, dataset_id, config=None):
    """slices are ordered dictionaries containing slice_id and raw_xy only."""
    cfg=config or PreregistrationConfig(); start=time.perf_counter()
    if not slices: raise ValueError('Empty series')
    ids=[s['slice_id'] for s in slices]
    if len(set(ids)) != len(ids): raise ValueError('Duplicate slice ids')
    samples=[];sample_labels=[]
    for s in slices:
        xy=np.asarray(s['raw_xy']); valid=np.isfinite(xy).all(1);finite=xy[valid]
        idx=sample_indices(len(finite),cfg.sample_size,cfg.seed,s['slice_id']);samples.append(finite[idx])
        sample_labels.append(np.asarray(s['labels'],dtype=str)[valid][idx] if s.get('labels') is not None else None)
    middle=len(slices)//2; transforms={middle:np.eye(3)}; metadata={middle:{'status':'reference','warnings':[],'parent':None,'path_length':0}}
    reference=samples[middle]
    if len(reference)<cfg.min_points or np.linalg.matrix_rank(reference-reference.mean(0))<2:
        metadata[middle].update(status='failed',warnings=['invalid_reference_geometry'])
    edges=[]
    for indices in [range(middle+1,len(slices)), range(middle-1,-1,-1)]:
        for i in indices:
            parent=i-1 if i>middle else i+1
            try:
                pair,info=register_pair(samples[i],samples[parent],cfg,sample_labels[i],sample_labels[parent])
            except Exception as exc:
                pair=np.eye(3); info={'status':'failed','warnings':['fit_failed'],'error':f'{type(exc).__name__}: {exc}'}
            transforms[i]=transforms[parent]@pair
            warnings=list(info['warnings'])
            if metadata[parent]['status'] in ('failed','blocked'):
                info['status']='blocked'; warnings.append('failed_reference_chain')
            elif metadata[parent]['warnings']:
                warnings.append('uncertain_reference_chain')
                if info['status']=='estimated': info['status']='warning'
            metadata[i]={**info,'warnings':warnings,'parent':ids[parent], 'path_length':abs(i-middle)}
            edges.append({'source':ids[i],'target':ids[parent],'pair_matrix':pair.tolist(),**info})
    # Two-step overlap check: advisory drift/shape-change flag; no extra fitting.
    for i in range(len(slices)):
        if i==middle: continue
        j=i-2 if i>middle+1 else i+2 if i<middle-1 else None
        if j is not None and len(samples[i]) and len(samples[j]):
            a,b=transform_points(samples[i],transforms[i]),transform_points(samples[j],transforms[j])
            extent=max(float(np.linalg.norm(np.ptp(b,axis=0))),1e-12)
            value=_residual(a,b)/extent; metadata[i]['two_step_residual']=value
            if value>cfg.residual_warning:
                metadata[i]['warnings'].append('two_step_drift_or_shape_change')
                if metadata[i]['status']=='estimated': metadata[i]['status']='warning'
    # Drift checks happen after fitting; propagate their uncertainty along the
    # already-composed path as well as the pair-fit warnings above.
    for indices in [range(middle+1,len(slices)),range(middle-1,-1,-1)]:
        for i in indices:
            parent=i-1 if i>middle else i+1
            if metadata[parent]['warnings'] and 'uncertain_reference_chain' not in metadata[i]['warnings']:
                metadata[i]['warnings'].append('uncertain_reference_chain')
                if metadata[i]['status']=='estimated': metadata[i]['status']='warning'
    parameters=asdict(cfg)
    if not cfg.annotation_weight:parameters.pop('annotation_weight')  # retain existing geometry-only frame IDs
    frame_contract={'ids':ids,'config':parameters,'coordinate_hashes':[hashlib.sha256(np.asarray(s['raw_xy'],dtype='<f8').tobytes()).hexdigest() for s in slices]}
    if cfg.annotation_weight:frame_contract['annotation_hashes']=[hashlib.sha256(json.dumps(np.asarray(s['labels'],dtype=str).tolist()).encode()).hexdigest() if s.get('labels') is not None else None for s in slices]
    frame_id='display:'+str(dataset_id)+':'+hashlib.sha256(json.dumps(frame_contract,sort_keys=True).encode()).hexdigest()[:16]
    return {'schema_version':'1.0','dataset_id':dataset_id,'frame_id':frame_id,
        'reference_slice':ids[middle],'coordinate_convention':'column homogeneous: display = T @ raw; T_child = T_parent @ T_child_to_parent',
        'input_coordinate_key':'spatial','output_coordinate_key':'display_preregistered_xy',
        'method':'rigid multistart trimmed ICP (SE2, scipy KDTree + Kabsch)'+('; annotation-assisted candidate ranking' if cfg.annotation_weight else ''),
        'parameters':parameters,'qc_dependency':'none; all QC computed from raw inputs',
        'edges':edges,'slices':{ids[i]:{**metadata[i],'matrix':transforms[i].tolist(),
            'inverse_matrix':np.linalg.inv(transforms[i]).tolist(),'n_points':len(slices[i]['raw_xy']),
            'fit_points':len(samples[i])} for i in range(len(slices))},
        'elapsed_seconds':time.perf_counter()-start,
        'warning':'Same display ROI does not imply identical anatomy. Cumulative drift remains possible; this is not downstream alignment.'}
