#!/usr/bin/env python3
"""Config-driven, immutable review of regional tissue coverage and body envelopes."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pyvista as pv
import spateo as st
from matplotlib.colors import to_rgba
from coverage_envelope import make_grid, coverage, regional_candidate, smooth_envelope, inside, repair_supported_gaps


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(1048576),b''):h.update(chunk)
    return h.hexdigest()


def measure_save(mesh,points,name,color,path):
    mesh.cell_data['mesh_label']=np.full(mesh.n_cells,name)
    mesh.cell_data['mesh_label_rgba']=np.tile(to_rgba(color),(mesh.n_cells,1))
    st.tdr.save_model(mesh,str(path),binary=True);loaded=st.tdr.read_model(str(path))
    assert np.array_equal(mesh.points,loaded.points) and np.array_equal(mesh.faces,loaded.faces)
    assert np.array_equal(mesh['mesh_label_rgba'],loaded['mesh_label_rgba'])
    assert loaded.n_open_edges==0
    hit=inside(points,loaded)
    f=st.tdr.model_morphology(loaded,pc=pv.PolyData(points[hit]))
    f.update(Cells_inside=int(hit.sum()),Cell_density_unrounded=int(hit.sum())/loaded.volume,
             Sphericity=float(np.pi**(1/3)*(6*loaded.volume)**(2/3)/loaded.area))
    return {'vtk':str(path),'sha256':digest(path),'features':f,'source_selected_cells':len(points),
            'inside_fraction':float(hit.mean()),'open_edges':0,'roundtrip_passed':True,
            'components_kept':len(loaded.split_bodies()),'components_removed':0,
            'removed_volume_fraction':0,'components':[]}


def preview(body,meshes,colors,center,basis,path):
    width=np.ptp((body.points-center)@basis[0])
    cams=[[(0,-width*.65,width*1.6),(0,0,0),(0,1,0)],[(0,0,width*2),(0,0,0),(0,1,0)],
          [(0,-width*2,0),(0,0,0),(0,0,1)],[(width*2,0,0),(0,0,0),(0,0,1)]]
    p=pv.Plotter(shape=(2,2),off_screen=True,window_size=(1400,1000),border=False)
    def oriented(mesh):
        m=mesh.copy();m.points=(mesh.points-center)@basis.T;m.point_data.clear();m.cell_data.clear();return m
    b=oriented(body)
    for i,cam in enumerate(cams):
        p.subplot(i//2,i%2);p.set_background('white');p.add_text(['Isometric','XY','XZ','YZ'][i],font_size=14,color='#222222')
        p.add_mesh(b,color='#BFC7CD',opacity=.2 if meshes else 1,smooth_shading=True)
        for m,color in zip(meshes,colors):p.add_mesh(oriented(m),color=color,smooth_shading=True)
        p.camera_position=cam;p.enable_parallel_projection();p.camera.parallel_scale=width*.4
    p.show(screenshot=str(path))


def run(config):
    out=Path(config['output_dir']).resolve()
    if out.exists():raise FileExistsError(out)
    pointpath=Path(config['point_cloud']).resolve();pc=st.tdr.read_model(str(pointpath))
    if 'obs_index' not in pc.point_data or len(np.unique(pc['obs_index']))!=pc.n_points:
        raise ValueError('Unique obs_index required')
    labels=np.asarray(pc[config['label_key']]).astype(str)
    sources={str(pointpath):digest(pointpath)}
    oldbody=st.tdr.read_model(config['previous_body'])
    sources[config['previous_body']]=digest(config['previous_body'])
    names=[s['name'] for s in config['tissues']]
    if len(set(n.casefold() for n in names))!=len(names) or any(not n or '/' in n or '\\' in n or n.casefold() in ['body','.','..'] for n in names):
        raise ValueError('Unsafe or duplicate tissue names')
    oldmeshes={s['name']:st.tdr.read_model(s['vtk']) for s in config['tissues']}
    for s in config['tissues']:sources[s['vtk']]=digest(s['vtk'])
    grid=make_grid(pc.points,[oldbody,*oldmeshes.values()],long_axis=config.get('long_axis',220))
    center=np.asarray(pc.points).mean(0);_,_,basis=np.linalg.svd(pc.points-center,full_matrices=False)
    if np.linalg.det(basis)<0:basis[-1]*=-1
    out.mkdir(parents=True)
    report={'status':'review_candidate','config':config,'source_hashes':sources,'source_cells_changed':False,
            'grid':{'origin':grid.origin.tolist(),'shape':list(grid.shape),'spacing':grid.spacing},
            'basis_rows':basis.tolist(),'center':center.tolist(),'meshes':{},'warnings':[]}
    (out/'config.json').write_text(json.dumps(config,indent=2))
    selected={};colors=[]
    for spec in config['tissues']:
        name=spec['name'];mask=np.isin(labels,spec['values']);pts=pc.points[mask]
        if not mask.any():raise ValueError(f'No cells for {name}')
        old=oldmeshes[name];before,oldhit=coverage(pts,old,basis[0],center)
        chosen=old;after=before;attempts=[];params=None;eligible=[];target_met=False
        if spec.get('refine',False):
            for q in config.get('quantiles',[.02,.01,.005]):
                mesh,settings=regional_candidate(pts,old,grid,basis[0],center,quantile=q,
                    sigma=config.get('tissue_sigma',3.2),min_component_cells=config.get('min_component_cells',20))
                local_repairs=[]
                for _ in range(config.get('local_repair_passes',1)):
                    mesh,repair=repair_supported_gaps(pts,mesh,grid,sigma=config.get('tissue_sigma',3.2),
                                                    min_component_cells=config.get('min_component_cells',20))
                    local_repairs.append(repair)
                settings['local_repairs']=local_repairs
                cov,hit=coverage(pts,mesh,basis[0],center)
                ratio=mesh.volume/old.volume
                accepted=(cov['global']>=config.get('global_target',.985)
                          and (cov['min_axial'] is None or cov['min_axial']>=config.get('regional_target',.9))
                          and (cov['min_spatial_block'] is None or cov['min_spatial_block']>=config.get('block_target',.8))
                          and ratio<=config.get('max_volume_ratio',2.0))
                attempts.append({'quantile':q,'coverage':cov,'volume_ratio':ratio,'accepted':bool(accepted)})
                st.tdr.save_model(mesh,str(out/f'{name}.q{q}.candidate.vtk'),binary=True)
                if accepted:chosen=mesh;after=cov;params=settings;target_met=True;break
                # A soft regional target must not silently discard a clear improvement.
                # Retain minimally expansive Pareto-improving candidates with warnings.
                improves=(cov['global']>before['global']+.005 and ratio<=config.get('max_volume_ratio',2.)
                          and all(cov[k] is None or before[k] is None or cov[k]>=before[k] for k in ['min_axial','min_spatial_block']))
                if improves:eligible.append((ratio,mesh,cov,settings))
            if params is None and eligible:
                _,chosen,after,params=min(eligible,key=lambda item:item[0])
                report['warnings'].append(f'{name}: improved all coverage summaries within volume guard, but regional targets remain unmet; review required.')
            if params is None:report['warnings'].append(f'{name}: no candidate met volume/improvement guards; previous geometry retained.')
        selected[name]=chosen;colors.append(spec['color'])
        record=measure_save(chosen,pts,name,spec['color'],out/f'{name}.vtk')
        newhit=inside(pts,chosen)
        record.update(label=spec.get('label',name),components_before=len(old.split_bodies()),
            coverage_review={'before':before,'after':after,'attempts':attempts,'parameters':params,
                'volume_change_fraction':chosen.volume/old.volume-1,
                'previous_geometry_retained':params is None,'requested_targets_met':target_met if spec.get('refine',False) else None,
                'newly_outside_obs_ids':np.asarray(pc['obs_index'])[mask][oldhit&~newhit].astype(str).tolist(),
                'remaining_outside_obs_ids':np.asarray(pc['obs_index'])[mask][~newhit].astype(str).tolist()},
            display_note='Regional and spatial-block source-cell coverage reviewed. Original cells are unchanged. Sparse supported components are retained; coverage is not anatomical validation.')
        report['meshes'][name]=record
        print(name,round(before['global'],4),'->',round(after['global'],4),'min axial',after['min_axial'],'volume',chosen.volume/old.volume,flush=True)
    body,bodyqc=smooth_envelope(pc.points,list(selected.values()),grid,
                              sigma=config.get('body_sdf_sigma',2.5),clearance_voxels=config.get('body_clearance_voxels',1.8))
    note='Smooth signed-distance display envelope of point-supported body and retained tissue union, with a numerical containment certificate. Clearance and filled body cavities are display choices, not anatomical growth. Tissue surfaces can overlap and are not a space-filling segmentation.'
    bodyrecord=measure_save(body,pc.points,'body','#BFC7CD',out/'body.vtk')
    bodyrecord.update(label='Body · smooth envelope',components_before=len(oldbody.split_bodies()),display_note=note,
                      envelope_review=bodyqc,volume_change_fraction=body.volume/oldbody.volume-1)
    report['meshes']['body']=bodyrecord
    report['body_qc']=bodyqc
    report['body_qc']['previous_body_volume']=float(oldbody.volume)
    report['body_qc']['new_body_volume']=float(body.volume)
    for path,hash_ in sources.items():assert digest(path)==hash_
    preview(body,[],[],center,basis,out/'body.four_views.png')
    preview(body,list(selected.values()),colors,center,basis,out/'all_tissues.four_views.png')
    for spec in config['tissues']:
        if spec.get('refine',False):preview(body,[selected[spec['name']]],[spec['color']],center,basis,out/f'{spec["name"]}.four_views.png')
    (out/'manifest.json').write_text(json.dumps(report,indent=2))
    print('COMPLETE',out,'union/body',bodyqc['union_over_body_voxel_fraction'],flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',required=True,type=Path)
    args=parser.parse_args();run(json.loads(args.config.read_text()))
