"""Attach display-only registration to an existing QC run and replay exact ROI."""
from __future__ import annotations
import json
import time
from pathlib import Path

def preregister_outputs(input_dir, dataset_id, config=None, display_points=1800):
    import numpy as np
    import pandas as pd
    from anndata import read_h5ad
    from spateo.preprocessing.slice_quality import _resolve_slice_labels, _resolve_xy, _sha256, _choose_layer, _matrix_row_sum
    from spateo.preprocessing.slice_preregistration import PreregistrationConfig, register_series, sample_indices, transform_points
    root=Path(input_dir); cfg=config or PreregistrationConfig(); started=time.perf_counter()
    destination=root/'display_preregistration'
    destination.mkdir(exist_ok=False)  # immutable fit; reruns require a new QC output directory
    payload_path=root/'slice_quality_display_payload.json'
    payload=json.loads(payload_path.read_text()); metrics=pd.read_csv(root/'slice_quality_metrics.csv',dtype={'slice_id':str,'source_slice_id':str})
    sources={s['path']:s for s in payload['provenance']['sources']}; slices=[]
    # Read each file once, following exactly the detector's resolved labels and coordinates.
    for source_file,group in metrics.groupby('source_file',sort=False):
        contract=sources[source_file]
        before=_sha256(source_file)
        if before!=contract['sha256']: raise ValueError('Source changed since QC: '+source_file)
        adata=read_h5ad(source_file,backed='r')
        try:
            coord=contract['coordinate_source']
            key=coord.split(':',1)[1].split('[')[0] if coord.startswith('obsm:') else None
            xkey,ykey=coord.split(':',1)[1].split(',') if coord.startswith('obs:') else (None,None)
            x,y,_=_resolve_xy(adata,key,xkey,ykey)
            slice_key=contract['slice_source'].split(':',1)[1] if contract['slice_source'].startswith('obs:') else 'auto'
            labels,_,_=_resolve_slice_labels(adata,slice_key,key,Path(source_file).stem)
            _,matrix,_=_choose_layer(adata,contract['count_layer']); totals=_matrix_row_sum(matrix)
            for row in group.itertuples():
                idx=np.flatnonzero(labels.astype(str)==str(row.source_slice_id))
                if len(idx)!=row.n_locations: raise ValueError('QC/registration slice identity mismatch: '+row.slice_id)
                item={'slice_id':row.slice_id,'raw_xy':np.column_stack([x[idx],y[idx]]),
                    'obs_names':np.asarray(adata.obs_names[idx],dtype=str),'source_row':idx,
                    'source_file':source_file,'counts':totals[idx], 'coordinate_source':coord}
                if cfg.annotation_weight and contract.get('celltype_key'):
                    item['labels']=np.asarray(adata.obs[contract['celltype_key']].iloc[idx],dtype=str)
                for k in ['cell_id','CellID','source_obs_name']:
                    if k in adata.obs: item[k]=np.asarray(adata.obs[k].iloc[idx],dtype=str)
                slices.append(item)
        finally: adata.file.close()
        if _sha256(source_file)!=before: raise ValueError('Source mutated during read: '+source_file)
    by_id={s['slice_id']:s for s in slices}; slices=[by_id[s] for s in metrics.sort_values('slice_index').slice_id]
    result=register_series(slices,dataset_id,cfg)
    result['source_checksums']={p:s['sha256'] for p,s in sources.items()}
    if cfg.annotation_weight:result['annotation_sources']={p:s.get('celltype_key') for p,s in sources.items()}
    result['input_coordinate_sources']=sorted({s['coordinate_source'] for s in slices})
    result['input_coordinate_key'] = (result['input_coordinate_sources'][0].split(':',1)[1].split('[')[0] if len(result['input_coordinate_sources'])==1 else 'see input_coordinate_sources')
    result['display_points_per_slice']=display_points
    result['files']={}; payload['preregistration']=result
    for i,s in enumerate(slices):
        sid=s['slice_id']; info=result['slices'][sid]; t=np.asarray(info['matrix']); raw=s['raw_xy']; display=transform_points(raw,t)
        filename=f'slice_{i:04d}.npz'
        np.savez_compressed(destination/filename,raw_xy=raw,display_xy=display,obs_names=s['obs_names'],
            source_row=s['source_row'],**{k:s[k] for k in ['cell_id','CellID','source_obs_name','labels'] if k in s})
        info.update(source_file=s['source_file'],coordinate_source=s['coordinate_source'],cache_file=filename)
        result['files'][filename]=_sha256(destination/filename)
        idx=sample_indices(len(raw),display_points,cfg.seed,sid)
        finite=np.isfinite(raw[idx]).all(1); idx=idx[finite]
        payload['point_samples'][sid]={'x':raw[idx,0].tolist(),'y':raw[idx,1].tolist(),
            'display_x':display[idx,0].tolist(),'display_y':display[idx,1].tolist(),
            'counts':([None]*len(idx) if not sources[s['source_file']].get('expression_capture_available',True) else s['counts'][idx].tolist()),'depth':[0.5]*len(idx),
            'obs_names':s['obs_names'][idx].tolist(),'source_row':s['source_row'][idx].tolist(),
            'n_total':len(raw),'n_displayed':len(idx),'coordinate_source':s['coordinate_source'],
            'source_file':s['source_file'],'frame_id':result['frame_id']}
    result['total_elapsed_seconds']=time.perf_counter()-started
    (destination/'manifest.json').write_text(json.dumps(result,indent=2,ensure_ascii=False))
    (root/'slice_quality_display_payload.raw.json').write_bytes(payload_path.read_bytes())
    payload_path.write_text(json.dumps(payload,ensure_ascii=False))
    return result


def export_roi(input_dir, roi_path, output_csv):
    import csv
    import numpy as np
    from spateo.preprocessing.slice_quality import _sha256
    from spateo.preprocessing.slice_preregistration import roi_mask, inverse_roi_polygon
    root=Path(input_dir)/'display_preregistration'; manifest=json.loads((root/'manifest.json').read_text())
    roi=json.loads(Path(roi_path).read_text())
    if roi['frame_id']!=manifest['frame_id']: raise ValueError('ROI coordinate frame mismatch')
    if roi.get('coordinate_mode','display')!='display': raise ValueError('Only common display-frame ROI may be replayed')
    selected=[str(s) for s in (roi.get('slice_ids') or list(manifest['slices']))]; bounds=roi['bounds']; counts={}; polygons={}
    with Path(output_csv).open('x',newline='') as f:
        writer=csv.writer(f);writer.writerow(['slice_id','source_file','source_row','obs_name','raw_x','raw_y','display_x','display_y','frame_id','preregistration_status'])
        for sid in selected:
            info=manifest['slices'][sid]; path=root/info['cache_file']
            if _sha256(path)!=manifest['files'][info['cache_file']]: raise ValueError('Coordinate cache checksum mismatch')
            with np.load(path,allow_pickle=False) as data:
                raw=data['raw_xy'];display=data['display_xy']
                source_rows=data['source_row'];obs_names=data['obs_names']
                mask=roi_mask(raw,info['matrix'],bounds); counts[sid]=int(mask.sum())
                polygons[sid]=inverse_roi_polygon(info['matrix'],bounds).tolist()
                for j in np.flatnonzero(mask):
                    writer.writerow([sid,info['source_file'],int(source_rows[j]),obs_names[j],*raw[j],*display[j],manifest['frame_id'],info['status']])
    summary={'counts':counts,'raw_polygons':polygons,'frame_id':manifest['frame_id'],'bounds':bounds,
        'warning':manifest['warning'],'selection':'all input points, exact display rectangle test'}
    Path(str(output_csv)+'.json').write_text(json.dumps(summary,indent=2));return summary


def attach_preregistration_html(input_dir, output_html):
    root=Path(input_dir); payload_path=root/'slice_quality_display_payload.json'
    if not payload_path.exists(): return
    payload=json.loads(payload_path.read_text())
    if not payload.get('preregistration'): return
    path=Path(output_html); page=path.read_text()
    extension='<script id="preregistrationPayload" type="application/json">'+json.dumps({'preregistration':payload['preregistration']},ensure_ascii=False).replace('</','<\\/')+'</script>\n<script>'+Path(__file__).with_name('slice_quality_roi.js').read_text()+'</script>'
    path.write_text(page.replace('</body>',extension+'</body>'))
