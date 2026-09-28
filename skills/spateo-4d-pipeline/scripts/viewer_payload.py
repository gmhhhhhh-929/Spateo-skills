"""Stage 4: portable capability-based viewer; scientific values are never display-scaled."""
from pathlib import Path
import json


def build_payload(config, manifest):
    import numpy as np
    import anndata as ad
    from pipeline_runtime import restore, valid_outputs
    stages=manifest['stages']; p=config['dashboard']; key=config['alignment']['aligned_key']
    def active(stage):
        return stages[stage]['status'] in ('completed','reused','imported')
    def load(stage):
        record=stages[stage]
        if not valid_outputs(record): raise ValueError('Viewer checkpoint failed hash verification: '+stage)
        pair=[ad.read_h5ad(record['outputs'][k+'_h5ad']['path']) for k in ('stage1','stage2')]
        for d in pair: d.uns=restore(dict(d.uns))
        return pair
    def cloud(a, cap, before=False):
        idx=np.arange(a.n_obs)
        if cap and len(idx)>cap: idx=np.sort(np.random.default_rng(config['runtime']['seed']).choice(idx,cap,replace=False))
        annotation=config['subset']['annotation_key']
        d={'ids':a.obs_names[idx].tolist(),'xyz':np.asarray(a.obsm[key])[idx].tolist(),
           'groups':a.obs[annotation].astype(str).iloc[idx].tolist() if annotation in a.obs else ['cells']*len(idx),
           'total':a.n_obs,'features':{}}
        if before: d['before']=np.asarray(a.obsm[config['alignment']['spatial_key']])[idx].tolist()
        for name in ['speed',*config['metrics']['selected']]:
            if name in a.obs: d['features'][name]=a.obs[name].to_numpy()[idx].astype(float).tolist()
        mk=config['mapping']['key']
        if 'V_'+mk in a.obsm:
            v=np.asarray(a.obsm['V_'+mk])[idx]
            d.update(rawVectors=v.tolist(),endpoints=np.asarray(a.obsm['X_'+mk])[idx].tolist())
            d['features']['displacement']=np.linalg.norm(v,axis=1).tolist()
        vf=a.uns.get(config['morphofield']['key'],{})
        if vf.get('implementation') == 'spateo-native-rbf':
            np.testing.assert_allclose(vf['X'],a.obsm[key])
            v=np.asarray(vf['V'])[idx]; d['fieldVectors']=v.tolist();d['features']['speed']=np.linalg.norm(v,axis=1).tolist()
        return d
    base=load('alignment')
    payload={'schema':'spateo-4d-viewer/v1','runId':manifest['run_id'],'runStatus':'completed',
        'labels':config['labels'],'unit':config['inputs']['coordinate_unit'],'subset':config['subset'],
        'stages':{k:('completed' if k=='dashboard' else v['status']) for k,v in stages.items()},
        'alignment':[cloud(base[0],p['max_points'],True),cloud(base[1],p['max_target_points'],True)],
        'analysis':None,'trajectories':None,'glm':{},'tables':{},'qc':{},
        'limits':{'vectors':p['max_vectors'],'trajectories':p['max_trajectories']},
        'defaultFeature':p['default_feature'],'modelTime':config['trajectory']['t_end'],'modelDirection':config['trajectory']['direction']}
    for stage in ('alignment','mapping','morphofield','trajectory','metrics'):
        for label in ('qc','mapping_qc','field_qc','trajectory_qc','glm_summary'):
            out=stages[stage].get('outputs',{}).get(label)
            if out: payload['qc'][label]=json.loads(Path(out['path']).read_text(encoding='utf-8'))
    latest=next((s for s in ('metrics','morphofield','mapping') if active(s)),None)
    if latest:
        a,b=load(latest); payload['analysis']=[cloud(a,p['max_points']),cloud(b,p['max_target_points'])]
    if active('trajectory'):
        a,_=load('trajectory'); fate=a.uns[config['trajectory']['key']]
        ids=np.linspace(0,a.n_obs-1,min(a.n_obs,p['max_trajectories'])).astype(int)
        payload['trajectories']={'total':a.n_obs,'ids':a.obs_names[ids].tolist(),
            'paths':[np.asarray(fate['prediction'][str(i)] if str(i) in fate['prediction'] else fate['prediction'][i]).tolist() for i in ids],
            'times':np.asarray(next(iter(fate['t'].values()))).tolist()}
    if active('metrics'):
        outputs=stages['metrics']['outputs']
        if 'glm_curves' in outputs: payload['glm']=json.loads(Path(outputs['glm_curves']['path']).read_text(encoding='utf-8'))
        for k,v in outputs.items():
            if k.startswith('glm_degs_') and not k.endswith('_selected'):
                payload['tables'][k.removeprefix('glm_degs_')]=Path(v['path']).read_text(encoding='utf-8')
    # Coordinate-only diagnostics can also be added to historical verified runs.
    # They do not alter the registration or imply anatomical correctness.
    from registration_qc import pair_qc
    source, target = payload['alignment']
    payload['registrationReview'] = {
        name: pair_qc(source[name], target[name], source['groups'], target['groups'])
        for name in ('before', 'xyz')
    }
    payload['registrationReview']['scope'] = 'exported points; inspect counts if display sampling is enabled'
    payload['capabilities']=['alignment']
    if payload['analysis']: payload['capabilities'].append('mapping')
    if latest in ('metrics','morphofield'):payload['capabilities'].append('field')
    if payload['trajectories']:payload['capabilities'].append('trajectory')
    if active('metrics'):payload['capabilities'].append('features')
    return payload


def render_html(payload):
    from plotly.offline import get_plotlyjs
    template=Path(__file__).resolve().parents[1]/'assets/viewer.html'
    data=json.dumps(payload,ensure_ascii=False,separators=(',',':'),allow_nan=False).replace('<','\\u003c')
    return template.read_text(encoding='utf-8').replace('__PLOTLY__',get_plotlyjs()).replace('__PAYLOAD__',data)


def attach_comparison(payload, label, directory):
    """Attach a hashed benchmark after joining IDs and verifying original coordinates."""
    import hashlib
    import numpy as np
    from registration_qc import pair_qc
    directory=Path(directory)
    record=json.loads((directory/'comparison.json').read_text())
    path=directory/'coordinates.npz'
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if record.get('status')!='completed' or digest!=record.get('coordinates_sha256'):
        raise ValueError('Comparison checkpoint is incomplete or has changed')
    xyz=[]
    with np.load(path,allow_pickle=False) as data:
        for k,cloud in enumerate(payload['alignment'],1):
            ids=data[f'stage{k}_ids'].tolist()
            if len(ids)!=len(set(ids)) or len(cloud['ids'])!=len(set(cloud['ids'])):
                raise ValueError('Duplicate cell identity in comparison')
            lookup={cell:i for i,cell in enumerate(ids)}
            if not set(cloud['ids']) <= set(ids): raise ValueError('Comparison cell IDs do not match')
            idx=np.array([lookup[cell] for cell in cloud['ids']])
            np.testing.assert_allclose(data[f'stage{k}_before'][idx],cloud['before'],rtol=0,atol=1e-5)
            if data[f'stage{k}_labels'][idx].tolist()!=cloud['groups']:
                raise ValueError('Comparison annotations do not match')
            values=data[f'stage{k}_benchmark_rigid'][idx]
            if not np.isfinite(values).all(): raise ValueError('Nonfinite comparison coordinates')
            xyz.append(values.tolist())
    comparisons=payload.setdefault('registrationComparisons',[])
    comparisons.append({'id':'comparison_'+str(len(comparisons)), 'label':label, 'coordinates':xyz,
        'qc':{'xyz':pair_qc(*xyz,*[c['groups'] for c in payload['alignment']])},
        'provenance':{'coordinates_sha256':digest,'inputs':record['inputs'],
                      'parameters':record['parameters'],'spateo_root':record['spateo_root']}})
    return payload


def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--payload',required=True);p.add_argument('--output',required=True)
    p.add_argument('--comparison',nargs=2,action='append',metavar=('LABEL','BENCHMARK_DIR'),default=[])
    args=p.parse_args(); path=Path(args.output)
    if path.exists():raise FileExistsError(path)
    payload=json.loads(Path(args.payload).read_text(encoding='utf-8'))
    for label,directory in args.comparison: attach_comparison(payload,label,directory)
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(render_html(payload), encoding='utf-8')

if __name__=='__main__':main()
