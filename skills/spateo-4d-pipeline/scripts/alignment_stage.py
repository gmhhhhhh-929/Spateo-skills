"""Stage 1: identity-checked point clouds, shared-frame registration and H5AD export."""
from pathlib import Path
from pipeline_runtime import validate_input, save_pair, write_json, restore


def attach_pointcloud(data, filename, key):
    import numpy as np
    import spateo as st
    pc = st.tdr.read_model(filename)
    if 'obs_index' not in pc.point_data:
        raise ValueError('VTK must contain obs_index; row order is not an identity contract')
    ids = np.asarray(pc.point_data['obs_index']).astype(str)
    if len(ids) != len(set(ids)) or set(ids) != set(data.obs_names):
        raise ValueError('VTK obs_index must be a unique exact match to paired H5AD obs_names')
    order = {cell: i for i, cell in enumerate(ids)}
    xyz = np.asarray(pc.points)[[order[cell] for cell in data.obs_names]]
    if xyz.shape != (data.n_obs, 3) or not np.isfinite(xyz).all():
        raise ValueError('Invalid VTK point coordinates')
    if key in data.obsm:
        preserved = key + '_before_vtk'
        if preserved in data.obsm: raise ValueError('Refusing to overwrite preserved pre-VTK coordinates')
        data.obsm[preserved] = data.obsm[key].copy()
    data.obsm[key] = xyz.copy()


def save_cloud(data, key, path, annotation):
    import numpy as np
    import spateo as st
    kwargs = {'groupby': annotation} if annotation in data.obs else {}
    pc, _ = st.tdr.construct_pc(data, spatial_key=key, **kwargs)
    pc.point_data['obs_index'] = data.obs_names.to_numpy(dtype=str)
    if annotation in data.obs: pc.point_data[annotation] = data.obs[annotation].astype(str).to_numpy()
    for name, values in data.obsm.items():
        if name.startswith(('V_', 'X_')) or name == 'velocity':
            if values.shape == (data.n_obs,3): pc.point_data[name] = np.asarray(values)
    st.tdr.save_model(pc, str(path))
    restored = st.tdr.read_model(str(path))
    np.testing.assert_allclose(restored.points, data.obsm[key])
    np.testing.assert_array_equal(restored.point_data['obs_index'].astype(str), data.obs_names)
    return path


def alignment(config, directory):
    import anndata as ad
    import numpy as np
    import spateo as st
    from scipy.spatial import cKDTree
    p, w = config['alignment'], config['workflow']
    a, b = (ad.read_h5ad(config['inputs'][k]) for k in ('stage1', 'stage2'))
    for i, data in enumerate((a, b), 1):
        data.uns = restore(dict(data.uns))
        pc = config['inputs'].get('pointcloud' + str(i))
        if pc: attach_pointcloud(data, pc, p['spatial_key'] if w['entry']=='alignment' else p['aligned_key'])
        validate_input(data, config)
    if w['entry'] == 'alignment':
        for data in (a, b):
            st.pp.normalize_total(data, layer=p['counts_layer'], out_layer='normalized', target_sum=p['target_sum'], inplace=True)
            st.pp.log1p_layer(data, layer='normalized', out_layer=p['log_layer'], set_X=False, inplace=True)
        common = [g for g in a.var_names if g in set(b.var_names)]
        if not common: raise ValueError('No common genes')
        # Explicit reference subsets avoid copying both full expression matrices
        # merely to select reference cells. Native sampling semantics are kept.
        from spateo._native import sample
        count = min(p['n_sampling'], a.n_obs, b.n_obs)
        references = [data[sample(arr=np.asarray(data.obs_names), n=count,
                                  method=p['sampling_method'], X=data.obsm[p['spatial_key']])].copy()
                      for data in (a,b)]
        annotation = config['subset']['annotation_key']
        if p['use_annotation'] and not all(annotation in d.obs for d in (a,b)):
            raise ValueError('Annotation-informed alignment requires labels in both inputs')
        rep_layer = [p['log_layer'], annotation] if p['use_annotation'] else p['log_layer']
        rep_field = ['layer','obs'] if p['use_annotation'] else 'layer'
        aligned, _, _, _ = st.align.morpho_align_ref(
            models=[a, b], models_ref=references, rep_layer=rep_layer, rep_field=rep_field, spatial_key=p['spatial_key'],
            key_added=p['aligned_key'], genes=common, mode=p['mode'],
            n_sampling=min(p['n_sampling'], a.n_obs, b.n_obs), sampling_method=p['sampling_method'],
            max_iter=p['max_iter'], nonrigid_start_iter=p['nonrigid_start_iter'],
            nn_init=p['nn_init'], iter_key_added=None,
            device=config['runtime']['device'], verbose=False)
    else:
        aligned = [a, b]
        for data in aligned:
            meta=data.uns.get('spateo_4d_frame',{})
            if meta and (meta.get('id') != w['frame_id'] or meta.get('unit') != config['inputs']['coordinate_unit']):
                raise ValueError('Imported frame metadata disagrees with the declared frame_id or unit')
        for data in aligned:
            if p['aligned_key'] not in data.obsm:
                raise ValueError('Imported pair lacks the declared aligned_key')
        if w['entry'] == 'field':
            vf = a.uns.get(config['morphofield']['key'], {})
            if vf.get('method') != 'sparsevfc' or vf.get('implementation') != 'spateo-native-rbf':
                raise ValueError('Field entry requires a compatible native fitted field')
            np.testing.assert_allclose(vf['X'], a.obsm[p['aligned_key']])
            if 'normalized' not in a.layers: raise ValueError('Field entry needs normalized expression for GLMs')
    qc = {'entry': w['entry'], 'coordinate_unit': config['inputs']['coordinate_unit'], 'frame_id': w['frame_id'] or directory.parent.parent.name,
          'registration_mode': p['mode'], 'nearest_neighbor_is_not_anatomical_validation': True}
    outputs = {}
    for name, original, data in zip(('stage1','stage2'), (a,b), aligned):
        np.testing.assert_array_equal(data.obsm[p['spatial_key']], original.obsm[p['spatial_key']])
        coords = np.asarray(data.obsm[p['aligned_key']])
        if coords.shape != (data.n_obs, 3) or not np.isfinite(coords).all(): raise ValueError('Invalid aligned coordinates')
        data.uns['spateo_4d_frame'] = {'id': qc['frame_id'], 'unit': qc['coordinate_unit'], 'key': p['aligned_key']}
        qc[name] = {'n_cells': data.n_obs, 'min': coords.min(0).tolist(), 'max': coords.max(0).tolist()}
        outputs[name+'_pointcloud'] = save_cloud(data, p['aligned_key'], directory/(name+'.vtk'), config['subset']['annotation_key'])
    for key, title in [(p['spatial_key'],'before'), (p['aligned_key'],'after')]:
        x,y = [np.asarray(d.obsm[key]) for d in aligned]
        da,db=cKDTree(y).query(x)[0],cKDTree(x).query(y)[0]
        qc[title]={'symmetric_nn_mean': float((da.mean()+db.mean())/2), 'source_to_target_p95': float(np.quantile(da,.95))}
    from registration_qc import pair_qc
    annotation = config['subset']['annotation_key']
    if all(annotation in d.obs for d in aligned):
        qc['annotation_qc'] = {
            title: pair_qc(aligned[0].obsm[key], aligned[1].obsm[key],
                           aligned[0].obs[annotation], aligned[1].obs[annotation])
            for key, title in [(p['spatial_key'], 'before'), (p['aligned_key'], 'after')]
        }
    qc['anatomical_status'] = 'requires_review'
    qc['analysis_coordinate_key'] = p['aligned_key']
    qc['nn_init'] = p['nn_init']
    qc['use_annotation'] = p['use_annotation']
    write_json(directory/'qc.json', qc)
    return {**save_pair(*aligned, directory), **outputs, 'qc': directory/'qc.json'}
