from pathlib import Path
from pipeline_runtime import validate_input, save_pair, write_json, restore
from alignment_stage import save_cloud

def mapping(config, a, b, directory):
    import numpy as np
    import pandas as pd
    import spateo as st

    s, p, al = config["subset"], config["mapping"], config["alignment"]
    # AnnData expects category-ordered color arrays when slicing categoricals.
    # Real reconstruction inputs may instead store label->color dictionaries.
    for data in (a,b):
        for name,value in list(data.uns.items()):
            if name.endswith('_colors') and isinstance(value,dict) and name[:-7] in data.obs:
                column=data.obs[name[:-7]]
                if hasattr(column,'cat'):
                    data.uns.setdefault('spateo_4d_palettes',{})[name]=value.copy()
                    data.uns[name]=np.asarray([value.get(str(k),'#999999') for k in column.cat.categories],dtype=str)
    if s["group"] is None:
        raise ValueError("Select an observed annotation group for same-type mapping; run groups independently")
    if s["group"] is not None:
        pair = []
        for data in (a, b):
            labels = data.obs[s["annotation_key"]]
            if labels.isna().any():
                raise ValueError("Missing subset annotation")
            part = data[labels.astype(str) == str(s["group"])].copy()
            if not part.n_obs:
                raise ValueError("Requested subset has no cells in one stage")
            pair.append(part)
        a, b = pair
    common = [g for g in a.var_names if g in b.var_names]
    totals_a = dict(
        zip(a.var_names, np.asarray(a.layers[al["counts_layer"]].sum(0)).ravel())
    )
    totals_b = dict(
        zip(b.var_names, np.asarray(b.layers[al["counts_layer"]].sum(0)).ravel())
    )
    common = [g for g in common if totals_a[g] > 0 and totals_b[g] > 0]
    if not common:
        raise ValueError("No expressed shared genes in selected subsets")
    a, b = a[:, common].copy(), b[:, common].copy()
    if a.n_obs * b.n_obs > p["max_pairs"]:
        raise ValueError(
            "Mapping exceeds mapping.max_pairs; choose an explicit subset or increase the reviewed memory budget"
        )
    for data in (a, b):
        if np.any(np.asarray(data.layers[al["counts_layer"]].sum(1)).ravel() <= 0):
            raise ValueError("Gene harmonization left zero-library cells")
        st.pp.normalize_total(
            data,
            layer=al["counts_layer"],
            out_layer="normalized",
            target_sum=p["target_sum"],
            inplace=True,
        )
        st.pp.log1p_layer(
            data,
            layer="normalized",
            out_layer=al["log_layer"],
            set_X=True,
            inplace=True,
        )
        if s["annotation_key"] not in data.obs:
            data.obs[s["annotation_key"]] = "all cells (display group)"
    import warnings
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        _, pi = st.tdr.cell_directions(
            adataA=a,
            adataB=b,
            layer=al["log_layer"],
            spatial_key=al["aligned_key"],
            key_added=p["key"],
            alpha=p["alpha"],
            numItermax=p["numItermax"],
            numItermaxEmd=p["numItermaxEmd"],
            device=config["runtime"]["device"],
            inplace=True,
        )
    messages=[str(w.message) for w in caught]
    write_json(directory/'mapping_warnings.json',messages)
    if any('numItermax reached' in message for message in messages):
        raise RuntimeError('OT solver reached its iteration limit; increase numItermaxEmd before accepting mapping')
    row_error=float(np.max(np.abs(pi.sum(1)-1/a.n_obs)))
    col_error=float(np.max(np.abs(pi.sum(0)-1/b.n_obs)))
    if row_error > 1e-4/a.n_obs or col_error > 1e-4/b.n_obs:
        raise ValueError('Transport marginal constraints failed')
    write_json(directory/'mapping_qc.json', {'source_cells':a.n_obs,'target_cells':b.n_obs,
        'shared_genes':len(common),'row_marginal_max_error':row_error,'column_marginal_max_error':col_error,
        'total_mass':float(pi.sum()),'warnings':messages})
    x, v = a.obsm["X_" + p["key"]], a.obsm["V_" + p["key"]]
    if not np.isfinite(v).all() or not np.isfinite(pi).all():
        raise ValueError("Mapping produced nonfinite results")
    np.testing.assert_allclose(x - a.obsm[al["aligned_key"]], v)
    np.savez_compressed(
        directory / "transport.npz",
        pi=pi,
        source_ids=a.obs_names.to_numpy(dtype=str),
        target_ids=b.obs_names.to_numpy(dtype=str),
    )
    table = pd.DataFrame(
        {
            "cell_id": a.obs_names,
            "group": a.obs[s["annotation_key"]].astype(str).to_numpy(),
            "displacement": np.linalg.norm(v, axis=1),
        }
    )
    table.to_csv(directory / "mapping.csv", index=False)
    summary = table.groupby("group")["displacement"].agg(
        n="size",
        mean="mean",
        median="median",
        p05=lambda x: x.quantile(0.05),
        p95=lambda x: x.quantile(0.95),
    )
    summary.to_csv(directory / "mapping_summary.csv")
    return {
        **save_pair(a, b, directory),
        "stage1_pointcloud": save_cloud(a, al["aligned_key"], directory/"stage1.vtk", s["annotation_key"]),
        "stage2_pointcloud": save_cloud(b, al["aligned_key"], directory/"stage2.vtk", s["annotation_key"]),
        "mapping_qc": directory/"mapping_qc.json",
        "transport": directory / "transport.npz",
        "mapping_table": directory / "mapping.csv",
        "mapping_summary": directory / "mapping_summary.csv",
    }



def morphofield(config, a, b, directory):
    import numpy as np
    import spateo as st

    p = config["morphofield"].copy()
    key = p.pop("key")
    p["M"] = min(p["M"], a.n_obs)
    xyz = a.obsm[config["alignment"]["aligned_key"]]
    st.tdr.morphofield_sparsevfc(
        a,
        spatial_key=config["alignment"]["aligned_key"],
        V_key="V_" + config["mapping"]["key"],
        key_added=key,
        NX=xyz.copy(),
        **p,
    )
    if not np.isfinite(a.uns[key]["V"]).all():
        raise ValueError("Vector-field fit produced nonfinite values")
    raw = a.obsm["V_" + config["mapping"]["key"]]
    fitted = np.asarray(a.uns[key]["V"])
    cosine = np.sum(raw*fitted,1) / np.maximum(np.linalg.norm(raw,axis=1)*np.linalg.norm(fitted,axis=1), 1e-12)
    qc = {"implementation": a.uns[key].get("implementation"), "iterations": int(a.uns[key]["iteration"]),
          "raw_magnitude_quantiles": np.quantile(np.linalg.norm(raw,axis=1),[0,.5,.95,1]).tolist(),
          "fitted_magnitude_quantiles": np.quantile(np.linalg.norm(fitted,axis=1),[0,.5,.95,1]).tolist(),
          "median_cosine": float(np.median(cosine)), "rmse": float(np.sqrt(np.mean((raw-fitted)**2)))}
    a.obsm["velocity"] = fitted.copy()
    cloud = save_cloud(a, config["alignment"]["aligned_key"], directory/"field.vtk", config["subset"]["annotation_key"])
    write_json(directory/"field_qc.json", qc)
    return {**save_pair(a,b,directory), "field_qc": directory/"field_qc.json", "field_pointcloud": cloud}



def trajectory(config, a, b, directory):
    import numpy as np
    import spateo as st

    p = config["trajectory"]
    st.tdr.morphopath(
        a,
        vf_key=config["morphofield"]["key"],
        key_added=p["key"],
        t_end=p["t_end"],
        interpolation_num=p["interpolation_num"],
        direction=p["direction"],
        cores=1,
    )
    for values in a.uns[p["key"]]["prediction"].values():
        if not np.isfinite(values).all():
            raise ValueError("Trajectory contains nonfinite coordinates")
    fate = a.uns[p["key"]]
    paths = np.stack([np.asarray(fate['prediction'][i]) for i in range(a.n_obs)])
    if paths.shape[2] != 3: raise ValueError('Unexpected native trajectory orientation; expected cells x time x XYZ')
    times = np.stack([np.asarray(fate['t'][i]) for i in range(a.n_obs)])
    np.savez_compressed(directory/'trajectories.npz', coordinates=paths, times=times, cell_ids=a.obs_names.to_numpy(dtype=str))
    xyz = a.obsm[config['alignment']['aligned_key']]
    margin = .1*np.maximum(np.ptp(xyz,axis=0), 1e-9)
    outside = ((paths < xyz.min(0)-margin)|(paths > xyz.max(0)+margin)).any(2)
    write_json(directory/'trajectory_qc.json', {'n_paths': len(paths), 'model_time_end': p['t_end'],
        'fraction_points_outside_source_box_plus_10_percent': float(outside.mean()),
        'model_time_is_not_days': True})
    return {**save_pair(a,b,directory), 'trajectories': directory/'trajectories.npz', 'trajectory_qc': directory/'trajectory_qc.json'}


