"""Stage 3: native field features and gene association, with actual fitted curves."""
from pathlib import Path
from pipeline_runtime import save_pair, write_json


def metrics(config, a, b, directory):
    import numpy as np
    import pandas as pd
    import spateo as st
    from scipy import sparse
    from scipy.stats import spearmanr
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p = config['metrics']
    st.tdr.morphofield_velocity(a, vf_key=config['morphofield']['key'], key_added='velocity')
    a.obs['speed'] = np.linalg.norm(a.obsm['velocity'],axis=1)
    extras, summaries, curves = {}, {}, {}
    for i, axis in enumerate('xyz'): a.obs['velocity_'+axis] = a.obsm['velocity'][:,i]
    if any(k.startswith('jacobian_') for k in p['selected']):
        st.tdr.morphofield_jacobian(a, vf_key=config['morphofield']['key'], key_added='jacobian')
        J=np.asarray(a.uns['jacobian'])
        if J.shape != (3,3,a.n_obs) or not np.isfinite(J).all(): raise ValueError('Invalid Jacobian tensor')
        a.obs['jacobian_frobenius']=np.sqrt((J**2).sum(axis=(0,1)))
        for i,axis in enumerate('xyz'):
            for j,other in enumerate('xyz'): a.obs['jacobian_'+axis+other]=J[i,j,:]
        np.savez_compressed(directory/'jacobian.npz',jacobian=J,cell_ids=a.obs_names.to_numpy(dtype=str))
        extras['jacobian']=directory/'jacobian.npz'

    for metric in p['selected']:
        if metric != 'speed' and not metric.startswith(('jacobian_', 'velocity_')):
            getattr(st.tdr, 'morphofield_'+metric)(a, vf_key=config['morphofield']['key'], key_added=metric)
        if not np.isfinite(a.obs[metric].to_numpy()).all(): raise ValueError('Nonfinite '+metric)
        if metric not in p['glm_metrics']: continue
        genes = p['glm_genes']
        if genes == '*':
            counts = a.layers[config['alignment']['counts_layer']]
            genes = a.var_names[np.asarray((counts>0).sum(0)).ravel() >= p['glm_min_cells']].tolist()
        missing = set(genes)-set(a.var_names)
        if missing: raise ValueError('Unknown GLM genes: '+str(sorted(missing)))
        if a.obs[metric].nunique()<4: raise ValueError('Spline GLM needs four distinct feature values: '+metric)
        key='glm_degs_'+metric
        st.tl.glm_degs(a, layer='normalized', genes=genes, key_added=key,
            fullModelFormulaStr=f'~cr({metric}, df=3)', reducedModelFormulaStr='~1',
            qval_threshold=None, llf_threshold=None)
        result=a.uns[key]['glm_result'].copy()
        # Spateo removes failed genes from its table; restore them for a complete test ledger.
        failed=[g for g in genes if g not in result.index]
        if failed:
            failures=pd.DataFrame({'status':'fail','family':'NB2','pval':1.,'qval':1.},index=failed)
            result=pd.concat([result,failures])
        result['selected']=(result['status']=='ok') & (result['qval']<=p['qval_threshold'])
        if p['llf_threshold'] is not None: result['selected'] &= result['log-likelihood']<=p['llf_threshold']
        correlations=a.uns[key]['correlation']
        result['spearman_rho']=[float(spearmanr(correlations[g][metric], correlations[g]['expression']).statistic) if g in correlations and correlations[g]['expression'].nunique()>1 else np.nan for g in result.index]
        result.to_csv(directory/(key+'.csv'),index_label='gene')
        result[result.selected].to_csv(directory/(key+'.selected.csv'),index_label='gene')
        extras[key]=directory/(key+'.csv')
        extras[key+'_selected']=directory/(key+'.selected.csv')
        top=result[result.status=='ok'].head(p['glm_top_plots']).index.tolist()
        curves[metric]={}
        for j,g in enumerate(top):
            df=correlations[g].sort_values(metric)
            thin=np.linspace(0,len(df)-1,min(1500,len(df))).astype(int)
            shown=df.iloc[thin]
            curves[metric][g]={c: shown[c].astype(float).tolist() for c in [metric,'expression','mu','ci_lower','ci_upper']}
            curves[metric][g].update(qval=float(result.loc[g,'qval']),selected=bool(result.loc[g,'selected']))
            fig,ax=plt.subplots(figsize=(6,4))
            ax.scatter(df[metric],df.expression,s=3,alpha=.18,color='#64748b',rasterized=True)
            ax.plot(df[metric],df.mu,color='#0d9488',lw=2)
            ax.fill_between(df[metric].to_numpy(),df.ci_lower.to_numpy(),df.ci_upper.to_numpy(),alpha=.18,color='#0d9488')
            ax.set(xlabel=metric,ylabel='Size-factor normalized expression',title=f'{g} · q={result.loc[g,"qval"]:.3g}')
            fig.tight_layout(); path=directory/f'{metric}_gene_{j:02d}.png';fig.savefig(path,dpi=150);plt.close(fig)
            extras[f'{metric}_plot_{j}']=path
        summaries[metric]={'tested':len(genes),'failed':len(failed),'selected':int(result.selected.sum()),
            'formula':f'~cr({metric}, df=3)','layer':'normalized','NB2_dispersion_alpha':1.,
            'FDR_family':'all requested genes within this feature','not_causal':True}
        a.uns[key]={'glm_result':result, 'correlation':{g:correlations[g] for g in top}}
    a.obs[p['selected']].to_csv(directory/'metrics.csv',index_label='cell_id')
    write_json(directory/'glm_summary.json',summaries);write_json(directory/'glm_curves.json',curves)
    return {**save_pair(a,b,directory),**extras,'metrics_table':directory/'metrics.csv',
        'glm_summary':directory/'glm_summary.json','glm_curves':directory/'glm_curves.json'}


def gp(config, a, b, directory):
    import numpy as np
    import spateo as st

    p = config["gp"]
    missing = set(p["genes"]) - set(a.var_names)
    if missing:
        raise ValueError("Unknown GP genes: " + str(sorted(missing)))
    if set(p["genes"]) & set(a.obs.columns):
        raise ValueError("GP gene IDs conflict with observation field names")
    result = st.tdr.gp_interpolation(
        a,
        target_points=a.obsm[config["alignment"]["aligned_key"]].copy(),
        keys=p["genes"],
        spatial_key=config["alignment"]["aligned_key"],
        layer="normalized",
        training_iter=p["training_iter"],
        method=p["method"],
        inducing_num=min(p["inducing_num"], a.n_obs),
        device=config["runtime"]["device"],
        verbose=False,
    )
    if not np.isfinite(np.asarray(result.X)).all():
        raise ValueError("GP produced nonfinite expression")
    result.write_h5ad(directory / "gene_interpolation.h5ad")
    return {"gene_interpolation": directory / "gene_interpolation.h5ad"}

