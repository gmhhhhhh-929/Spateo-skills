"""Verify the actual full-cell server run without moving expression data off server."""
import sys,json
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from pipeline_runtime import read_json,load_pair,valid_outputs,write_json
m=read_json(sys.argv[1]);out=Path(sys.argv[2]);assert m['status']=='completed'
assert all(valid_outputs(v) for v in m['stages'].values() if v['status']=='completed')
a,b=load_pair(m['stages']['alignment']['outputs'])
for k,d in zip(('stage1','stage2'),(a,b)):
 source=ad.read_h5ad(m['config']['inputs'][k]);assert d.obs_names.equals(source.obs_names)
 np.testing.assert_array_equal(d.obsm['spatial_3d'],source.obsm['spatial_3d'])
 assert (d.layers['counts']-source.X).nnz==0
 del source
mapped,target=load_pair(m['stages']['mapping']['outputs'])
assert mapped.obs.anno.astype(str).eq('CNS').all() and target.obs.anno.astype(str).eq('CNS').all()
np.testing.assert_allclose(mapped.obsm['X_cells_mapping']-mapped.obsm['aligned'],mapped.obsm['V_cells_mapping'])
tr=np.load(m['stages']['trajectory']['outputs']['trajectories']['path'])
np.testing.assert_array_equal(tr['cell_ids'],mapped.obs_names)
np.testing.assert_allclose(tr['coordinates'][:,0,:],mapped.obsm['aligned'])
assert np.isfinite(tr['coordinates']).all() and np.all(np.diff(tr['times'],axis=1)>0)
features,_=load_pair(m['stages']['metrics']['outputs']);assert features.obs_names.equals(mapped.obs_names)
np.testing.assert_allclose(np.trace(features.uns['jacobian'],axis1=0,axis2=1),features.obs.divergence)
assert np.isfinite(features.obs[m['config']['metrics']['selected']].to_numpy()).all()
for feature in m['config']['metrics']['glm_metrics']:
 table=pd.read_csv(m['stages']['metrics']['outputs']['glm_degs_'+feature]['path'])
 assert len(table)==48 and table.status.eq('ok').all() and table.qval.between(0,1).all()
result={'passed':True,'run_id':m['run_id'],'checks':['all output SHA256 hashes','full-cell identity preservation','unchanged original coordinates','unchanged original counts','same-annotation mapping','mapped endpoint/vector equality','trajectory seed identity and starting position','trajectory finiteness and monotonic model time','Jacobian trace equals divergence','finite scalar metrics','336 successful GLM fits with valid q values'], 'full_cells':[a.n_obs,b.n_obs],'mapped_cells':[mapped.n_obs,target.n_obs],'trajectory_shape':list(tr['coordinates'].shape),'glm_tests':336}
write_json(out,result);print(json.dumps(result))
