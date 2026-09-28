"""Fit legacy SparseVFC on the exact corrected native mapping input."""
import argparse,json,hashlib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--max-iter',type=int,default=500);args=p.parse_args()
import anndata as ad,numpy as np,spateo as st
out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
a=ad.read_h5ad(args.input);a.uns.clear();X=np.asarray(a.obsm['aligned'],float);V=np.asarray(a.obsm['V_cells_mapping'],float)
st.tdr.morphofield_sparsevfc(a,spatial_key='aligned',V_key='V_cells_mapping',key_added='compare',NX=X.copy(),M=100,lambda_=0.02,restart_num=1,restart_seed=[0],MaxIter=args.max_iter,inplace=True)
f=a.uns['compare'];fitted=np.asarray(f['V'])
assert fitted.shape==V.shape and np.isfinite(fitted).all()
np.savez_compressed(out/'field.npz',cell_ids=a.obs_names.to_numpy(dtype=str),X=X,V_raw=V,V_fit=fitted)
info={'spateo':st.__file__,'input_sha256':hashlib.sha256(Path(args.input).read_bytes()).hexdigest(),'parameters':{'M':100,'lambda_':.02,'MaxIter':args.max_iter,'restart_num':1,'restart_seed':[0]},'beta':float(f['beta']),'iteration':int(f['iteration']),'raw_magnitude_quantiles':np.quantile(np.linalg.norm(V,axis=1),[0,.5,.95,1]).tolist(),'fitted_magnitude_quantiles':np.quantile(np.linalg.norm(fitted,axis=1),[0,.5,.95,1]).tolist(),'median_cosine_to_raw':float(np.median(np.sum(V*fitted,axis=1)/np.maximum(1e-12,np.linalg.norm(V,axis=1)*np.linalg.norm(fitted,axis=1))))}
(out/'summary.json').write_text(json.dumps(info,indent=2))
print(info)
