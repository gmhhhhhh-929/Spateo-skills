import json
from pathlib import Path
import numpy as np
import anndata as ad
root=Path(__file__).resolve().parents[3]/'refinement-20260928'
legacy=np.load(root/'legacy_field_500/field.npz')
m=json.loads((root/'analysis/runs/local_anatomy_refined_v1/manifest.json').read_text())
a=ad.read_h5ad(m['stages']['morphofield']['outputs']['stage1_h5ad']['path'])
f=a.uns['VecFld_morpho'];new=np.asarray(f['V']);old=legacy['V_fit']
np.testing.assert_array_equal(a.obs_names.to_numpy(dtype=str),legacy['cell_ids'])
np.testing.assert_allclose(a.obsm['aligned'],legacy['X'],rtol=0,atol=0)
np.testing.assert_allclose(a.obsm['V_cells_mapping'],legacy['V_raw'],rtol=0,atol=0)
norm=lambda v:np.linalg.norm(v,axis=1)
r={'same_input_ids_coordinates_displacements':True,'native_beta':float(f['beta']),
   'native_median_magnitude':float(np.median(norm(new))),'legacy_median_magnitude':float(np.median(norm(old))),
   'median_pointwise_magnitude_ratio':float(np.median(norm(new)/norm(old))),
   'median_cosine_old_to_native':float(np.median(np.sum(new*old,axis=1)/(norm(new)*norm(old)))),
   'magnitude_ratio_quantiles':np.quantile(norm(new)/norm(old),[0,.05,.5,.95,1]).tolist(),
   'legacy_solver':'Dynamo SparseVFC Gaussian-uniform mixture EM; velocity-weighted controls; lambda*sigma2*K regularization',
   'native_solver':'spateo-native-rbf; spatial TRN controls; Huber IRLS; lambda*I ridge',
   'interpretation':'This dataset-specific comparison is not a universal conversion factor. Legacy direction-only acceptance does not verify vector magnitudes. Neither result is a calibrated physical velocity.'}
(root/'field_comparison.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
