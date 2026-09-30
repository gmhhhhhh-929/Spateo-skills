"""A coordinate-identifiable case verifies native G_init, IDs and conservation."""
import copy
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import pipeline_runtime as rt
from morphogenesis_stage import mapping


def test_spatial_initialization_preserves_correspondence_under_row_permutation(tmp_path):
    rng = np.random.default_rng(814)
    xyz = rng.normal(size=(20, 3)) * [5, 20, 3]
    a = ad.AnnData(sparse.csr_matrix(np.ones((20, 8))),
                   obs=pd.DataFrame({'anno': pd.Categorical(['CNS'] * 20)},
                                    index=[f'source-{i}' for i in range(20)]))
    a.layers['counts'] = a.X.copy()
    a.obsm['aligned'] = xyz.copy()
    b = a[rng.permutation(20)].copy()
    b.obs_names = ['target-' + i.removeprefix('source-') for i in b.obs_names]
    c = copy.deepcopy(rt.DEFAULTS)
    c['subset']['group'] = 'CNS'
    c['mapping']['initialization'] = 'aligned_spatial'
    outputs = mapping(c, a, b, tmp_path)
    saved = ad.read_h5ad(outputs['stage1_h5ad'])
    np.testing.assert_allclose(saved.obsm['X_cells_mapping'], xyz, atol=1e-6)
    np.testing.assert_allclose(saved.obsm['V_cells_mapping'], 0, atol=1e-6)
    for key in ['initial_transport', 'transport']:
        archive = np.load(outputs[key])
        np.testing.assert_allclose(archive['pi'].sum(0), 1/20, atol=1e-7)
        np.testing.assert_allclose(archive['pi'].sum(1), 1/20, atol=1e-7)
        np.testing.assert_array_equal(archive['source_ids'], a.obs_names)
        np.testing.assert_array_equal(archive['target_ids'], b.obs_names)
    qc = json.loads(outputs['mapping_qc'].read_text())
    assert qc['initialization']['method'] == 'aligned_spatial'
    assert qc['initialization']['transport_cost'] == 0


def test_mapping_initialization_is_validated(tmp_path):
    c = copy.deepcopy(rt.DEFAULTS)
    for name in ['a.h5ad', 'b.h5ad']:
        (tmp_path / name).touch()
    c['inputs'].update(stage1=str(tmp_path/'a.h5ad'), stage2=str(tmp_path/'b.h5ad'), coordinate_unit='um')
    c['mapping']['initialization'] = 'unrecognized'
    p = tmp_path/'bad.json'
    p.write_text(json.dumps(c))
    with pytest.raises(ValueError, match='mapping.initialization'):
        rt.canonical_config(p)
