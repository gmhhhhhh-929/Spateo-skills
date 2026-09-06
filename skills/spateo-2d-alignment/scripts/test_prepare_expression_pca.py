"""Synthetic shared-PCA and access-contract tests; no biological input files."""
import argparse
from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

import prepare_expression_pca as prep


class TestPrepareExpressionPCA(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.slice_dir = self.root/'perturbed_slices'
        self.slice_dir.mkdir()
        self.source = self.root/'expression_source.h5ad'
        self.rng = np.random.default_rng(103)
        self.ids = np.array([f'cell-{i:03d}' for i in range(48)])
        self.genes = np.array([f'gene-{i:02d}' for i in range(12)])
        self.counts = self.rng.poisson(3,size=(48,12)).astype(np.float64)
        self.counts[:14,:4] += 7
        self.counts[14:28,4:8] += 8
        self.counts[28:42,8:] += 9
        self.wanted = self.ids[:42]
        self.write_source(self.counts)
        self.input_paths = []
        for i in range(3):
            ids = self.ids[i*14:(i+1)*14]
            obs = pd.DataFrame({'z':np.full(14,0.7+i*1.3)},index=ids)
            obj = ad.AnnData(X=sparse.csr_matrix(np.ones((14,2))),obs=obs)
            obj.obsm['spatial'] = self.rng.normal(size=(14,2))+i
            obj.obsm['X_pca'] = np.tile([1.,0.],(14,1))
            obj.obsm['forbidden_original_xyz'] = self.rng.normal(size=(14,3))
            path = self.slice_dir/f'input_{i}.h5ad'
            obj.write_h5ad(path)
            self.input_paths.append(path)

    def tearDown(self):
        self.temp.cleanup()

    def write_source(self,values,order=None,storage='csr',path=None):
        order = np.arange(48) if order is None else order
        values = values[order]
        matrix = sparse.csr_matrix(values) if storage=='csr' else sparse.csc_matrix(values) if storage=='csc' else values
        # Source X has real gene dimension as required by AnnData. Deliberately
        # invalid old PCA/reference arrays test that no source obsm is read.
        obj = ad.AnnData(X=sparse.csr_matrix(np.ones((48,12))),obs=pd.DataFrame(index=self.ids[order]),
                         var=pd.DataFrame(index=self.genes))
        obj.layers['counts_X'] = matrix
        obj.layers['counts'] = sparse.csr_matrix(values+100)
        obj.obsm['X_pca'] = np.tile([1.,0.],(48,1))
        obj.obsm['spatial_3d'] = np.arange(48*3).reshape(48,3).astype(float)
        obj.obs['annotation_should_not_be_read'] = ['hidden']*48
        obj.write_h5ad(self.source if path is None else path)

    def args(self,name='result',**overrides):
        arguments = dict(slice_dir=str(self.slice_dir),expression_source=str(self.source),
                         expression_layer='auto',matrix_state='counts',output_dir=str(self.root/name),
                         n_hvg=2000,n_pcs=50,seed=71,annotation_key=None,physical_z_key='physical_z')
        arguments.update(overrides)
        return argparse.Namespace(**arguments)

    def execute(self,args):
        with redirect_stdout(io.StringIO()):
            return prep.prepare(args)

    def scores_by_id(self,output):
        results = {}
        for path in sorted((Path(output)/'slice_h5ad').glob('*.h5ad')):
            obj = ad.read_h5ad(path)
            results.update(zip(obj.obs_names,np.asarray(obj.obsm['X_pca'])))
        return np.vstack([results[cell] for cell in self.wanted])

    def test_shared_joint_basis_no_annotation_and_dimension_cap(self):
        args = self.args()
        manifest = self.execute(args)
        self.assertEqual(manifest['schema_version'],'expression-pca-v1')
        self.assertTrue(manifest['shared_basis'])
        self.assertFalse(manifest['fitted_on_spatial'])
        self.assertFalse(manifest['fitted_on_annotation'])
        self.assertEqual(manifest['effective_n_pcs'],11)
        self.assertEqual(manifest['effective_n_hvg'],12)
        self.assertEqual(manifest['expression_source']['selected_layer'],'layers/counts_X')
        with np.load(Path(args.output_dir)/'basis.npz') as basis:
            np.testing.assert_array_equal(basis['fit_cell_ids'],self.wanted)
            means,loadings = basis['gene_means'],basis['loadings']
            for record in manifest['slices']:
                obj = ad.read_h5ad(Path(args.output_dir)/record['path'])
                self.assertEqual(set(obj.obsm),{'spatial','X_pca'})
                self.assertNotIn('anno',obj.obs)
                reconstructed = obj.X@loadings-(means@loadings)[None,:]
                np.testing.assert_allclose(reconstructed,obj.obsm['X_pca'],atol=2e-6)
                self.assertEqual(record['basis_sha256'],manifest['basis_sha256'])
                self.assertEqual(record['output_h5ad_sha256'],prep.hash_file(Path(args.output_dir)/record['path']))
                self.assertEqual(record['features_sha256'],prep.hash_numeric(obj.obsm['X_pca']))
                self.assertEqual(record['cell_ids_sha256'],prep.hash_ids(obj.obs_names))
        scores = self.scores_by_id(args.output_dir)
        np.testing.assert_allclose(scores.mean(axis=0),0,atol=2e-7)
        self.assertGreater(np.linalg.norm(scores[:14].mean(axis=0)),0.1)

    def test_source_coordinates_old_pca_and_slice_reference_never_read(self):
        actual_getitem = h5py.Dataset.__getitem__
        actual_hash = prep.hash_file
        source = self.source.resolve()
        slices = set(p.resolve() for p in self.input_paths)
        def guarded_getitem(node,selection):
            filename = Path(node.file.filename).resolve()
            if filename==source:
                self.assertTrue(node.name.startswith(('/layers/counts_X/','/obs/_index','/var/_index')),node.name)
            elif filename in slices:
                self.assertIn(node.name,('/obs/_index','/obs/z','/obsm/spatial'))
            return actual_getitem(node,selection)
        def guarded_hash(path):
            self.assertNotIn(Path(path).resolve(),slices|{source})
            return actual_hash(path)
        with mock.patch.object(h5py.Dataset,'__getitem__',new=guarded_getitem),mock.patch.object(prep,'hash_file',side_effect=guarded_hash):
            manifest = self.execute(self.args())
        self.assertFalse(manifest['source_coordinates_read'])
        self.assertFalse(manifest['preexisting_pca_read'])
        self.assertTrue(all('/obsm/' not in r['dataset'] for r in manifest['read_access_ledger'] if r['role']=='expression_source'))

    def test_source_row_order_mapping_and_superset_cells(self):
        first = self.args('first')
        self.execute(first)
        reordered = self.root/'reordered_expression.h5ad'
        changed = self.counts.copy()
        changed[42:] *= 100000
        self.write_source(changed,order=self.rng.permutation(48),path=reordered)
        second = self.args('second',expression_source=str(reordered))
        self.execute(second)
        np.testing.assert_array_equal(self.scores_by_id(first.output_dir),self.scores_by_id(second.output_dir))

    def test_dense_and_csc_storage(self):
        reference = self.args('reference')
        self.execute(reference)
        expected = self.scores_by_id(reference.output_dir)
        for storage in ('dense','csc'):
            path = self.root/(storage+'.h5ad')
            self.write_source(self.counts,storage=storage,path=path)
            args = self.args(storage,expression_source=str(path))
            self.execute(args)
            np.testing.assert_allclose(self.scores_by_id(args.output_dir),expected,atol=2e-5)

    def test_inputs_unchanged(self):
        paths = [self.source]+self.input_paths
        before = {path:hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
        self.execute(self.args())
        after = {path:hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
        self.assertEqual(before,after)

    def test_noncounts_not_guessed_and_declared_log1p_not_logged_twice(self):
        logvalues = np.log1p(self.counts/self.counts.sum(axis=1)[:,None]*10000)
        self.write_source(logvalues)
        with self.assertRaisesRegex(ValueError,'integer-like'):
            self.execute(self.args('wrong_counts'))
        self.assertFalse((self.root/'wrong_counts').exists())
        args = self.args('declared_log',matrix_state='log1p')
        manifest = self.execute(args)
        self.assertFalse(manifest['pca']['normalization_applied'])
        self.assertFalse(manifest['pca']['log1p_applied'])
        for record in manifest['slices']:
            obj = ad.read_h5ad(Path(args.output_dir)/record['path'])
            rows = [int(cell.split('-')[1]) for cell in obj.obs_names]
            cols = [int(gene.split('-')[1]) for gene in obj.var_names]
            np.testing.assert_allclose(obj.X.toarray(),logvalues[np.ix_(rows,cols)],atol=1e-6)

    def test_constant_onehot_nonfinite_and_zero_cells_rejected(self):
        invalid = [np.ones((48,12)),np.eye(12)[np.arange(48)%12]]
        for values in invalid:
            self.write_source(values)
            with self.assertRaises(ValueError):
                self.execute(self.args())
        for bad in ('nan','zero','negative'):
            values = self.counts.copy()
            if bad=='nan':values[0,0]=np.nan
            elif bad=='zero':values[0]=0
            else:values[0,0]=-1
            self.write_source(values)
            with self.assertRaises(ValueError):
                self.execute(self.args())

    def test_optional_annotation_is_copied_but_never_fit(self):
        for path in self.input_paths:
            with h5py.File(path,'r+') as handle:
                handle['obs'].create_dataset('arbitrary_annotation',data=np.array(['group']*14,dtype=h5py.string_dtype()))
        first = self.args('without')
        self.execute(first)
        second = self.args('with',annotation_key='arbitrary_annotation')
        manifest = self.execute(second)
        np.testing.assert_array_equal(self.scores_by_id(first.output_dir),self.scores_by_id(second.output_dir))
        for record in manifest['slices']:
            obj = ad.read_h5ad(Path(second.output_dir)/record['path'])
            self.assertIn('arbitrary_annotation',obj.obs)
            self.assertFalse(manifest['fitted_on_annotation'])

    def test_tiny_cell_dimension_cap_and_zero_norm_pca_rejection(self):
        values = sparse.csr_matrix([[2.,3.,5.,1.,6.],[5.,2.,3.,4.,2.],[3.,5.,2.,7.,1.]])
        result = prep.normalize_select_pca(values,np.array(['a','b','c','d','e']),'counts',2000,50,71)
        self.assertEqual(result[-1]['effective_n_pcs'],2)
        self.assertEqual(result[3].shape,(3,2))
        # The first row is exactly the global centroid in a rank-two log-space
        # matrix. A zero PCA norm would make cosine alignment ill-defined.
        values = sparse.csr_matrix([[2.,2.,2.],[3.,2.,1.],[1.,3.,2.],[2.,1.,3.]])
        with self.assertRaisesRegex(ValueError,'zero-norm'):
            prep.normalize_select_pca(values,np.array(['a','b','c']),'log1p',2000,50,71)

    def test_no_z_no_annotation_shared_pca_and_both_validators(self):
        reference = self.args('with_z_reference')
        self.execute(reference)
        renamed = []
        for path,number in zip(self.input_paths,[9,3,20]):
            with h5py.File(path,'r+') as handle:
                del handle['obs']['z']
            new_path = path.with_name(f'ordered_SL{number}.h5ad')
            path.rename(new_path)
            renamed.append(new_path)
        self.input_paths = renamed
        args = self.args('without_z')
        manifest = self.execute(args)
        self.assertFalse(manifest['physical_z_available'])
        self.assertEqual(manifest['slice_order_source'],'numeric_SL')
        self.assertEqual([r['input_SL_number'] for r in manifest['slices']],[3,9,20])
        for record in manifest['slices']:
            obj = ad.read_h5ad(Path(args.output_dir)/record['path'])
            self.assertNotIn('z',obj.obs)
            self.assertNotIn('physical_z',obj.obs)
            self.assertNotIn('anno',obj.obs)
            self.assertIsNone(record['physical_z'])
        np.testing.assert_array_equal(self.scores_by_id(reference.output_dir),self.scores_by_id(args.output_dir))
        package = Path(__file__).resolve().parents[1]
        for pipeline,style in [('pairwise-rigid','pairwise'),('continuity-guided','continuity')]:
            result = subprocess.run([sys.executable,str(package/'pipelines'/pipeline/'validate_inputs.py'),
                '--slice-dir',str(Path(args.output_dir)/'slice_h5ad'),'--representation','expression-pca',
                '--annotation-qc','off','--filename-style',style],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr+'\n'+result.stdout)
            validation = json.loads(result.stdout)
            self.assertEqual(validation['physical_z_order'],'metadata_absent_SL_order_only')

    def test_mixed_z_availability_rejected(self):
        with h5py.File(self.input_paths[0],'r+') as handle:
            del handle['obs']['z']
        with self.assertRaisesRegex(ValueError,'Mixed physical z availability'):
            self.execute(self.args())
        self.assertFalse((self.root/'result').exists())

    def test_all_z_aliases_and_explicit_key_must_agree_exactly(self):
        cases = [
            ('consistent_aliases','physical_z',None,None),
            ('consistent_explicit','custom_z',None,None),
            ('conflicting_alias','physical_z','z','offset'),
            ('conflicting_explicit','custom_z','custom_z','offset'),
            ('exact_not_approximate','physical_z','slice_z','nextafter'),
            ('nonfinite_secondary','physical_z','z','nonfinite'),
            ('nonconstant_secondary','physical_z','slice_z','nonconstant'),
        ]
        for name,explicit_key,changed_key,change in cases:
            with self.subTest(name=name):
                for i,path in enumerate(self.input_paths):
                    with h5py.File(path,'r+') as handle:
                        for key in ('physical_z','z','slice_z','custom_z'):
                            values = np.full(14,0.7+i*1.3)
                            if i==0 and key==changed_key:
                                if change=='offset':values += 99
                                elif change=='nextafter':values = np.nextafter(values,np.inf)
                                elif change=='nonfinite':values[0] = np.nan
                                else:values[0] += 1
                            if key in handle['obs']:del handle['obs'][key]
                            handle['obs'].create_dataset(key,data=values)
                args = self.args(name,physical_z_key=explicit_key)
                if change is not None:
                    message = 'finite constant' if change in ('nonfinite','nonconstant') else 'Conflicting physical z columns'
                    with self.assertRaisesRegex(ValueError,message):
                        self.execute(args)
                    self.assertFalse(Path(args.output_dir).exists())
                else:
                    manifest = self.execute(args)
                    np.testing.assert_array_equal([r['physical_z'] for r in manifest['slices']],
                                                  [0.7+i*1.3 for i in range(3)])
                    expected_keys = {'/obs/physical_z','/obs/z','/obs/slice_z'}
                    if explicit_key=='custom_z':expected_keys.add('/obs/custom_z')
                    for path in self.input_paths:
                        reads = [r for r in manifest['read_access_ledger']
                                 if r['file']==str(path.resolve()) and r['dataset'] in expected_keys]
                        self.assertEqual({r['dataset'] for r in reads},expected_keys)
                        self.assertTrue(all(r['read_calls']==1 for r in reads))
                    for record in manifest['slices']:
                        obj = ad.read_h5ad(Path(args.output_dir)/record['path'])
                        np.testing.assert_array_equal(obj.obs['physical_z'],obj.obs['z'])
                        np.testing.assert_array_equal(obj.obs['z'],np.full(len(obj),record['physical_z']))

    def test_no_z_requires_unambiguous_numeric_sl_order(self):
        for path in self.input_paths:
            with h5py.File(path,'r+') as handle:
                del handle['obs']['z']
        with self.assertRaisesRegex(ValueError,'ordering requires one numeric SL'):
            self.execute(self.args())
        for i,path in enumerate(self.input_paths):
            path.rename(path.with_name(f'prefix{i}_SL7.h5ad'))
        with self.assertRaisesRegex(ValueError,'numeric SL identifiers must be unique'):
            self.execute(self.args())


if __name__=='__main__':
    unittest.main(verbosity=2)
