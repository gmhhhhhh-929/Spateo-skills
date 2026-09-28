"""Run one pinned Spateo backend for a controlled two-timepoint comparison.

Run in separate processes/environments for old and native Spateo. Fixed reference
IDs isolate alignment from the replaced sampler. Never import both backends into
one interpreter, or reuse a comparison's coordinates for downstream inference
without a separately accepted analysis run.
"""
import argparse
import hashlib
import importlib.metadata
import json
import random
import sys
import time
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--source-root', help='Pinned checkout containing spateo/')
    parser.add_argument('--reference-policy', choices=['fixed', 'library'], default='fixed')
    parser.add_argument('--sample-size', type=int, default=2000)
    parser.add_argument('--annotation-aware', action='store_true')
    parser.add_argument('--initialization', choices=['nn', 'identity'], default='nn')
    args = parser.parse_args()
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    config_path = Path(args.config).resolve()
    cfg = json.loads(config_path.read_text())
    if args.source_root:
        sys.path.insert(0, str(Path(args.source_root).resolve()))
    import numpy as np
    import scipy.sparse as sparse
    import anndata as ad
    import spateo as st
    import torch
    from registration_qc import pair_qc
    seed = cfg['runtime']['seed']
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    root = Path(st.__file__).resolve().parent
    if args.source_root and root != Path(args.source_root).resolve() / 'spateo':
        raise RuntimeError('Imported Spateo does not match the requested source checkout')
    record = {'config_sha256': sha(config_path), 'spateo_root': str(root),
              'parameters': vars(args), 'seed': seed, 'status': 'running',
              'source_hashes': {str(p.relative_to(root)): sha(p) for p in sorted((root/'alignment').rglob('*.py'))},
              'versions': {p: importlib.metadata.version(p) for p in ('numpy','scipy','anndata','torch','POT')}}
    def write_record():
        (out/'comparison.json').write_text(json.dumps(record, indent=2))
    write_record()
    started = time.monotonic()
    pair, refs, reference_ids, input_provenance = [], [], [], []
    alignment = cfg['alignment']
    for key in ('stage1','stage2'):
        path = Path(cfg['inputs'][key])
        if not path.is_absolute(): path = config_path.parent/path
        original = ad.read_h5ad(path, backed='r')
        original.uns.clear()
        # Full cells need coordinates for BA_transform, not duplicated expression.
        # Only the exact selected reference cells enter the scientific optimizer.
        full = ad.AnnData(obs=original.obs.copy(), var=original.var.copy(),
                          obsm={alignment['spatial_key']:np.asarray(original.obsm[alignment['spatial_key']]).copy()})
        n = min(args.sample_size, original.n_obs)
        fixed_idx = np.random.default_rng(19491001).choice(original.n_obs,n,replace=False)
        fixed_ids = original.obs_names[fixed_idx].tolist()
        reference_ids.append(fixed_ids)
        if args.reference_policy=='library':
            from spateo.alignment.utils import downsampling
            ids = downsampling([full],n_sampling=n,sampling_method='random',spatial_key=alignment['spatial_key'])[0].obs_names.tolist()
        else:
            ids = fixed_ids
        a = original[ids].to_memory()
        counts = a.layers.get(alignment['counts_layer'])
        if counts is None:
            if not alignment.get('x_is_counts'): raise ValueError('A verified count layer is required')
            counts = a.X
        x = sparse.csr_matrix(counts, dtype=np.float64)
        if x.data.size and (np.min(x.data)<0 or not np.isfinite(x.data).all()):
            raise ValueError('Invalid expression input')
        totals = np.asarray(x.sum(axis=1)).ravel()
        if np.any(totals<=0): raise ValueError('Zero-library reference cells')
        target = alignment.get('target_sum')
        if target is None:
            raw = original.layers.get(alignment['counts_layer'],original.X)
            all_totals=np.concatenate([np.asarray(raw[i:i+10000].sum(axis=1)).ravel() for i in range(0,original.n_obs,10000)])
            target=float(np.median(all_totals[all_totals>0]))
        normalized = sparse.diags(target/totals) @ x
        logged = normalized.tocsr(); logged.data = np.log1p(logged.data)
        a.layers['benchmark_log1p'] = logged
        refs.append(a); pair.append(full)
        input_provenance.append({'path':str(path),'sha256':sha(path),'shape':list(original.shape),'normalization_target':target})
        original.file.close()
    (out/'fixed_reference_ids.json').write_text(json.dumps(reference_ids))
    gene_set=set(pair[1].var_names); genes=[g for g in pair[0].var_names if g in gene_set]
    annotation=cfg['subset']['annotation_key']
    layers=['benchmark_log1p',annotation] if args.annotation_aware else 'benchmark_log1p'
    fields=['layer','obs'] if args.annotation_aware else 'layer'
    extra={}
    if args.initialization=='identity': extra['nn_init']=False
    record.update(inputs=input_provenance,genes=len(genes),initialization_layer='X',fixed_reference_sha256=sha(out/'fixed_reference_ids.json'))
    write_record()
    aligned, aligned_refs, _, _ = st.align.morpho_align_ref(
        models=pair,models_ref=refs,n_sampling=args.sample_size,sampling_method='random',
        rep_layer=layers,rep_field=fields,genes=genes,spatial_key=alignment['spatial_key'],
        key_added='benchmark',mode='SN-S',iter_key_added=None,max_iter=alignment['max_iter'],
        nonrigid_start_iter=alignment['nonrigid_start_iter'],
        device=cfg['runtime']['device'],verbose=False,**extra)
    arrays={}
    for i,a in enumerate(aligned):
        for name in ('benchmark','benchmark_rigid','benchmark_nonrigid'):
            arrays[f'stage{i+1}_{name}']=np.asarray(a.obsm[name])
        arrays[f'stage{i+1}_before']=np.asarray(a.obsm[alignment['spatial_key']])
        arrays[f'stage{i+1}_ids']=a.obs_names.to_numpy(dtype=str)
        arrays[f'stage{i+1}_labels']=a.obs[annotation].to_numpy(dtype=str)
    np.savez_compressed(out/'coordinates.npz',**arrays)
    record['actual_reference_ids']=[a.obs_names.tolist() for a in aligned_refs]
    record['qc']={name:pair_qc(aligned[0].obsm[name],aligned[1].obsm[name],aligned[0].obs[annotation],aligned[1].obs[annotation])
                  for name in ('benchmark_rigid','benchmark_nonrigid')}
    record.update(status='completed',elapsed_seconds=time.monotonic()-started,coordinates_sha256=sha(out/'coordinates.npz'))
    write_record()


if __name__=='__main__':
    main()
