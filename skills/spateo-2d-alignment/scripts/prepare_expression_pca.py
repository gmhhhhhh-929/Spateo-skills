#!/usr/bin/env python3
"""Build one shared expression PCA basis without reading reference coordinates.

Only explicitly whitelisted HDF5 datasets are read. The expression source is
never loaded as AnnData and never hashed as a whole file. AnnData is used only
to write new sanitized output objects assembled from allowed arrays.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import time

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, svds

SCHEMA_VERSION = 'expression-pca-v1'
FEATURE_HASH_ENCODING = 'sha256(np.ascontiguousarray(X_pca,dtype="<f4").tobytes(order="C"))'
ID_HASH_ENCODING = 'sha256(sorted IDs as compact ensure_ascii=False UTF-8 JSON)'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def decode(value):
    return value.decode('utf-8') if isinstance(value,(bytes,np.bytes_)) else str(value)


def strings(values):
    return np.array([decode(v) for v in np.asarray(values).reshape(-1)],dtype=str)


def hash_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def hash_ids(ids, sort=True):
    values = strings(ids).tolist()
    if sort:
        values = sorted(values)
    return hashlib.sha256(json.dumps(values,ensure_ascii=False,separators=(',',':')).encode()).hexdigest()


def hash_numeric(values,dtype='<f4'):
    return hashlib.sha256(np.ascontiguousarray(values,dtype=dtype).tobytes(order='C')).hexdigest()


def hash_csr(matrix):
    value = matrix.tocsr(copy=True)
    value.sum_duplicates()
    value.sort_indices()
    value.eliminate_zeros()
    digest = hashlib.sha256()
    for a,dtype in ((np.array(value.shape),'<i8'),(value.indptr,'<i8'),
                    (value.indices,'<i8'),(value.data,'<f8')):
        digest.update(np.ascontiguousarray(a,dtype=dtype).tobytes())
    return digest.hexdigest()


class ReadLedger:
    """Aggregate dataset access, including the exact permitted HDF5 paths."""
    def __init__(self):
        self.entries = {}

    def read(self,node,selection=(),role='expression_source'):
        require(isinstance(node,h5py.Dataset),'Expected HDF5 dataset: '+node.name)
        value = node[selection]
        key = (role,str(Path(node.file.filename).resolve()),node.name)
        entry = self.entries.setdefault(key,{'role':role,'file':key[1],'dataset':node.name,
                                             'read_calls':0,'elements_read':0,'numeric_bytes_read':0})
        entry['read_calls'] += 1
        entry['elements_read'] += int(np.asarray(value).size)
        entry['numeric_bytes_read'] += int(np.asarray(value).nbytes)
        return value

    def records(self):
        return [self.entries[k] for k in sorted(self.entries)]


def read_column(node,ledger,role):
    if isinstance(node,h5py.Dataset):
        values = ledger.read(node,role=role)
        return strings(values) if np.asarray(values).dtype.kind in 'OSU' else np.asarray(values)
    require(isinstance(node,h5py.Group) and 'codes' in node and 'categories' in node,
            'Unsupported column encoding: '+node.name)
    codes = np.asarray(ledger.read(node['codes'],role=role),dtype=int)
    categories = read_column(node['categories'],ledger,role)
    require(np.all(codes>=0) and np.all(codes<len(categories)),'Missing categorical values: '+node.name)
    return categories[codes]


def read_index(handle,group_name,ledger,role):
    group = handle[group_name]
    name = decode(group.attrs.get('_index','_index'))
    require(name in group,'Missing '+group_name+' IDs')
    ids = strings(read_column(group[name],ledger,role))
    require(len(set(ids))==len(ids),'Duplicate '+group_name+' IDs')
    require(all(len(v)>0 for v in ids),'Empty '+group_name+' IDs')
    return ids


def inspect_slices(slice_dir,annotation_key,physical_z_key,ledger):
    paths = sorted(Path(slice_dir).resolve().glob('*.h5ad'))
    require(bool(paths),'No H5AD slices found')
    records = []
    seen = set()
    for path in paths:
        before = path.stat()
        with h5py.File(path,'r') as handle:
            ids = read_index(handle,'obs',ledger,'slice_metadata')
            require(not seen.intersection(ids),'A cell ID occurs in multiple slices')
            seen.update(ids)
            require('spatial' in handle['obsm'],'Slice lacks obsm/spatial: '+str(path))
            xy = np.asarray(ledger.read(handle['obsm']['spatial'],role='slice_spatial'),dtype=np.float64)
            require(xy.shape==(len(ids),2) and np.isfinite(xy).all(),'Slice spatial must be finite N x 2')
            keys = [key for key in dict.fromkeys([physical_z_key,'physical_z','z','slice_z'])
                    if key is not None and key in handle['obs']]
            z_key = keys[0] if keys else None
            z_values = None
            for key in keys:
                values = np.asarray(read_column(handle['obs'][key],ledger,'slice_metadata'),dtype=np.float64)
                require(values.shape==(len(ids),) and np.isfinite(values).all() and len(np.unique(values))==1,
                        'Physical z must be one finite constant per slice: '+key+' in '+str(path))
                if z_values is None:
                    z_values = values
                else:
                    require(np.array_equal(z_values,values),
                            'Conflicting physical z columns: '+z_key+' and '+key+' in '+str(path))
            annotation = None
            if annotation_key is not None and annotation_key in handle['obs']:
                annotation = strings(read_column(handle['obs'][annotation_key],ledger,'slice_annotation'))
                require(len(annotation)==len(ids),'Annotation length mismatch')
        content = {'cell_ids_ordered_sha256':hash_ids(ids,False),'spatial_sha256':hash_numeric(xy,'<f8'),
                   'physical_z_sha256':hash_numeric(z_values,'<f8') if z_values is not None else None,
                   'annotation_sha256':hash_ids(annotation,False) if annotation is not None else None}
        content_hash = hashlib.sha256(json.dumps(content,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        records.append({'input_path':str(path),'input_stat':{'size':before.st_size,'mtime_ns':before.st_mtime_ns},
                        'ids':ids,'xy':xy,'physical_z':float(z_values[0]) if z_values is not None else None,'annotation':annotation,
                        'input_z_key':z_key,'allowed_input_content':content,'input_allowed_content_sha256':content_hash})
    available = [r['physical_z'] is not None for r in records]
    require(all(available) or not any(available),
            'Mixed physical z availability is not allowed: either every slice provides z or all omit it')
    if all(available):
        records.sort(key=lambda r:(r['physical_z'],r['input_path']))
        require(len(set(r['physical_z'] for r in records))==len(records),'Slices must have distinct increasing physical z')
        for record in records:
            record['slice_order_source'] = 'physical_z'
    else:
        for record in records:
            matches = re.findall(r'(?:^|[_-])SL(\d+)(?=[_.-]|$)',Path(record['input_path']).name,flags=re.IGNORECASE)
            require(len(matches)==1,'Missing physical z: slice ordering requires one numeric SL identifier per filename')
            record['input_SL_number'] = int(matches[0])
            record['slice_order_source'] = 'numeric_SL'
        require(len(set(r['input_SL_number'] for r in records))==len(records),
                'Missing physical z: numeric SL identifiers must be unique')
        records.sort(key=lambda r:r['input_SL_number'])
    return records


def select_expression_node(handle,choice):
    choices = ('counts_X','counts','X') if choice=='auto' else (choice,)
    for name in choices:
        if name=='X' and 'X' in handle:
            return handle['X'],'X'
        if name!='X' and 'layers' in handle and name in handle['layers']:
            return handle['layers'][name],'layers/'+name
    raise ValueError('Requested expression matrix is absent; auto tries layers/counts_X, layers/counts, then X')


def read_wanted_matrix(node,source_rows,n_source_cells,n_genes,ledger):
    """Read CSR selected rows, dense selected rows, or stream selected CSC rows."""
    order = np.argsort(source_rows,kind='stable')
    ordered_rows = np.asarray(source_rows,dtype=np.int64)[order]
    blocks = []
    if isinstance(node,h5py.Dataset):
        require(node.shape==(n_source_cells,n_genes),'Expression matrix shape does not match IDs')
        for start in range(0,len(ordered_rows),256):
            selected = ordered_rows[start:start+256]
            blocks.append(sparse.csr_matrix(np.asarray(ledger.read(node,(selected,slice(None))))))
        packed = sparse.vstack(blocks,format='csr')
        return packed[np.argsort(order)].tocsr()
    require(isinstance(node,h5py.Group),'Unsupported expression storage')
    encoding = decode(node.attrs.get('encoding-type',''))
    shape = tuple(int(v) for v in node.attrs.get('shape',()))
    require(shape==(n_source_cells,n_genes),'Sparse expression shape does not match IDs')
    require(encoding in ('csr_matrix','csc_matrix'),'Expression must be dense, CSR or CSC')
    indptr = np.asarray(ledger.read(node['indptr']),dtype=np.int64)
    if encoding=='csr_matrix':
        require(len(indptr)==n_source_cells+1,'Invalid CSR indptr')
        # Adjacent selected rows are read together; chunks cap temporary memory.
        cursor = 0
        while cursor<len(ordered_rows):
            end = cursor+1
            while end<len(ordered_rows) and end-cursor<2048 and ordered_rows[end]==ordered_rows[end-1]+1:
                end += 1
            first,last = int(ordered_rows[cursor]),int(ordered_rows[end-1])+1
            left,right = int(indptr[first]),int(indptr[last])
            data = np.asarray(ledger.read(node['data'],slice(left,right)))
            indices = np.asarray(ledger.read(node['indices'],slice(left,right)),dtype=np.int64)
            require(np.all((indices>=0)&(indices<n_genes)),'Invalid CSR gene index')
            blocks.append(sparse.csr_matrix((data,indices,indptr[first:last+1]-left),shape=(last-first,n_genes)))
            cursor = end
        packed = sparse.vstack(blocks,format='csr')
        return packed[np.argsort(order)].tocsr()
    require(len(indptr)==n_genes+1,'Invalid CSC indptr')
    target = np.full(n_source_cells,-1,dtype=np.int64)
    target[np.asarray(source_rows,dtype=np.int64)] = np.arange(len(source_rows))
    rows,columns,values = [],[],[]
    for start in range(0,n_genes,128):
        end = min(start+128,n_genes)
        left,right = int(indptr[start]),int(indptr[end])
        source_indices = np.asarray(ledger.read(node['indices'],slice(left,right)),dtype=np.int64)
        data = np.asarray(ledger.read(node['data'],slice(left,right)))
        require(np.all((source_indices>=0)&(source_indices<n_source_cells)),'Invalid CSC row index')
        mapped = target[source_indices]
        keep = mapped>=0
        rows.append(mapped[keep])
        columns.append(np.repeat(np.arange(start,end),np.diff(indptr[start:end+1]))[keep])
        values.append(data[keep])
    return sparse.coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(columns))),
                             shape=(len(source_rows),n_genes)).tocsr()


def load_expression(source_path,wanted_ids,layer,ledger):
    source_path = Path(source_path).resolve()
    before = source_path.stat()
    with h5py.File(source_path,'r') as handle:
        source_ids = read_index(handle,'obs',ledger,'expression_source')
        genes = read_index(handle,'var',ledger,'expression_source')
        mapping = {cell:i for i,cell in enumerate(source_ids)}
        missing = [cell for cell in wanted_ids if cell not in mapping]
        require(not missing,'Expression source lacks '+str(len(missing))+' required cell IDs')
        source_rows = np.array([mapping[cell] for cell in wanted_ids],dtype=np.int64)
        node,selected_layer = select_expression_node(handle,layer)
        matrix = read_wanted_matrix(node,source_rows,len(source_ids),len(genes),ledger)
    matrix = matrix.astype(np.float64).tocsr()
    matrix.sum_duplicates()
    matrix.sort_indices()
    matrix.eliminate_zeros()
    audit = {'path':str(source_path),'selected_layer':selected_layer,'n_source_cells':len(source_ids),
        'n_wanted_cells':len(wanted_ids),'n_source_genes':len(genes),'wanted_expression_matrix_sha256':hash_csr(matrix),
        'expression_matrix_hash_encoding':'canonical CSR shape/indptr/indices <i8 and values <f8',
        'source_cell_ids_sha256':hash_ids(source_ids),'source_cell_ids_ordered_sha256':hash_ids(source_ids,False),
        'source_gene_ids_sha256':hash_ids(genes),'source_gene_ids_ordered_sha256':hash_ids(genes,False),
        'input_stat':{'size':before.st_size,'mtime_ns':before.st_mtime_ns},
        'source_file_sha256':None,'source_full_file_hash_performed':False,
        'source_file_hash_omission_reason':'Whole-file hashing would read unapproved spatial/reference bytes'}
    return matrix,genes,audit


def validate_matrix(matrix,state):
    require(matrix.ndim==2 and min(matrix.shape)>=3,'Expression needs at least three cells and genes')
    require(np.isfinite(matrix.data).all(),'Expression contains nonfinite values')
    require(np.all(matrix.data>=0),'Expression must be nonnegative')
    if state=='counts':
        require(np.all(np.abs(matrix.data-np.rint(matrix.data))<=1e-6),
                'Counts mode requires integer-like values; use --matrix-state log1p only for declared log1p data; no inverse transform is guessed')
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    require(np.all(np.isfinite(totals)) and np.all(totals>0),'Zero-total or nonfinite-total cells are not valid PCA/cosine inputs')
    strict_onehot = np.all(np.diff(matrix.indptr)==1) and np.all(matrix.data==1)
    require(not strict_onehot,'Strict one-hot expression is not an expression PCA input')
    return totals


def normalize_select_pca(matrix,genes,state,n_hvg,n_pcs,seed):
    """Fit only expression values; deliberately no coordinate/annotation args."""
    require(n_hvg>=3 and n_pcs>=2,'Request at least three genes and two PCs')
    totals = validate_matrix(matrix,state)
    normalized = matrix.copy().astype(np.float64)
    if state=='counts':
        for start in range(0,normalized.shape[0],4096):
            end = min(start+4096,normalized.shape[0])
            left,right = normalized.indptr[start],normalized.indptr[end]
            normalized.data[left:right] *= np.repeat(10000.0/totals[start:end],np.diff(normalized.indptr[start:end+1]))
        np.log1p(normalized.data,out=normalized.data)
    # Declared log1p data are used as-is: neither library-size renormalization
    # nor a second logarithm is applied.
    means_all = np.asarray(normalized.mean(axis=0)).ravel()
    second = np.zeros(normalized.shape[1],dtype=np.float64)
    for start in range(0,normalized.shape[0],4096):
        squared = normalized[start:start+4096].copy()
        np.square(squared.data,out=squared.data)
        second += np.asarray(squared.sum(axis=0)).ravel()
    variance = np.maximum(second/normalized.shape[0]-means_all**2,0)
    variable = np.flatnonzero(variance>1e-12)
    require(len(variable)>=3,'Fewer than three variable genes; constant expression is not valid PCA input')
    # Gene IDs break exact variance ties; final basis columns are in gene-ID
    # order, making mapping independent of input expression row order.
    ranked = sorted(variable,key=lambda i:(-variance[i],str(genes[i])))
    selected = np.array(sorted(ranked[:min(n_hvg,len(ranked))],key=lambda i:str(genes[i])),dtype=int)
    selected_matrix = normalized[:,selected].tocsr()
    effective_pcs = min(n_pcs,selected_matrix.shape[1]-1,selected_matrix.shape[0]-1)
    require(effective_pcs>=2,'Effective PCA dimension must be at least two')
    means = np.asarray(selected_matrix.mean(axis=0)).ravel()
    def matvec(v):
        v = np.asarray(v).reshape(-1)
        return np.asarray(selected_matrix@v).ravel()-means@v
    def rmatvec(v):
        v = np.asarray(v).reshape(-1)
        return np.asarray(selected_matrix.T@v).ravel()-means*v.sum()
    def matmat(v):
        return selected_matrix@v-(means@v)[None,:]
    def rmatmat(v):
        return selected_matrix.T@v-means[:,None]*v.sum(axis=0)[None,:]
    operator = LinearOperator(selected_matrix.shape,dtype=np.float64,matvec=matvec,rmatvec=rmatvec,
                              matmat=matmat,rmatmat=rmatmat)
    _,singular,vt = svds(operator,k=effective_pcs,which='LM',tol=1e-6,maxiter=3000,random_state=seed)
    descending = np.argsort(singular)[::-1]
    singular,loadings = singular[descending],vt[descending].T
    rank_threshold = max(float(singular[0]),1.0)*1e-10
    numerical_rank = int(np.sum(singular>rank_threshold))
    require(numerical_rank>=2,'Expression has fewer than two nonconstant PCA directions')
    if numerical_rank<effective_pcs:
        effective_pcs = numerical_rank
        singular,loadings = singular[:effective_pcs],loadings[:,:effective_pcs]
    for i in range(effective_pcs):
        if loadings[np.argmax(abs(loadings[:,i])),i]<0:
            loadings[:,i] *= -1
    scores = np.asarray(selected_matrix@loadings-(means@loadings)[None,:],dtype=np.float32)
    require(np.isfinite(scores).all(),'PCA features contain nonfinite values')
    require(np.all(np.linalg.norm(scores.astype(np.float64),axis=1)>1e-12),'PCA has zero-norm cells')
    require(np.all(np.var(scores.astype(np.float64),axis=0)>1e-12),'PCA contains constant columns')
    require(not (np.isin(scores,[0,1]).all() and np.all(scores.sum(axis=1)==1)),
            'PCA unexpectedly resembles one-hot features')
    audit = {'requested_n_hvg':n_hvg,'effective_n_hvg':len(selected),'requested_n_pcs':n_pcs,
        'effective_n_pcs':effective_pcs,'dimension_cap':min(n_pcs,len(selected)-1,len(scores)-1),
        'numerical_rank_retained':numerical_rank,'matrix_state':state,'normalization':
        'total-count normalize to 10000 across all genes, then natural log1p' if state=='counts' else
        'declared log1p values used unchanged; no normalization or extra logarithm',
        'hvg_method':'top population variance of lognormalized expression; gene-ID tie breaking',
        'algorithm':'scipy.sparse.linalg.svds, centered sparse LinearOperator','seed':seed,
        'centered':True,'unit_variance_scaled':False,'fitted_on_spatial':False,'fitted_on_annotation':False,
        'normalization_applied':state=='counts','log1p_applied':state=='counts',
        'input_total_min':float(totals.min()),'input_total_max':float(totals.max()),
        'selected_gene_ids_sha256':hash_ids(genes[selected],False),'joint_features_sha256':hash_numeric(scores),
        'explained_variance':(singular**2/(len(scores)-1)).tolist(),
        'finite':True,'nonzero_norm_all_cells':True,'strict_onehot':False}
    return selected_matrix,genes[selected],selected,scores,means,loadings,singular,audit


def verify_stat(path,expected):
    current = Path(path).stat()
    require(current.st_size==expected['size'] and current.st_mtime_ns==expected['mtime_ns'],
            'An input changed during preparation: '+str(path))


def prepare(args):
    import anndata as ad
    started = time.time()
    output = Path(args.output_dir).resolve()
    require(not output.exists(),'Output must be a new directory')
    source_path = Path(args.expression_source).resolve()
    slice_paths = list(Path(args.slice_dir).resolve().glob('*.h5ad'))
    require(source_path not in [p.resolve() for p in slice_paths],
            'Expression source must be separate from the slice XY inputs')
    ledger = ReadLedger()
    slices = inspect_slices(args.slice_dir,args.annotation_key,args.physical_z_key,ledger)
    wanted_ids = np.array(sorted(cell for r in slices for cell in r['ids']),dtype=str)
    expression,genes,source_audit = load_expression(source_path,wanted_ids,args.expression_layer,ledger)
    matrix,selected_genes,selected,scores,means,loadings,singular,pca = normalize_select_pca(
        expression,genes,args.matrix_state,args.n_hvg,args.n_pcs,args.seed)
    output.mkdir(parents=True)
    destination = output/'slice_h5ad'
    destination.mkdir()
    np.savez_compressed(output/'basis.npz',gene_ids=selected_genes,source_gene_indices=selected,
                        gene_means=means,loadings=loadings,singular_values=singular,
                        explained_variance=singular**2/(len(scores)-1),fit_cell_ids=wanted_ids)
    basis_hash = hash_file(output/'basis.npz')
    index = {cell:i for i,cell in enumerate(wanted_ids)}
    outputs = []
    for i,record in enumerate(slices):
        ids = record['ids']
        rows = np.array([index[cell] for cell in ids],dtype=int)
        # The CS01/YEXPRESSION tokens satisfy the bundled legacy pairwise
        # runner's filename grammar; they do not encode biological identity.
        slice_id = f'CS01_SL{i:04d}_YEXPRESSION'
        name = slice_id+'.Spatial.h5ad'
        obs = pd.DataFrame(index=pd.Index(ids,name=None))
        obs['cell_id'] = ids
        obs['slice_id'] = slice_id
        obs['slice_index'] = np.full(len(ids),i,dtype=np.int32)
        if record['physical_z'] is not None:
            obs['physical_z'] = record['physical_z']
            obs['z'] = record['physical_z']
        annotation_key = args.annotation_key if record['annotation'] is not None else None
        if annotation_key:
            require(annotation_key not in obs.columns,'Annotation key conflicts with preserved identity/z metadata')
            obs[annotation_key] = record['annotation']
        obj = ad.AnnData(X=matrix[rows].astype(np.float32).tocsr(),obs=obs,
                         var=pd.DataFrame(index=pd.Index(selected_genes,name=None)))
        obj.obsm['spatial'] = record['xy'].copy()
        obj.obsm['X_pca'] = np.ascontiguousarray(scores[rows],dtype=np.float32)
        obj.uns['expression_pca'] = {'schema_version':SCHEMA_VERSION,'mode':'expression-pca',
            'shared_basis':True,'basis_sha256':basis_hash,'feature_key':'X_pca',
            'fitted_on_spatial':False,'fitted_on_annotation':False,'manifest':'../expression_pca_manifest.json'}
        path = destination/name
        obj.write_h5ad(path,compression='gzip')
        # Only newly generated sanitized outputs may be read and fully hashed.
        with h5py.File(path,'r') as verify:
            require(set(verify['obsm'])=={'spatial','X_pca'},'Unexpected output coordinate/feature key')
            require(np.array_equal(verify['obsm']['X_pca'][()],scores[rows]),'Output PCA readback mismatch')
            require(np.array_equal(verify['obsm']['spatial'][()],record['xy']),'Output XY readback mismatch')
        outputs.append({'slice_id':slice_id,'slice_index':i,'name':name,'path':'slice_h5ad/'+name,
            'output_h5ad_sha256':hash_file(path),'n_cells':len(ids),'cell_ids_sha256':hash_ids(ids),
            'features_sha256':hash_numeric(scores[rows]),'feature_hash_encoding':FEATURE_HASH_ENCODING,
            'basis_sha256':basis_hash,'feature_key':'X_pca','spatial_key':'spatial',
            'physical_z':record['physical_z'],'physical_z_available':record['physical_z'] is not None,
            'slice_order_source':record['slice_order_source'],'input_SL_number':record.get('input_SL_number'),
            'annotation_key':annotation_key,
            'input_path':record['input_path'],'input_z_key':record['input_z_key'],
            'input_file_sha256':None,'input_full_file_hash_performed':False,
            'input_allowed_content_sha256':record['input_allowed_content_sha256'],
            'input_allowed_content':record['allowed_input_content'],'input_stat':record['input_stat']})
    verify_stat(source_path,source_audit['input_stat'])
    for record in slices:
        verify_stat(record['input_path'],record['input_stat'])
    manifest = {'schema_version':SCHEMA_VERSION,'status':'frozen_before_alignment',
        'closed_at_unix':time.time(),'mode':'expression-pca','feature_mode':'expression-pca',
        'shared_basis':True,'basis_scope':'all and only requested slice cell IDs, jointly across every slice',
        'physical_z_available':slices[0]['physical_z'] is not None,'slice_order_source':slices[0]['slice_order_source'],
        'fitted_on_spatial':False,'fitted_on_annotation':False,'annotation':False,
        'fit_used_annotations':False,'fit_used_spatial':False,'basis_file':'basis.npz','basis_sha256':basis_hash,
        'feature_key':'X_pca','spatial_key':'spatial','feature_hash_encoding':FEATURE_HASH_ENCODING,
        'cell_ids_hash_encoding':ID_HASH_ENCODING,'n_cells':len(wanted_ids),'n_slices':len(slices),
        'cell_ids_sha256':hash_ids(wanted_ids),'n_hvg':pca['effective_n_hvg'],'n_pcs':pca['effective_n_pcs'],
        'requested_n_hvg':args.n_hvg,'requested_n_pcs':args.n_pcs,'effective_n_hvg':pca['effective_n_hvg'],
        'effective_n_pcs':pca['effective_n_pcs'],'pca':pca,'expression_source':source_audit,'slices':outputs,
        'read_access_ledger':ledger.records(),'source_coordinates_read':False,'source_annotations_read':False,
        'preexisting_pca_read':False,'input_files_modified':False,
        'input_hash_policy':'Source and slice H5ADs are never whole-file hashed; hash only approved datasets, because other keys may contain references. New sanitized outputs are fully hashed.',
        'implementation_sha256':hash_file(__file__),'runtime_seconds':time.time()-started,
        'runtime':{'python':sys.version,'numpy':np.__version__,'scipy':__import__('scipy').__version__,
                   'h5py':h5py.__version__,'anndata':ad.__version__}}
    with (output/'expression_pca_manifest.json').open('x') as handle:
        json.dump(manifest,handle,ensure_ascii=False,indent=2,allow_nan=False)
    print(json.dumps({'output':str(output),'n_cells':len(wanted_ids),'n_slices':len(slices),
                      'effective_n_hvg':pca['effective_n_hvg'],'effective_n_pcs':pca['effective_n_pcs'],
                      'basis_sha256':basis_hash},ensure_ascii=False))
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--slice-dir',required=True)
    parser.add_argument('--expression-source',required=True)
    parser.add_argument('--expression-layer',choices=('auto','counts_X','counts','X'),default='auto')
    parser.add_argument('--matrix-state',choices=('counts','log1p'),default='counts')
    parser.add_argument('--output-dir',required=True)
    parser.add_argument('--n-hvg',type=int,default=2000)
    parser.add_argument('--n-pcs',type=int,default=50)
    parser.add_argument('--seed',type=int,default=7021)
    parser.add_argument('--annotation-key',default=None,help='Optional slice obs field to copy; never used in PCA')
    parser.add_argument('--physical-z-key',default='physical_z',help='Slice obs Z field; fallbacks physical_z, z, slice_z')
    prepare(parser.parse_args())


if __name__=='__main__':
    main()
