# Current IO source API

Source commit: `615644f88613bea8ceb2e2df1e2391d16de55ec1`. Automatic APIs return SpatialReadResult.

## spateo/io/general/_serialization.py

```python
def save(file, path)
def load(path, backend=None)
```

## spateo/io/general/_tabular.py

```python
def read_csv(**kwargs)
```

## spateo/io/single/_formats.py

```python
def read_h5ad(filename, **kwargs)
def read_10x_h5(filename, *, genome: str | None=None, gex_only: bool=True, backup_url: str | None=None)
def read_10x_mtx(path, *, var_names: Literal['gene_symbols', 'gene_ids']='gene_symbols', make_unique: bool=True, gex_only: bool=True, prefix: str | None=None, compressed: bool=True)
```

## spateo/io/single/_read.py

```python
def read(path, backend='python', **kwargs)
```

## spateo/io/spatial/_atera.py

```python
def read_atera(path: PathLike, *, library_id: Optional[str]=None, load_image: bool=True, image_key: str='dapi', image_max_dim: int=4096, load_boundaries: bool=True, load_nucleus_boundaries: bool=True, load_cell_groups: bool=True, cell_groups_csv: Optional[PathLike]=None, load_he_image: bool=False, he_image: Optional[PathLike]=None, he_alignment_csv: Optional[PathLike]=None, he_max_dim: int=4096, cache_file: Optional[PathLike]=None)
```

## spateo/io/spatial/_image.py

```python
def read_image(adata: AnnData, filename: str, scale_factor: float, slice: Optional[str]=None, img_layer: Optional[str]=None)
```

## spateo/io/spatial/_merfish.py

```python
def read_merfish(path: Union[str, Path], *, counts_file: str, meta_file: str, load_images: bool=True, z_layer: Optional[int]=3, load_boundaries: bool=True)
```

## spateo/io/spatial/_nanostring.py

```python
def read_nanostring(path: Union[str, Path], *, counts_file: str, meta_file: str, fov_file: Optional[str]=None)
```

## spateo/io/spatial/_provenance.py

```python
def spatial_file_manifest(path: PathLike, *, max_files: int=10000)
```

## spateo/io/spatial/_seqfish.py

```python
def read_seqfish(path: Union[str, Path], *, counts_file: str, meta_file: str, load_images: bool=True, load_labels: bool=True)
```

## spateo/io/spatial/_seqscope.py

```python
def read_seqscope(matrix_dir: PathLike, positions_path: PathLike, binsize: Optional[int]=1, add_props: bool=True)
```

## spateo/io/spatial/_slideseq.py

```python
def read_slideseq(path: Union[str, Path], *, counts_file: str='MappedDGEForR.csv', bead_file: str='BeadLocationsForR.csv', load_images: bool=True)
```

## spateo/io/spatial/_starmap_plus.py

```python
def read_starmap_plus(path: Union[str, Path], *, counts_file: str, meta_file: str, spatial_file: str, reorient_xy: bool=False, dtype: str='float32')
```

## spateo/io/spatial/_stereoseq.py

```python
def read_bgi_as_dataframe(path: Union[str, PathLike], label_column: Optional[str]=None)
def read_stereoseq(path, *, bin_size=None, chemistry=None, load_images=True, max_memory_bytes=1024 ** 3)
def read_bgi_agg(path: Union[str, PathLike], stain_path: Optional[Union[str, PathLike]]=None, binsize: int=1, gene_agg: Optional[Dict[str, Union[List[str], Callable[[str], bool]]]]=None, prealigned: bool=False, label_column: Optional[str]=None, version: Literal['stereo']='stereo')
def read_bgi(path: Union[str, PathLike], binsize: Optional[int]=None, segmentation_adata: Optional[AnnData]=None, labels_layer: Optional[str]=None, labels: Optional[Union[np.ndarray, str, PathLike]]=None, seg_binsize: int=1, label_column: Optional[str]=None, add_props: bool=True, version: Literal['stereo']='stereo')
```

## spateo/io/spatial/_visium.py

```python
def read_visium(path: Union[str, 'PathLike[str]'], genome: Optional[str]=None, *, count_file: str='filtered_feature_bc_matrix.h5', library_id: Optional[str]=None, load_images: bool=True, load_qc_images: bool=False, source_image_path: Optional[Union[str, 'PathLike[str]']]=None, hires_image_path: str='spatial/tissue_hires_image.png', lowres_image_path: str='spatial/tissue_lowres_image.png', scalefactors_path: str='spatial/scalefactors_json.json')
```

## spateo/io/spatial/_visium_hd.py

```python
def read_visium_hd_bin(path: Union[str, Path], sample: Optional[str]=None, binsize: int=16, count_h5_path: str='filtered_feature_bc_matrix.h5', count_mtx_dir: str='filtered_feature_bc_matrix', tissue_positions_path: str='spatial/tissue_positions.parquet', hires_image_path: str='spatial/tissue_hires_image.png', lowres_image_path: str='spatial/tissue_lowres_image.png', scalefactors_path: str='spatial/scalefactors_json.json')
def read_visium_hd_seg(path: Union[str, Path], sample: Optional[str]=None, cell_segmentations_path: str='graphclust_annotated_cell_segmentations.geojson', count_h5_path: str='filtered_feature_cell_matrix.h5', hires_image_path: str='spatial/tissue_hires_image.png', lowres_image_path: str='spatial/tissue_lowres_image.png', scalefactors_path: str='spatial/scalefactors_json.json')
def read_visium_hd(path: Union[str, Path], data_type: Literal['bin', 'cellseg']='bin', sample: Optional[str]=None, binsize: int=16, count_h5_path: str='filtered_feature_bc_matrix.h5', count_mtx_dir: str='filtered_feature_bc_matrix', tissue_positions_path: str='spatial/tissue_positions.parquet', cell_segmentations_path: str='graphclust_annotated_cell_segmentations.geojson', cell_matrix_h5_path: str='filtered_feature_cell_matrix.h5', hires_image_path: str='spatial/tissue_hires_image.png', lowres_image_path: str='spatial/tissue_lowres_image.png', scalefactors_path: str='spatial/scalefactors_json.json')
```

## spateo/io/spatial/_xenium.py

```python
def read_xenium(path: Union[str, Path], *, library_id: Optional[str]=None, load_image: bool=True, image_key: str='morphology_focus', image_max_dim: int=4096, load_boundaries: bool=True, cache_file: Optional[Union[str, Path]]=None)
```

## spateo/io/spatial/auto/_automatic.py

```python
def read_spatial(path: Union[str, Path], *, technology: Optional[str]=None, load: bool=True, load_images: bool=True, max_memory_bytes: int=_DEFAULT_MEMORY, max_files: int=10000, max_depth: int=4, stereoseq_bin_size: Optional[int]=None, stereoseq_chemistry: Optional[str]=None)
```

## spateo/io/spatial/auto/_contracts.py

```python
def read_core(candidate: Candidate, budget)
```

## spateo/io/spatial/auto/_stereo.py

```python
def read_gem(path, budget, bin_size=None, chemistry=None)
def read_gef(path, budget, bin_size=None, chemistry=None)
def read_stereo_core(path, budget, bin_size=None, chemistry=None)
```
