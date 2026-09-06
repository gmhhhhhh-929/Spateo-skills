# 固定提交的精确 API 签名

源码：`d6aa68addc475dd0b56f69cebe7823b1f79933a9`；签名从 AST 提取，未导入 optional reader 依赖。此表用于核对参数；语义与平台陷阱见 [platforms.md](platforms.md)。

## spateo/io/general/_serialization.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/general/_serialization.py)

```python
save(file, path)
```

```python
load(path, backend=None)
```

## spateo/io/general/_tabular.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/general/_tabular.py)

```python
read_csv(**kwargs)
```

## spateo/io/single/_formats.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/single/_formats.py)

```python
read_h5ad(filename, **kwargs)
```

```python
read_10x_h5(filename, *, genome: str | None=None, gex_only: bool=True, backup_url: str | None=None) -> AnnData
```

```python
read_10x_mtx(path, *, var_names: Literal['gene_symbols', 'gene_ids']='gene_symbols', make_unique: bool=True, gex_only: bool=True, prefix: str | None=None, compressed: bool=True) -> AnnData
```

## spateo/io/single/_read.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/single/_read.py)

```python
read(path, backend='python', **kwargs)
```

## spateo/io/spatial/_atera.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_atera.py)

```python
read_atera(path: PathLike, *, library_id: Optional[str]=None, load_image: bool=True, image_key: str='dapi', image_max_dim: int=4096, load_boundaries: bool=True, load_nucleus_boundaries: bool=True, load_cell_groups: bool=True, cell_groups_csv: Optional[PathLike]=None, load_he_image: bool=False, he_image: Optional[PathLike]=None, he_alignment_csv: Optional[PathLike]=None, he_max_dim: int=4096, cache_file: Optional[PathLike]=None) -> AnnData
```

## spateo/io/spatial/_geometry.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_geometry.py)

```python
alpha_shape(x: np.ndarray, y: np.ndarray, alpha: float=1, buffer: float=1, vectorize: bool=True) -> Tuple[Union[MultiPolygon, Polygon], List]
```

```python
get_concave_hull(path: str, binsize: int=20, min_agg_umi: Optional[int]=None, alpha: float=1.0, buffer: Optional[float]=None) -> Tuple[Union[MultiPolygon, Polygon], List]
```

## spateo/io/spatial/_image.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_image.py)

```python
add_image_layer(adata: AnnData, img: np.ndarray, scale_factor: float, slice: Optional[str]=None, img_layer: Optional[str]=None) -> AnnData
```

```python
read_image(adata: AnnData, filename: str, scale_factor: float, slice: Optional[str]=None, img_layer: Optional[str]=None) -> AnnData
```

## spateo/io/spatial/_merfish.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_merfish.py)

```python
read_merfish(path: Union[str, Path], *, counts_file: str, meta_file: str, load_images: bool=True, z_layer: Optional[int]=3, load_boundaries: bool=True) -> AnnData
```

## spateo/io/spatial/_nanostring.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_nanostring.py)

```python
read_nanostring(path: Union[str, Path], *, counts_file: str, meta_file: str, fov_file: Optional[str]=None) -> AnnData
```

## spateo/io/spatial/_provenance.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_provenance.py)

```python
spatial_file_manifest(path: PathLike, *, max_files: int=10000) -> dict[str, Any]
```

## spateo/io/spatial/_seqfish.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_seqfish.py)

```python
read_seqfish(path: Union[str, Path], *, counts_file: str, meta_file: str, load_images: bool=True, load_labels: bool=True) -> AnnData
```

## spateo/io/spatial/_seqscope.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_seqscope.py)

```python
read_seqscope(matrix_dir: PathLike, positions_path: PathLike, binsize: Optional[int]=1, add_props: bool=True) -> AnnData
```

## spateo/io/spatial/_slideseq.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_slideseq.py)

```python
read_slideseq(path: Union[str, Path], *, counts_file: str='MappedDGEForR.csv', bead_file: str='BeadLocationsForR.csv', load_images: bool=True) -> AnnData
```

## spateo/io/spatial/_starmap_plus.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_starmap_plus.py)

```python
read_starmap_plus(path: Union[str, Path], *, counts_file: str, meta_file: str, spatial_file: str, reorient_xy: bool=False, dtype: str='float32') -> AnnData
```

## spateo/io/spatial/_stereoseq.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_stereoseq.py)

```python
read_bgi_as_dataframe(path: Union[str, PathLike], label_column: Optional[str]=None) -> pd.DataFrame
```

```python
read_bgi_agg(path: Union[str, PathLike], stain_path: Optional[Union[str, PathLike]]=None, binsize: int=1, gene_agg: Optional[Dict[str, Union[List[str], Callable[[str], bool]]]]=None, prealigned: bool=False, label_column: Optional[str]=None, version: Literal['stereo']='stereo') -> AnnData
```

```python
read_bgi(path: Union[str, PathLike], binsize: Optional[int]=None, segmentation_adata: Optional[AnnData]=None, labels_layer: Optional[str]=None, labels: Optional[Union[np.ndarray, str, PathLike]]=None, seg_binsize: int=1, label_column: Optional[str]=None, add_props: bool=True, version: Literal['stereo']='stereo') -> AnnData
```

## spateo/io/spatial/_utils.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_utils.py)

```python
bin_indices(coords: np.ndarray, coord_min: float, binsize: int=50) -> int
```

```python
centroids(bin_indices: np.ndarray, coord_min: float=0, binsize: int=50) -> float
```

```python
contour_to_geo(contour)
```

```python
get_points_props(data: pd.DataFrame) -> pd.DataFrame
```

```python
get_label_props(labels: np.ndarray) -> pd.DataFrame
```

```python
get_bin_props(data: pd.DataFrame, binsize: int) -> pd.DataFrame
```

```python
in_concave_hull(p: np.ndarray, concave_hull: Union[Polygon, MultiPolygon]) -> np.ndarray
```

```python
in_convex_hull(p: np.ndarray, convex_hull: Union[Delaunay, np.ndarray]) -> np.ndarray
```

```python
bin_matrix(X: Union[np.ndarray, spmatrix], binsize: int) -> Union[np.ndarray, csr_matrix]
```

```python
get_coords_labels(labels: np.ndarray) -> pd.DataFrame
```

## spateo/io/spatial/_visium.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_visium.py)

```python
read_visium(path: Union[str, 'PathLike[str]'], genome: Optional[str]=None, *, count_file: str='filtered_feature_bc_matrix.h5', library_id: Optional[str]=None, load_images: bool=True, load_qc_images: bool=False, source_image_path: Optional[Union[str, 'PathLike[str]']]=None, hires_image_path: str='spatial/tissue_hires_image.png', lowres_image_path: str='spatial/tissue_lowres_image.png', scalefactors_path: str='spatial/scalefactors_json.json') -> AnnData
```

## spateo/io/spatial/_visium_hd.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_visium_hd.py)

```python
read_visium_hd_bin(path: Union[str, Path], sample: Optional[str]=None, binsize: int=16, count_h5_path: str='filtered_feature_bc_matrix.h5', count_mtx_dir: str='filtered_feature_bc_matrix', tissue_positions_path: str='spatial/tissue_positions.parquet', hires_image_path: str='spatial/tissue_hires_image.png', lowres_image_path: str='spatial/tissue_lowres_image.png', scalefactors_path: str='spatial/scalefactors_json.json') -> AnnData
```

```python
read_visium_hd_seg(path: Union[str, Path], sample: Optional[str]=None, cell_segmentations_path: str='graphclust_annotated_cell_segmentations.geojson', count_h5_path: str='filtered_feature_cell_matrix.h5', hires_image_path: str='spatial/tissue_hires_image.png', lowres_image_path: str='spatial/tissue_lowres_image.png', scalefactors_path: str='spatial/scalefactors_json.json') -> AnnData
```

```python
read_visium_hd(path: Union[str, Path], data_type: Literal['bin', 'cellseg']='bin', sample: Optional[str]=None, binsize: int=16, count_h5_path: str='filtered_feature_bc_matrix.h5', count_mtx_dir: str='filtered_feature_bc_matrix', tissue_positions_path: str='spatial/tissue_positions.parquet', cell_segmentations_path: str='graphclust_annotated_cell_segmentations.geojson', cell_matrix_h5_path: str='filtered_feature_cell_matrix.h5', hires_image_path: str='spatial/tissue_hires_image.png', lowres_image_path: str='spatial/tissue_lowres_image.png', scalefactors_path: str='spatial/scalefactors_json.json') -> AnnData
```

```python
write_visium_hd_cellseg(adata: AnnData, path: Union[str, Path], sample: Optional[str]=None)
```

## spateo/io/spatial/_xenium.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/_xenium.py)

```python
read_xenium(path: Union[str, Path], *, library_id: Optional[str]=None, load_image: bool=True, image_key: str='morphology_focus', image_max_dim: int=4096, load_boundaries: bool=True, cache_file: Optional[Union[str, Path]]=None) -> AnnData
```

## spateo/io/spatial/auto/__init__.py

[固定提交源码](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/spateo/io/spatial/auto/__init__.py)

```python
detect_spatial_technologies(path: PathLike, *, technology: Optional[str]=None, min_confidence: float=0.0) -> List[SpatialReadMatch]
```

```python
detect_spatial_technology(path: PathLike, *, technology: Optional[str]=None, min_confidence: float=0.5, strict: bool=True) -> SpatialReadMatch
```

```python
read_auto_spatial(path: PathLike, *, technology: Optional[str]=None, min_confidence: float=0.5, strict: bool=True, return_match: bool=False, **reader_kwargs: Any) -> Any
```

`SpatialReadMatch` 另有 `.load_reader()` 与 `.read(**overrides)`；`read_spatial_auto = read_auto_spatial`。

`spateo/data_io.py` 本身没有 reader 实现函数：由 anndata 兼容导出。源文件 SHA256 见 [source_manifest.json](source_manifest.json)。
