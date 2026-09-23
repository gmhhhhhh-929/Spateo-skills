# Point-cloud API contract

This skill is verified against `gmhhhhhh-929/spateo-release` commit
`615644f88613bea8ceb2e2df1e2391d16de55ec1`. Hashes are recorded in
[source_manifest.json](source_manifest.json). Reinspect the signatures and behavior before using another revision.

## Public calls

```python
pc, plot_cmap = st.tdr.construct_pc(
    adata,
    layer="X",
    spatial_key="spatial",
    groupby=None,
    key_added="groups",
    mask=None,
    colormap="rainbow",
    alphamap=1.0,
)
st.tdr.save_model(pc, "point_cloud.vtk", binary=True)
pc = st.tdr.read_model("point_cloud.vtk")
```

`construct_pc` returns a `pyvista.PolyData` and the colormap that a plotting call should use. The point model contains `point_data[key_added]` and `point_data["obs_index"]`. Categorical labels also materialize `point_data[f"{key_added}_rgba"]`; continuous values retain their scalar array and return the requested colormap instead.

The pinned implementation begins by copying the AnnData. Budget for that additional memory on large datasets; preview sampling does not reduce construction memory or the full VTK cardinality.

`save_model` accepts a string path and silently overwrites an existing file, so the wrapper creates only a new or empty output directory. A point cloud is a single `PolyData` and is saved as lowercase `.vtk`. Spateo requires `.vtm` for a `MultiBlock`; do not promise `.vtk` for a future multi-block artifact without first converting it to a single dataset and reviewing the semantic loss.

## Coloring modes

- `groupby=None`: uniform categorical label `same`.
- `groupby="obs_key"`: categorical or numeric values from `adata.obs`.
- `groupby="gene"`: numeric expression from `adata.X` or the selected layer.
- `groupby=("gene_a", "gene_b")`: row-wise sum of the selected genes. The current implementation fails for a Python list despite its docstring; use a tuple.
- External values are joined to an in-memory AnnData copy by exact observation ID and then passed through the obs route. The source H5AD is never rewritten.

For categorical data, use a complete label-to-color dictionary when stable biological colors matter. A named colormap assigns colors to sorted unique labels, so adding or renaming a category can change colors. A color list follows that same sorted category order. Numeric obs values are continuous; convert them to strings before invocation if they represent categories.

The current source treats an alpha list as one alpha per point, not one per category. The wrapper therefore accepts only a scalar or a label-to-alpha dictionary covering every visible category. Continuous values do not encode `alphamap`; set their opacity in the plotting call. The pinned implementation replaces masked category names with `mask` and lets a scalar `alphamap` overwrite their zero alpha. The wrapper instead builds the original categories without the source `mask` argument, supplies neutral internal color/alpha values for hidden categories omitted from dictionaries, then sets requested categories' RGBA alpha to zero. The VTK therefore retains the biological label while the manifest records requested and effective display mappings.

## Required validation

- `adata.obsm[spatial_key]` is numeric, finite, non-empty, and exactly `(n_obs, 3)`.
- Full coordinate rank is required by default. `--allow-planar` is an explicit exception for an intentionally planar dataset represented with three columns.
- `obs_names` are unique. Gene coloring additionally requires unique `var_names`, present genes, and a present layer.
- Labels are complete and finite. An external table has exactly one row for every observation ID and no extra IDs.
- Palette and alpha dictionaries cover every visible category. The literal label `mask` requires an explicit mask declaration because Spateo reserves it.

The wrapper saves with `st.tdr.save_model`, reloads with `st.tdr.read_model`, and compares point count, coordinates, point-data names and values before reporting success.
