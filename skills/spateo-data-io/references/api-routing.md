# API routing

Use `import spateo as st`. Current source is pinned in source_manifest.json.

| Input | Entry point | Contract |
| --- | --- | --- |
| Supported spatial bundle / GEM / GEF | `st.io.read_spatial(path, ...)` | SpatialReadResult with named entries and diagnostics. |
| Existing H5AD, preserve metadata | `anndata.read_h5ad(path)` | AnnData with its stored metadata; the skill CLI uses this for H5AD conversion. |
| Spateo direct H5AD reader | `st.io.read_h5ad(path)` | Adds Spateo UMI/pp metadata; do not use this to infer that an AGG is UMI data. |
| 10x H5 | `st.io.read_10x_h5(path, genome=None, gex_only=True)` | Features-by-barcodes storage becomes observations-by-genes; direct default filters non-Gene-Expression features. |
| 10x MTX | `st.io.read_10x_mtx(path, var_names='gene_symbols', make_unique=True, compressed=True)` | MTX, feature and barcode files; uncompressed v3 bundles need `compressed=False`. Preserve stable feature IDs if symbols change. |
| CSV / TSV | `st.io.read(str(path))` or `st.io.read_csv(filepath_or_buffer=..., sep=...)` | DataFrame, not spatial AnnData; CSV helper accepts keyword arguments. |
| Known platform, custom filenames | `st.io.read_<platform>(...)` | AnnData; direct readers retain their own format and loading semantics. Do not pass the unified automatic reader options indiscriminately. |
| Stereo-seq | `st.io.read_stereoseq(path, ...)` | Native GEM/GEF core reader; inspect bin/cell representation and recorded units. |
| Legacy BGI aggregate | `st.io.read_bgi_agg(...)` | Spatial image grid; not expression-matrix input for QC/alignment. |

`st.read_h5ad` is an AnnData compatibility export and differs from `st.io.read_h5ad`. `st.io.save/load` persist Python objects, not H5AD; use `adata.write_h5ad` for the pipeline. Read only trusted pickle files.

Automatic aliases return the result container; there is no `(adata, match)` return mode. Old flat technology modules and scored detector calls must not be copied from historical examples. Use the maintained `spateo/io/{general,single,spatial}` exports.
