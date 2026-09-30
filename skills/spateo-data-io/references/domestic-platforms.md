# Domestic spatial-platform outputs

These readers consume alignment/counting/decoding outputs. They do not align FASTQ/BCL, decode microscopy, identify cells from stains, or prove a platform is cheaper. Vendor names are retained separately from shared storage formats. Read the exact installed signatures before passing direct-reader arguments.

## Supported native contracts

| Technology / producer | Files and required fields | Reader and automatic evidence | Coordinate interpretation |
| --- | --- | --- | --- |
| Stereo-seq / STOmics, BGI | GEM/GEM2 text with gene, x, y, counts; supported bin/cell GEF schemas | `read_stereoseq`; native GEM fields or GEF schema. `read_spatial` is automatic. | Preserve source bin/cell and available offset/resolution metadata. V1/V2 chemistry cannot be inferred from the container schema. |
| SeekSpace / SeekGene | Standard MEX trio (`matrix.mtx[.gz]`, `features.tsv[.gz]`, `barcodes.tsv[.gz]`) plus `cell_locations.tsv[.gz]` with `Cell_Barcode`, `X`, `Y` | `read_seekspace`; distinctive cell-location fields and complete MEX contract | Source chip/image pixels; no automatic conversion by a presumed 0.2653 µm factor. |
| BMKMANU S1000 / Biomarker | **Aggregated** MEX trio plus plural `barcodes_pos.tsv[.gz]`: three headerless columns barcode, pos_w, pos_h | `read_bmkmanu`; aggregated position companion and validated barcode join | Native display XY. Physical scale is not established solely by these files. |
| Salus STS | Standard MEX trio plus `spatial.txt[.gz]`: headerless whitespace-delimited barcode, x, y | `read_salus`; position companion and full ID/shape/value contract | Native source pixels; preserve supplied orientation. |
| Singleron CeleScope space export | `filtered_feature_bc_matrix.h5` and/or `raw_feature_bc_matrix.h5`, plus `spatial/positions_list.csv` with six Visium-style position columns; H5 chemistry description `Spatial3` | `read_singleron`; distinctive chemistry metadata and complete position contract for auto detection | `x=pxl_col_in_fullres`, `y=pxl_row_in_fullres`, full-resolution image pixels. |

All four new matrix readers align positions to the matrix's barcode IDs. Stable feature IDs are the variable index; gene-name metadata is preserved. Repeated gene symbols are not made into new identities by adding suffixes. A generic MEX/H5 without distinctive source evidence cannot uniquely identify a vendor.

Images and scale metadata are optional: the common asset loader can retain compatible single-frame raster images and recognized scale JSON within the resource budget. A vendor's complete output directory may also contain masks, boundaries, alignments and web reports; their existence does not mean this core reader imports them into AnnData. Review the actual asset status in `uns['spateo_io']` and do not imply automatic count-image registration.

**BMKMANU raw export is a different contract.** Singular `barcode_pos.tsv` can contain five columns: barcode, subarea_col_index, subarea_row_index, col_index, row_index. Do not reinterpret these as three-column XY. Use the producer's BSTMatrix aggregation/export procedure to obtain the supported plural `barcodes_pos.tsv.gz` and matching MEX files, preserving the selected bin/aggregation level.

## Primary evidence and acquisition

- SeekGene code writes the MEX bundle and `cell_locations.tsv.gz`: [SeekSpaceTools report.py, pinned source](https://github.com/seekgene/SeekSpaceTools/blob/2434074a53117ad8e56136566e07d91b09ed1c9c/src/seekspacetools/run/report.py#L377).
- Biomarker's output guide: [BMKMANU spatial export documentation](https://www.biomarker.com.cn/archives/28318). Confirm the actual delivered export version; raw and aggregated position files are not interchangeable.
- Salus's paper-linked repository: [SalusSTS](https://github.com/xuzaoxu/SalusSTS). Obtain expression and spatial companions for the same sample from the repository/study links.
- Singleron's maintained pipeline: [CeleScope](https://github.com/singleron-RD/CeleScope). This support applies to the verified space export contract, not every unrelated Singleron assay.
- Stereo-seq official example source: [MOSTA downloads](https://db.cngb.org/stomics/mosta/download/). These datasets establish the original Stereo-seq example; they are not V2 datasets.

The source repository's [domestic-format documentation](https://github.com/gmhhhhhh-929/spateo-release/blob/d884216b2f030bcb5616b01a8f950d31e1e6d833/docs/technicals/domestic_spatial_io_zh.md) and [validation record](https://github.com/gmhhhhhh-929/spateo-release/blob/d884216b2f030bcb5616b01a8f950d31e1e6d833/docs/technicals/domestic_spatial_io_validation_20260930.json) record pinned evidence and verified public examples in more detail. The public BMKMANU example is [GEO GSM8816652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8816652): native MEX, coordinates and an image, with no H5AD input. Link users to a dataset-specific official download page when one is known; a vendor homepage is guidance, not a guessed data URL. Synthetic contract tests do not establish compatibility with every commercial delivery version.

## Automatic and explicit examples

```python
import spateo as st

# A complete native output directory supplies the evidence.
result = st.io.read_spatial('/data/vendor/sample-output', load_images=False)
print(result.report)

# Leave all matrices pending until a chosen entry is requested.
result = st.io.read_spatial('/data/collection', lazy=True, load_images=False)
for key, entry in result.datasets.items():
    print(key, entry.technology, entry.status)

# Explicit readers are useful after the export contract is known.
adata = st.io.read_seekspace('/data/seekspace/sample')
# When one directory has several representations, request all outcomes:
result = st.io.read_singleron('/data/singleron/outs', return_result=True)
```

The four new explicit readers accept `load_images`, `max_memory_bytes` and `return_result`; they do not accept arbitrary filename or reader kwargs. Their default returns a unique successful AnnData. Use `return_result=True` to inspect all named outcomes when several samples or raw/filtered representations are present. Use the common `read_spatial(..., lazy=True)` entry point for lazy loading.

Do not manually select the platform merely to force an unrelated layout past discovery. If a producer delivers only a processed H5AD, treat it as a separately labelled H5AD input, not a demonstration of the native downstream-file reader.
