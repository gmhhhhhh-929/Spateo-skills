---
name: spateo-data-io
description: Read spatial-platform outputs, including domestic platforms, with source-verified Spateo IO; use automatic or lazy loading, explain recovery actions, and export validated AnnData before slice QC or alignment.
---

# Spateo data IO

Stage 2: environment → **IO** → slice quality (including its viewer) → 2D alignment → 3D pipeline → 4D pipeline. The current 3D implementation begins with a validated point-cloud VTK and can continue to full-body or annotation-specific surface meshes.

Use the source-verified API at `gmhhhhhh-929/spateo-release`, pinned in [source_manifest.json](references/source_manifest.json). Read [API routing](references/api-routing.md) for explicit readers, [automatic outcomes](references/auto-and-errors.md) for lazy loading and recovery, and the relevant row in [platform contracts](references/platforms.md). Exact signatures and IO hashes are in [source-api.md](references/source-api.md). Check signatures again for another revision.

## Choose the route

- For a spatial directory or Stereo-seq GEM/GEF, use `st.io.read_spatial(path)`. It returns **SpatialReadResult**, not AnnData. All samples and representations remain named entries; it does not select a highest-scoring platform.
- For discovery without core loading, use `st.io.read_spatial(path, load=False, load_images=False)` and inspect `result.report`. There is no public `detect_spatial_technologies` API in this revision.
- For on-demand materialization, use `st.io.read_spatial(path, lazy=True, load_images=False)`. Reports and iteration stay metadata-only; access `result.adata` for a unique complete scope, or call a selected `entry.materialize()`. This defers full loading; it is not backed or out-of-core AnnData.
- For H5AD, use `anndata.read_h5ad` to preserve existing type/metadata. For 10x expression-only H5/MTX or a known nonstandard layout, use an explicit reader and inspect its signature. Expression matrices alone do not supply spatial coordinates.
- Generic CSV/TSV reading returns a DataFrame. Join expression, metadata and coordinates by verified IDs before constructing AnnData; equal row counts do not establish identity.
- For known SeekSpace, BMKMANU, Salus STS or CeleScope space bundles, each `st.io.read_<tech>` invokes its independent platform module. Automatic reading calls the same platform core. See [domestic reader routing](references/domestic-platforms.md#independent-readers-and-shared-core) when tracing implementation or comparing both routes.

## Inspect outcomes, then export

```python
import spateo as st

result = st.io.read_spatial('/data/sample/outs', load_images=False)
print(result.report)
# This property is valid only for exactly one dataset and a complete, successful scope.
if result.status == 'ok' and len(result.datasets) == 1:
    adata = result.adata
# For collections inspect each entry instead; do not concatenate independent specimens.
for key, entry in result.datasets.items():
    print(key, entry.technology, entry.representation, entry.status, entry.to_dict()['diagnostics'])
```

`ready` means the core input contract passed. Image warnings, chemistry declarations, physical units and biological suitability still need review. Keep failed, unresolved and deferred entries visible. Inspect each diagnostic's `recovery` actions in `result.report` or `entry.to_dict()` before retrying. Lazy materialization and `entry.load(max_memory_bytes=...)` use the same resolved adapter and can remain deferred. Avoid repeatedly materializing unchanged failures; consult [retry behavior](references/auto-and-errors.md). A memory budget is not an OS memory limit. Automatic aliases `read_auto_spatial` and `read_spatial_auto` have this same result contract; old `strict`, `min_confidence`, `return_match` and arbitrary reader kwargs are invalid.

Run helpers from this skill directory:

```bash
python scripts/spateo_io.py discover /data/sample/outs --technology visium
python scripts/spateo_io.py read /data/collection --lazy --dataset-key EXACT_DISCOVERED_KEY \
  --output-dir /output/one-sample
python scripts/spateo_io.py read /data/sample/outs --output-dir /output/io-run \
  --matrix-semantics counts --coordinate-unit fullres_pixel
python scripts/spateo_io.py convert /data/sample.h5ad /output/copy.h5ad --reader read_h5ad
python scripts/validate_anndata.py /output/copy.h5ad --require-spatial \
  --matrix-semantics counts --coordinate-unit um
```

`read` writes a complete diagnostic report plus independently named H5AD files for ready entries. `--dataset-key` restricts exports; adding `--lazy` also limits materialization to those keys. It does not change the full-scope status. Exit code 2 means incomplete/failed scope or export, even if some H5ADs were saved. Outputs never overwrite existing paths. Each successful file is reread and compared across X, identities, layers, coordinates, raw and metadata/images before publication. A serialization failure stays an export failure; do not delete metadata to hide it. A saved report is an audit, not a resumable loader; create new handles from the native source in a later process.

For repeated gene symbols, inspect the stable ID metadata and explicitly use `--feature-id-column gene_ids` (or the actual unique column). This changes var_names while retaining their original values in `var['source_var_name']`; it is recorded in the output manifest. The wrapper refuses ambiguous IDs or an unhandled `.raw` feature migration. It never invents suffixes to disguise duplicate identities.

## Handoff contract

These inputs are standard downstream files from vendor alignment/counting or decoding pipelines, not FASTQ/BCL processing. A shared MEX/H5 layout does not establish the platform; require the source contract or a supported explicit declaration. Do not rename a generic matrix to manufacture platform evidence. Read [domestic platforms](references/domestic-platforms.md) for verified formats, evidence and known support limits; do not assume lower cost or a vendor's unspecified delivery format.

Confirm observation unit (cell/spot/bead/bin), feature IDs, source expression semantics, complete ID coverage, finite XY/XYZ, slice/FOV identity, units and frame, images and vendor transformations. Native reports inventory paths and sizes; they are not content hashes. The wrapper hashes successfully read explicit single-file inputs and output H5ADs; discovery never scans an entire file to compute a hash. Directory-content hashes are explicitly incomplete. Record hashes of important source bundle files when reproducibility requires them.

For serial slices, deliver per-specimen data and actual slice order to `spatial-slice-quality-qc`; its nested viewer belongs to the same QC stage. Follow [alignment-contract.md](references/alignment-contract.md) for later expression PCA or annotation modes. Do not invent annotations, physical z, coordinate conversions or registration. Preserve raw counts separately from derived expression. BGI AGG is an XY image grid, not a cell-by-gene matrix.

## Verification

Use [validation.md](references/validation.md) for reproducible source and skill tests and their limits. Real API synthetic tests establish execution/contracts, not biological accuracy or universal vendor-format support.
