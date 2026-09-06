# Input contract

Use one H5AD per slice and one specimen per invocation. Slice cell counts can differ. `obs_names` must be nonempty, unique within each slice, and globally unique across the specimen. Preserve a map if an upstream import must namespace duplicate raw IDs.

The H5AD contains `obsm['spatial']` as finite N×2 unaligned XY and `obsm['X_pca']` as finite N×D features with nonzero row norms. The continuity CLI can select other key names. The pairwise runner fixes these two names. Preflight rejects extra `obsm` entries, including reference/aligned spatial arrays, and explicit reference-coordinate columns in `obs` before loading the selected arrays. Keep original data unchanged; create a sanitized input copy instead.

## Slice names and physical order

- Pairwise names must match `CS<digits>_SL<digits>_Y<word-characters>.Spatial.h5ad`, for example `CS01_SL001_YSAMPLE.Spatial.h5ad` and `CS01_SL002_YSAMPLE.Spatial.h5ad`. All files belong to the same CS stage; SL numbers must be distinct. Use a flat directory: preflight rejects any nested H5AD because the preserved pairwise loader recursively discovers files.
- Continuity names must contain `SL<digits>`, for example `SL001.h5ad`; the same pairwise naming convention works for both.
- Both original algorithms sort by numeric SL, not by directory enumeration or `obs['slice_ID']`. Assign SL numbers in ascending physical-z order when exporting the per-slice copies. Retain the original slice identifier as metadata.
- Optional `obs['z']` or `obs['physical_z']` must be finite and constant within each slice. If both exist they must agree. If present, preflight requires physical z to increase strictly in SL order. If absent, it explicitly reports that only SL order was checked. No coordinates are rescaled or reordered silently.
- `slice_ID`, `slice_id`, and `physical_z` are useful metadata but are not required inputs to either original algorithm. They do not override filename order. Store the slice-to-z table for downstream 3D assembly.

## Feature meanings

| Representation | Contract |
| --- | --- |
| `expression-pca` | One joint PCA basis across all the specimen's slices, fitted from documented expression/count preprocessing and a common ordered gene space. The basis and feature hashes must match `expression_pca_manifest.json`. Independent fits per slice or adjacent pair, and annotation one-hot stored under `X_pca`, are invalid. Annotation columns are optional. |
| `annotation-onehot` | Exact 0/1 one-hot rows with a single active column. One dictionary across all slices; the same `obs[annotation_key]` label always has the same vector and distinct labels have distinct vectors. Labels must be nonblank and nonmissing. UTF-8 cell IDs and labels are supported. Choose the annotation field requested for the task. |
| `spatial-only` | Every cell has the identical 30-dimensional vector of ones. Pairwise spatial-only mode also constructs this representation inside the original loader; the preflight makes the input contract explicit. |

The pairwise native loader creates temporary IDs `<filename-slice-id>:<row-index>`. The packaged CLI saves the input ID map and restores original IDs in `aligned_coordinates.csv.gz`, checking native raw XY and slice membership. Continuity keeps original `obs_names` and writes numeric SL identifiers in its native coordinate table.

## Preflight alone

```bash
python /path/to/spateo-2d-alignment/pipelines/pairwise-rigid/validate_inputs.py \
  --slice-dir /path/to/blind-slices --representation annotation-onehot --annotation-key anno
```

For the continuity directory's equivalent checker, use `--filename-style continuity`. Both public run CLIs perform preflight automatically. Expression mode finds `expression_pca_manifest.json` in the slice directory's parent, or accepts `--pca-provenance`. It verifies a shared fitted basis and per-slice feature/identity/file hashes. Preflight cannot reconstruct undocumented PCA provenance; use `scripts/prepare_expression_pca.py` and the authorized expression source when it is missing or contradictory.

The PCA preparation script accepts multiple input slice H5ADs and a single expression source spanning their cell IDs. It fits only the requested cells together, maps by ID rather than row position, writes new files and preserves XY and available physical Z. The expression source may contain additional cells; those do not participate in the PCA fit. For separately stored expression sources, first construct one common-gene expression table with globally unique IDs and a recorded gene intersection; do not join independent PC scores.
