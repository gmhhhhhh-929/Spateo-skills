# Shared expression PCA without annotation or reference coordinates

Use `scripts/prepare_expression_pca.py` when alignment should use expression
PCA and the existing per-slice `X_pca` is absent, unverified or annotation
one-hot. Annotation is optional. Every slice uses one joint expression basis
fitted to all and only the requested cells from a single biological specimen.
This interface does not pool species or developmental stages.

```bash
python scripts/prepare_expression_pca.py \
  --slice-dir /path/to/unaligned_or_perturbed_slices \
  --expression-source /path/to/original_expression.h5ad \
  --expression-layer auto --matrix-state counts \
  --output-dir /path/to/new_expression_preparation \
  --n-hvg 2000 --n-pcs 50 --seed 7021
```

The output directory must be new. Outputs are `slice_h5ad/`, `basis.npz` and
`expression_pca_manifest.json`. Use the new `slice_h5ad/` as alignment input.
Validators discover the manifest beside that directory, or accept its path
through `--pca-provenance`. Each output H5AD contains only `obsm['spatial']`
and `obsm['X_pca']`; no old feature, reference or 3D coordinate array is copied.
`X` stores the selected lognormalized expression as CSR. Cell IDs are preserved
in their original per-slice order, with slice metadata and physical z when it
was supplied.

When available, input slice z comes from `obs['physical_z']`, `obs['z']`, or
`obs['slice_z']`; `--physical-z-key` supplies another explicit column. Every
present candidate column is checked once, including all standard aliases and
the explicit column. Each must contain one finite constant per slice, and all
present columns must agree exactly, without a numerical tolerance. Conflicting
columns are rejected rather than silently replaced in the output. Z must also
be distinct between slices. All slices must either supply z or omit it; mixed
availability is rejected.

When all slices omit z, every input filename must contain one unique numeric
`SL` identifier, such as `sample_SL7.h5ad`. The numerical SL values determine
slice order, while no physical spacing is invented. The manifest records
`physical_z_available: false` and `slice_order_source: numeric_SL`; output
objects have slice indices but no `z` or `physical_z` obs column. PCA fitting
does not depend on either z or slice numbering. Missing, repeated or ambiguous
SL identifiers fail with an ordering error when no physical z is available.

Output files use `CS01_SL0000_YEXPRESSION.Spatial.h5ad` and
successive ordinals to satisfy the bundled legacy pairwise filename grammar.
`CS01` and `YEXPRESSION` are compatibility tokens, not biological identities.

By default no annotation field is read. `--annotation-key YOUR_COLUMN` can
copy an available slice annotation into the new output; it still does not
enter feature selection, normalization or PCA. Missing annotation is supported.
The expression source is never used to supply annotation.

## Relation to official Spateo joint PCA

The [official Spateo alignment tutorial](https://spateo-release.readthedocs.io/en/latest/tutorials/notebooks/3_alignment/1.%20Basic%20usage%20of%20Spateo%20alignment%20for%202D%20slices.html)
uses a shared expression representation between slices. The official
[`group_pca` implementation](https://spateo-release.readthedocs.io/en/latest/_modules/spateo/alignment/utils.html#group_pca)
accepts a list of slices: it concatenates them, computes one PCA and splits
the scores back into the individual objects. This naturally extends to more
than two slices. A conceptual call, after compatible preprocessing and common
gene handling, is:

```python
import spateo as st

# slices: a list of expression-bearing AnnData slices with compatible genes
st.align.group_pca(slices, pca_key="X_pca")
```

Our audited preparer follows the same **one shared basis across slices**
principle. It does not execute that official function and is not numerically
equivalent to it. The official implementation uses concatenation with the
default shared-gene join and batch-aware HVG handling. Normalization occurs
before that call: the linked tutorial normalizes each slice to its median
total counts and applies log1p. This preparer explicitly uses total-count
normalization to 10,000 and population-variance gene selection. It exposes and
hashes the common basis and uses sparse centered PCA while keeping reference
coordinates inaccessible. Do not describe its results as the direct output
of official `group_pca`. The implementation comparison was also checked against
source revision `d6aa68addc475dd0b56f69cebe7823b1f79933a9`.

## Matrix state and feature computation

`--expression-layer auto` checks `layers['counts_X']`, then `layers['counts']`,
then `X`. An explicit `counts_X`, `counts` or `X` selects just that matrix.
Failure of counts validation does not trigger fallback to another matrix.

* `--matrix-state counts`, the default, requires finite, nonnegative,
  integer-like counts and positive cell totals. Counts are normalized across
  all source genes to 10,000 per requested cell, then transformed once with
  natural `log1p`.
* `--matrix-state log1p` declares that the selected matrix already contains
  appropriate log1p expression. Values are used as-is, without library-size
  renormalization or another logarithm. The program never guesses an inverse
  transform. The user is responsible for correctly declaring the matrix state.

Gene selection ranks population variance over all requested cells, using gene
IDs to break ties. The selected matrix is centered through a sparse
`LinearOperator`; `scipy.sparse.linalg.svds` fits the shared basis once. No dense
cells-by-genes centered array, annotation, spatial coordinate, per-slice PCA
fit or per-gene unit-variance scaling is used.

Tiny datasets cap the requested PCs to
`min(requested_pcs, selected_genes - 1, requested_cells - 1)` and record both
requested and effective dimensions. Numerically zero singular directions are
removed, and at least two meaningful directions must remain. Zero-total cells,
strict one-hot input, constant expression, nonfinite PCA, zero-norm PCA cells
and constant PCA columns are rejected. Unverified preexisting `X_pca` is never
read or reused.

## Access and provenance contract

The original expression H5AD is opened read-only with `h5py`. Reads are limited
to its observation IDs, gene IDs and selected expression matrix. Dense, CSR
and CSC matrix storage are supported. CSR and dense inputs read the requested
rows; CSC inputs stream columns and retain only requested rows for fitting.
All source cell IDs may be indexed to match the requested cells, while extra
source cells never enter normalization, HVG selection or PCA.

Slice H5AD reads are limited to IDs, `obsm['spatial']` as finite N×2, physical z
and any explicitly requested annotation. Source coordinates, original 3D
coordinates, source annotations and old per-slice PCA are not read. The
manifest records actual HDF5 dataset accesses in `read_access_ledger`.

Original source and input slice files are **not** hashed as whole files,
because that would read unapproved reference bytes. Instead, allowed matrix,
ID, gene, XY, z and optional annotation contents receive deterministic hashes,
and file size/mtime are recorded and checked for changes. New sanitized H5ADs
and `basis.npz` receive full-file SHA-256 hashes. The manifest states these
different hash scopes explicitly rather than presenting a dataset hash as a
whole-file hash.

`basis.npz` contains selected gene IDs, gene means, loadings, singular values,
explained variance and canonical fit-cell IDs. Its whole-file hash is repeated
in every slice record. Feature hashes use contiguous little-endian float32
scores in row-major order; identity hashes use sorted compact UTF-8 JSON.
Output H5AD hashes bind identities, row order, features and spatial coordinates
together. The manifest records `mode: expression-pca`, `shared_basis: true`,
`fitted_on_spatial: false` and `fitted_on_annotation: false`.

Run the synthetic contract tests with:

```bash
cd scripts
python -m unittest -v test_prepare_expression_pca.py
```

They verify a genuinely joint basis, ID mapping under reordered/superset
expression sources, no-annotation operation, no-z/no-annotation outputs through
both validators, mixed-z and ambiguous-order rejection, consistent-z aliases,
conflicting explicit/alias z fields and invalid secondary-z rejection, dense/CSR/CSC support,
bytewise unchanged inputs, explicit log1p handling and invalid-input rejection. Dataset
read traps additionally reject attempted access to source coordinates, source
annotations, old PCA, reference arrays or whole original-file hashing.
