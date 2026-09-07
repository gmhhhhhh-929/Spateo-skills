---
name: spateo-2d-alignment
description: Prepare shared expression PCA across multiple slices and align serial spatial-transcriptomics data with Spateo rigid alignment and gated continuity repair, with or without tissue annotations. Use for ordered slice H5AD inputs and preserved cell identities; not developmental timepoint mapping or morphogenesis.
---

# Spateo 2D slice alignment

Choose one of the two bundled pipelines. Keep the user's choice and explicit parameters; the shipped defaults are configurable starting points.

| Pipeline | Representation | Use |
| --- | --- | --- |
| `pipelines/pairwise-rigid` | Shared expression PCA, shared annotation one-hot, or identical ones30 vectors | Adjacent full-data SN-S rigid alignment, including a spatial-only control |
| `pipelines/continuity-guided` | Shared expression PCA, shared annotation one-hot, or identical ones30 vectors | Full-data SN-S alignment followed by gated continuity repairs; annotations are optional in expression mode |

Read [input-contract.md](references/input-contract.md) before preparing or accepting inputs. Verify `X_pca` content, rather than inferring its meaning from its name. Expression PCs must come from a common gene space and one joint fitted basis for the specimen; annotation one-hot needs one shared dictionary across every slice. Numerical checks alone cannot establish expression PCA provenance.

Choose `expression-pca` when only expression is available or annotations are unreliable. Choose `annotation-onehot` when the user wants trusted tissue labels to define correspondence. Keep these as explicit alternatives: do not silently replace the chosen representation or use unreliable labels for repair. Expression mode defaults to annotation QC off and can run with no annotation column. `--annotation-qc provided` explicitly enables the selected labels for continuity QC.

For expression inputs, read [expression-pca.md](references/expression-pca.md). The preparation script fits **one PCA across all slices of the specimen**, then assigns scores back by original cell ID. It does not fit separate slice or adjacent-pair bases. Example:

```bash
python /path/to/spateo-2d-alignment/scripts/prepare_expression_pca.py \
  --slice-dir /path/to/unaligned-slices \
  --expression-source /path/to/original-expression.h5ad \
  --matrix-state counts --output-dir /path/to/new-expression-inputs
```

Use the generated `slice_h5ad` directory in either pipeline. Its neighboring `expression_pca_manifest.json` authenticates the common gene space, basis and each slice's features. The source reader opens expression and identifiers only. It preserves the supplied unaligned XY and never reads source reference coordinates or an existing ambiguous `X_pca`. Already log-normalized expression requires the explicit `log1p` state and correct layer; do not normalize or log twice. Input adaptation preserves every cell; it does not silently filter zero-count cells.

Use an existing compatible Spateo interpreter and run `scripts/check_runtime.py --device <cpu-or-index>`. See [runtime.md](references/runtime.md) for dependencies and the limits of this check. For remote execution, use the known host, interpreter, input directory and a new output path; ask only for missing execution context that prevents proceeding.

Prepare new blind H5AD copies when inputs need adaptation. Preserve cell IDs, existing unaligned XY values and biological slice metadata. Keep reference/aligned XY arrays out of the inference inputs. The bundled preflight requires exactly the selected spatial and representation keys in `obsm`; annotation values belong in `obs`. Order filenames by physical z when available. The legacy runner requires names such as `CS01_SL001_YSAMPLE.Spatial.h5ad`; generic names without `SL` are not accepted.

For pairwise alignment, read [pairwise-rigid.md](references/pairwise-rigid.md), then run:

```bash
python /path/to/spateo-2d-alignment/pipelines/pairwise-rigid/run.py \
  --slice-dir /path/to/blind-slices \
  --output-dir /path/to/new-pairwise-run \
  --representation expression-pca --device 0
```

Use `--representation annotation-onehot` with the intended `--annotation-key`, or `--representation spatial-only` with constant ones30 inputs. `--dry-run` prints the exact invocation without reading data or running alignment.

For continuity repair, read [continuity-guided.md](references/continuity-guided.md). With shared expression PCA and no annotation dependency:

```bash
python /path/to/spateo-2d-alignment/pipelines/continuity-guided/run.py \
  --stage specimen --slice-dir /path/to/new-expression-inputs/slice_h5ad \
  --output-dir /path/to/new-continuity-run \
  --representation expression-pca --annotation-qc off \
  --profile generalized --postprocess auto --device 0
```

Use `--representation annotation-onehot --annotation-key anno` for annotation correspondence. Generalized continuity uses expression-neighbor initialization off and anonymous multi-label terminal acceptance; it has no ROD preference. The compatibility default is `legacy`, preserving the previous initialization and tissue-priority behavior. The generalized profile remains opt-in because frozen validation found average Drosophila gains alongside Planarian and individual-sample regressions. See [validation.md](references/validation.md); do not silently choose a profile using reference accuracy.

Both pipelines write into new output directories and retain native outputs. Check the frozen manifest, coordinate hash, all-cell identity coverage and slice mapping before delivering coordinates. Review continuity QC and accepted/rejected operations; an absent or rejected repair does not itself indicate execution failure. Insufficient annotation evidence skips annotation-dependent proposals; it is not permission to invent labels. Use biological z from input metadata for downstream 3D assembly: the pairwise runner's display z is not a measured coordinate.

If reference-based accuracy is requested, freeze candidates before evaluation and use the user's convention. Keep reference coordinates and scores out of inference. For authorized method development, designate development and test specimens before scoring, record all tried configurations, and freeze one policy before viewing new test scores. New perturbations of the same specimen test perturbation robustness, not independent biological generalization. This package includes no benchmark data.

Source hashes, dependency extraction and migration changes are in [source_migration.json](provenance/source_migration.json). The small [smoke test](scripts/smoke_test.py) validates packaging and synthetic input/geometry behavior; it does not repeat a GPU accuracy benchmark.

## Review and workflow handoff

For inspection, use [spatial-before-after-viewer](../spatial-before-after-viewer/SKILL.md) and [spatial-pointcloud-viewer](../spatial-pointcloud-viewer/SKILL.md). For pair QC, component review, sampling and recorded replay, use [lightweight-spatial-alignment-workflow](../lightweight-spatial-alignment-workflow/SKILL.md). These companion tools restore the fixed pre-zebrafish support release; they do not change either alignment pipeline's inference code.
