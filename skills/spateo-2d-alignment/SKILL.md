---
name: spateo-2d-alignment
description: Align serial 2D spatial-transcriptomics slices with pairwise Spateo SN-S rigid alignment or annotation-guided continuity repair. Use for ordered slice H5AD inputs and preserved cell identities; not for developmental timepoint mapping or 3D morphogenesis.
---

# Spateo 2D slice alignment

Choose one of the two bundled pipelines. Keep the user's choice and explicit parameters; the shipped defaults are configurable starting points.

| Pipeline | Representation | Use |
| --- | --- | --- |
| `pipelines/pairwise-rigid` | Shared expression PCA, shared annotation one-hot, or identical ones30 vectors | Adjacent full-data SN-S rigid alignment, including a spatial-only control |
| `pipelines/continuity-guided` | Shared annotation one-hot | Full-data SN-S alignment followed by the original gated continuity repairs |

Read [input-contract.md](references/input-contract.md) before preparing or accepting inputs. Verify `X_pca` content, rather than inferring its meaning from its name. Expression PCs must come from a common gene space and one joint fitted basis for the specimen; annotation one-hot needs one shared dictionary across every slice. Numerical checks alone cannot establish expression PCA provenance.

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

For continuity repair, read [continuity-guided.md](references/continuity-guided.md), then run:

```bash
python /path/to/spateo-2d-alignment/pipelines/continuity-guided/run.py \
  --stage specimen --slice-dir /path/to/blind-slices \
  --output-dir /path/to/new-continuity-run \
  --annotation-key anno --postprocess auto --repair-policy auto --device 0
```

Both pipelines write into new output directories and retain native outputs. Check the frozen manifest, coordinate hash, all-cell identity coverage and slice mapping before delivering coordinates. Review continuity QC and accepted/rejected operations; an absent or rejected repair does not itself indicate execution failure. Use biological z from input metadata for downstream 3D assembly: the pairwise runner's display z is not a measured coordinate.

If reference-based accuracy is requested, freeze the candidate before accessing reference coordinates and use the user's evaluation convention. Do not tune the alignment against reference accuracy. This package includes no benchmark data or historical accuracy outputs.

Source hashes, dependency extraction and migration changes are in [source_migration.json](provenance/source_migration.json). The small [smoke test](scripts/smoke_test.py) validates packaging and synthetic input/geometry behavior; it does not repeat a GPU accuracy benchmark.
