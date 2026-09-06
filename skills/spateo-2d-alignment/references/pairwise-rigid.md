# Pairwise rigid mode

`pipelines/pairwise-rigid/run.py` is the portable CLI. The two files under its `runners/` directory are unchanged source copies: the normal runner supports either expression PCA or annotation features, and the spatial-only runner replaces features with a constant vector.

The CLI uses an explicit `full-data-sns` preset: all cells, SN-S, no second stage, no reference downsampling, no initial coordinate override, cosine dissimilarity, 300 iterations and batch size 800. `--seed`, `--max-iter`, `--batch-size`, `--device` and `--z-display-spacing` are configurable. The seed default is 20260817; it is a reproducibility default, not a requirement for a new study. The bootstrap seeds Python, NumPy and Torch because the source runner's `RNG_SEED` alone controls its loader sampling.

The source implementation retains its scientific kwargs, including `beta=0.05`, `K=50`, `lambdaVF=1`, `partial_robust_level=1`, `sigma2_init_scale=2`, disabled expression nearest-neighbor initialization, and chunked computation. Spatial-only preserves `sigma2_end=None`. This migration does not alter these equations or parameter implementations.

The raw runner's own defaults differ from the public preset and remain unchanged:

| Raw environment option | Original default | Public CLI preset |
| --- | --- | --- |
| `USE_MORPHO_ALIGN_REF` | `1` | `0` |
| `USE_SPATEO_INTERNAL_DOWNSAMPLING` | `1` | `0` |
| `RNG_SEED` | `20260612` | Configurable, default `20260817` |
| `MAX_CELLS_PER_SLICE` | `0` | `0` |
| `STAGE1_MODE` / `STAGE2_MODE` | `SN-S` / `none` | Same |
| `BATCH_SIZE` / `STAGE1_MAX_ITER` | `800` / `300` | Configurable, same defaults |
| `Z_DISPLAY_SPACING` | `40` | Configurable, default `40` |

An expert can invoke the preserved raw runners with an explicitly reviewed environment to use their original advanced options. Their bare imports require `DATASET_ROOT` and `OUTDIR`; they do not implement `--help`. Use the public CLI for routine execution and immutable output paths.

The public output includes `invocation.json`, `input_cell_id_map.csv.gz`, `run.log`, `native/`, `aligned_coordinates.csv.gz` with original IDs, and `frozen_manifest.json`. A failed run remains as an incomplete directory; inspect the log and use a new output path for a retry. No evaluator or reference data is invoked by the CLI.

Native `z` is SL-relative and native `z_display` is viewer spacing. Neither replaces measured physical z. The public final coordinate table intentionally contains only original cell IDs, slice IDs and aligned XY.
