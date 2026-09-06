# Continuity-guided mode

`pipelines/continuity-guided/run.py` validates the exact shared annotation one-hot contract, then invokes `engine.py`. The engine preserves the selected source Spateo-first pipeline's scientific functions, parameters, gates and fallbacks. `_core.py` and `_continuity_scoring.py` contain only its reachable I/O, geometry, plotting and scoring dependencies; they are not additional alignment entrypoints.

The native defaults are 300 Spateo iterations, batch size 800, seed 20260817, SN-S, cosine dissimilarity and `beta=0.05`. All cells participate in the baseline and exported coordinates. `--sample-cap 1800` limits deterministic annotation-aware scoring/plotting samples; it does not downsample the final cell table. CLI arguments preserve their original configurable defaults.

`--postprocess auto --repair-policy auto` applies only repairs whose blind gates pass: internal interface translation, short terminal block repair, priority-annotation continuity correction and sparse terminal singleton rescue. `--postprocess none` returns the original full-data Spateo baseline without those repairs. `--priority-annotation` is configurable; the inherited default `ROD` is used only where supported, with the source's deterministic fallback for other annotations. Do not force that label into a different specimen. The inherited `learned_probe` option is retained for source reproducibility and is not the general workflow recommendation.

The source prohibits reflections and rejects `spatial_3d` inputs. The packaged preflight additionally rejects any extra spatial/reference `obsm` before the engine opens H5AD objects. Keep reference data outside the supplied input directory.

Outputs include `aligned_coordinates.csv.gz`, `slice_transforms.csv`, `spateo_pair_summary.csv`, `interface_operations.csv`, `blind_qc_summary.json`, two blind QC plots and `frozen_manifest.json`. Inspect operations and QC rather than assuming every proposed repair is accepted. Failed optional repair proposals retain their previous coordinates according to the source implementation.

The frozen manifest hashes the engine and outputs. Module-level source and migrated hashes are separately recorded in the skill's provenance. Run `scripts/verify_bundle.py` before changing code or trusting a copied package; a changed source helper requires new provenance and validation.
