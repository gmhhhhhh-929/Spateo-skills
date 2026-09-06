# Continuity-guided mode

`pipelines/continuity-guided/run.py` validates the selected feature contract, then invokes `engine.py`. Choose `--representation annotation-onehot`, `expression-pca`, or `spatial-only`. Expression mode verifies the joint-basis manifest, supports absent annotations and defaults to `--annotation-qc off`. `_core.py` and `_continuity_scoring.py` provide I/O, geometry, plotting and scoring helpers; they are not additional alignment entrypoints.

The defaults are 300 Spateo iterations, batch size 800, seed 20260817, SN-S, cosine dissimilarity and `beta=0.05`. All cells participate in the baseline and exported coordinates. `--sample-cap 1800` limits deterministic scoring/plotting samples; it does not downsample the final cell table. Annotation-aware sampling is used only when annotations are enabled.

`--profile generalized` turns off expression-neighbor initialization and replaces the named tissue preference with anonymous, supported-label proposals. A terminal proposal needs at least three supported shared labels; its rigid motion must be bounded, improve whole-tissue continuity, improve independent labels in aggregate, and avoid excessive deterioration in any checked label. Proposals that fail retain the previous coordinates. The same rules apply to every specimen. These gates reduce unsupported edits; blind continuity is not proof of anatomical correctness.

`--profile legacy` is the compatibility default. It retains the earlier `nn_init=True` and `ROD` priority default for reproduction. `--profile generalized` is opt-in: frozen validation found average Drosophila improvements but Planarian regressions and severe individual failures. It has not established consistently better accuracy across specimens. See [validation.md](validation.md). The legacy `learned_probe` bypass is not accepted by the generalized profile.

`--postprocess auto --repair-policy auto` runs gated internal interface translation, short terminal-block repair and sparse terminal-singleton rescue. `--postprocess none` returns the full-data Spateo initial alignment. With annotation QC off, the runtime uses an internal geometry-only marker and skips annotation-dependent proposals; this marker is not a tissue assignment. Supplied expression PCs still define Spateo correspondence. Enable `--annotation-qc provided --annotation-key KEY` only when the labels are suitable for this extra role; missing labels then fail validation rather than silently mixing labeled and unlabeled slices.

When annotation QC is off, legacy edge summaries may report `shared_annotations=1` or `semantic_labels=1` for the internal whole-tissue geometry marker. These counts are not biological annotation evidence. Use the explicit `blind_contract` and `annotation_qc` fields to interpret the scores.

The short terminal-block trigger and other inherited geometric gates remain limited heuristics. Large smooth drift cannot generally be separated from genuine anatomy using continuity alone. Preserve uncertain outputs and report the limitation rather than forcing a visually smooth correction.

The source prohibits reflections and rejects `spatial_3d` inputs. The packaged preflight additionally rejects any extra spatial/reference `obsm` before the engine opens H5AD objects. Keep reference data outside the supplied input directory.

Outputs include `aligned_coordinates.csv.gz`, `slice_transforms.csv`, `spateo_pair_summary.csv`, `interface_operations.csv`, `blind_qc_summary.json`, two blind QC plots and `frozen_manifest.json`. Inspect operations and QC rather than assuming every proposed repair is accepted. Failed optional repair proposals retain their previous coordinates according to the source implementation.

The frozen manifest records the selected feature mode, annotation-QC behavior, effective profile parameters, engine and output hashes. Source and updated module hashes are recorded in the skill's provenance. Run `scripts/verify_bundle.py` before trusting a copied package; changed helpers require updated provenance and validation. Original reference coordinates are used only by a separate evaluator after result freezing.
