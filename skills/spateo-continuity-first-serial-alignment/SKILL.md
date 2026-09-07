---
name: spateo-continuity-first-serial-alignment
description: Run blind serial Spateo alignment with gated continuity repairs using the packaged continuity-guided pipeline. Use for sparse terminal slices, displaced terminal blocks, or internal interfaces, with annotation one-hot or verified shared expression PCA.
---

# Spateo Continuity-First Serial Alignment

This is the compatibility entrypoint for the source skill of the same name.
Install it beside `spateo-2d-alignment`, or use the complete repository checkout.
Read [the canonical skill](../spateo-2d-alignment/SKILL.md) and its
[continuity contract](../spateo-2d-alignment/references/continuity-guided.md).
This entrypoint uses the already published pre-zebrafish implementation;
it does not ship another pipeline or change its defaults.

Run from this skill directory with sanitized inputs and a new output directory:

```bash
python ../spateo-2d-alignment/pipelines/continuity-guided/run.py \
  --stage specimen --slice-dir /path/to/blind-slices \
  --output-dir /path/to/new-run \
  --representation annotation-onehot --annotation-key anno \
  --profile legacy --postprocess auto --device 0
```

For expression-only data, prepare one shared PCA basis across the specimen's
slices using the canonical skill's preparer, then select
`--representation expression-pca --annotation-qc off`. Do not reuse independently
fitted per-slice PCA or silently reinterpret annotation one-hot as expression.

Preserve source data and stable cell IDs. Keep reference coordinates out of
inference. Freeze coordinates, inputs, parameters and QC before reference-based
evaluation. Review accepted and rejected repairs and inspect full-point
[before/after](../spatial-before-after-viewer/SKILL.md) and
[3D pointcloud](../spatial-pointcloud-viewer/SKILL.md) results.
Use the user's explicitly selected evaluation protocol; this alias does not
add an evaluator or a new accuracy convention.

Legacy remains the compatibility default. The existing generalized profile is
opt-in and has mixed validation results; it is not evidence of universal gains.
No zebrafish image priors, initialization experiments or later scoring changes
are included.
