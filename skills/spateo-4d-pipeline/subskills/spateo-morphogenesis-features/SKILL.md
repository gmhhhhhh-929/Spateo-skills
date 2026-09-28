---
name: spateo-morphogenesis-features
description: Compute native Spateo morphogenesis velocity, acceleration, curvature, curl, divergence, torsion and Jacobian features, test associated genes with glm_degs, and export complete test tables and actual fitted relationship curves.
---

# 3 · Features and associated genes

Implementation: [features_stage.py](../../scripts/features_stage.py). Read the [native feature caveats](../../references/protocol-migration.md) before interpretation.

Enter from the field checkpoint of subskill 2 or set `workflow.entry="field"` with native field/source and target H5ADs plus the declared frame. A trajectory alone without the underlying field is insufficient for differential geometry. Feature computation depends on the fitted field, so changing only the trajectory horizon does not invalidate features.

Supported scalar selections: `speed`, `velocity_x/y/z`, `acceleration`, `curl`, `divergence`, `curvature`, `torsion`, `jacobian_frobenius`, and each `jacobian_xx` … `jacobian_zz`. Velocity is saved as a 3-vector; speed is its magnitude. Jacobian is saved as the complete `3 × 3 × cells` tensor and cell-ID NPZ whenever a Jacobian scalar is requested. Explicit components keep their signs; curl/torsion scalars use the native magnitude convention. Never give a tensor directly to a scalar spline formula.

For each requested `metrics.glm_metrics`:

1. Use `layers['normalized']` (size-factor normalized counts), not log1p X. `glm_genes` is an explicit gene list or `"*"` for all genes passing `glm_min_cells`; an empty list disables gene testing.
2. Call `st.tl.glm_degs` with `fullModelFormulaStr=f"~cr({metric}, df=3)"` and reduced `~1`. Preserve the complete tested-gene universe, successful and failed statuses, p/q-values, likelihood and a separate selection flag. BH correction applies within each feature's tested genes; comparisons across features require a declared additional correction.
3. Save full and selected CSVs, Spearman rho as a descriptive association measure, test-count summaries and the actual fitted `mu`/confidence intervals from Spateo. A spline association need not be monotonic or have large Spearman rho; neither test demonstrates causality.
4. Save top fitted curves and PNGs. Explicitly label nonsignificant top-ranked genes if there are no hits. Do not silently drop fit failures or copy notebook organism-specific gene IDs into another dataset.

```bash
python scripts/run_features.py --config /project/field.json --project /project/analysis
```

For all-gene screens, review memory/runtime: Spateo retains per-cell fit tables. `glm_top_plots` limits retained display curves, not the tested gene universe. Native NB2 currently uses dispersion alpha=1. Spatial autocorrelation and one specimen per timepoint limit biological inference.

Optional expression GP remains a separate disabled-by-default branch in the shared runner; it is spatial interpolation, not inferred continuous developmental time. GO enrichment needs a verified organism gene mapping and database and is not a prerequisite or an automatically completed output.

Interpret FDR selection separately from association strength. Report full-cell Spearman rho and the tested gene universe; a significant nonlinear spline may have weak monotonic correlation. Feature-colored paths inherit the seed-cell value unless features are explicitly evaluated along the path. For plotting and interpretation, follow the [viewer contract](../../references/viewer-contract.md).
