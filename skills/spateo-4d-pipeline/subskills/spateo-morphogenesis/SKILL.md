---
name: spateo-morphogenesis
description: Run native cross-timepoint cell mapping, SparseVFC, trajectories, geometric metrics and optional GLM or gene GP from tracked 3D alignment outputs.
---

# spateo-morphogenesis

Read [protocol migration](../../references/protocol-migration.md) and the parent [config contract](../../references/config-contract.md). Use the parent runner with a valid alignment parent manifest:

```bash
python scripts/run_4d_pipeline.py --config /project/config.json --project /project/analysis \
  --parent-manifest /project/analysis/runs/aligned/manifest.json --stop-after gp
```

The command above is run from the main pipeline directory. The compatibility filename `scripts/morphogenesis_10dpa_14dpa_cns.py` also accepts the same v3 arguments. It does not trust arbitrary stale H5AD prerequisites.

Select the requested biological group in both stages; intersect expressed genes in source order and normalize with native `st.pp`. Mapping writes X_<key>, V_<key>, a transport matrix with ordered IDs and displacement CSV. SparseVFC learns the field at source coordinates; no mesh is mandatory. Native trajectories, metrics and GP each have their own checkpoint/output. GLMs use normalized counts with a correctly interpolated formula; genes and thresholds must be explicit. GP genes default empty and GP is off unless requested.

Inspect vectors, finite values, metric variability and extrapolation. Do not interpret a mapping as observed lineage, model time as elapsed biological time, or spatial expression GP as temporal interpolation. Optional static mesh plots are outside the default runner; missing plots do not justify altering scientific coordinates.
