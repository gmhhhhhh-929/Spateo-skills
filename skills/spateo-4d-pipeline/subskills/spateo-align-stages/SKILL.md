---
name: spateo-align-stages
description: Align two 3D timepoints using native Spateo morpho_align_ref with explicit count and coordinate contracts.
---

# spateo-align-stages

Read the parent [config contract](../../references/config-contract.md). Run the shared engine from the parent skill directory:

```bash
python scripts/run_4d_pipeline.py --config /project/config.json --project /project/analysis --stop-after alignment
```

This produces a partial run with two aligned H5ADs and QC JSON. Preserve original XYZ and counts. Configure the same units and choose SN-S (rigid) or SN-N (nonrigid) deliberately. Alignment does not require an annotation; subset annotations are checked later. Review shared genes, sample sizes and actual aligned bounds. The source returns models, reference models, transport matrices and reference transport matrices. The final tuple member is not sigma2.

Continue through a child run using the partial run's `--parent-manifest`; validated alignment checkpoints are reused. The compatibility filename `scripts/alignment_10dpa_14dpa.py` now accepts the v3 runner arguments, not the old notebook flags.
