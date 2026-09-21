---
name: spateo-4d-pipeline
description: Run tracked two-timepoint 3D alignment and native Spateo morphogenesis, including cell mapping, SparseVFC, trajectories, geometric metrics, optional GLMs/GP and an interactive dashboard. Use for complete 4D analyses or child runs with verified checkpoint reuse.
---

# Spateo 4D pipeline

Stage 6 of the collection: environment → IO → slice quality with viewer → 2D alignment → reserved 3D reconstruction/backbone/interpolation → **4D**. Start here with two already reconstructed 3D AnnData objects. The reserved 3D skill is not a prerequisite implementation to invent.

Use current native Spateo at commit `615644f88613bea8ceb2e2df1e2391d16de55ec1`. Dynamo is not required. The supported scientific sequence comes from the user's [protocol notebooks](https://github.com/gmhhhhhh-929/Spateo-protocol-files/tree/b11ae99fbdc4ae46d41880e9306ab7e5c2751ac5/code/04_alignment_and_morphogenesis). Read [protocol migration](references/protocol-migration.md) before reproducing their settings or plots.

## Inputs and decisions

Require finite `(n_cells, 3)` coordinates with actual 3D support, matching units, unique IDs, and a verified raw-count layer. The runner does not silently treat normalized X as counts: `alignment.x_is_counts=true` is an explicit assertion, followed by numeric checks. Preserve original coordinates under their original key. Cell annotations are needed only when selecting a biological subset; a missing annotation is never inferred.

Confirm source→target temporal order, count layer, coordinate unit, subset, rigid/nonrigid alignment mode, mapping budget and requested GLM/GP genes. Mapping can allocate dense pairwise costs; `mapping.max_pairs` bounds the reviewed cell-pair count but is not a RAM guarantee. Use the provided CPU defaults as starting settings, not validated biological parameters.

## Configure and execute

Copy [assets/config.template.json](assets/config.template.json), then read [config-contract.md](references/config-contract.md) for supported fields and migration from old configs. Paths resolve relative to that JSON and are saved absolute. Unknown settings fail instead of being ignored.

```bash
python scripts/run_4d_pipeline.py --config /project/config.json --project /project/analysis --dry-run
python scripts/run_4d_pipeline.py --config /project/config.json --project /project/analysis
```

The dry run imports the actual environment and hashes inputs/implementation without running scientific stages or creating a run. Check the resolved settings and cost before compute. Each execution creates a new run and records separate H5AD checkpoints, CSV/NPZ outputs, logs, hashes and a manifest. Failures mark remaining stages blocked; optional disabled stages are skipped. Completed parents are never modified.

## Ordered workflow and subskills

1. [spateo-align-stages](subskills/spateo-align-stages/SKILL.md): native normalization/log layers → `morpho_align_ref` → aligned H5AD pair and coordinate QC.
2. [spateo-morphogenesis](subskills/spateo-morphogenesis/SKILL.md): select subset/shared expressed genes → normalization → `cell_directions` → SparseVFC → trajectories and geometric metrics; optional GLMs and spatial GP expression interpolation.
3. [spateo-render-dashboard](subskills/spateo-render-dashboard/SKILL.md): show source/target, mapped endpoints, displacement/metrics, vectors and available GLM tables. Display sampling does not change scientific outputs.

[spateo-manage-runs](subskills/spateo-manage-runs/SKILL.md) owns checkpoint verification and stage state; [spateo-refine-analysis](subskills/spateo-refine-analysis/SKILL.md) handles child configs and dependency-aware reruns. Install this entire directory including subskills and scripts.

## Review and handoff

Compare alignment in the exact coordinate key used for mapping. Check subset counts, shared genes, mapped endpoint/vector consistency, finite fields/metrics and trajectory bounds. Mapping is an inferred correspondence, not lineage tracing. Integration time is model time; `t_end=1` does not claim one day. GP interpolates expression over space, not continuous developmental time.

Return the manifest, config, run/parent IDs, stage states, checkpoint/table paths and dashboard. The trajectory, metric and GP branches have separate artifacts, all indexed by the manifest. State test/biological limits from [validation.md](references/validation.md). Do not report skipped GLM/GP or the reserved 3D stage as completed analysis.
