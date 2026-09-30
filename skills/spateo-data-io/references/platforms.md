# Platform-specific decisions

Read only the relevant row, then inspect its source signature in source-api.md.

For SeekSpace, BMKMANU S1000, Salus STS and Singleron space, read the exact native-file contracts and primary evidence in [domestic-platforms.md](domestic-platforms.md). Shared matrix storage alone does not identify a platform.

| Platform | Core layout / key decisions |
| --- | --- |
| Visium | Feature H5 or MTX plus `spatial/tissue_positions*.csv`/Parquet; XY is full-resolution pixel column/row, joined by barcode. Keep tissue flags and scalefactors. Simultaneous alternate matrix encodings need contract resolution, not an assumption of equivalence. |
| Visium HD | Named square-bin resolutions and segmented-cell representations are separate entries. Cell segmentation uses supported GeoJSON geometry and matching cell IDs. Do not combine different resolutions as independent cells. |
| Xenium / Atera | Cell-feature H5 plus cells metadata; vendor identity distinguishes Atera. Atera is a preview layout. Preserve platform coordinate units; H&E/vendor affines are metadata, not estimated registration. |
| MERFISH | Matched cell-by-gene and cell-metadata tables per region. Join explicit cell IDs. Keep native string IDs (including leading zeros), raw int64 values and supplied Z; keep regions separate unless their frames are documented. |
| seqFISH | Matched counts and cell-coordinate groups; supported label identifiers are IDs. Preserve section identity and native XY/XYZ. Matrix orientation requires matching IDs; no positional fallback or invalid-to-zero conversion. |
| NanoString / CosMx | Expression and metadata tables with FOV/cell identity; local pixel coordinates across FOVs do not imply one global frame. |
| Slide-seq | DGE and bead-location tables; preserve barcode identity and exact raw counts during streamed loading; both routes use the native streaming parser. |
| STARmap PLUS | Expression and spatial metadata tables; a processed expression filename does not establish raw-count semantics. Valid TYPE declaration rows are schema metadata, not observations. The direct default keeps source numeric precision; an explicit lossy cast is rejected. |
| Stereo-seq GEM | Native bin/cell resolution is preserved unless explicit `stereoseq_bin_size` is supplied. Counts are validated without legacy uint16 truncation; inspect headers, offset/resolution and observation representation. |
| Stereo-seq GEF | Native bin and cell schemas; default uses smallest stored bin resolution, an explicit size selects an existing resolution. CellBin remains cells. Schema versions do not determine chemistry. |
| Open-ST | Existing H5AD can use the H5AD reader; no dedicated native-output automatic reader is validated here. |
| Seq-Scope | Use explicit `read_seqscope`; not part of automatic technology discovery. Supply matrix and position files, choose barcode or bin mode explicitly. |

Stereo chemistry accepts a user-declared V1/V2 label or unknown; never infer it from GEF schema version. Binning changes the observation unit and must be recorded. `read_stereoseq` currently calls automatic reading and is not an independent parser. That native GEM/GEF route and legacy `read_bgi`/`read_bgi_agg` have different contracts. Legacy read_bgi chooses one of binsize, bin columns, segmentation labels or label-column modes; review its signature and counts range. It does not perform segmentation. Do not turn a mask or AGG image into a gene matrix.

The automatic core emphasizes complete expression/ID/coordinate validation. Direct readers expose richer platform assets and custom filenames with their own memory options (`load_image` versus `load_images`, boundary/label/transcript flags). Inspect actual available arguments; the unified reader does not forward arbitrary kwargs. Missing/duplicate IDs or nonfinite coordinates require a recorded correction upstream, not random coordinates or positional joins.

## Comparing automatic and direct results

Check each route against native source values and IDs, then compare equivalent observations/features. Agreement between two routes that share a parser is not an independent accuracy estimate.

| Convention | How to compare or act |
| --- | --- |
| 10x features | Direct H5/MTX and platform readers default to Gene Expression-only features; automatic reading retains all source features. Align verified feature IDs before comparing the corresponding matrix. Duplicate gene symbols are retained by platform readers; the generic MTX API separately defaults to `make_unique=True`. |
| Required metadata | Direct readers now reject missing/duplicate IDs and invalid coordinates instead of silently dropping matrix cells, keeping the first duplicate or joining by row position. Restore matching source companions; do not weaken validation to reproduce the former partial output. |
| Raw values | Native raw-count readers reject negative, fractional or nonfinite values and unsafe aggregation/casts. Keep int64 precision for large counts. STARmap processed expression is a separate supported representation that permits finite signed/fractional values. |
| Stereo/Seq-Scope coordinates | A bin's origin and center are different coordinates. Native Stereo-seq bin output and legacy `read_bgi(add_props=True)` geometry use different conventions; compare at the same bin size by bin identity and the documented origin/center transform. Seq-Scope is direct-only and its binned output uses bin centers; do not count it as an automatic detection success. |
| CosMx frames | `obsm['spatial']` contains FOV-local pixels; available source global columns live in `obsm['spatial_fov']`. Direct `fov_file` metadata is retained, not implicitly added as a second offset. Keep FOV identity in cell IDs. |
| STARmap orientation | Direct `reorient_xy=False` preserves source XY(Z). With `True`, the direct reader applies `x'=max(y)-y`, `y'=max(x)-x`, leaving Z unchanged. Automatic reading preserves the source frame. |
| Images | Automatic loading has a cumulative 32 MiB decoded-image budget per entry, additionally bounded by remaining memory; raising `max_memory_bytes` does not raise that image ceiling. Large or multiframe images can remain by path with asset status while core data is ready. Legacy direct readers have platform-specific image options and may eagerly load a larger image. Compare available source pixels and asset statuses, not identical `uns` keys or image counts. |

An invalid present HD H5 now raises its actual core error; the direct HD bin reader falls back to MEX only when H5 is absent. Optional missing/corrupt images do not discard an otherwise valid core, and a missing lowres no longer hides a valid hires. Legacy direct readers do not acquire automatic lazy/resume or memory-budget parameters from these accuracy fixes.
