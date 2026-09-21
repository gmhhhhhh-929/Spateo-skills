# Platform-specific decisions

Read only the relevant row, then inspect its source signature in source-api.md.

| Platform | Core layout / key decisions |
| --- | --- |
| Visium | Feature H5 or MTX plus `spatial/tissue_positions*.csv`; XY is full-resolution pixel column/row, joined by barcode. Keep tissue flags and scalefactors. Simultaneous alternate matrix encodings need contract resolution, not an assumption of equivalence. |
| Visium HD | Named square-bin resolutions and segmented-cell representations are separate entries. Cell segmentation uses supported GeoJSON geometry and matching cell IDs. Do not combine different resolutions as independent cells. |
| Xenium / Atera | Cell-feature H5 plus cells metadata; vendor identity distinguishes Atera. Atera is a preview layout. Preserve platform coordinate units; H&E/vendor affines are metadata, not estimated registration. |
| MERFISH | Matched cell-by-gene and cell-metadata tables per region. Join explicit cell IDs. Keep regions separate unless their frames are documented. |
| seqFISH | Matched counts and cell-coordinate groups; supported label identifiers are IDs. Preserve section identity and native XY/XYZ. |
| NanoString / CosMx | Expression and metadata tables with FOV/cell identity; local pixel coordinates across FOVs do not imply one global frame. |
| Slide-seq | DGE and bead-location tables; preserve barcode identity and counts during streamed loading. |
| STARmap PLUS | Expression and spatial metadata tables; a processed expression filename does not establish raw-count semantics. Valid TYPE declaration rows are schema metadata, not observations. |
| Stereo-seq GEM | Native bin/cell resolution is preserved unless explicit `stereoseq_bin_size` is supplied. Counts are validated without legacy uint16 truncation; inspect headers, offset/resolution and observation representation. |
| Stereo-seq GEF | Native bin and cell schemas; default uses smallest stored bin resolution, an explicit size selects an existing resolution. CellBin remains cells. Schema versions do not determine chemistry. |
| Seq-Scope | Use explicit `read_seqscope`; not part of automatic technology discovery. Supply matrix and position files, choose barcode or bin mode explicitly. |

Stereo chemistry accepts a user-declared V1/V2 label or unknown; never infer it from GEF schema version. Binning changes the observation unit and must be recorded. The newer `read_stereoseq`/automatic core and legacy `read_bgi`/`read_bgi_agg` have different contracts. Legacy read_bgi chooses one of binsize, bin columns, segmentation labels or label-column modes; review its signature and counts range. It does not perform segmentation. Do not turn a mask or AGG image into a gene matrix.

The automatic core emphasizes complete expression/ID/coordinate validation. Direct readers expose richer platform assets and custom filenames with their own memory options (`load_image` versus `load_images`, boundary/label/transcript flags). Inspect actual available arguments; the unified reader does not forward arbitrary kwargs. Missing/duplicate IDs or nonfinite coordinates require a recorded correction upstream, not random coordinates or positional joins.
