# Handoff to slice QC and alignment

For each biological specimen, provide observation-by-gene AnnData, stable cell and feature identifiers, raw-count semantics and their layer, original finite XY or XYZ, physical units/frame, slice/FOV identity and explicitly ordered slices. Keep specimens and timepoints independent during serial-slice QC/2D alignment. Do not fabricate physical z from numeric filenames.

Stage 3 runs `spatial-slice-quality-qc`, including `subskills/spatial-slice-quality-viewer`. Pass its traceable keep/exclude review to stage 4 (`spateo-2d-alignment`). Retaining a slice is not biological certification.

IO does not create PCA or alignment transforms. Downstream expression-PCA mode requires one verified shared basis fit from the appropriate specimen's measured expression, with matching cell IDs, feature order and provenance manifest. Matching dimensionality alone is insufficient. Annotation one-hot requires observed labels and a shared category order; it is not counts or expression PCA. Spatial-only constant vectors belong to an explicitly chosen downstream mode.

The read-only validator can audit numeric representation constraints but cannot prove PCA origin or registration accuracy. It never repairs an input. Generic IO preserves XYZ; `--require-2d` rejects three columns when an explicit XY handoff is requested. Derive XY under a new key with recorded meaning rather than dropping z silently.

The 3D stage can now construct and review a point-cloud VTK after the aligned sections have been combined into an H5AD with verified XYZ coordinates. It does not yet implement surface, voxel, reconstructed-cell, backbone, or spatial-interpolation models. Already reconstructed, validated 3D H5AD inputs may enter the 4D pipeline directly; do not mark pending 3D phases as completed.
