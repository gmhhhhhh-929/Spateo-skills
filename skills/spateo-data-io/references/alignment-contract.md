# Handoff to slice QC and alignment

For each biological specimen, provide observation-by-gene AnnData, stable cell and feature identifiers, raw-count semantics and their layer, original finite XY or XYZ, physical units/frame, slice/FOV identity and explicitly ordered slices. Keep specimens and timepoints independent during serial-slice QC/2D alignment. Do not fabricate physical z from numeric filenames.

Stage 3 runs `spatial-slice-quality-qc`, including `subskills/spatial-slice-quality-viewer`. Pass its traceable keep/exclude review to stage 4 (`spateo-2d-alignment`). Retaining a slice is not biological certification.

IO does not create PCA or alignment transforms. Downstream expression-PCA mode requires one verified shared basis fit from the appropriate specimen's measured expression, with matching cell IDs, feature order and provenance manifest. Matching dimensionality alone is insufficient. Annotation one-hot requires observed labels and a shared category order; it is not counts or expression PCA. Spatial-only constant vectors belong to an explicitly chosen downstream mode.

The read-only validator can audit numeric representation constraints but cannot prove PCA origin or registration accuracy. It never repairs an input. Generic IO preserves XYZ; `--require-2d` rejects three columns when an explicit XY handoff is requested. Derive XY under a new key with recorded meaning rather than dropping z silently.

The reserved 3D stage will eventually consume aligned serial sections for reconstruction, backbone analysis and interpolation. Already reconstructed, validated 3D inputs may enter the implemented 4D pipeline directly; mark the unimplemented 3D workflow as external rather than completed by these skills.
