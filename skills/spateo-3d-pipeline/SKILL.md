---
name: spateo-3d-pipeline
description: Identify the reserved 3D reconstruction, backbone analysis and gene-interpolation stage in the Spateo workflow; implementation is pending and this entrypoint does not execute analysis.
---

# Spateo 3D pipeline — reserved

Stage 5: consumes reviewed 2D-aligned serial sections and will provide 3D model reconstruction → backbone analysis and spatial gene interpolation. Its outputs will feed the two-timepoint 4D pipeline.

Implementation is intentionally pending. Do not invent an executable workflow, dependency list, output contract or passing test result for this stage. Preserve this position in workflow plans and state that the implementation is unavailable. Already reconstructed and validated 3D H5AD inputs can go directly to the implemented 4D stage, with their external reconstruction provenance recorded.

The user's current request reserves this stage even though protocol reference notebooks contain potential future material. Implement it only in a separate authorized task with agreed code and validation scope.
