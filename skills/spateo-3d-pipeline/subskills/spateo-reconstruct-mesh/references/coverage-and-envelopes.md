# Coverage-aware surface review and smooth display envelopes

Use this route when a tissue surface misses visible cells, sparsely sampled ends
disappear, cleanup removes plausible thin structures, or a convex body looks
faceted. It supplements the first-pass builder; it is not a universal replacement
for the user's preferred reconstruction or an anatomical segmentation algorithm.

## Diagnose before changing geometry

Compare the exact annotation selection and coordinate frame used by mesh and
points. Display all selected original cells for coverage review, not a sampled or
already-filtered subset. Marker radius can visually protrude even when its center
is inside; count centers geometrically. Separately inspect opaque surfaces to
distinguish transparency/depth-sorting artifacts from real spikes or holes.

Measure coverage on the final saved surface, after smoothing/component filtering,
not only on its pre-extraction voxel mask. Record global containment plus spatial
coverage in equal-width long-axis bins and transverse/3D blocks, with cell counts.
Use PCA only to define a reproducible audit axis, not an anatomical axis or an
alignment transformation. Tiny bins are noisy: report their counts and define
minimum support explicitly. The supplied audit uses ten axial bins, 6×6×3 spatial
blocks, and a 30-cell reporting threshold as configurable review conventions.

A 95% overall score can coexist with almost zero coverage at sparse endpoints.
Conversely, six missed cells in one sparse bin should not trigger inflation of an
entire organ merely to satisfy a percentage. Keep absolute counts, point IDs,
regional scores and volume changes together. Compare within the same coordinate
frame and the same original selection; do not compare filtered-cell and raw-cell
denominators as though they were identical.

## Distinguish the constraints

- **Own-cell coverage:** a tissue mesh covers its selected source-cell centers.
- **Containment:** tissue triangles remain inside the display body.
- **Union occupancy:** tissue union volume divided by body volume, on one shared
  physical grid. Count overlaps once; report overlap and uncovered fractions too.

Independent lineage meshes need not partition the body. Labels can intermingle,
samples can be sparse, real cavities can exist, and an intentionally expanded
display shell adds empty margin. Never sum tissue volumes to claim body coverage,
inflate every tissue to fill the shell, or present a space-filling nearest-label
segmentation as measured anatomy. A mutually exclusive partition is a different
model requiring an explicit user choice and its own validation.

## Repair the relevant cause

1. For systematic boundary loss, inspect the density isovalue, final smoothing
   shrinkage and physical bandwidth. Uniform scaling is not a general solution.
2. For density heterogeneity, compare a smooth local threshold field with the
   global threshold. `regional_candidate` calibrates overlapping axial windows;
   independent spatial-block QC catches its transverse blind spots.
3. `repair_supported_gaps` optionally lowers thresholds only in bounded patches
   around missed cells that have same-label neighbors. Its support radius,
   neighbor count, threshold fraction and antialiasing scale are explicit choices,
   not biological constants. Unsupported cells remain in the original point cloud
   and coverage report; they are not silently relabeled as noise.
4. Preserve supported disconnected regions. A fixed percentage of the largest
   component can erase genuine small anatomy in a heterogeneous tissue. Record
   component cell support and protect previously accepted substantial components.
5. For attached bulbs or spikes, compare thickness-aware opening and mild surface
   smoothing. Audit severed necks, component splits, lost-cell IDs and volume loss.
   Stronger opening can cut valid bridges; aggressive Taubin filtering can create
   pointed overshoot. Do not assume more smoothing always improves a mesh.

Select the least expansive candidate meeting reviewed coverage targets. Topology,
source preservation and requested body containment are hard gates. Coverage
targets are review goals: if none meet every goal, a bounded candidate improving
global, axial and block summaries can be returned **with explicit warnings**;
otherwise retain the prior mesh. Report unmet targets and rejected candidates.
Never silently weaken a threshold or report a fallback as passing. Parameters
learned from one lineage/stage should be retested, not hard-coded as universal.

## Smooth body construction

For a display envelope, combine point-supported body occupancy with the retained
tissue union on a shared isotropic grid. Smooth its signed-distance volume field,
extract an isosurface, add the smallest reviewed positive clearance that passes
containment, and recompute morphology. Filling cavities here is a display choice
for the **body**, not permission to fill tissue cavities. Do not force planarians
into ellipsoids or treat display-envelope expansion as biological growth.

Convex hull halfspaces are valid only for convex shells. For a non-convex smooth
shell, checking vertices alone is insufficient: triangles can bridge across a
concavity. The supplied `triangle_containment` checks signed distance at triangle
centroids against enclosing radii, recursively subdividing unresolved triangles.
This uses the distance function's 1-Lipschitz bound. It requires a correctly
oriented, closed, non-self-intersecting shell and is a numerical certificate,
not exact arithmetic. Outside/contacting or unresolved triangles fail. Recheck
after serialization, coordinate rounding and display decimation; do not assume a
full-resolution certificate automatically covers a browser LOD.

## Reusable runner

`scripts/review_meshes.py --config /absolute/path/review.json` uses
`scripts/coverage_envelope.py`. It requires NumPy, SciPy, scikit-image, PyVista
with `voxelize_binary_mask(reference_volume=...)` support (tested with 0.46.5),
Matplotlib and the validated Spateo runtime.

Minimal configuration (paths must be real, output directory must not exist):

```json
{
  "point_cloud": "/results/sample.pc.vtk",
  "label_key": "lineage",
  "previous_body": "/results/previous/body.vtk",
  "output_dir": "/results/coverage_review_v2",
  "long_axis": 220,
  "global_target": 0.985,
  "regional_target": 0.9,
  "block_target": 0.8,
  "max_volume_ratio": 2.0,
  "local_repair_passes": 1,
  "body_sdf_sigma": 3.0,
  "body_clearance_voxels": 1.0,
  "tissues": [
    {"name": "neural", "label": "Neural", "values": ["neural"],
     "vtk": "/results/previous/neural.vtk", "color": "#D98CB3", "refine": true}
  ]
}
```

List every tissue that must be enclosed, including unchanged ones with
`"refine": false`. Coordinate units are inherited, never inferred. The runner
saves candidates separately, selected labeled VTKs, raw-cell coverage and IDs,
source hashes, Spateo readback checks, recomputed morphology, union/overlap voxel
estimates and four-view previews. It refuses overwrite. A failed run may leave
diagnostic files; never consume a run without its completed manifest and previews.
Large grids fail the memory bound rather than allocating without limit.

Inspect the per-tissue views and point overlays, not only the combined view.
Return unresolved low-support regions and large volume changes to the user for
review. Keep technical detail in audit/export artifacts when the presentation UI
should stay uncluttered. Hiding diagnostics must not erase provenance.
