# Registration review and version comparisons

Use this procedure when a tissue is visibly displaced, when a symmetric specimen may have inverted anatomical axes, or when a user asks whether a Spateo update changed alignment.

1. Check input and output cell IDs and original coordinates. Recover the full rigid transform from original/aligned coordinates if needed; inspect determinant, singular values and residual. A determinant near +1 only excludes reflection. A proper 180-degree rotation can still reverse dorsoventral orientation and be anatomically wrong.
2. Report whole-specimen and per-annotation bidirectional distances separately. Source-to-target and target-to-source tails can differ with tissue growth or incomplete coverage. `registration_qc.py` also reports nearest-label agreement; this is density- and annotation-dependent, not a universal acceptance threshold.
3. Review dorsal/ventral and anterior/posterior evidence in the original 3D reconstruction. For planarian labels `dorsal epidermal`/`ventral epidermal`, opposite mean-Z ordering is an explicit clue. Z is not universally an anatomical axis and global label centroids may be spatially biased; inspect lateral views and other tissues before concluding.
4. Keep biological size differences visible. A nonrigid overlap can remove the deformation being studied. A reference notebook may plot a nonrigid coordinate key but map cells with a rigid key. Record both rather than assuming those images show the same scientific frame.

## Controlled runtime comparison

`../scripts/compare_registration.py` runs one backend per process and writes immutable coordinate/QC/provenance outputs. Run on the remote compute host with the verified H5AD inputs; no raw matrices need transfer. A pinned source checkout must import from the requested root, and its environment must satisfy its dependencies. Do not replace the user's working native environment to run legacy Dynamo code.

First compare old/native with **fixed identical reference cell IDs**, common normalization arithmetic, common genes, device/dtype/seed, iterations and initialization. Then compare `--reference-policy library` to isolate automatic sampling. The old/native source audit must include actual imported files, not just version strings. Failed/missing runs are not evidence of numerical agreement.

```bash
python compare_registration.py --config pair.json --source-root /pinned/old-checkout --output /runs/old-fixed
python compare_registration.py --config pair.json --source-root /pinned/native-checkout --output /runs/native-fixed
```

Additional same-backend candidates can use `--annotation-aware` and `--initialization identity` to test native Spateo's expression-plus-label representation and initialization sensitivity. They are diagnostic candidates, not automatically accepted corrections. The helper exports both rigid and nonrigid coordinates. Compare by exact cell IDs, annotate the reference sample and run's input hashes, and review all relevant anatomical axes. `allow_flip=False` does not prohibit a 180-degree proper rotation.

The comparison helper has to be executed in the actual backend environments before its output can support a version claim. Source equality is a code audit only. If remote access is unavailable, finish coordinate diagnostics and viewer work, preserve a ready comparison command, and clearly report that runtime comparison and affected downstream reruns remain pending.

Once a corrected frame is accepted, create a new immutable native analysis run. Recompute mapping, field, trajectories and features because all depend on registration. A display-only repair cannot fix invalid scientific coordinates.

## What `nn_init` actually does

In this pinned Spateo implementation, `_coarse_rigid_alignment` aggregates coordinates and expression into voxels, searches nearest candidates in the chosen **expression/embedding representation** (KL distance for a matrix layer; Euclidean for embeddings), and estimates a robust coarse rotation and translation. This operation is not PCA; PCA is only relevant if a user explicitly supplies a PCA representation. The accepted inlier pairs also guide later rigid updates when `nn_init=True`, so switching it off removes both the coarse transform and this continuing guidance. Coordinate normalization and the main variational registration still run. `allow_flip=False` excludes the optional reflection branch, but cannot exclude a proper 180° rotation that reverses dorsal/ventral anatomy.

For already roughly co-oriented reconstructed specimens, compare with `nn_init=False`. For arbitrary input poses, blindly disabling initialization can instead cause another local optimum. Keep the general default and select a dataset-specific setting using anatomy, labels and before/after evidence; do not assert that the reference workflow always produces a correct pose.
