# Validation and limits

The shared expression-PCA interface is implemented and functionally validated.
It does not establish that expression input has the same accuracy as annotation
one-hot input, or that either continuity profile works best for every specimen.

The compatibility default remains `--profile legacy`. The optional
`--profile generalized` removes the fixed tissue-name priority and disables
nearest-neighbor initialization. A single policy was frozen before the final
tests; no per-specimen selection from reference scores was used.

## Frozen annotation-input benchmark

All reported accuracy compares recovered coordinates with ID-matched original
3D coordinates. A whole specimen receives one proper XY similarity transform
(positive uniform scale, rotation and translation), with physical Z unchanged.
A2 and A5 count cells with XYZ error strictly below 2% and 5% of the original
whole-specimen XY bounding-box diagonal. There is no per-slice or per-group
reference fitting. Means weight biological specimens equally.

| Original perturbations | Legacy A2 / A5 (%) | Generalized A2 / A5 (%) |
| --- | ---: | ---: |
| Drosophila, 9 specimens | 54.605 / 92.978 | 57.999 / 94.532 |
| Planarian, 8 specimens | 91.074 / 99.947 | 89.169 / 99.839 |

These 17 specimens were previously evaluated, so this comparison is
retrospective, including its designated development and holdout partitions.

Independent perturbations tested 9 historical holdout specimens with two new
seeds each and 10 additional Drosophila stages with one seed each. Seeds were
averaged within a biological specimen before averaging specimens. The table
below uses only cases where both methods completed.

| Cohort | Paired biological specimens / cases | Change in A2 / A5 (percentage points) |
| --- | ---: | ---: |
| Additional Drosophila stages | 9 / 9 | +1.110 / +6.333 |
| Historical Drosophila, new perturbations | 5 / 10 | +0.626 / +7.950 |
| Historical Planarian, new perturbations | 4 / 8 | −3.246 / −0.229 |

Additional-stage completion was 9/10 for legacy and 10/10 for generalized.
The legacy failure was an SVD failure in coarse NN initialization on a
single-label input. The corresponding generalized output had A2/A5 of only
0.177/1.368%; completing a run did not mean accurate recovery. Failure accuracy
was NA, never zero. Means over different success sets are not paired effects.

The largest individual new-perturbation Drosophila A2 regression was 55.848
percentage points. These results do **not** support universally better
stability, improved Planarian accuracy, or choosing generalized automatically.
Two additional conservative refiners produced no accepted coordinate changes
in their 24 development runs and are not included in the selected profile.

## Expression and packaging checks

Joint PCA preparation and full 300-iteration, no-annotation continuity
alignment completed for one 14-slice Planarian specimen (120,473 cells) and one
16-slice Drosophila specimen (28,103 cells), using 2,000 genes and 50 shared PCs
per specimen. Every cell and slice assignment was retained; input XY and Z
were unchanged. Saved rigid transforms replayed output XY within 5.1e-11.
No annotation or original reference coordinates entered those runs. This is
functional validation, not a new expression-mode A2/A5 benchmark.

Synthetic tests cover joint-basis reconstruction, ID redistribution, absent
annotation, absent Z with explicit slice order, sparse and dense expression,
matrix-state declarations, input tampering and invalid features. They do not
replace real accuracy evaluation. See
[enhancement_validation.json](../provenance/enhancement_validation.json) for
validation scopes and [source_migration.json](../provenance/source_migration.json)
for implementation hashes. The package contains no benchmark data or reference
coordinates.
