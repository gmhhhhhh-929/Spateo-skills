---
name: spatial-slice-quality-qc
description: Quality-control ordered spatial-transcriptomics slices before alignment, using bilateral tissue-loss and capture-loss evidence, multi-domain review, display-only rigid preregistration, and traceable retain/exclude outputs. Use with the companion spatial-slice-quality-viewer for the current multi-panel report.
---

# Spateo Referee QC

Use this skill on one biological series or an explicit collection of independent
series. Preserve source H5AD, coordinate keys, cell identities and prior outputs.
Read [methods and complete workflow](references/methods.md) before interpreting
thresholds, modifying the policy or preparing scientific methods text.

For questions about individual workflow controls or ambiguous input/output, use
[the node-by-node Chinese guide](references/NODE_EXPLANATIONS.md) or
[the clickable diagram](references/workflow_explained.html). Node IDs link the
explanations to the editable SVG; regenerate these together after diagram edits.

## Current decision contract

Use `policies/tissue_loss_v3.json` for new v3 runs. The v3 resolver adds a
separate bilateral loss route to the retained joint-review v2 score route.
Point loss must be supported by area or density loss; capture loss requires
both count and gene deficits with explicit capture-evidence availability.
The primary window and at least `min(2, available windows)` must confirm the
same loss route. Do not require a second evidence domain or an aggregate-score
threshold for this separate route. See methods for exact definitions.

Anatomical protection requires positive evidence of a smooth local transition
in both point count and area, normal molecular/continuity evidence, and no
independent loss evidence. Normal expression alone is not proof of intact tissue.
The legacy joint-review route keeps stage-1 thresholds K=0.129 and E=0.700.
High-score direct exclusion requires detector `exclude`, two-sided neighborhood
and no anatomical protection; otherwise the row enters internal review.
The review boundary is 0.540. **Both bands execute evidence checks.**
Use the packaged policy JSON for the actual domain/severity/window requirements;
never substitute an older keep-only band or copy a tier from another policy.
A lower score means weaker aggregate anomaly, so its exclusion evidence gates
must be no weaker than those of the higher score band. The runtime validates
this relation. Say "higher/lower exclusion-evidence threshold", rather than an
ambiguous "stricter QC".

This v3 policy is experimental; historical v2 metric-stress results do not
validate v3. Apply it using `--application-scope experimental_policy`.
That scope clears certified calls and labels the audit and viewer accordingly.
Do not set independent-validation flags or reuse historical certification for it.
Public actions are retain/exclude; the engine and audit use `keep` as the retain
alias. Internal review and each failed condition remain in the audit. Retain
means evidence did not justify exclusion, not proven health. Endpoints,
whole-series degradation and missing reference measurements can remain missed.

## Runtime and inputs

Use a compatible Spateo environment (tested with Python 3.10, NumPy 1.26.4,
SciPy 1.13.1, pandas 2.2.3 and AnnData 0.10.9). The bundled `runtime/` contains
only the two current Referee modules. The existing adapter loads these in memory
alongside installed Spateo; it never edits the installed package. A source
checkout can be selected explicitly with `--spateo-source`.

Resolve the slice column or discrete z, actual order, original x/y and expression
layer explicitly. One-H5AD-per-slice directories use natural filename order;
use `--manifest` with path/order when it differs. Multiple specimens must not
share one artificial neighbor sequence. Verify count semantics: annotation
one-hot is not measured expression; normalized X is not raw captured counts.
Missing expression stays missing and cannot support capture evidence.
For v3, rescan raw inputs or explicitly rescore a complete metric cache carrying
the original geometry and measured-capture availability contract. A legacy audit
without `loss_evidence_version=bilateral-v3` is not a v3 recalculation. Do not
invent availability from finite normalized proxies. `SliceQCConfig` enables
`tissue_loss_enabled` by default; set it false only for a labeled legacy replay.

## Run the current workflow

From this skill directory, using new output directories:

```bash
python scripts/run_slice_quality_qc.py scan \
  --input /path/to/input.h5ad --output-dir /path/to/qc_run \
  --slice-key slices --spatial-key spatial --layer counts \
  --window 3 --write-display-payload --no-report
```

Use `scan-batch` for an explicit dataset manifest. `--window auto` compares 3/5/7
when supported; it changes neighborhood selection, not the locked policy.

When shared-coordinate ROI comparison is required:

```bash
python scripts/run_slice_quality_qc.py preregister \
  --input-dir /path/to/qc_run --dataset-id specimen_A
```

The display fit is rigid 2D: seed 13, at most 1000 fit points, 24 angles,
50 iterations and the closest 80% correspondences. Compose adjacent transforms
into the middle-slice reference. No scale, reflection, nonrigid deformation or
smoothing. Save matrices, reference chains, inverses, point caches and failures.
Optional observed-annotation ranking uses `--annotation-weight 0.5`; missing
support and ambiguous fits remain visible. QC always uses original coordinates.

Publish and render through the companion viewer skill:

```bash
python subskills/spatial-slice-quality-viewer/scripts/build_viewer.py \
  --input-dir /path/to/qc_run --output-dir /path/to/new_report \
  --policy policies/tissue_loss_v3.json --allow-unvalidated-policy \
  --application-scope experimental_policy \
  --language en
```

This runs the existing publication function with complete binary resolution,
then the existing detailed renderer. For a batch, supply several input
folders with distinct basenames. To publish tables without rendering, use
`publish --complete-binary --policy ... --allow-unvalidated-policy --application-scope experimental_policy`.
The explicit `--allow-unvalidated-policy` flag consents to experimental use,
not certification. Its default is false; do not silently bypass a validation
failure. Alternatively omit `--policy` to render an existing complete v3 audit
without changing its decisions. Primary and confirming rows must carry the
same `loss_evidence_config` as the supplied policy; mismatches require rescoring.

## Verify and hand off

Check every original slice appears exactly once; every internal review has a
matched tier, complete window diagnostics and a recorded outcome. Confirm both
score bands are enabled in the actual policy. Compare any new exclusions with
the prior audit; preserve score, source and coordinate provenance.
Also verify `tissue_loss_exclusion_gate`, route, bilateral deficit fractions,
primary-window confirmation, confirming/available windows and failed reasons.
Do not count initial review candidates as final detections. For known injected
defects, detection is final exclude; report real-data exclusion counts separately
from sensitivity or accuracy, which require ground-truth labels.

The viewer must show the current multi-panel evidence page and only final
retain/exclude states. It must retain preregistration warnings and experimental
policy scope. Do not describe a completed fit as anatomical correspondence.

For an exported display ROI, replay on full cached points:

```bash
python scripts/run_slice_quality_qc.py export-roi \
  --input-dir /path/to/qc_run --roi /path/to/roi.json \
  --output-csv /path/to/new_selected_points.csv
```

Verify frame_id and cache hashes. Exact inverse selection is a polygon; a raw
axis-aligned bounding box is not equivalent. Sample counts in the browser are
not full-data counts. Never delete source slices or apply alignment policy
files automatically.

Historical experiment launchers, dataset-specific rescue scripts and the old
separate low-band audit prototype are intentionally absent. Publication now
executes both review bands through the same scientific evaluator. The retained
`slice_quality_report.py` is a dependency of the collection/diagnostic interfaces;
`slice_quality_visualization.py` generates the current final detailed report.
