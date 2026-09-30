---
name: spatial-slice-quality-viewer
description: Generate the current Spateo Referee retain/exclude report with overview, statistics, multi-metric and bilateral-loss evidence, all slices, methods, and shared-frame ROI. Supports one QC run or a dataset collection using the companion QC skill's canonical renderer.
---

# Spateo Referee Viewer

Install as part of the complete `spatial-slice-quality-qc` directory. This skill owns the orchestration;
the QC skill owns the canonical renderer and scientific runtime. Do not copy
templates into a second viewer or convert review labels in presentation code.

## Build

```bash
python scripts/build_viewer.py \
  --input-dir /path/to/qc_run --output-dir /path/to/new_report \
  --policy ../../policies/tissue_loss_v3.json \
  --allow-unvalidated-policy \
  --application-scope experimental_policy --language en
```

The v3 policy evaluates bilateral tissue/capture losses in parallel with the
legacy low- and high-score evidence checks. Historical v2 validation is not v3
validation; the experimental policy has no independent biological certification.
The report must display route, actual deficit fractions, primary/confirming
window evidence, missing-evidence limitations and retain/exclude public labels.
The audit and engine retain their compatible `keep` alias; do not relabel review.
The current gate values, score branches and source functions are defined in
[the QC methods](../../references/methods.md).
Do not describe the low band as keep-only, or mark the new policy independently
certified. Preserve a historical input audit only when the user requests that
policy; omitting `--policy` renders its existing complete decisions unchanged.
The explicit development flag does not grant biological certification. V3
requires both primary and confirming evidence configuration to match its policy;
rendering cannot manufacture missing evidence or repair a configuration mismatch.

Accept metrics, manifest and display payload. Without `--policy`, require a
complete binary audit matching every input slice and its original order.
The bundled Referee runtime is selected automatically; `--spateo-source` is an
explicit checkout override. `--qc-skill` selects a non-sibling installation.
Output must be new. Several distinct input directory names create a collection.

For an input with preregistration, reuse verified raw/display coordinate caches
and original capture values. Never replace them with another coordinate key to
make the layout prettier. `--source-h5ad` is available for a single dataset when
source point-level measurements are needed; honor the manifest contract.

## Verify

The detailed report retains Overview, Statistics, Slice evidence, All slices,
and Methods. Evidence contains fixed-bandwidth KDE, captured counts or detected
genes when available, connected components, and linked 3/5-slice windows.
Unavailable values stay unavailable; sample-only calculations are disclosed.
Methods must show actual per-band gate values, not stale hardcoded thresholds.

Verify final calls against the generated audit; low-band review is fully checked
when the joint policy is supplied. With v3, verify the additional tissue-loss
gate against the audit even for low-score slices; it does not require a second
anomaly domain. Review remains internal. Confirm experimental
policy and input-validation limitations appear in each report and collection.

Check raw/display switching, shared coordinate extents, focal-slice ROI dragging
in all evidence panels, neighboring selections, sample/total counts, fit warnings,
frame_id and inverse-polygon export. The same rectangle is not proof of the same
anatomy. Full-point replay uses the QC skill's `export-roi` against the original
run recorded in `viewer_workflow.json`.

Use browser inspection when available. Distinguish real browser verification from
DOM/canvas substitutes or static checks. Keep the input/audit/script hashes and
source run path in `viewer_workflow.json`; no source H5AD is edited.
