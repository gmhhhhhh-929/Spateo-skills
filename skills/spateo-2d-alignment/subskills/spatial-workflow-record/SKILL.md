---
name: spatial-workflow-record
description: Track, audit, and summarize spatial transcriptomics alignment workflows as chronological records, v0.2 state pointers, compact dashboards, clean coordinate links, provenance validators, and operation ledgers.
---

# Spatial Workflow Record

## Package paths

Run the relative commands below from this skill directory. Keep this subskill inside the complete `spateo-2d-alignment` directory.
The shared support tools are the fixed pre-zebrafish snapshot, located in
`../../pipelines/pairwise-rigid`. Verify entrypoints with its `skill.lock.yaml`.
This packaging restoration does not change alignment algorithms or defaults.


## Scope

Use this skill to make spatial alignment work traceable. It organizes evidence
and decisions; it must not compute transforms, edit h5ad files, overwrite
coordinate caches, or mark a candidate final without explicit user decision.

## v0.2 Record Model

Each operation step should contain:

```text
step.yaml
record.json
decision.md
provenance/
  skill_run.json
  validator_report.json
  input_manifest.json
recipes/
qc/
review/
light/
full/
exports/
logs/
```

`record.json` remains the machine-readable operation record. `step.yaml` is the
short human summary. `decision.md` records manual choices or why a candidate was
rejected/not needed.

## State Pointers

Dashboards and handoffs must read `states/*.yaml` before scanning steps:

- `slice_level_baseline.yaml`;
- `current_review.yaml`;
- `current_accepted.yaml`;
- `final.yaml`.

Each state should link to clean coordinates, audit points, review HTML, recipe
chain, validation report, row count, sha256, and validator status.

## Dashboard Guidance

The dashboard should default to a compact chronological table, not an image
gallery. The top section shows current states and highlights the clean
coordinate CSV. The operations table shows one row per step with concise action
text, status, validator status, and links to records, manifests, metrics,
recipes, logs, and viewers.

Do not display bulky audit CSVs as the latest coordinate file. Link them as
debug/audit artifacts only.
