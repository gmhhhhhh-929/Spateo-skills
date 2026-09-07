# Skills completeness audit — 2026-09-07

The previous publication had only three discoverable skill entries. The source
collection's review, QC, sampling, component, replay and record skills were
excluded when packaging the two primary execution paths. This was incomplete
for the requested end-to-end alignment workflow.

## Result

- All **14 SKILL.md entries** in the user-designated server `skills/` directory
  are now represented; the repository has **17 entries** including the three
  existing environment, Data IO and alignment skills.
- `spatial-before-after-viewer` and `spatial-pointcloud-viewer` each include their
  executable renderer and a script-hash record. They can be installed separately.
- Other companion workflows resolve shared support scripts in
  `skills/spateo-2d-alignment/pipelines/pairwise-rigid`. Install those skills
  together with `spateo-2d-alignment`, preserving sibling directory names.
- `spateo-continuity-first-serial-alignment` is a compatibility entry to the
  published `continuity-guided` implementation. Its old duplicate scripts and
  rejected experimental engine are deliberately not reintroduced.
- The two pipeline directories remain `pairwise-rigid` and `continuity-guided`.
  ROI, viewers and workflow helpers are supporting tools, not additional
  pipelines.

The complete per-skill inventory, excluded duplicate files, source timestamps,
source hashes, packaged hashes and 21 preserved published Python hashes are in
[companion_migration.json](skills/spateo-2d-alignment/provenance/companion_migration.json).
The older `source_migration.json` records the original narrower publication;
its historical exclusion list is superseded for the restored files by this new
companion manifest.

## Pre-zebrafish boundary

The baseline is personal commit `89bbd1a20414908feb380b2e28730607349a5fce`
and organization commit `c44091bbfad102fc76844467eba6f4b4fe6942da`. Their trees
were identical before this restoration. All 21 previously published Python
files remain byte-for-byte unchanged.

The server's fixed v0.2.3 support release matches every file in the original
archive retained during publication on September 6. Source skill documents
carry September 2–3 modification times; several also match the installed local
copies. The server source has no Git history, so timestamps are supporting
evidence, not a claimed historical Git revision. The downloaded snapshot and
its hashes are retained with the local audit.

No ZESTA/zebrafish image references, image-guided code, later initialization or
endpoint experiments, modified accuracy scoring, biological datasets or result
viewers are uploaded. The existing shared-PCA and optional generalized profile
from the published pre-zebrafish baseline are retained.

## Verification

Executed checks:

- 17 skill frontmatters passed the skill-creator validator.
- Fixed-release entrypoint and restored runtime hashes verified; existing
  published Python code unchanged.
- All 9 restored workflow tests passed: lock resolution, hash mismatch,
  spatial-only provenance, balanced sampling, viewer metadata/policy, clean
  export and identity mapping, and workflow state/dashboard generation.
- Both standalone viewers actually rendered the same synthetic 120-cell,
  3-slice input without sampling. Summary point/slice counts matched; the
  before/after contact sheet was opened and inspected. Embedded JavaScript
  passed `node --check` for both generated HTML files.
- Python syntax, local Markdown links and the two-pipeline directory contract
  passed the collection validator.

One **test-only** change quotes the fixture's runner path so checkouts with
spaces work. The corresponding release-manifest test hash was updated; all
runtime entrypoint hashes remain original. These are packaging/function checks,
not a new biological accuracy benchmark or browser-interaction certification.

From the repository root:

```bash
python scripts/validate_collection.py
python skills/spateo-2d-alignment/pipelines/pairwise-rigid/tests/test_lightweight_workflow.py
```

Use a Python environment with the relevant dependencies from
`requirements-local.txt` for the runtime tests and tools. Spateo alignment needs
the separate validated Spateo environment described by the environment skill.
