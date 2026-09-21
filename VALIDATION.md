# Validation and scope

## Native IO / 4D update — 2026-09-21

Source: `gmhhhhhh-929/spateo-release@615644f88613bea8ceb2e2df1e2391d16de55ec1`. Protocol reference: `gmhhhhhh-929/Spateo-protocol-files@b11ae99fbdc4ae46d41880e9306ab7e5c2751ac5`. The collection now has six top-level entrypoints (including reserved 3D) and 26 total skills.

| Check | Executed result | Scope |
| --- | --- | --- |
| Native source IO / Stereo-seq / native runtime | 96 passed, including real Visium | Real public implementations, synthetic platform/error contracts, native runtime and no-Dynamo source checks. |
| New skill behavior | 19 passed | 7 IO tests; 12 4D cases including parametrized dependency invalidation. Scientific calls run without mocks, with a guard rejecting Dynamo imports. |
| Existing QC + relocated viewer | 16 passed | 15 preserved QC tests and a new scan→nested-viewer CLI→complete-report integration test. |
| Skill structure | 26 entrypoints passed skill-creator validation | YAML, names and unfinished scaffolds. |
| Collection packaging | Passed | 26 skills, 6 top-level, links, Python syntax, 25 locked alignment entrypoints, preserved snapshots and viewer hashes. |
| Frozen 2D bundle | 16 source file hashes matched | No new 2D scientific algorithm/accuracy claim. |
| Browser QA | Passed | Offline Plotly rendered; flow, target and overlay views and curl coloring worked; target=40/80 and overlay=80/80; no browser console errors. |

The source suite initially skipped the optional real-Visium test because its environment variable was unset. It was then run with the locally available public adult-mouse-brain Visium data and passed; the final IO smoke entrypoint run passed all 103 cases (96 source + 7 skill IO) with no skips. The actual IO CLI exported **2,702 spots × 32,285 features** and reread all persisted values. Its first attempt correctly rejected 40 duplicated gene symbols; explicit `--feature-id-column gene_ids` succeeded and retained original symbols in `var['source_var_name']`. No raw data are included in this repository.

The 4D synthetic run executes reference alignment, cell directions, SparseVFC, trajectories, all five geometric metrics, two successful GLM fits and one GP interpolation, then writes an offline dashboard. Tests preserve original coordinates/counts, verify mapped endpoints and transport IDs, read all checkpoints, check failed/partial run states, enforce parent immutability and verify content-hash-based reuse. The test inputs use 40 cells per stage and 12 genes, CPU, a 5-iteration alignment with nonrigid_start_iter=1, 20 field iterations and 2 GP training iterations. These are execution smoke settings, not scientific analysis defaults.

Issues found and resolved: retired confidence-based IO calls; incorrect auto-reader return assumptions; required AnnData 0.10 serialization of native preprocessing history (reversible tagged mappings); logger streams surviving closed stage logs; too-short reference-alignment iteration schedules; source/nonrigid coordinate-key inconsistency; notebook GLM formula interpolation and normalized-expression semantics; missing mapping summaries; target-view count and GLM likelihood display. Source algorithms themselves were not patched.

Warnings were retained: Torch JIT deprecation, a Numba duplicate-compilation warning, a SciPy import deprecation and statsmodels' default negative-binomial dispersion. Successful execution does not establish biological validity, lineage tracing, benchmark accuracy, parameter adequacy or numerical equivalence with old notebooks. Full protocol datasets, GPU runs, publication mesh plots and the reserved 3D pipeline were not tested or implemented in this update.

Reproduce from this repository with a compatible native environment:

```bash
PYTHONPATH=/path/to/spateo-release python -m pytest tests/native_skills -q
PYTHONPATH=/path/to/spateo-release python -m pytest tests/referee -q
python scripts/validate_collection.py
python skills/spateo-2d-alignment/scripts/verify_bundle.py
# Optional real-data test: set SPATEO_VISIUM_DATA to a valid Visium bundle.
PYTHONPATH=/path/to/spateo-release python skills/spateo-data-io/scripts/smoke_source.py --source-root /path/to/spateo-release
```

Runtime versions, source hashes and test counts are recorded in [validation/native-io-4d.json](validation/native-io-4d.json).

## Historical validation below

The remaining sections preserve earlier publication evidence and counts. They do not describe the current IO API, current total skill count, or newly executed 4D tests.

## Source basis

Environment and Data IO source: `gmhhhhhh-929/spateo-release` commit `d6aa68addc475dd0b56f69cebe7823b1f79933a9`. The Data IO source manifest records 32 reviewed implementation/test file hashes. Alignment snapshots and subsequent changes are identified separately in the [migration manifest](skills/spateo-2d-alignment/provenance/source_migration.json).

## Current alignment validation

The current bundle verifies **16 code files** through `scripts/verify_bundle.py`. This includes the preserved pairwise runners, continuity helpers and entrypoints, the frozen anonymous policy, validators and expression preparation/tests. The earlier five-file verification described the initial imported snapshot, not the full current bundle. Changes to representation handling and optional continuity behavior are recorded explicitly; current continuity code is not claimed to be byte-identical to the original release.

Executed synthetic checks comprise 13 expression-preparation tests, 10 expression-mode/provenance checks, eight existing packaging/input/geometry smoke checks, and 18 anonymous-policy/runner checks. Coverage includes joint-basis reconstruction, cell-ID mapping, no annotation, no physical z with explicit numeric slice order, dense/CSR/CSC expression, declared preprocessing state, manifest tampering, invalid representations, input immutability, exact fallback and supported-label acceptance gates. Fake Spateo calls in mode tests verify that supplied features reach the backend unchanged; they are distinct from the real runs below.

Two full real **no-annotation expression** preparations and continuity runs completed, explicitly selecting `--profile generalized`, 300 iterations and batch size 800:

| Specimen | Cells | Slices | Shared expression basis |
| --- | ---: | ---: | --- |
| One Planarian development specimen | 120,473 | 14 | 2,000 genes, 50 PCs |
| One Drosophila development specimen | 28,103 | 16 | 2,000 genes, 50 PCs |

Every cell and slice assignment was retained. Input XY and physical z were unchanged, and saved rigid transforms replayed output XY within `5.1e-11`. No annotation or original reference coordinates entered these runs. These GPU executions establish functionality on two specimens; they do not measure expression-mode A2/A5 accuracy.

Two additional real **annotation packaging equivalence** runs used the same specimens, generalized profile, 300 iterations, batch size 800 and seed 20260817. All recovered per-cell XY values were exactly equal to the separately frozen `universal_nn_off` outputs. Transform differences between CSV and JSON serialization were at most `5.684341886080802e-14`. Cell identities, slice mapping, finite coordinates, proper SE(2) transforms and unchanged input files/XY all passed. No reference data or accuracy evaluation was used in this comparison. This establishes packaging equivalence for these two development specimens, not every possible input.

The [enhancement validation record](skills/spateo-2d-alignment/provenance/enhancement_validation.json) identifies runtime code hashes and subsequent changes. After these real runs, the CLI default was restored to legacy and the pipeline identifier changed to its functional name; computation for the explicitly selected generalized profile was unchanged. PCA preparation subsequently gained all-slices-without-z ordering support and explicit rejection of conflicting Z columns, covered by 13 synthetic tests. The PCA numerical algorithm did not change; the real runtime snapshots and their inputs remain immutable.

## Accuracy interpretation and default profile

**Legacy remains the compatibility default; generalized is opt-in.** The frozen annotation-input study used one whole-specimen proper XY similarity fit, preserved physical z, and counted XYZ errors strictly below 2%/5% of the original whole-specimen XY diagonal for A2/A5. It did not refit individual slices or groups. Biological specimens were weighted equally, with repeated perturbation seeds averaged within specimens first.

The original 17-specimen comparison is retrospective. Independent perturbations additionally covered nine historical holdout specimens with two new seeds each and ten additional Drosophila stages with one seed each. Generalized improved average Drosophila results but regressed on Planarian means and had severe individual regressions. One legacy run failed while generalized completed that case with very low accuracy. Failures remain NA rather than being scored as zero; completion and success-set averages do not establish accurate recovery or paired improvement. See [the full aggregate results and limitations](skills/spateo-2d-alignment/references/validation.md).

These findings do not support a universal accuracy or stability improvement, automatic profile selection, or an expression-mode accuracy claim.

## Initial publication checks retained as historical evidence

The following environment and Data IO results describe the original publication only; current IO results supersede them. **This alignment extension did not rerun their installation or reader suites.**

- All three `SKILL.md` entrypoints passed the skill creator's frontmatter and naming checks at initial publication; current skill validation also passed for the alignment extension.
- The environment verifier passed nine contract tests on Python 3.9 and Python 3.12, covering installation identity, portable paths, missing dependencies, JSON output and structured failures. The original two scientific smoke functions retained the source AST.
- Data IO hashes matched the pinned checkout. Thirteen synthetic API/CLI cases and thirteen unchanged source IO tests passed using real source readers. They checked matrix orientation, ID joins, axes, images, H5AD roundtrips, ambiguous detection and invalid alignment representations. See [Data IO validation](skills/spateo-data-io/references/validation.md) for cases and package versions.
- Independent copied-directory checks exercised portable paths, IO, alignment geometry and CLI behavior. They caught and resolved XYZ preservation, zero-norm features, nested unreviewed slice files, empty labels and UTF-8 metadata hashing.
- Three bounded CPU runs passed in an existing Spateo environment: pairwise expression-PCA, pairwise spatial-only and continuity-guided, each with two slices of 80 synthetic cells, three iterations and batch size 40. All IDs were retained and coordinates were finite. These preceded the final initial-publication preflight/UTF-8 fixes, subsequently covered by synthetic checks. The [historical CPU report](skills/spateo-2d-alignment/provenance/cpu_smoke.json) records the sequence and environment. Those original runs did not certify the separately pinned source checkout or measure biological accuracy.

## Reproduce focused checks

Run from the repository root, using the interpreter appropriate to each skill:

```bash
python skills/setup-spateo-environment/scripts/test_verify_environment.py
python skills/spateo-2d-alignment/scripts/verify_bundle.py
python skills/spateo-2d-alignment/scripts/check_runtime.py --device cpu
python skills/spateo-2d-alignment/scripts/smoke_test.py
python skills/spateo-2d-alignment/scripts/tests_expression_mode.py
python -m unittest discover -s skills/spateo-2d-alignment/scripts \
  -p test_prepare_expression_pca.py -v
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$SPATEO_CHECKOUT" \
  python skills/spateo-data-io/scripts/smoke_source.py --source-root "$SPATEO_CHECKOUT"
```

The environment contract tests require `packaging`; scientific and Data IO checks require their documented libraries. A minimal YAML-validation environment is not a Spateo runtime.

## Remaining limits

The pinned source package supports Python 3.10–3.12. Initial local Data IO reader tests used an existing Python 3.9 scientific environment; they were partial reader verification, not a supported full Spateo installation certification. The pinned source alignment API was checked statically. The real alignment runs above used the existing execution environment and do not certify a fresh installation of that separate source revision.

A fresh complete Spateo installation, full environment scientific/3D smoke and all supported platform datasets were not rerun for this alignment extension. The real expression and annotation checks are scoped as stated above; the initial small CPU smoke is historical evidence, not a description of the current validation's entire extent. Run the supplied runtime checks in the intended environment and retain their results for each analysis.
