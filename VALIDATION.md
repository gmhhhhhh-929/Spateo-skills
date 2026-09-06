# Validation and scope

This publication packages three skills and two alignment pipelines. It does not include biological datasets or historical benchmark results.

## Source basis

Environment and Data IO source: `gmhhhhhh-929/spateo-release` commit `d6aa68addc475dd0b56f69cebe7823b1f79933a9`. The Data IO source manifest records 32 reviewed implementation/test file hashes. Alignment sources are identified separately in the [migration manifest](skills/spateo-2d-alignment/provenance/source_migration.json).

## Executed checks

- All three `SKILL.md` entrypoints passed the skill creator's `quick_validate.py` frontmatter and naming checks.
- The environment verifier passed its nine contract tests on Python 3.9 and Python 3.12. They cover installation identity, portable paths, missing dependencies, JSON output and structured failure handling. The original two scientific smoke functions retain the source AST.
- Data IO source hashes match the pinned checkout. Thirteen synthetic API/CLI cases and thirteen unchanged source IO tests passed. They call the real source readers and check matrix orientation, ID joins, coordinate axes, images, H5AD roundtrips, ambiguous detection and rejection of invalid alignment representations. See [Data IO validation](skills/spateo-data-io/references/validation.md) for the exact executed cases, source tests and package versions.
- The five imported alignment source files passed `scripts/verify_bundle.py`. The two pairwise runners are byte-identical to their source snapshots. Continuity dependencies retain the scientific function bodies documented in the migration manifest; import names, the removal of unused standalone entrypoints, and UTF-8 metadata hashing for non-ASCII cell IDs/annotations are recorded. The UTF-8 change affects provenance serialization, not alignment coordinates or scoring.

- Independent tests copied the skills to a separate directory and passed the thirteen IO cases, eight alignment packaging/input/geometry checks, command-line help/dry-run, and coordinate identity export. This caught and resolved XYZ preservation, zero-norm PCA, unreviewed nested slice files, empty annotation labels and non-ASCII metadata hashing.
- Three bounded CPU executions passed in an existing Spateo environment: pairwise expression-PCA, pairwise spatial-only and continuity-guided. Each used two slices of 80 synthetic cells, three iterations and batch size 40. All 160 cell IDs were retained uniquely, coordinates were finite, and outputs were frozen. The pairwise runs also verified that inherited runner settings could not override the explicit mode and iteration arguments. The environment reported Spateo 1.1.2, Torch 2.5.1+cu121, NumPy 1.26.4, SciPy 1.13.1 and AnnData 0.10.9; these runs used CPU. These CPU runs preceded the final preflight and UTF-8 hashing fixes, which were covered by the subsequent local synthetic and independent checks. The [CPU report](skills/spateo-2d-alignment/provenance/cpu_smoke.json) records this sequence. They do not certify the separately pinned source checkout or measure biological accuracy.

## Reproduce focused checks

Run from the repository root, using the interpreter appropriate to each skill:

```bash
python skills/setup-spateo-environment/scripts/test_verify_environment.py
python skills/spateo-2d-alignment/scripts/verify_bundle.py
python skills/spateo-2d-alignment/scripts/check_runtime.py --device cpu
python skills/spateo-2d-alignment/scripts/smoke_test.py
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$SPATEO_CHECKOUT" \
  python skills/spateo-data-io/scripts/smoke_source.py --source-root "$SPATEO_CHECKOUT"
```

The environment contract tests require `packaging`; the scientific and Data IO checks require the libraries documented by their skills. A minimal YAML-validation environment is not a Spateo runtime.

## Limits

The source package supports Python 3.10–3.12. The local Data IO reader tests used an existing Python 3.9 scientific environment; they are partial reader verification, not a supported full Spateo installation certification. The intended source revision's alignment API was checked statically; that alone does not establish numerical or GPU compatibility.

A complete fresh Spateo installation, full environment scientific/3D smoke, all platform datasets and the large biological GPU accuracy benchmark were not rerun for this packaging task. Run the supplied runtime checks in the execution environment and retain the results for each analysis.
