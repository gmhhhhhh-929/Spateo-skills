# Validation

Run against the pinned Spateo checkout with its existing Python 3.10–3.12 environment:

```bash
PYTHONPATH=/path/to/spateo-release python scripts/smoke_source.py --source-root /path/to/spateo-release
```

This verifies the pinned commit and hashes of every Python file under `spateo/io`, then executes the source IO suite, native-runtime tests and the skill's behavioral CLI tests. The skill checks full and partial collection outcomes, lazy selected loading, discovery without full-file hashing, native domestic auto/direct routes, stable feature IDs, reversed coordinate joins, missing-file recovery, CLI exports, H5AD roundtrips and rejected overwrites. Readers are real; fault-injection tests separately simulate manifest serialization/publication failures and verify no false successful output is left behind.

The domestic independent-reader tests disable automatic entry points and parsing, trace both routes to the same `_tech.py.read_core`, and recursively inspect imports. Additional legacy-platform tests disable the automatic entry point while exercising direct readers. This does not imply that every public direct API is independent: `read_stereoseq` remains a convenience wrapper around automatic reading, and Seq-Scope has no automatic route. Expected counts, coordinates, IDs and images are specified independently of the parser; matching two readers alone is not sufficient evidence of correctness.

After reviewing and committing an intentional source update, regenerate the source API and manifest:

```bash
python scripts/refresh_source_manifest.py --source-root /path/to/spateo-release
python scripts/refresh_source_manifest.py --source-root /path/to/spateo-release --check
```

Generation refuses uncommitted IO changes so that the named commit describes the hashed files. Updating the manifest is not a substitute for tests or API review. The CLI records its expected source commit without claiming to have verified the running installation; use the smoke runner to establish that match.

For a focused CLI run:

```bash
PYTHONPATH=/path/to/spateo-release python -m pytest scripts/test_spateo_io.py -q
```

The tests use tiny synthetic native files, including two observations and duplicate gene symbols, rather than preassembled H5ADs for the new domestic input readers. Inspect the actual command output and recorded run reports for counts; a previous run's test count is not a claim about the current checkout.

Source tests cover multiple synthetic platform contracts and errors. They do not establish all vendor-version compatibility, biological annotations, large-data performance, or alignment accuracy. A successful numeric audit does not identify normalized data as raw counts. An optional-asset warning can coexist with a ready core matrix.

## Current verified revision, 2026-09-30

Pinned source `22b91af5888930838a1d4a46057404db86574ffa`: all **48 IO Python file hashes** match. The complete source IO, native-runtime and skill CLI smoke finished with **431 passed, 1 skipped** in 24.99 seconds. The source IO scope is **405 passed, 1 skipped**, including **168 new checks** (73 10x, 61 imaging-platform, 34 other native-reader checks). The skipped optional source test needs a supplied real-Visium path; this round separately completed the full real Visium audit below.

All 18 changed/new source Python files passed compilation, isort/Black and whitespace checks. Repository-wide `make check` still reports three unchanged baseline isort failures, documented in the source audit; it is not a passing whole-repository check. The smoke log retains 99 expected fixture/runtime warnings. Skill frontmatter and whitespace checks also passed.

## Accuracy audit and limits

The [non-domestic reader validation](https://github.com/gmhhhhhh-929/spateo-release/blob/22b91af5888930838a1d4a46057404db86574ffa/docs/technicals/non_domestic_reader_validation_zh.md) and its [machine-readable record](https://github.com/gmhhhhhh-929/spateo-release/blob/22b91af5888930838a1d4a46057404db86574ffa/docs/technicals/non_domestic_reader_validation_20260930.json) describe per-platform coverage, fixes, successful and rejected cases, and intentional route differences. Report actual source-test and smoke-run totals for the pinned revision after running the commands above; these are case pass rates, not a population-wide accuracy estimate.

Independent native fixtures cover H5/MEX, CSV/TSV/gzip/Parquet, GeoJSON and GEM/GEF layouts with shuffled IDs, duplicate gene symbols, leading-zero identifiers, large counts, fractional XY/XYZ, optional images, invalid inputs and H5AD roundtrips. Tests compare expected source values independently of the parser. Processed STARmap values, 10x Gene Expression filtering, bin center/origin conventions, FOV frames, direct-only Seq-Scope and the Stereo wrapper are recorded separately rather than forcing identical outputs.

Fresh complete-source audits in this round use the native Visium `V1_Adult_Mouse_Brain` and Slide-seq `Puck_180413_7` datasets. Both automatic and direct matrices, identities and coordinates match their complete original sources, with four successful H5AD roundtrips. Visium hires/lowres pixels and scales agree exactly. Slide-seq's original 6030 × 6030 bead image is loaded by the direct reader; automatic reading retains its path as `deferred_resource` because its decoded 36,360,900 bytes exceed 32 MiB. This is intended asset budgeting, not a core count failure. Other platforms in this round have native-schema fixture coverage rather than new public production-data audits; there is no fresh MOSTA/ARTISTA full read in this round.

The earlier [domestic reader validation](https://github.com/gmhhhhhh-929/spateo-release/blob/82002ba0910a0fa29874f1da92e49a76bc9186a7/docs/technicals/domestic_reader_refactor_validation_20260930.json) separately records 237 source IO passes, one skip, 38 independent-reader checks, and the full native BMKMANU audit. It is historical evidence for that revision, not the current test total.
