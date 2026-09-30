# Validation

Run against the pinned Spateo checkout with its existing Python 3.10–3.12 environment:

```bash
PYTHONPATH=/path/to/spateo-release python scripts/smoke_source.py --source-root /path/to/spateo-release
```

This verifies the pinned commit and hashes of every Python file under `spateo/io`, then executes the source IO suite, native-runtime tests and the skill's behavioral CLI tests. The skill checks full and partial collection outcomes, lazy selected loading, discovery without full-file hashing, native domestic auto/direct routes, stable feature IDs, reversed coordinate joins, missing-file recovery, CLI exports, H5AD roundtrips and rejected overwrites. Readers are real; fault-injection tests separately simulate manifest serialization/publication failures and verify no false successful output is left behind.

The source's independent-reader tests also disable automatic entry points and parsing before calling each direct reader, trace both routes to the same `_tech.py.read_core`, and recursively inspect the dependency graph for imports back into the automatic layer. Expected counts, coordinates, IDs and images are specified independently of the parser; matching two readers alone is not sufficient evidence of correctness.

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

Pinned source `82002ba0910a0fa29874f1da92e49a76bc9186a7`: all **47 IO Python file hashes** match. The complete source IO, native-runtime and skill CLI smoke finished with **263 passed, 1 skipped** in 13.67 seconds in the validation environment. The source IO scope is **237 passed, 1 skipped**; the 38 new independent-reader checks all passed. They comprise 9 cases per technology, one additional Singleron raw/filtered representation check and one recursive dependency guard. The skipped real-Visium test needs an externally supplied dataset path.

The [refactor validation record](https://github.com/gmhhhhhh-929/spateo-release/blob/82002ba0910a0fa29874f1da92e49a76bc9186a7/docs/technicals/domestic_reader_refactor_validation_20260930.json) records a fresh real-native BMKMANU audit: direct and automatic reading agree exactly on counts, `obs`, `var`, coordinates and PNG images; all 25,239,573 original matrix records and the H5AD roundtrip match. The recorded reader is `spateo.io.spatial._bmkmanu.read_core`. Metadata equality is checked before H5AD serialization, which can categorize string columns in memory. The other three technologies use independently specified native-schema fixtures, not real production datasets in this round.

All 25 changed source Python files passed compilation, isort/Black and whitespace checks. The repository-wide `make check` still finds three unchanged baseline isort issues, recorded separately; this is not claimed as a passing whole-repository check. Runtime warnings during the combined smoke are retained in its log and do not change its test results.

## Previous baseline, 2026-09-30

Pinned source `d884216b2f030bcb5616b01a8f950d31e1e6d833`: all 34 IO Python file hashes match. The full skill smoke (source IO, native runtime, and CLI behavior) completed with **225 passed, 1 skipped**. The focused CLI scope contains 19 passing checks: 12 packaged here and 7 existing repository checks. `quick_validate.py` and whitespace validation also passed.

The [source validation record](https://github.com/gmhhhhhh-929/spateo-release/blob/d884216b2f030bcb5616b01a8f950d31e1e6d833/docs/technicals/domestic_spatial_io_validation_20260930.json) separately records 199 passing source IO tests, one skip and 76 new domestic/lazy/recovery cases. Its public native BMKMANU example (GSM8816652) has 50,970 observations × 54,752 features; all 25,239,573 source matrix entries, 39,824,941 total counts, coordinates, axis IDs and H5AD roundtrip were checked exactly. This is one real public dataset and two successful reads of it, not a cross-platform biological accuracy estimate. The other new platforms currently have verified-schema fixture tests.
