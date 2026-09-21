# Validation

Run against the pinned Spateo checkout with its existing Python 3.10–3.12 environment:

```bash
PYTHONPATH=/path/to/spateo-release python scripts/smoke_source.py --source-root /path/to/spateo-release
```

This executes native automatic-reading, Stereo-seq and native-runtime source tests plus the skill's behavioral IO tests. The skill checks full and partial collection outcomes, deferred discovery, CLI exports, H5AD roundtrip preservation, rejected overwrites, and failed data contracts. No reader is mocked. See the repository `VALIDATION.md` and machine-readable reports for the actual run and limits.

Source tests cover multiple synthetic platform contracts and errors. They do not establish all vendor-version compatibility, biological annotations, large-data performance, or alignment accuracy. A successful numeric audit does not identify normalized data as raw counts. An optional-asset warning can coexist with a ready core matrix.
