# Validation and limits

Run the tests from the repository root in a compatible native Spateo environment:

```bash
PYTHONPATH=/path/to/spateo-release python -m pytest tests/native_skills -q
```

Tests use synthetic 3D H5AD inputs and real public Spateo calls for cross-stage alignment, mapping, SparseVFC, trajectories, all five geometric metrics, GLM and GP. They verify output identity, original coordinate preservation, separate persisted checkpoints, input/implementation/output hashes, partial-run resume, failed-run state and dependency-aware child reuse. A test import guard rejects any Dynamo import. No biological accuracy claim follows from these smoke tests.

The repository VALIDATION.md records executed counts, runtime and any warnings. Publication mesh rendering, GPU runs, full-size protocol data, numerical equivalence to the historical notebooks, and biological validity require separate experiments. Offline dashboard generation is tested independently of browser/GPU capabilities.
