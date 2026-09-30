# Automatic, lazy and recovery outcomes

`read_spatial(path, ...)` returns `SpatialReadResult`. It identifies source contracts and selects a unique valid reader; it does not rank platforms by confidence. See [source-api.md](source-api.md) for the pinned signature. The skill CLI defaults optional image loading off.

## Choose the loading mode

```python
import spateo as st

# Inspect candidates, contracts and recovery advice without loading the core.
inspection = st.io.read_spatial('/data/sample', load=False, load_images=False)
print(inspection.report)

# Create on-demand handles. No complete core matrix has been validated yet.
pending = st.io.read_spatial('/data/sample', lazy=True, load_images=False)
print(pending.report)  # metadata-only; normally pending / deferred
adata = pending.adata  # loads only when the scope has one unambiguous entry
print(pending.report)  # fresh state after the attempted materialization
```

Do not combine `lazy=True` and `load=False`: use lazy handles or a discovery-only inspection. Lazy loading defers allocation; it does not turn text/MTX/GEM into an out-of-core AnnData. If the materialized matrix exceeds the budget, no partial matrix is returned. For an existing H5AD that must stay disk backed, use `anndata.read_h5ad(path, backed='r')` and close its file handle when finished.

For a collection, inspect entries without loading them, then choose exact discovered keys:

```python
result = st.io.read_spatial('/data/collection', lazy=True, load_images=False)
for key, entry in result.datasets.items():
    print(key, entry.technology, entry.representation, entry.status)
key = 'EXACT_KEY_FROM_REPORT'
entry = result.datasets[key]
entry.load(max_memory_bytes=4 * 1024**3)
if entry.status == 'ready':
    adata = entry.adata
else:
    print(entry.to_dict()['diagnostics'])
```

`entry.materialize()` returns AnnData or raises when unavailable. `entry.load()` returns the updated entry so callers can handle a deferred or failed result explicitly. Repeated attempts with an unchanged budget do not silently repeat an expensive failed materialization. To intentionally retry a recoverable failure, use `retry=True`; if source files changed, run `read_spatial` again to rediscover them. Do not reuse a stale handle after downloading or replacing files. Loading does not silently select another technology.

## Status and scope

| Outcome | Meaning and action |
| --- | --- |
| entry `ready` | Core ID/value/coordinate validation passed. Review optional assets, provenance and semantics. |
| entry `deferred` | Core loading is pending, was not requested, or exceeded resources. Read the diagnostic and its recovery actions. |
| entry `unresolved` | More than one valid reading interpretation remains. Use a precise input directory, or a verified technology declaration supported by the API. |
| entry `failed` | A required file, dependency or data contract failed. Preserve the report and correct the specific problem before retrying. |
| result `ok` | Every named input is ready and no scope-level error remains. Multiple ready inputs still cannot use `.adata`. |
| result `pending` | All entries are deferred and no scope error remains. This does not mean all values have been checked. |
| result `partial` | Mixed ready/deferred/failure or a scope error. A saved subset is not complete ingestion. |
| result `failed` | No usable/deferred result. Missing or unsupported paths return diagnostics. |

`result.report`, `.status`, key lookup and iteration do not materialize lazy data. `result.adata` is the deliberate loading trigger only in lazy mode and only for a unique complete scope. `result.write_report(path)` saves current diagnostics; the JSON does not serialize loader closures and cannot resume reading by itself. Recreate a result from the source in a new process.

Inventory depth/file limits and unfollowed symlinks can leave discovery incomplete. Resolve the scope before claiming complete ingestion. The resource budget includes already loaded collection entries; it is not an OS-enforced memory limit. Initial file probes estimate only the structures they can inspect cheaply and do not prove whole-run memory sufficiency.

## Recover without inventing data

Diagnostics in `result.report` or `entry.to_dict()` contain structured `recovery` actions with a message, relevant paths and official help links. The raw `entry.diagnostics` list is the underlying event history and does not add those report-only suggestions. They are advice (`automatic: false`), not downloaded files or completed repairs.

- **Missing matrix/ID/coordinate files:** identify the exact required companion files and retain the diagnostic. Obtain the complete output directory from the sequencing provider or the verified study's official download page. A platform documentation URL is not a sample-specific download URL. Do not invent one or reuse another sample's positions.
- **Download/access problems:** use the repository's documented HTTPS, browser, CLI or authenticated portal workflow. For a verified public URL, a resumable command such as `curl --fail --location --continue-at - --output filename 'VERIFIED_URL'` can be suggested; verify expected size/checksum and archive contents afterward. If the server lacks range support, restart into a new file. Do not copy session credentials into commands or reports. Native reading itself performs no network downloads.
- **Archive or nesting:** inspect/unpack the authorized archive to a new directory, keep related filenames together, and point discovery at the extracted sample root. Increase `max_depth`/`max_files` only after examining why scope was incomplete.
- **Memory:** omit optional images, select a specific sample/representation, or allocate a budget supported by the machine. Stereo-seq binning is a scientific output-resolution choice; do not automatically coarsen counts to make an import pass.
- **Dependencies:** use the existing Spateo environment and install only the missing dependency appropriate to the affected reader. Do not replace the entire environment for one optional backend.
- **ID, coordinate or count violations:** report the failing fields/identities. Correct the source export using verified mapping information; do not use row order, invented coordinates, fill values or duplicated labels as a repair.
- **Optional assets:** a ready core may retain deferred/missing image or geometry warnings. Explain what is and is not loaded and preserve the source path. A TIFF preview is not proof of count-image registration.

Automatic core readers retain all features and validate ID joins. Direct platform readers may use different filtering or metadata conventions; use [the comparison rules](platforms.md#comparing-automatic-and-direct-results) before interpreting a shape, coordinate or image difference as an error. Native provenance records inventory, validation and resolution reasoning with no confidence score. The inventory is not a SHA256 manifest, and external assets are not made portable merely by moving an H5AD.
