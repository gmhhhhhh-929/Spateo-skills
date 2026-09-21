# Automatic reading and explicit outcomes

`read_spatial(path, *, technology=None, load=True, load_images=True, max_memory_bytes=1073741824, max_files=10000, max_depth=4, stereoseq_bin_size=None, stereoseq_chemistry=None)` returns SpatialReadResult. The skill CLI defaults optional image loading off.

| Outcome | Meaning and action |
| --- | --- |
| entry `ready` | Core ID/value/coordinate validation passed; inspect optional-asset diagnostics and semantics before downstream use. |
| entry `deferred` | Loading was not requested or resources were insufficient. `entry.load(max_memory_bytes=...)` retries its resolved adapter, preserving diagnostics. |
| entry `unresolved` | Alternatives cannot be uniquely resolved by format contracts. Select a precise input directory or resolve source layout; no confidence override exists. |
| entry `failed` | Inspect diagnostics for missing dependencies, format/ID violations or read errors; preserve the report. |
| result `ok` | Every named input is ready and no scope-level error remains. Multiple ready inputs still cannot use `.adata`. |
| result `pending` | All entries deferred and no scope error. |
| result `partial` | Mixed ready/deferred/failure or a scope error. Never equate a saved subset with complete success. |
| result `failed` | No usable/deferred result. Missing or unsupported paths return diagnostics. |

`result.adata` requires exactly one entry and overall `ok`. `result.report` is regenerated after deferred loads; `result.write_report(path)` explicitly saves diagnostics and requires an existing parent directory. Reading itself writes no files.

Inventory depth/file limits and unfollowed symlinks can leave discovery incomplete. Resolve the scope before claiming complete ingestion. Increasing the allocation budget does not guarantee RAM availability; the budget includes already loaded collection entries. Changed sources between discovery and load fail the contract; rediscover instead of reusing stale loaders.

Automatic core readers retain all features and validate ID joins. Direct platform readers may filter or use different metadata fallbacks. Optional images/geometry are not all loaded by the core reader. Native provenance contains inventory, validation and resolution reasoning, with no probability or confidence score. The inventory is not a SHA256 manifest and external assets are not made portable by moving an H5AD.
