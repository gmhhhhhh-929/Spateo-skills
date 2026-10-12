# Environment audit evidence

The shared `scripts/audit_environment.py` entrypoint serves discovery, directory documentation, runtime reporting and saved-report comparison. Use its `--help` and subcommand help for exact supported arguments. It is an audit tool, not an installer.

## Documents and provenance

Save each audit into its own new or empty output directory outside the inspected prefix and source checkout; existing nonempty outputs are not overwritten. Markdown is the human-readable document; JSON retains structured evidence. Link the generated files in the user handoff. A cloud Page, Word document or PDF is unnecessary unless the user requests that format.

| Command | Outputs |
| --- | --- |
| `discover` | `environment-discovery.md`, `environment-discovery.json` |
| `tree` | `environment-tree.md`, `environment-tree.json`, `packages.json`, `checks/inventory.json` and its log |
| `report` | `environment-report.md`, `environment-report.json`, tree outputs, `packages.json`, per-probe JSON/logs and health logs in `checks/` |
| `compare` | `environment-comparison.md`, `environment-comparison.json` |

`report --baseline` adds comparison files to that new report directory. The package snapshot records installed distributions; it is not a resolved installation lockfile or a guarantee of cross-platform reproducibility.

Keep the interpreter's absolute executable and prefix, platform/architecture, run time, selected profile, installed distributions, source identity and requested checkout visible. Distinguish imported module location, editable-install metadata and Git revision; they answer different identity questions. A requested checkout must not silently become the imported source by changing `sys.path`.

Machine paths and host identifiers belong in that user's report. Do not commit real audit outputs, SSH details or private environment paths as reusable/public fixtures. Use synthetic fixtures for repository tests. Prior setup evidence is necessary to claim that an environment was created using this skill; otherwise label its origin unconfirmed.

## Capability states

| State | Supported statement | Limit |
| --- | --- | --- |
| `smoke_passed` | The named tiny synthetic operation completed and passed its checks in this interpreter. | Does not establish full pipeline behavior, performance or biological accuracy. |
| `api_only` | The required import/API/signature checks passed. | The downstream computation was not exercised. |
| `blocked` | A required check failed, timed out or found a missing prerequisite. | Preserve the failing evidence; identify whether it is environmental or probe-specific. |
| `not_tested` | This profile omitted the check, or no applicable evidence was obtained. | Do not call it ready or count it as a pass. |

Record a probe's actual outcome separately from a downstream capability's interpretation. If a dependency check or selected checkout identity fails, overall health must reflect it. Preserve independent successful computations rather than replacing every row with the first error. A metadata-only report never establishes runtime readiness, even when all listed packages are present.

The capability matrix should cover the installed profile's relevant downstream families: data IO, preprocessing, slice quality, 2D alignment, 3D point clouds, 3D meshes, native morphogenesis/trajectories and features. Describe the operation or API actually tested instead of inferring support from a package name. The `full` profile adds computational mesh/GPU evidence; its name does not mean that every scientific stage, device backend, display or optional dependency was validated.

## Overall health and command exits

| Report overall | Meaning |
| --- | --- |
| `blocked` | Inventory, dependency consistency, source-profile verification or audit implementation consistency failed. Independent capability evidence is still retained. |
| `metadata_only` | Metadata health checks passed, but runtime capabilities were not tested. |
| `partial` | Health checks passed, but at least one requested capability probe was blocked. |
| `checked_with_limits` | Requested health/capability checks passed at their stated evidence levels; omitted or full-pipeline behavior remains untested. |

Report exit `0` covers `metadata_only` and `checked_with_limits`; exit `1` covers `blocked` and `partial`, with diagnostic reports still saved. Exit `2` denotes a command/file error. Tree exit `1` means inventory or directory listing errors. Discovery and comparison success only confirm that the requested inventory/comparison was written; inspect recorded search errors and coverage before drawing conclusions.

## Interpretation and handoff

Explain overall health first, then which tasks have computational evidence, which have only API evidence, and what failed or remains untested. Mention relevant input contracts separately: data formats, count semantics, IDs, finite coordinates, slice order, shared representations and annotations still need validation when a real task starts.

For failures, report the first useful error with a link to the complete raw output and a concrete recommendation. Do not hide timeouts or convert them to unsupported success. Optional hardware absence limits the corresponding GPU capability; it does not by itself prove CPU workflows are unusable. A headless mesh computation does not establish interactive rendering.

Compare reports with their profiles and evidence coverage visible. Preserve missing optional evidence as unknown within supported schemas; reject unsupported schemas. Distinguish newly measured evidence from a true change in an already-tested operation. A failed inventory makes package additions/removals and same-environment identity unknown (`null`); it must not be treated as an empty installed environment. Package/source drift is descriptive evidence, not a command to upgrade or downgrade.
