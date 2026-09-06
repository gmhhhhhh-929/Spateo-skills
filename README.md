# Spateo Skills

Reusable agent skills for preparing a Spateo environment, reading data into AnnData, and aligning serial 2D tissue slices.

[中文说明](README.zh-CN.md)

## Skills

| Stage | Skill | What it provides |
| --- | --- | --- |
| 1 | [setup-spateo-environment](skills/setup-spateo-environment/SKILL.md) | Installation guidance, environment diagnosis, and capability verification for a separate Spateo source checkout or an installed package. |
| 2 | [spateo-data-io](skills/spateo-data-io/SKILL.md) | Source-backed reader selection, spatial platform detection, AnnData validation, and coordinate/metadata preservation. |
| 3 | [spateo-2d-alignment](skills/spateo-2d-alignment/SKILL.md) | Two packaged alignment pipelines, input validation, portable execution, and recorded source provenance. |

Each skill is a self-contained directory with a `SKILL.md` entrypoint. Supporting scripts and references live inside that directory, so a skill can be copied into an agent's skill search path without copying the entire repository.

## Alignment pipelines

Both pipelines belong to the **2D alignment skill**; they are execution choices within the same stage.

| Pipeline | Use | Source snapshot |
| --- | --- | --- |
| `pairwise-rigid` | Adjacent rigid Spateo alignment with expression PCA, annotation one-hot, or spatial-only input. | Previously identified as v0.2.3. |
| `continuity-guided` | Annotation-aware serial alignment with automatic continuity checks and repair. | Previously identified as v3.7.0. |

The folder names describe their behavior. Historical version identifiers are retained in provenance records to identify the exact source implementations. The [migration manifest](skills/spateo-2d-alignment/provenance/source_migration.json) records source hashes and packaging changes.

## Use

Clone this repository and let your agent read the relevant `SKILL.md`. With Codex, individual skills can also be installed by copying their directories into `~/.codex/skills/`.

```text
Use $setup-spateo-environment to prepare and verify a Spateo environment.
Use $spateo-data-io to inspect this dataset and convert it to AnnData.
Use $spateo-2d-alignment to align these ordered tissue slices with continuity-guided.
```

Spateo itself is installed from a **separate source checkout** or a compatible environment; this repository contains skills and alignment pipelines, not the complete Spateo library. The environment skill documents the required checkout and installation commands.

Data IO preserves the meaning of counts, annotations, spatial coordinates, units, and cell identifiers. It does not establish registration or silently manufacture an `X_pca`. Follow the alignment skill's input contract before selecting a representation and running either pipeline.

## Sources and verification

The environment skill and Data IO API review are based on [gmhhhhhh-929/spateo-release](https://github.com/gmhhhhhh-929/spateo-release/tree/d6aa68addc475dd0b56f69cebe7823b1f79933a9), commit `d6aa68addc475dd0b56f69cebe7823b1f79933a9`.

In that source tree, `spateo/data_io.py` is the AnnData compatibility entrypoint; maintained reader implementations are in `spateo/io/`. The Data IO skill follows those implementations and distinguishes automatic detection from explicit reader calls.

See [validation](VALIDATION.md) for executed checks and their limits, and [licensing and source notices](NOTICE.md) for the status of imported material. This initial collection contains the first three workflow stages.
