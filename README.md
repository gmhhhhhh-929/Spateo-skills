# Spateo Skills

Reusable agent skills for preparing a Spateo environment, reading data into AnnData, and aligning serial 2D tissue slices.

[中文说明](README.zh-CN.md)

## Skills

| Stage | Skill | What it provides |
| --- | --- | --- |
| 1 | [setup-spateo-environment](skills/setup-spateo-environment/SKILL.md) | Installation guidance, environment diagnosis, and capability verification for a separate Spateo source checkout or an installed package. |
| 2 | [spateo-data-io](skills/spateo-data-io/SKILL.md) | Source-backed reader selection, spatial platform detection, AnnData validation, and coordinate/metadata preservation. |
| 3 | [spateo-2d-alignment](skills/spateo-2d-alignment/SKILL.md) | Two alignment pipelines, shared expression-PCA preparation, optional annotation, strict input validation, and recorded implementation provenance. |

The collection contains **17 discoverable skills**: the three primary entries above and all 14 entries from the requested server skills directory. The companion workflows share the fixed support scripts in `spateo-2d-alignment/pipelines/pairwise-rigid`; install the collection together. Both viewer skills also bundle their own renderer and can be installed independently. Commands in companion skills run from the relevant skill directory.

## Alignment review and workflow companions

| Skill / 技能 | Purpose / 用途 |
| --- | --- |
| [detect-spatial-components](skills/detect-spatial-components/SKILL.md) | Connected-component detection · 组织连通分量 |
| [lightweight-spatial-alignment-workflow](skills/lightweight-spatial-alignment-workflow/SKILL.md) | Workflow orchestration · 工作流 |
| [remote-workflow-intake](skills/remote-workflow-intake/SKILL.md) | Remote execution context · 远程任务环境 |
| [spateo-continuity-first-serial-alignment](skills/spateo-continuity-first-serial-alignment/SKILL.md) | Alias to continuity-guided · 连续性配准入口 |
| [spateo-pairwise-qc](skills/spateo-pairwise-qc/SKILL.md) | Pairwise numeric and visual QC · 成对质控 |
| [spateo-pairwise-run](skills/spateo-pairwise-run/SKILL.md) | Pairwise job planning · 成对任务计划 |
| [spateo-roi-refine](skills/spateo-roi-refine/SKILL.md) | ROI/drop refinement · 局部修正 |
| [spatial-alignment-compose](skills/spatial-alignment-compose/SKILL.md) | Recipe composition and replay · 变换组合与重放 |
| [spatial-balanced-sample](skills/spatial-balanced-sample/SKILL.md) | Component-balanced sampling · 平衡采样 |
| [spatial-before-after-viewer](skills/spatial-before-after-viewer/SKILL.md) | Before/after and displacement viewer · 配准前后及位移查看 |
| [spatial-component-align](skills/spatial-component-align/SKILL.md) | Component candidates · 分量配准候选 |
| [spatial-pointcloud-viewer](skills/spatial-pointcloud-viewer/SKILL.md) | Full-points 3D, slice/pair/component viewer · 全点三维查看 |
| [spatial-triad-component-align](skills/spatial-triad-component-align/SKILL.md) | Three-slice QC and candidates · 三切片检查 |
| [spatial-workflow-record](skills/spatial-workflow-record/SKILL.md) | Records, state pointers and dashboard · 记录与看板 |

This restoration keeps the published alignment core at personal commit `89bbd1a` / organization commit `c44091b` unchanged. It includes no zebrafish image priors, later experimental alignment refinements, or later scoring changes. See the [completeness audit](COMPLETENESS.md).

## Alignment pipelines

Both pipelines belong to the **2D alignment skill**. Both accept verified shared expression PCA, annotation one-hot, or spatial-only input.

| Pipeline | Use | Source snapshot |
| --- | --- | --- |
| `pairwise-rigid` | Adjacent rigid Spateo alignment. | Previously identified as v0.2.3. |
| `continuity-guided` | Serial Spateo alignment with automatic continuity checks and optional repairs. | Derived from the v3.7.0 snapshot, with explicit representation and profile controls. |

For continuity-guided, **`--profile legacy` remains the default**. The optional `--profile generalized` disables nearest-neighbor initialization and replaces the named-tissue priority with anonymous supported-label proposals. Frozen validation found average Drosophila improvements, Planarian regressions and severe individual failures. It does not establish a generally better profile. See [accuracy and validation limits](skills/spateo-2d-alignment/references/validation.md).

Historical version identifiers, source hashes and subsequent implementation changes are recorded in the [migration manifest](skills/spateo-2d-alignment/provenance/source_migration.json).

## Use

Clone this repository and let your agent read the relevant `SKILL.md`. With Codex, copy the contents of `skills/` together into your configured skill directory; keep companion dependencies as siblings. The two viewers may also be copied independently.

```text
Use $setup-spateo-environment to prepare and verify a Spateo environment.
Use $spateo-data-io to inspect this dataset and convert it to AnnData.
Use $spateo-2d-alignment to align these ordered tissue slices using shared expression PCA without annotation.
```

Spateo itself is installed from a **separate source checkout** or a compatible environment; this repository contains skills and alignment pipelines, not the complete Spateo library. The environment skill documents the required checkout and installation commands.

Data IO preserves the meaning of counts, annotations, spatial coordinates, units, and cell identifiers. It does not establish registration or silently manufacture an `X_pca`. Follow the alignment input contract and select a representation explicitly.

For expression input, [prepare shared PCA](skills/spateo-2d-alignment/references/expression-pca.md) jointly over all requested slices of **one biological specimen**. The preparer matches cells to the expression source by ID and writes new sanitized slices, `basis.npz`, and `expression_pca_manifest.json`. It never pools species or developmental stages. Both validators require the shared basis and matching per-slice feature, identity and file hashes; matching feature dimensions alone is insufficient. An explicit `--pca-provenance` can select the manifest.

Annotation is not required for expression PCA. In continuity-guided expression or spatial-only mode, annotation QC defaults to `off`, including when labels are present. Use `--annotation-qc provided --annotation-key KEY` only when those labels should participate in QC; missing labels then fail validation. Annotation one-hot mode requires provided labels. With QC off, annotation-dependent proposals are skipped and geometric continuity checks remain available. Physical z is preserved when supplied; otherwise unique numeric SL filenames define order without inventing spacing.

## Sources and verification

The environment skill and Data IO API review are based on [gmhhhhhh-929/spateo-release](https://github.com/gmhhhhhh-929/spateo-release/tree/d6aa68addc475dd0b56f69cebe7823b1f79933a9), commit `d6aa68addc475dd0b56f69cebe7823b1f79933a9`.

In that source tree, `spateo/data_io.py` is the AnnData compatibility entrypoint; maintained reader implementations are in `spateo/io/`. The Data IO skill follows those implementations and distinguishes automatic detection from explicit reader calls.

The current alignment extension passed synthetic contract checks, two full no-annotation expression runs, and two real annotation packaging comparisons with identical recovered XY coordinates. Expression runs establish functionality, not expression-mode accuracy. See [validation](VALIDATION.md) for the executed checks and their limits, and [licensing and source notices](NOTICE.md) for the status of imported material. The collection contains the first three workflow stages and distributes no biological datasets or reference coordinates.
