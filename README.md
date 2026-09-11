# Spateo Skills

Reusable agent skills for preparing a Spateo environment, reading data into AnnData, screening slice quality, reviewing evidence, and aligning serial 2D tissue slices.

[中文说明](README.zh-CN.md)

## Skills

| Stage | Skill | What it provides |
| --- | --- | --- |
| 1 | [setup-spateo-environment](skills/setup-spateo-environment/SKILL.md) | Installation guidance, environment diagnosis, and capability verification for a separate Spateo source checkout or an installed package. |
| 2 | [spateo-data-io](skills/spateo-data-io/SKILL.md) | Source-backed reader selection, spatial platform detection, AnnData validation, and coordinate/metadata preservation. |
| 3 | [spatial-slice-quality-qc](skills/spatial-slice-quality-qc/SKILL.md) | Four-domain QC, full review of both score bands and display-only rigid preregistration. |
| 4 | [spatial-slice-quality-viewer](skills/spatial-slice-quality-viewer/SKILL.md) | Current binary multi-panel reports, collections and traceable ROI. |
| 5 | [spateo-2d-alignment](skills/spateo-2d-alignment/SKILL.md) | Two alignment pipelines, shared expression-PCA preparation, optional annotation, strict input validation, and recorded implementation provenance. |

The top level exposes **five skill entrypoints**. All 14 alignment companion
entries live inside `spateo-2d-alignment/subskills/`, including the viewers,
QC, sampling, ROI and replay tools. See the
[alignment subskill catalog](skills/spateo-2d-alignment/references/companion-skills.md).

```text
skills/
├── setup-spateo-environment/
├── spateo-data-io/
├── spatial-slice-quality-qc/
├── spatial-slice-quality-viewer/
└── spateo-2d-alignment/
    ├── SKILL.md
    ├── pipelines/
    │   ├── pairwise-rigid/
    │   └── continuity-guided/
    ├── subskills/                 # 14 alignment workflow subskills
    │   ├── spatial-before-after-viewer/
    │   ├── spatial-pointcloud-viewer/
    │   └── …
    ├── scripts/
    ├── references/
    └── provenance/
```

Install the complete alignment directory; its main `SKILL.md` routes to the
nested workflows. No zebrafish-stage updates are included. See the
[completeness audit](COMPLETENESS.md).

## Alignment pipelines

Both pipelines belong to the **2D alignment skill**. Both accept verified shared expression PCA, annotation one-hot, or spatial-only input.

| Pipeline | Use | Source snapshot |
| --- | --- | --- |
| `pairwise-rigid` | Adjacent rigid Spateo alignment. | Previously identified as v0.2.3. |
| `continuity-guided` | Serial Spateo alignment with automatic continuity checks and optional repairs. | Derived from the v3.7.0 snapshot, with explicit representation and profile controls. |

For continuity-guided, **`--profile legacy` remains the default**. The optional `--profile generalized` disables nearest-neighbor initialization and replaces the named-tissue priority with anonymous supported-label proposals. Frozen validation found average Drosophila improvements, Planarian regressions and severe individual failures. It does not establish a generally better profile. See [accuracy and validation limits](skills/spateo-2d-alignment/references/validation.md).

Historical version identifiers, source hashes and subsequent implementation changes are recorded in the [migration manifest](skills/spateo-2d-alignment/provenance/source_migration.json).

## Use

Clone this repository and let your agent read the relevant `SKILL.md`. Copy the desired top-level skill directory into your configured skill directory. For alignment, copy all of `spateo-2d-alignment/`, including `subskills/`; use its main entrypoint to access the nested workflows.

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

The current alignment extension passed synthetic contract checks, two full no-annotation expression runs, and two real annotation packaging comparisons with identical recovered XY coordinates. Expression runs establish functionality, not expression-mode accuracy. See [validation](VALIDATION.md) for the executed checks and their limits, and [licensing and source notices](NOTICE.md) for the status of imported material. The collection contains five top-level skills and distributes no biological datasets or reference coordinates.

## Spateo Referee

Install both complete QC and viewer directories as siblings. Joint-review v2 checks both review bands, with maximum-domain requirements 0.65 for the lower band and 0.60 for the higher band, retaining the other safeguards. It has metric-stress evidence only: use `experimental_policy`, not biological certification. See [methods, thresholds and validation](skills/spatial-slice-quality-qc/references/methods.md) and the [complete workflow SVG](skills/spatial-slice-quality-qc/references/workflow.svg). The package retains current runtime dependencies; historical experiment launchers and biological data are not distributed.
