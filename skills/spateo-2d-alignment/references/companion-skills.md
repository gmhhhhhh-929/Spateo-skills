# Alignment subskills / 配准子技能

Install or copy the entire `spateo-2d-alignment` directory. Start from its
`SKILL.md`, then read only the subskill needed for the current operation.
Nested entries are explicit workflow references; automatic recursive discovery
is not required. Run a subskill's relative commands from its own directory.

安装时保留整个 `spateo-2d-alignment` 目录；统一从主 `SKILL.md` 进入，
再按当前任务读取对应子技能。子技能的相对命令从各自目录执行。

| Skill / 技能 | Purpose / 用途 |
| --- | --- |
| [detect-spatial-components](../subskills/detect-spatial-components/SKILL.md) | Connected-component detection · 组织连通分量 |
| [lightweight-spatial-alignment-workflow](../subskills/lightweight-spatial-alignment-workflow/SKILL.md) | Workflow orchestration · 工作流 |
| [remote-workflow-intake](../subskills/remote-workflow-intake/SKILL.md) | Remote execution context · 远程任务环境 |
| [spateo-continuity-first-serial-alignment](../subskills/spateo-continuity-first-serial-alignment/SKILL.md) | Alias to continuity-guided · 连续性配准入口 |
| [spateo-pairwise-qc](../subskills/spateo-pairwise-qc/SKILL.md) | Pairwise numeric and visual QC · 成对质控 |
| [spateo-pairwise-run](../subskills/spateo-pairwise-run/SKILL.md) | Pairwise job planning · 成对任务计划 |
| [spateo-roi-refine](../subskills/spateo-roi-refine/SKILL.md) | ROI/drop refinement · 局部修正 |
| [spatial-alignment-compose](../subskills/spatial-alignment-compose/SKILL.md) | Recipe composition and replay · 变换组合与重放 |
| [spatial-balanced-sample](../subskills/spatial-balanced-sample/SKILL.md) | Component-balanced sampling · 平衡采样 |
| [spatial-before-after-viewer](../subskills/spatial-before-after-viewer/SKILL.md) | Before/after and displacement viewer · 配准前后及位移查看 |
| [spatial-component-align](../subskills/spatial-component-align/SKILL.md) | Component candidates · 分量配准候选 |
| [spatial-pointcloud-viewer](../subskills/spatial-pointcloud-viewer/SKILL.md) | Full-points 3D, slice/pair/component viewer · 全点三维查看 |
| [spatial-triad-component-align](../subskills/spatial-triad-component-align/SKILL.md) | Three-slice QC and candidates · 三切片检查 |
| [spatial-workflow-record](../subskills/spatial-workflow-record/SKILL.md) | Records, state pointers and dashboard · 记录与看板 |

The two viewer subskills bundle their own scripts and may optionally be copied
independently. All other workflow subskills should stay inside this package so
their shared pipeline paths remain valid. Moving directories changes no
alignment, rendering, or evaluation algorithms.
