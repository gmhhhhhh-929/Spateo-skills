# Spateo Skills

本仓库整理空间转录组分析的前三个技能：环境配置、Data IO（转为AnnData）、2D切片配准。

[English](README.md)

| 阶段 | 技能入口 | 内容 |
| --- | --- | --- |
| 1 | [setup-spateo-environment](skills/setup-spateo-environment/SKILL.md) | 安装、诊断和核验独立Spateo环境，区分Spateo源码仓库与本Skills仓库。 |
| 2 | [spateo-data-io](skills/spateo-data-io/SKILL.md) | 根据实际源码选择读取接口，处理平台识别与歧义，验证AnnData并保留坐标、单位、注释和细胞ID。 |
| 3 | [spateo-2d-alignment](skills/spateo-2d-alignment/SKILL.md) | 两套配准pipeline、输入预检、可移植运行入口及源文件溯源。 |

每个技能都自包含在自己的目录中，以`SKILL.md`为入口；脚本与详细参考也位于该目录。可让代理直接读取入口，也可将单个技能目录复制到Codex的`~/.codex/skills/`后使用。

## 两套配准pipeline

两者共同属于第3个技能，目录按功能命名：

| 名称 | 功能 | 来源 |
| --- | --- | --- |
| `pairwise-rigid` | 相邻切片的Spateo刚性配准，可使用expression PCA、注释one-hot或spatial-only输入。 | 原v0.2.3实现。 |
| `continuity-guided` | 注释感知的连续切片配准，包含自动连续性检查及修复。 | 原v3.7.0实现。 |

原版本号保留在[迁移记录](skills/spateo-2d-alignment/provenance/source_migration.json)中，用于核对源文件与算法实现。记录同时列出移植后的哈希和为独立打包进行的修改。

## 使用顺序

1. 使用`$setup-spateo-environment`准备环境。Spateo库需要另行克隆或安装；不要在本Skills仓库根目录执行Spateo的`pip install -e .`。
2. 使用`$spateo-data-io`读取数据并检查AnnData。保留表达计数、注释、原始空间坐标、单位和细胞ID的原有含义。
3. 使用`$spateo-2d-alignment`检查切片输入，明确表征类型与切片顺序，再选择pipeline执行。Data IO阶段不等同于完成配准，也不会自动把注释one-hot认作expression PCA。

## 来源与验证

环境技能及Data IO接口审查基于[用户的Spateo源码仓库](https://github.com/gmhhhhhh-929/spateo-release/tree/d6aa68addc475dd0b56f69cebe7823b1f79933a9)，固定提交为`d6aa68addc475dd0b56f69cebe7823b1f79933a9`。

当前源码的`spateo/data_io.py`是AnnData兼容入口，维护中的读取实现位于`spateo/io/`。Data IO技能以这些实际实现为依据，并区分自动识别与显式调用。

已执行的检查及测试边界见[VALIDATION.md](VALIDATION.md)，导入代码的许可与来源见[NOTICE.md](NOTICE.md)。本轮先交付前三个阶段。
