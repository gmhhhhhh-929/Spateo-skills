# Spateo Skills

本仓库整理空间转录组分析的前三个技能：环境配置、Data IO（转为 AnnData）、2D 切片配准。

[English](README.md)

| 阶段 | 技能入口 | 内容 |
| --- | --- | --- |
| 1 | [setup-spateo-environment](skills/setup-spateo-environment/SKILL.md) | 安装、诊断和核验独立 Spateo 环境，区分 Spateo 源码仓库与本 Skills 仓库。 |
| 2 | [spateo-data-io](skills/spateo-data-io/SKILL.md) | 根据实际源码选择读取接口，处理平台识别与歧义，验证 AnnData 并保留坐标、单位、注释和细胞 ID。 |
| 3 | [spateo-2d-alignment](skills/spateo-2d-alignment/SKILL.md) | 两套配准 pipeline、共享表达 PCA 准备、可选注释、严格输入预检及实现溯源。 |

每个技能都自包含在自己的目录中，以 `SKILL.md` 为入口；脚本与详细参考也位于该目录。可让代理直接读取入口，也可将单个技能目录复制到 Codex 的 `~/.codex/skills/` 后使用。

## 两套配准 pipeline

两者共同属于第 3 个技能，均支持经过审核的共享 expression PCA、注释 one-hot 和 spatial-only 输入。

| 名称 | 功能 | 来源 |
| --- | --- | --- |
| `pairwise-rigid` | 相邻切片的 Spateo 刚性配准。 | 原 v0.2.3 实现。 |
| `continuity-guided` | 连续切片配准，包含自动连续性检查及可选修复。 | 基于原 v3.7.0 快照，增加显式表征类型与 profile 控制。 |

continuity-guided 的兼容默认值仍是 **`--profile legacy`**。显式选择 `--profile generalized` 可关闭最近邻初始化，并用匿名、有支持度的标签候选替代固定组织名称优先级。冻结验证中，果蝇平均指标有所提升，涡虫平均指标下降，部分样本出现明显失败；这些结果不能证明 generalized 普遍更准确或更稳定。详见[准确率与验证边界](skills/spateo-2d-alignment/references/validation.md)。

原版本号、源文件哈希及后续实现修改保留在[迁移记录](skills/spateo-2d-alignment/provenance/source_migration.json)中。

## 使用顺序

1. 使用 `$setup-spateo-environment` 准备环境。Spateo 库需要另行克隆或安装；不要在本 Skills 仓库根目录执行 Spateo 的 `pip install -e .`。
2. 使用 `$spateo-data-io` 读取数据并检查 AnnData，保留表达计数、注释、空间坐标、单位和细胞 ID 的原有含义。
3. 使用 `$spateo-2d-alignment` 检查切片输入，明确表征类型与切片顺序，再选择 pipeline 执行。Data IO 不等同于配准，也不会把注释 one-hot 自动认作 expression PCA。

## 无注释的表达 PCA 输入

准备器按细胞 ID 对应原始表达矩阵，为**同一个生物样本的全部待配准切片**联合计算一次 PCA。它不混合不同物种或发育时期，也不分别为每张切片拟合不同基底。以下命令从仓库根目录运行，示例中的输入路径需换成实际文件，输出目录必须是新目录：

```bash
python skills/spateo-2d-alignment/scripts/prepare_expression_pca.py \
  --slice-dir inputs/specimen_slices \
  --expression-source inputs/specimen_expression.h5ad \
  --matrix-state counts --n-hvg 2000 --n-pcs 50 \
  --output-dir runs/specimen_pca

python skills/spateo-2d-alignment/pipelines/continuity-guided/run.py \
  --stage specimen --slice-dir runs/specimen_pca/slice_h5ad \
  --output-dir runs/specimen_alignment \
  --representation expression-pca --annotation-qc off --profile generalized --device 0
```

该示例显式使用已做完整功能验证的 generalized；它仍是可选方案，不能据此推断准确率优于 legacy。准备器默认不读取或复制注释。输出包括新的 `slice_h5ad/`、`basis.npz` 和 `expression_pca_manifest.json`；两套 pipeline 都审核共同基底、逐片特征、细胞身份及文件哈希，不能仅凭 PCA 列数相同认定输入合格。可用 `--pca-provenance` 显式指定 manifest。

continuity-guided 的 expression-pca 和 spatial-only 模式默认 `--annotation-qc off`，即使存在注释也不自动用于 QC。只有明确需要可靠标签参与检查时，才使用 `--annotation-qc provided --annotation-key KEY`；此时缺失注释会报错。annotation-onehot 模式必须提供注释。关闭注释 QC 后，依赖标签的修复会跳过，几何连续性检查仍可运行。

已有物理 z 会保留；全部切片均无 z 时，使用文件名中唯一的数字 SL 编号确定顺序，不虚构物理间距。表达预处理参数、矩阵状态声明和无 z 输入要求见[共享表达 PCA 说明](skills/spateo-2d-alignment/references/expression-pca.md)。

## 来源与验证

环境技能及 Data IO 接口审查基于[Spateo 源码仓库](https://github.com/gmhhhhhh-929/spateo-release/tree/d6aa68addc475dd0b56f69cebe7823b1f79933a9)，固定提交为 `d6aa68addc475dd0b56f69cebe7823b1f79933a9`。

该源码的 `spateo/data_io.py` 是 AnnData 兼容入口，维护中的读取实现位于 `spateo/io/`。Data IO 技能以这些实际实现为依据，并区分自动识别与显式调用。

本次配准扩展通过合成契约测试、两个真实样本的完整无注释 expression 运行，以及两个真实 annotation 样本的打包等价性验证；后者的逐细胞输出 XY 与冻结优化结果完全相同。表达模式验证的是可运行性，不代表已完成其 A2/A5 准确率评估。初次环境与 Data IO 验证仍作为已有证据，本次未重做这两个阶段。具体检查与边界见 [VALIDATION.md](VALIDATION.md)，许可与来源见 [NOTICE.md](NOTICE.md)。仓库不分发生物数据或参考坐标。
