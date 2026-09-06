# 从 Data IO 到 2D 配准的输入契约

这是本技能的交接约定，不是宣称所有 Spateo reader 会自动生成这些字段。先保留 reader 的原始 AnnData，再在独立派生对象上准备下游输入；坐标、图像和表征变换需可追溯。

| 项目 | 交接要求 |
| --- | --- |
| 观察单位 | 明确 cell/spot/bead/bin；`.X` 的行必须和这些观察一一对应，列为可识别基因。`__type='AGG'` 的 XY 图像矩阵不符合。 |
| 身份 | 原始 ID 唯一且无缺失。跨样本重复 barcode 用明确 sample+barcode 键并另存原 ID；记录所有过滤、去重、交集丢弃。 |
| 坐标 | 指定 obsm key，形状 `(n_obs,2)`，数值有限，XY 轴、单位、原点及镜像方向有说明。保留原坐标；已有 XYZ 不直接丢 Z，应显式派生 XY 并保留 Z。 |
| 切片 | `slice_id` 与数值顺序／真实 z 由实验设计或可信 metadata 提供，不能按字符串词典序或像素图层猜生物顺序。每片一个物理平面，避免混入局部原点重复的 FOV。 |
| counts | 稀疏非负整数基因计数可在原 X 或 layers['counts']；若源文件是 normalized/processed，明确其语义，不伪造 counts。注意 float dtype 可以存整数 counts，dtype 本身不能证明已归一化。 |
| 图像与几何 | 原图、scale、已给定 affine/WKT 继续保存；若坐标被转换，记录映射，不假设图像自动同步。 |
| 来源 | source commit/包版本、reader+参数、核心文件 SHA256、ID 集合摘要、坐标单位／框架、输出 SHA256。文件清单不是哈希。 |

### 交给本次发布的 2D runner

完整 canonical AnnData 继续保存图像、几何、额外坐标和 metadata；另写最小配准副本，`obsm` 只保留 `spatial` 与 `X_pca`，所有切片的 `obs_names` 全局唯一。pairwise-rigid 的原 runner 文件名严格为 `CS<stage数字>_SL<slice数字>_Y<word chip>.Spatial.h5ad`，例如 `CS01_SL001_YSAMPLE.Spatial.h5ad`；一个独立样本同一 stage，SL 唯一。其配准输入目录应为平面目录，不能混入嵌套 H5AD 或旧版副本，以免 runner 递归发现未审查文件。continuity-guided 接受文件名含 `SL<数字>` 的 `*.h5ad`。二者按 SL 数字排序，而不直接消费 `obs['slice_id']`、`z` 或 `physical_z` 来排序；在导出时用可信物理次序编码 SL，保留 z/physical_z。该 2D preflight 若检查 z/physical_z，则要求所有切片均有该 metadata 且随 SL **严格递增**；不改原坐标数值来强行满足排序。

pairwise-rigid 可用共享 expression PCA 或共享 annotation one-hot；continuity-guided 的标准输入是共享字典的 annotation one-hot 加 `obs[annotation_key]`。该特定 spatial-only 模式采用每细胞相同的 `ones(30)`。这些是下游 runner 约定，不能反过来裁掉 canonical IO 文件的 images、其他 obsm 或原始表达。

## 表征选择

**Expression PCA。** 先确认真实基因表达来源，保留 counts，按下游所需归一化/HVG/scaling 方案，在比较切片共享特征集合和同一个 PCA 基底内生成 `obsm['X_pca']` 或下游指定 key。独立每片 PCA 的列不能直接作为共同特征比较。记录处理步骤、基因顺序、components、随机种子与 fit/apply 范围；由 counts 计算 PCA 是后续派生，IO reader 并不会自动完成。数值上检查每行对应、有限、非全常数、非 annotation one-hot；本次 Spateo 2D cosine 输入还要求每行范数非零，验证器在 expression-pca 模式下据此拒绝零向量。欧氏距离单独使用 PCA 时零向量不必然非法，应另按实际算法契约验收。**任何数值检查都不能单独证明 PCA 来自表达**，必须核对来源 manifest。不要读取“正确配准坐标”来训练或选择表达表示。

**Annotation one-hot。** 只有明确选择注释模式时使用。确认指定 obs 列无缺失，定义所有相关切片共享且可复用的类别列表；存储列顺序。每行只有一个 1，其余为 0，且 argmax 反解严格等于该行注释；某片未出现的类别仍保留空列。不能把 one-hot 写成基因 X 或称为 expression PCA。可将它写到需要该布局的派生数据中的 X_pca，但 manifest 必须写 `annotation_onehot`。

**Spatial-only。** 只有下游明确要求时，在派生对象中生成每行相同、有限、非零常数向量，并检查全行相同；零向量会令某些 cosine 实现未定义。保留 expression 与 annotation 的原始信息，但算法输入只使用请求的空间模式。IO 不默认替换原 X_pca。

## 配准前执行检查

```bash
# 原始读入对象：只读结构、基因矩阵和空间检查。
python scripts/validate_anndata.py sample.h5ad \
  --matrix-semantics counts --spatial-key spatial --coordinate-unit um

# expression PCA：这里的 pass 仅表示数值结构合格，来源仍要查 manifest。
python scripts/validate_anndata.py slice.h5ad \
  --matrix-semantics counts --coordinate-unit um --alignment-mode expression-pca

# categories.json 是共享的有序类别字符串列表，不是每片重新排序。
python scripts/validate_anndata.py slice.h5ad \
  --matrix-semantics counts --coordinate-unit um \
  --alignment-mode annotation-onehot --annotation-key anno \
  --categories-json categories.json
```

验证器不自动修改 ID、坐标、X、layers 或 X_pca。普通 IO 模式接受有限 XY 或 XYZ；`--require-2d` 和任何 alignment-mode 才强制 XY。若基因矩阵在 layers，使用 `--matrix-key counts`；默认检查 X。没有可信类别顺序时 one-hot 验证会失败，要求补齐证据，而不是猜列含义。`--alignment-mode spatial-only` 检查已有非零常数表征，不在原数据生成它。

验收失败时，描述具体不满足项；可以在新的派生数据中修正已授权的问题，重新验证后交付。坐标分箱/centroid、图像加载或 vendor affine 的保存都不是完成了配准；配准由下游方法另外估计、冻结并评价。
