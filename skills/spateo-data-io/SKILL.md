---
name: spateo-data-io
description: 将平台原始输出、10x 矩阵或 H5AD 正确读成 Spateo AnnData，审查自动识别、细胞身份、坐标单位、图像及来源，并准备到 2D 配准的输入契约。适用于 Data IO 和数据转换；不执行分割、生成生物注释或估计配准变换。
---

# Spateo Data IO → AnnData

本技能依据用户源码 `gmhhhhhh-929/spateo-release` 的提交 **`d6aa68addc475dd0b56f69cebe7823b1f79933a9`** 编写。维护实现位于 `spateo/io/{general,single,spatial,spatial/auto}`；`spateo/data_io.py` 只兼容导出 AnnData 的常用类和读取函数。不要把旧版 `spateo.io.bgi`、`spateo.io.tenx` 或注册装饰器里的历史示例当作本版本 API。精确签名与源码摘要见 [source-api.md](references/source-api.md) 和 [source_manifest.json](references/source_manifest.json)。使用其他提交时先复核相应 reader 签名与行为。

## 选择输入路径

先确定一个数据集对应哪一张物理切片／区域、观察单位是细胞、spot、bead、bin 还是分子，原始文件如何按 ID 对应。保留原始文件；转换输出写到独立新目录。记录表达矩阵语义、坐标列与单位；不知道单位时标注未知，不能因为列名叫 `spatial` 就当作微米或已经配准。

- H5AD、10x H5/MTX、普通表格：读 [api-routing.md](references/api-routing.md) 的通用路由。`st.io.read()` 读表格返回 DataFrame；10x 表达矩阵不会自动产生空间坐标。
- 已知空间平台：读 [platforms.md](references/platforms.md) 中对应平台，用显式 reader 和输入文件。MERFISH、seqFISH、STARmap PLUS 的 `counts_file` 等必需参数不能依据过时示例省略。
- 未知目录：先执行 `st.io.detect_spatial_technologies(path)`，检查候选的 `technology/path/kwargs/evidence/confidence`，再读。[auto-and-errors.md](references/auto-and-errors.md) 说明歧义、内存开关、缓存与错误处理。一个目录有多个区域、HD 分辨率或技术时，逐项选择；`strict=False` 只是取排名第一，不证明选对数据。
- BGI 分子／read 表：先决定 `read_bgi_agg` 的 XY 网格或 `read_bgi` 的 bin／既有标签模式。网格 `.X` 的两个轴不是细胞与基因，不能直接交给表达 QC/PCA/2D 配准。此 reader 不替用户做细胞分割。

## 读取、验收、保存

在现有的支持环境中运行真实源码 API；此源码主包支持 Python **3.10–3.12**。按需检查 `anndata/numpy/scipy/pandas/h5py` 与目标 reader 的可选依赖，不默认安装完整 GPU 环境。不要将本技能脚本当成安装器。

```python
import spateo as st

# 先查看候选，显式选择单一数据集。
matches = st.io.detect_spatial_technologies("/data/sample/outs")
for m in matches:
    print(m.technology, m.path, dict(m.kwargs), m.confidence, m.evidence)

adata, match = st.io.read_auto_spatial(
    "/data/sample/outs", technology="visium", load_images=False,
    strict=True, return_match=True,
)
# 写之前核查 below 的身份、单位、有限性与矩阵语义。
adata.write_h5ad("/output/new-run/sample.h5ad")
```

验收必须包含：矩阵方向与来源、`obs_names` 唯一性、基因身份、读入前后 ID 覆盖、观察单位、所选 XY 形状与有限性、切片/FOV 身份、单位和坐标框架、图像与变换 metadata。Reader 可能丢弃无法匹配的细胞、保留 NaN、或在缺少 ID 时按行对应；不能以“成功返回 AnnData”代替核验。精确存储位置与平台例外见 [platforms.md](references/platforms.md)。

可用辅助脚本：

```bash
python scripts/spateo_io.py detect /data/sample/outs --technology visium
python scripts/spateo_io.py convert /data/sample/outs /output/new-run/sample.h5ad \
  --reader read_visium --kwargs-json /output/reader_kwargs.json
python scripts/validate_anndata.py /output/new-run/sample.h5ad \
  --matrix-semantics counts --spatial-key spatial --coordinate-unit fullres_pixel
```

CLI 是本技能提供的薄封装：参数直接传给实际 API、拒绝覆盖已有输出，不重排数据或修补数值；验证器只读输入。图像加载开关因平台而异，不能统一传 `load_images=False`。输入 SHA、版本、reader 参数与结果摘要保存在 CLI 的外置 manifest；原生 `spatial_file_manifest` 是路径／字节数清单，**不含内容哈希**。H5AD 保存后再读回并核验坐标、ID、矩阵及需保留的 images/metadata；遇到不可序列化内容时保留原值并明确外置，不静默删除。

## 向 2D 配准交接

读 [alignment-contract.md](references/alignment-contract.md)。交接的是按正确 ID 对应的细胞/spot/bin × 基因 AnnData、可解释的原始 XY、切片顺序，以及下游所选模式的表征。Data IO 不生成配准结果，也不证明样本已在公共坐标系。

- 基因 counts 保留在原始对象或明确 `layers['counts']`；表达归一化与 PCA 是后续派生步骤，不能把 one-hot 或 PCA 写成“原始基因 counts”。
- expression PCA 必须从原始表达计算并保留共同特征、共同 PCA 基底及拟合来源；仅凭数组形状无法证明来源。
- annotation one-hot 需要真实注释和所有相关切片共享的类别顺序；缺失标签不可静默作为普通类别。不要自行推断注释。
- spatial-only 的非零常数表征属于下游明确请求的模式，不是 Data IO 的默认产物。

辅助验证例：`python scripts/validate_anndata.py slice.h5ad --matrix-semantics counts --coordinate-unit um --alignment-mode annotation-onehot --annotation-key anno --categories-json categories.json`，类别文件是跨切片共享的有序字符串列表。表达 PCA 模式会报告“来源尚需人工/manifest 核验”，不会把数值检查当成无泄漏证明。

## 验证与限制

真实源码的可重复合成测试：`PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/path/to/source python scripts/smoke_source.py --source-root /path/to/source`。它检查实际 public API 的 10x/Visium/Atera/BGI/Seq-Scope、自动识别歧义、保存读回和验证器拒绝错误输入；不会 mock reader 或执行配准。测试范围与已知源码限制见 [validation.md](references/validation.md)。不宣称未测试平台已通过真实大数据验证。
