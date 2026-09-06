# 源码审查与可重复验证

## 依据与覆盖

本技能依据源码提交 `d6aa68addc475dd0b56f69cebe7823b1f79933a9`；[source_manifest.json](source_manifest.json) 列出所读 `spateo/data_io.py`、`spateo/io` 和 `tests/io` 共 32 个 Python 文件的 SHA256。实现签名从 AST 提取到 [source-api.md](source-api.md)。未使用旧 organization 技能作为 reader 权威，未复制大型原始数据或改动源码。

## 实际测试环境

此次审查使用本机已存在的 Python **3.9.6** 科学栈，真实导入该 checkout 的 `spateo.io`，没有 reader mock、没有替代实现、没有安装环境。这个结果只证明所测 **局部 IO** 在该现成环境运行成功，**不表示推荐或支持用 Python 3.9 安装 Spateo**。源码 `setup.py` 正式要求 `>=3.10,<3.13`，面向用户的安装与常规执行应选 Python 3.10–3.12。

包版本：anndata 0.10.9、NumPy 2.0.2、SciPy 1.13.1、pandas 2.3.3、h5py 3.14.0、Matplotlib 3.9.4、Pillow 11.3.0。机器化合成结果在 [smoke_result.json](smoke_result.json)，记录实际 Python、提交、版本和每项检查名。

## 重复合成测试

在技能目录中执行（`SOURCE` 指用户选定的、上述提交的源码 checkout）：

```bash
PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg PYTHONPATH="$SOURCE" \
  python scripts/smoke_source.py --source-root "$SOURCE"
```

脚本使用临时目录，完成后清理自己的合成输入，读写均不涉及原始实验数据。它核验源码 commit 和实际 import 路径，再调用 public reader；JSON 的 `status` 和 `n_tests` 是当前复核结果。覆盖：

- 10x H5 的 gene×cell→cell×gene 方向与 barcode 顺序，表达文件不自动生成坐标。
- Visium 按 ID 重排 metadata、XY 的 col/row 顺序、图片/scale/auto provenance、in_tissue 的原样保留及 H5AD 读回。
- 非有限坐标失败；通用 IO 保留 XYZ，而显式 2D 拒绝三列。
- one-hot 不冒充 expression PCA、共享类别顺序校验、PCA 零范数与 spatial-only 零向量拒绝。
- Atera preview 的真实检测、坐标身份与 scale 不重复应用。
- BGI AGG 网格与 bin/label×gene 区别、计数守恒、centroid、四选一参数、add_props=False 不生成坐标。
- Seq-Scope barcode 原样模式与 bin 聚合的身份/计数守恒。
- HD bin/cellseg 自动识别歧义及显式模式的路径选择。
- 通用 CSV 返回 DataFrame；CLI detect、真实转换、只读验证、输出拒绝覆盖及 XYZ 不截断的往返。

## 原始源码测试

此次还直接运行了以下未修改源码测试：

```bash
PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg PYTHONPATH="$SOURCE" \
  python -m pytest -q -p no:cacheprovider \
  "$SOURCE/tests/io/test_spatial_auto.py" \
  "$SOURCE/tests/io/test_io_layout.py" \
  "$SOURCE/tests/io/test_atera_visium.py"
```

结果：**13 passed**；14 条 Matplotlib/pyparsing 弃用警告，没有测试失败。`test_bgi.py` 依赖真实 fixture，`test_utils.py` 使用同一 fixture mixin，未冒充已运行；其源码已读。`test_image.py` 本提交为空文件。这不是整个 Spateo 测试套件通过的声明。

`skill-creator/scripts/quick_validate.py` 已通过本技能 frontmatter、名称和 scaffold 检查；三个脚本还进行了编译检查和真实 CLI 调用。该结构检查不替代科学输入契约验证。

## 验证边界

MERFISH、seqFISH、CosMx、Slide-seq、STARmap、HD cell segmentation 的大型真实数据／所有版本组合尚未由本次合成测试证明。相应路由来自固定源码静态审查，读入后仍要核对文件格式、ID 覆盖、内存和单位。source 内有已知边界（过时 docstring、BGI uint16/gene_agg 分箱、按行匹配、cache 失效、HD writer ID 规则），见 [auto-and-errors.md](auto-and-errors.md)。

本技能的验证器只报告数值与结构，不证明 PCA 的来源或生物注释正确性，也不进行分割、配准、准确率评价或读取参考配准坐标。CLI 在 manifest 中明确 `source_commit_verified_by_cli=False`：它记录实际导入模块文件与 SHA256，但只有合成测试／独立 provenance 检查验证指定 commit；用户换环境后应重新核对。目录输入只提供文件清单，核心输入内容 SHA256 需根据选定文件单独补充，不能把目录清单当完整内容哈希。
