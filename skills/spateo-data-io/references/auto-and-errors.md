# 自动识别、来源与错误处理

## 识别流程

`detect_spatial_technologies(path, technology=None, min_confidence=0.0)` 返回所有匹配的 `SpatialReadMatch`。每个包含 `technology,reader,path,kwargs,confidence,evidence`；`load_reader()` 加载 callable，`read(**overrides)` 用显式参数覆盖推断 kwargs 后调用。检测主要依赖文件名／目录结构，BGI 额外检查表头，**不验证 counts 文件内容**。缺失目录也可能表现为无候选，不能报告成成功读取。

`detect_spatial_technology(..., min_confidence=0.5, strict=True)` 没候选时抛 FileNotFoundError；strict 模式下最佳与其他候选 confidence 差值 ≤0.05 会抛 ValueError。候选按 confidence、证据数、路径倒序排序；confidence 是规则分数，不是统计概率。`strict=False` 直接取第一项，通常不适合一个目录内多个样本/分辨率。

`read_auto_spatial(..., return_match=False, **reader_kwargs)` 检测后读取，`return_match=True` 返回 `(adata, match)`。`read_spatial_auto` 是同一函数别名。识别支持 Atera、Xenium、Visium、HD bin/cellseg、Slide-seq、MERFISH、NanoString、seqFISH、STARmap PLUS、BGI；**不含 Seq-Scope**。技术名的常见 alias 包括 `visium_hd`（同时包含 bin/cellseg）、`cosmx/smi`（nanostring）、`merscope/vizgen`（merfish）、`stereoseq/stereo`（bgi）、`atera/wta`。未知 alias 抛 ValueError。

| 常见歧义 | 处理 |
| --- | --- |
| HD outs 同时有 cellseg 与 bin | 显式 `technology='visium_hd_cellseg'` 或 `'visium_hd_bin'`；若 bin 多分辨率，再选择精确 square_*um 路径。 |
| 多个 seqFISH section / MERFISH region | `detect_spatial_technologies` 列出所有组；逐候选的 path+kwargs 分别读取或显式 counts/meta 文件。仅 technology 不能消除同技术多组。 |
| 目录多个 BGI/GEM 文件 | 指向确切文件，核对 `binsize`。auto 推断 `binsize=1` 是空间 bin，不是细胞分割。 |
| 文件布局非标准，检测失败 | 选择实际存在的 direct reader，显式给文件名；不得改成随意最高分或猜数据来源。 |
| 一般 seqFISH/MERFISH counts/meta 文件存在但不匹配 | 核对文件分组与 ID。识别通过不替代数据连接验收。 |

对 BGI 使用 auto 且要按既有 label 读取时，必须覆盖推断的 `binsize=1`：`read_auto_spatial(gem, technology='bgi', binsize=None, label_column='cell_id')`。否则会违反四选一规则。

## Provenance 与缓存

auto 会调用 `record_spatial_io`，在 `uns['spateo_io']` 保存技术、绝对源路径、实际 reader、confidence、evidence、repr 形式参数与 manifest；Atera 标记 preview-xenium-v4，其他标记 stable 仅是代码中的标签，不保证所有版本/文件都通过测试。

`spatial_file_manifest(path,max_files=10000)` 递归文件清单用平行数组 `paths/sizes_bytes/roles` 和 `truncated`，路径相对 root；**不计算内容哈希，不复制图片**。一个文件输入只清点该文件。代码触及 max_files 即 truncated，不能把清单当完整哈希快照。direct reader 的记录结构不统一，使用外置 manifest 补充源码提交、重要输入 SHA256、参数、原始/输出 ID 集合、矩阵语义、坐标框架与单位。

不要重用选项不符的 Xenium/Atera cache_file。新输出用独立路径；保存前核验，保存后原样读回核验，外置 image/transcript 路径不会因为 H5AD 被搬走自动变成可用文件。

## 与该提交对应的陷阱

- **ID 丢弃或按行匹配：** Xenium/Atera/Slide-seq/HD cellseg/Seq-Scope 会取有匹配信息的子集；Visium/MERFISH 可保留无坐标行；seqFISH/STARmap 存在按行回退。记录丢失／重复／不匹配 ID 数量，不静默接受跨细胞配对。
- **单位不自动转换：** Xenium/Atera 的 local_px 回退、HD 的 array 坐标回退与一般平台标注都要实查来源。MERFISH vendor transform 和 Atera H&E affine 被保存，不是新估计的配准。
- **没有统一省内存开关：** Xenium/Atera 为 load_image（单数）；seqFISH 的 load_images 与 load_labels 分开；MERFISH 即使关闭图片仍可能整表读 transcripts；HD/NanoString 无 load_images 参数。按需估算图像、dense counts 与 transcript 内存，不能给无参数的 API 硬塞开关。
- **BGI uint16 计数：** 单条与聚合计数可能溢出，应与原始宽整数总量核验。本技能不修改源码 dtype。若不满足范围，报告限制并先提出可审查的转换/源码修复，不能把截断结果当 counts。
- **BGI `gene_agg` + binsize>1：** 当前 gene_agg 分支用原始 subset x/y 构建已分箱 shape，可能越界或错误定位。不要推荐这个组合；gene-specific 空间聚合应独立验证或先修复该处源码。
- **BGI 标签坐标：** seg_binsize>1 的 mask 扩展、origin 与实际 reads 坐标需额外小型已知几何测试；shape 警告不是通过。add_props=False 的 read_bgi 缺 spatial，不能用随机点补齐。
- **STARmap processed 表：** auto 优先 processed，不能称为原始 counts；必要时显式 counts_file 指 raw，仍检查数据语义。
- **通用读／保存边界：** `st.io.read_csv` 仅 kwargs；st.io.read 的 Path 表格及重复 sep 问题见 api-routing；`st.io.save` 是 pickle，不应作为跨语言 AnnData 文件。HD writer 会将数据截为 int32。
- **注释与实现不一致：** 以函数签名／函数体为准，包括 HD direct reader 不下钻 outs、STARmap 文件参数必填、read_h5ad 不重建 GeoDataFrame。不要复制不存在的旧 import 路径。

异常分流：缺文件→报告具体必需文件并重新定位；缺 optional backend→在现有环境中寻找该依赖或关闭该项允许的加载，不以虚假成功替代；格式/ID/坐标不符→保留原输入和失败诊断，输出清晰修正方案；OOME→停止大载入、使用平台支持的轻量选项或先准备有记录的子集，不把“重试成功”当完整数据保留证明。
