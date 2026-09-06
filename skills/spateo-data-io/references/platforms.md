# 平台输入、观察单位与空间保存

以下均描述固定提交实际实现。完整签名见 [source-api.md](source-api.md)，只读对应平台即可。`X` 是否原始 counts 还取决于所给文件：STARmap 的 processed expression 等不能仅凭 reader 名称称为 counts。

## 10x Visium 与 HD

| 入口 | 必需布局 | 观察与坐标 |
| --- | --- | --- |
| `read_visium(path, ...)` | `filtered_feature_bc_matrix.h5`（可由 count_file 修改）；`spatial/tissue_positions.parquet`、`tissue_positions.csv` 或 headerless `tissue_positions_list.csv`。允许 pipestance 根目录并下钻 outs。 | spot×gene；按 barcode 左连接。`spatial=[pxl_col_in_fullres,pxl_row_in_fullres]`，全分辨率图像像素；缺失坐标发警告但保留 NaN，须验收阻断配准。不会自动筛除 `in_tissue=0`。 |
| `read_visium_hd_bin(path, ...)` | **直接指向单一** `binned_outputs/square_016um` 等目录；H5 或 `filtered_feature_bc_matrix/`；positions 默认 parquet，可回退同名 CSV。 | bin×gene；坐标依次尝试 fullres col/row、pxl col/row、x/y、array col/row，单位随选中的字段。`binsize=16` 是 metadata，既不选择子目录，也不重新分箱。 |
| `read_visium_hd_seg(path, ...)` | 指向 segmented_outputs；`filtered_feature_cell_matrix.h5` 和 GeoJSON，默认 `graphclust_annotated_cell_segmentations.geojson`，有三个常见别名回退。需要 geopandas/shapely。 | cell×gene；`cell_id` 转为 `cellid_{9位补零ID}-1`，或用现成 `cellid`；取矩阵/多边形交集。`spatial` 来自 polygon centroid，单位沿用几何坐标，WKT 在 `obs['geometry']`。 |
| `read_visium_hd(path, data_type='bin'|'cellseg', ...)` | 仅分派到上面两个 reader。 | 不会自动从 `outs` 下钻具体 HD 子目录；需要自己选路径或用 auto 返回的 path。 |

Visium 图片在 `uns['spatial'][library_id]['images']['hires'/'lowres']`，scale JSON 在 `scalefactors`，fullres XY 本身不随 image downsample 改写。`load_images=False` 仍读 positions/scalefactors；`load_qc_images=True` 才额外载入诊断图。source image 参数仅记录路径，不估计变换。

HD reader 自动读取指定 image paths，无统一 `load_images` 参数；可传明确不存在的可选图像路径以跳过数组并在外置 manifest 记录该选择。HD `uns['spatial'][sample]` 保存图像与 scale；bin 另有 `binsize`。直接 HD reader 的来源记录不如 auto 完整，外置 manifest 必须补齐。

`write_visium_hd_cellseg(adata,path,sample=None)` 写 H5、GeoJSON、图像、scale 到新目录，需 `obs['geometry']` 与 `obsm['spatial']`。它将 `X` 转 int32，**不是无损保存任意浮点表达的通用出口**。先确认非负整数 counts；不得拿 PCA、one-hot、log expression 导出为 feature-cell counts。可回读的 ID 应采用 `cellid_000000001-1` 等格式：其他 ID 在 GeoJSON 端被另编序号但 H5 barcode 保持原值，可能失去交集。若确需转换，先在派生副本中显式映射并保存原 ID；缺失/无效 geometry 会跳过相应 GeoJSON 记录。

## Xenium 与 Atera

两者核心为 `cell_feature_matrix.h5` + `cells.parquet`/`cells.csv.gz`/`cells.csv`，输出 cell×gene，按 ID 连接并可能丢弃无 metadata 的矩阵行；Atera 还先对 metadata ID 去重。默认只保留 Gene Expression。optional boundaries 转为 WKT `obs['geometry']`；Xenium 不隐式加载巨大的 transcripts 表。

`read_xenium` 的 `spatial` 使用 `x_centroid,y_centroid`，标准 Xenium 为微米；也接受 `CenterX_local_px,CenterY_local_px` 但**没有做像素到微米换算**，遇到该回退不能沿用“单位微米”的 docstring。实验 JSON 作为 `uns['spatial'][library_id]['metadata']`；像素大小缺失时用 0.2125 µm/pixel。图像金字塔读取有 `image_max_dim=4096`，hires scalefactor 为 `downsample/pixel_size`。坐标不被改成图像像素。`load_image=False`、`load_boundaries=False` 可减小载入。

`read_atera` 对本提交支持的 Atera **preview-xenium-v4** 布局读取，不保证最终产品格式。auto 需 Atera 标记／多染色形态图证据，避免误判普通 Xenium；直接调用无标记只发警告。支持 `load_image/image_key`、`load_boundaries/load_nucleus_boundaries`、cell groups companion；核边界在 `obs['nucleus_geometry']`。H&E 文件可选，已有 3×3 alignment CSV 只保存为 `scalefactors['he_affine']`，H&E 图与 downsample 另存；此读入**不估计配准、不将该 affine 应用到细胞坐标**。大分子表仅进文件清单。

两者 `cache_file` 若存在就直接读取，**没有验证源文件/选项变化**；不同加载选项不得复用同一路径的缓存。默认 `cache_file=None`，必要时以来源哈希和选项区分新缓存。

## MERFISH / Vizgen MERSCOPE

`read_merfish(path, *, counts_file, meta_file, load_images=True, z_layer=3, load_boundaries=True)` 必须给配对 cell_by_gene 与 cell_metadata 文件；不提供 public 分子表→细胞矩阵自动聚合接口。源码有 private `_aggregate_transcripts_to_adata`，public reader 不调用它。只有 detected_transcripts 或 `.vzg` 不能代替所需 cell table。

counts 按 cell×gene 读为 int32 CSR；metadata 按 ID 左连接，`center_x/center_y` 别名规范化后写 `spatial`，通常是平台微米坐标，仍应核实源字段。缺失行可能导致 NaN。匹配分组的 detected_transcripts 若存在会整表读取，再最多保留 2,000,000 行副本于 `uns['spatial'][sample]['transcripts']`；`load_images=False` **不会跳过这张分子表**。大表的内存峰值仍取决于全表。

图像路径 `image_files[channel][z]` 即使不载入数组也保留；`z_layer=None` 跳过数组。默认只选择 z=3（实际不存在时看 warning/选择逻辑）；z 是图像层，不自动变为细胞切片标签。`transform_micron_to_mosaic_pixel` 保存 vendor transform，不改写 `spatial`；边界可在 `polygons` 和 `boundary_files`，不能假设总有 `obs['geometry']`。图像/多边形/metadata 混合对象应实际试写 H5AD。

## seqFISH

`read_seqfish(path, *, counts_file, meta_file, load_images=True, load_labels=True)` 读取 cell×gene 与坐标 metadata，可识别一目录中的 section 配对；配错 section 报错。counts 规范为 float32 CSR，但中间可能 densify，预估 cells×genes 内存。metadata 有 ID 交集则重索引；**无交集但行数一致会按行配对**，必须事先确认该对应，否则拒绝使用此结果。`spatial` 仅在 metadata 规范化出 center_x/center_y 时产生；成功读取可能没有空间矩阵。

`uns['spatial'][region]` 含 images、image_files、labels/label_files；图像失败路径单独记入 metadata。`load_images=False` 不等于 `load_labels=False`，省内存时显式同时关闭。坐标单位不在 reader 中统一转换，来源未知就声明未知。

## NanoString CosMx / SMI

`read_nanostring(path, *, counts_file, meta_file, fov_file=None)` 需要表达与 metadata 的 cell ID 和 FOV；拼接 `cellID_FOV` 成观察 ID，取交集。局部 XY 在 `obsm['spatial']`；已有 global XY 在 **`obsm['spatial_fov']`**。多 FOV 拼成一张物理切片时不能直接叠加局部原点；选择全局坐标并保留局部坐标，或先获得可信 FOV 布局。可选 fov_file 仅记录 FOV metadata，**不自动从偏移重建缺失全局坐标**。

`CellComposite/` 与 `CellLabels/` 中匹配 FOV 的图片默认会读取，无 `load_images` 关闭参数。`uns['spatial'][fov]['images']['hires'/'segmentation']` 存图；可从 metadata 或标签图提取 geometry WKT（局部像素框架）。读取 image 需要 image 后端、进度依赖 tqdm；如要使用全局 geometry，需另外明确坐标变换，不能只改 spatial key。

## Slide-seq / Slide-seqV2

`read_slideseq(path, *, counts_file='MappedDGEForR.csv', bead_file='BeadLocationsForR.csv', load_images=True)`：counts 文件为 gene×bead，输出 bead×gene；bead 文件需 barcode 和 XY。只保留交集，有缺失发警告；输出 `obs['barcode','xcoord','ycoord','sample']`、`spatial=[x,y]`。metadata 将单位标记为 pixel，但非标准输入需独立核实，而不是强行换算。图片在 `uns['spatial'][sample]['images']`；用 `load_images=False` 跳过。

## STARmap PLUS

`read_starmap_plus(path, *, counts_file, meta_file, spatial_file, reorient_xy=False, dtype='float32')` 的三个文件参数在签名里都必需。spatial_file 可显式传 None 启用内部发现；meta_file 指定一个缺失文件时内部允许不合并，但不要把不存在路径误报为已读取的 metadata。最好精确指定 counts/spatial 配对，auto 可能优先选 processed 矩阵，必须声明其表达语义。

表达 reader 依据首列 gene/features 标记或矩阵尺寸启发式判断 gene×observation 并转置，否则按 observation×gene 读取；确认实际方向及基因列，不能仅靠启发式。非数值表达会被转为 NaN 后填 0，故需预查输入表。若坐标没有可识别 ID 但行数相同，reader 可按行对应，验收前需外部证明。匹配的 XY(/Z) 全部写入 `obsm['spatial']`，因此 **有 Z 时可能为 `(n_obs,3)`**，对应原列同时在 obs；到 2D 必须显式派生 XY 并保留原 XYZ，不能把第三列当表征。不会自动生成有生物意义的切片顺序。`reorient_xy=True` 将 `(x,y)` 改为 `(max(y)-y,max(x)-x)`，默认 False；这是确定性坐标重排，不是注册/配准。平台 metadata 在 `uns['spatial'][group]['metadata']`。

## Stereo-seq / BGI 与 Seq-Scope

BGI 输入为 TSV / gz GEM，必需 `geneID,x,y` 和 **一个** `MIDCount/MIDCounts/UMICount/UMICounts`。`read_bgi_as_dataframe` 规范 total、spliced(EXONIC)、unspliced(INTRONIC)、label。坐标读为 uint32、计数读为 uint16，负数、非整数或大于 65,535 的单条计数必须预先检查，累加也应核对总量以发现溢出。

- `read_bgi_agg`：`.X` 是 XY 网格上的总 UMI，`uns['__type']='AGG'`。可能包含 stain、labels、spliced/unspliced、gene_agg layers；obs/var 名是网格索引，不是 cell/gene。保存 origin/bin metadata，后续已有分割可供 `read_bgi` 使用。`prealigned=True` 只按最小 RNA 坐标给 stain padding，不估计图像变换。`binsize>1` 与 `gene_agg` 联用存在未正确分箱原坐标的实现风险，见 errors。
- `read_bgi`：严格四选一 `binsize`、`segmentation_adata`、`labels` 或 `label_column`；segmentation_adata 与 labels_layer 必须成对且 AGG；labels 可是 numpy 阵列或 `.npy`。输出 bin/label×gene。`label_column` 中 0 是背景并被丢弃；不生成新的分割。默认 `add_props=True` 从 bin/标签几何生成 centroid/bbox/contour；`add_props=False` 在此实现中可能**不产生 spatial**。掩码对齐、seg_binsize/origin 是输入契约，不能依靠形状相同推定正确。
- BGI 空间 metadata 使用 `uns['spatial']['binsize','scale','scale_unit']` 的平面命名空间；不像 Visium 的 library→images。`version='stereo'` 通常 scale=0.5 µm（可由 ngs_tools 提供），数值 XY 保持原坐标单位，**没有乘 scale**。转换单位时只乘一次并记录。
- `read_seqscope(matrix_dir, positions_path, binsize=1, add_props=True)`：MTX + features/genes + barcodes，以及无表头 `barcode lane tile x y` 空白分隔位置表。`binsize=None` 保留 barcode；正整数分箱聚合（默认 1 仍可合并同位置条码）。positions 重复 barcode 报错，无匹配报错，部分无坐标则警告并丢弃。输出 spatial 为坐标或 bin centroid；`add_props=False` 仍生成 centroid。scale=1、unit=None；须查物理标定。**Seq-Scope 未列入 auto detector**，直接调用。

## 手动图像附件

`st.io.spatial.add_image_layer(adata,img,scale_factor,slice='sample',img_layer='hires')` 或 `st.io.read_image(adata,filename,scale_factor,...)` 修改该对象的 `uns['spatial'][slice]['images'][img_layer]` 和 `scalefactors[img_layer]`；这里 scale 的键不是 Visium 的 `tissue_hires_scalef`。提供非空字符串 slice/img_layer，避免默认 None 键影响 H5AD 序列化。read_image 使用 OpenCV 的原始通道顺序（彩色常为 BGR）；不做 RGB 转换、图像配准或 XY 变换。
