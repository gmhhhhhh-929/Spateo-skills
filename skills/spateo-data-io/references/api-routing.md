# API 路由与通用输入

权威是固定提交中的函数实现与导出表，而不是装饰器 `examples/related` 中的历史 `ov.*` 名称。所有下述 `st` 指 `import spateo as st`。

| 目标 | 本提交入口 | 输入与返回／重要差异 |
| --- | --- | --- |
| 通用读取 | `st.io.read(path, backend='python', **kwargs)` / `st.io.single.read` | `.h5ad` → AnnData；`.csv/.tsv/.txt/.gz` → pandas DataFrame；不会猜成空间 AnnData。其他后缀报 ValueError。 |
| H5AD | `st.io.read_h5ad(filename, **kwargs)` / `st.io.single.read_h5ad` | 原生 anndata 读取后设置 Spateo `uns['__type']='UMI'` 和 `uns['pp']`；因此原本 AGG 对象不宜用它重新判定类型。要保留原始类型用 `anndata.read_h5ad` 或 `st.io.read` 的 Python 路径。 |
| 10x H5 | `st.io.read_10x_h5(filename, *, genome=None, gex_only=True, backup_url=None)` | v3 `/matrix` 或 legacy genome group；磁盘 features×barcodes 转为 obs×var；保留 `gene_ids/feature_types` 等。`gex_only=True` 丢弃非 Gene Expression 特征。legacy 多 genome 必须选择。`backup_url` 仅保留参数，无下载实现。 |
| 10x MTX | `st.io.read_10x_mtx(path, *, var_names='gene_symbols', make_unique=True, gex_only=True, prefix=None, compressed=True)` | `matrix.mtx[.gz]` + `features.tsv[.gz]` / legacy `genes.tsv` + `barcodes.tsv[.gz]`。默认 v3 gzip；非压缩 v3 显式 `compressed=False`。基因符号去重时保留 `gene_ids`，后续跨片特征交集最好使用稳定基因 ID。 |
| CSV 工具 | `st.io.read_csv(**kwargs)` / `st.io.general.read_csv` | 仅关键字参数，例如 `filepath_or_buffer='metadata.csv', index_col=0`；返回 DataFrame，不是 AnnData。 |
| Python 对象 | `st.io.save(file, path)` / `st.io.load(path, backend=None)` | pickle/cloudpickle；不是 H5AD。仅加载可信 pickle。`save` 需要非空父目录，如 `./outputs/object.pkl`；裸文件名会触发 `os.makedirs('')`。 |
| AnnData 兼容入口 | `st.AnnData`, `st.concat`, `st.read_h5ad`, `st.read_csv` 等来自 `spateo/data_io.py` | 这些是 anndata re-export，与 `st.io.read_csv` 的 DataFrame 工具不同。`st.read` 在 anndata≥0.12 回退为 `read_h5ad`，不能作为跨后缀通用读取器。 |

`st.io.read` 的 Rust 分支依赖 `snapatac2`，返回后需要关闭文件且并非所有 Spateo 功能兼容；普通 AnnData 交接优先 Python。当前表格分支调用 `path.endswith`，实际传入 `Path` 的部分后缀可能报 AttributeError；转换为 `str(path)`。该分支自己传 `sep`，不要再用 kwargs 重复传入；需自定义分隔符时用 `st.io.read_csv(filepath_or_buffer=..., sep=...)`。

平台函数均从 `st.io` 和 `st.io.spatial` 导出。仅 `st.io.spatial` 额外公开 `read_bgi_as_dataframe`、`add_image_layer`、`bin_indices/centroids`、`get_*_props`、`bin_matrix/get_coords_labels` 等辅助函数；`alpha_shape/get_concave_hull` 顶层也导出。它们是分箱、形状或图像工具，不是配准算法。不存在 public `read_seqfish_plus`，即使某些装饰器 related 列表仍提到该名字。

普通细胞表格需显式构建：先按 ID 对齐 counts、metadata、XY，确认矩阵方向，使用 `AnnData(X=..., obs=..., var=...)`，再赋值 `obsm['spatial']`。同样长度不代表相同行序；连接后报告重复 ID、缺失 ID 和被删除行。不要给通用 `read_csv` 增加不存在的自动坐标识别能力。

保存使用 `adata.write_h5ad(new_path)`。普通 `st.io.single.read_h5ad` 不重建 GeoDataFrame；HD 源码注释中该说法不符合此函数实际实现，WKT 留在 `obs['geometry']`。如需 GeoDataFrame，应在使用时显式解析，保留 WKT 持久化字段。
