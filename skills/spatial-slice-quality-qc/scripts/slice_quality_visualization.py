"""Detailed, lightweight reporting for binary serial-slice QC results.

This module consumes the bounded display payload written by
``spateo.preprocessing.slice_quality.write_slice_quality_outputs``.  It never modifies the source
H5AD.  When the source H5AD is available on the compute host, the report can
read point-level count, detected-gene, and coordinate columns to add absolute
expression-capture and connected-component evidence.  It also adds fixed-
bandwidth two-dimensional KDE values, directional local-window explanations,
and publication-oriented summary tables to the final keep/exclude calls.
"""

from __future__ import annotations

import hashlib
import html
import json
import math
import re
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Union

import numpy as np
import pandas as pd
from scipy import ndimage, sparse
from scipy.sparse import csgraph
from scipy.spatial import cKDTree

_EPS = np.finfo(float).eps


_EVIDENCE_DEFINITIONS: tuple[dict[str, str], ...] = (
    {
        "key": "cell_density",
        "label": "细胞/位置密度",
        "label_en": "Cell/spot density",
        "domain": "density",
        "bad": "low",
    },
    {
        "key": "n_locations",
        "label": "细胞/位置数量",
        "label_en": "Number of cells/spots",
        "domain": "density",
        "bad": "low",
    },
    {
        "key": "largest_component_fraction",
        "label": "最大连通组织占比",
        "label_en": "Largest connected-tissue fraction",
        "domain": "density",
        "bad": "low",
    },
    {
        "key": "fragmentation",
        "label": "组织碎片化",
        "label_en": "Tissue fragmentation",
        "domain": "density",
        "bad": "high",
    },
    {
        "key": "hole_fraction",
        "label": "组织空洞比例",
        "label_en": "Tissue-hole fraction",
        "domain": "density",
        "bad": "high",
    },
    {
        "key": "knn_tail_ratio",
        "label": "稀疏邻域长尾",
        "label_en": "Sparse-neighborhood tail ratio",
        "domain": "density",
        "bad": "high",
    },
    {
        "key": "median_total_counts",
        "label": "每点捕获计数中位数",
        "label_en": "Median captured counts per spot",
        "domain": "expression",
        "bad": "low",
    },
    {
        "key": "median_n_genes",
        "label": "每点基因数中位数",
        "label_en": "Median detected genes per spot",
        "domain": "expression",
        "bad": "low",
    },
    {
        "key": "detected_gene_fraction",
        "label": "切片检出基因比例",
        "label_en": "Detected-gene fraction",
        "domain": "expression",
        "bad": "low",
    },
    {
        "key": "library_complexity",
        "label": "表达复杂度",
        "label_en": "Library complexity",
        "domain": "expression",
        "bad": "low",
    },
    {
        "key": "zero_fraction",
        "label": "零表达比例",
        "label_en": "Zero-expression fraction",
        "domain": "expression",
        "bad": "high",
    },
    {
        "key": "regional_low_depth_cluster_fraction",
        "label": "区域低捕获聚集比例",
        "label_en": "Regional low-capture cluster fraction",
        "domain": "expression",
        "bad": "high",
    },
    {
        "key": "local_low_depth_fraction",
        "label": "局部低捕获点比例",
        "label_en": "Local low-capture spot fraction",
        "domain": "expression",
        "bad": "high",
    },
    {
        "key": "median_pct_mito",
        "label": "线粒体比例中位数",
        "label_en": "Median mitochondrial fraction",
        "domain": "damage",
        "bad": "high",
    },
)


_REPORT_METRICS: tuple[dict[str, str], ...] = (
    {"key": "median_n_genes", "label": "Median genes / spot", "unit": "genes"},
    {
        "key": "median_total_counts",
        "label": "Median captured counts / spot",
        "unit": "captured counts",
    },
    {
        "key": "cell_density",
        "label": "Cell / spot density",
        "unit": "locations / coordinate-area",
    },
    {
        "key": "median_nn_distance",
        "label": "Observed spot spacing (resolution proxy)",
        "unit": "coordinate units",
    },
    {"key": "detected_genes", "label": "Total detected genes / slice", "unit": "genes"},
)


_EXCLUSION_EVIDENCE_GROUPS: tuple[dict[str, Any], ...] = (
    {
        "key": "expression_capture",
        "label": "表达捕获不足",
        "label_en": "Reduced expression capture",
        "metrics": (
            "median_total_counts",
            "median_n_genes",
            "detected_gene_fraction",
            "library_complexity",
            "zero_fraction",
        ),
        "panel": "expression",
    },
    {
        "key": "tissue_amount",
        "label": "组织量或密度不足",
        "label_en": "Reduced tissue amount or density",
        "metrics": ("cell_density", "n_locations"),
        "panel": "kde",
    },
    {
        "key": "fragmentation",
        "label": "组织碎片化或连通性下降",
        "label_en": "Fragmentation or reduced connectedness",
        "metrics": ("fragmentation", "largest_component_fraction"),
        "panel": "components",
    },
    {
        "key": "expression_continuity",
        "label": "跨切片表达不连续",
        "label_en": "Cross-slice expression discontinuity",
        "metrics": ("expression_profile", "celltype_composition"),
        "panel": "additional",
    },
    {
        "key": "regional_low_capture",
        "label": "区域低表达聚集",
        "label_en": "Regional low-capture cluster",
        "metrics": ("local_low_depth_fraction", "regional_low_depth_cluster_fraction"),
        "panel": "expression",
    },
    {
        "key": "damage",
        "label": "损伤或线粒体信号",
        "label_en": "Damage or mitochondrial signal",
        "metrics": ("median_pct_mito",),
        "panel": "additional",
    },
)


def _number(value: Any) -> Optional[float]:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _fmt(value: Any) -> str:
    number = _number(value)
    if number is None:
        return "NA"
    magnitude = abs(number)
    if magnitude and (magnitude < 0.01 or magnitude >= 1_000_000):
        return f"{number:.2e}"
    if magnitude < 10:
        return f"{number:.3f}".rstrip("0").rstrip(".")
    if magnitude < 100:
        return f"{number:.1f}".rstrip("0").rstrip(".")
    return f"{number:,.0f}"


def _relative_delta(actual: float, expected: float) -> Optional[float]:
    if abs(expected) <= _EPS:
        return None
    return actual / expected - 1.0


def _evidence_for_row(
    row: Mapping[str, Any], language: str = "zh"
) -> list[dict[str, Any]]:
    evidence: list[dict[str, Any]] = []
    for definition in _EVIDENCE_DEFINITIONS:
        key = definition["key"]
        actual = _number(row.get(key))
        expected = _number(row.get(f"{key}_expected"))
        anomaly = _number(row.get(f"{key}_anomaly")) or 0.0
        if actual is None or expected is None:
            continue
        delta = _relative_delta(actual, expected)
        adverse = actual < expected if definition["bad"] == "low" else actual > expected
        if anomaly < 0.28 and not (
            adverse and delta is not None and abs(delta) >= 0.25
        ):
            continue
        if language == "en":
            direction = (
                "lower"
                if actual < expected
                else "higher" if actual > expected else "similar"
            )
            if delta is None:
                delta_text = f"changed from {_fmt(expected)} to {_fmt(actual)} relative to the local window"
            else:
                delta_text = f"was {abs(delta) * 100:.1f}% {direction} than the local-window expectation"
            label = definition["label_en"]
            sentence = f"{label} {delta_text} (observed {_fmt(actual)}; expected {_fmt(expected)})"
        else:
            direction = (
                "低" if actual < expected else "高" if actual > expected else "相当"
            )
            if delta is None:
                delta_text = f"由 {_fmt(expected)} 变为 {_fmt(actual)}"
            else:
                delta_text = f"比局部窗口期望{direction} {abs(delta) * 100:.1f}%"
            label = definition["label"]
            sentence = (
                f"{label}{delta_text}（观测 {_fmt(actual)}；期望 {_fmt(expected)}）"
            )
        priority = anomaly
        if adverse and key == "cell_density" and delta is not None and delta <= -0.25:
            priority += 1.6
        elif adverse and key == "n_locations" and delta is not None and delta <= -0.35:
            priority += 1.2
        elif (
            adverse
            and key in {"fragmentation", "hole_fraction"}
            and delta is not None
            and delta >= 0.5
        ):
            priority += 0.9
        elif (
            adverse
            and key in {"median_total_counts", "median_n_genes"}
            and delta is not None
            and delta <= -0.25
        ):
            priority += 0.8
        evidence.append(
            {
                **definition,
                "label": label,
                "actual": actual,
                "expected": expected,
                "delta": delta,
                "direction": direction,
                "delta_text": delta_text,
                "anomaly": anomaly,
                "adverse": adverse,
                "priority": priority,
                "sentence": sentence,
            }
        )
    evidence.sort(key=lambda item: (item["adverse"], item["priority"]), reverse=True)
    return evidence


def _likely_causes(evidence: Sequence[Mapping[str, Any]], language: str = "zh") -> str:
    keys = {str(item["key"]) for item in evidence if item.get("adverse")}
    causes: list[str] = []
    translations = {
        "density": (
            "local tissue loss or tearing, spot loss, or under-segmentation"
            if language == "en"
            else "局部组织缺失/撕裂、捕获点丢失或分割漏检"
        ),
        "fragmentation": (
            "section fragmentation, detached tissue, or sparse segmentation"
            if language == "en"
            else "切片破碎、组织块脱落或稀疏分割"
        ),
        "holes": (
            "tissue loss, a dry-spot/capture failure, or a true anatomical cavity"
            if language == "en"
            else "组织空洞、dry-spot/捕获失败或真实腔室"
        ),
        "expression": (
            "reduced RNA capture, RNA degradation, or insufficient effective sequencing/molecule capture"
            if language == "en"
            else "RNA 捕获效率下降、降解或有效测序/分子捕获不足"
        ),
        "damage": (
            "cellular damage or a stress-associated expression shift"
            if language == "en"
            else "细胞损伤或应激相关表达偏移"
        ),
        "other": (
            "biological plane-of-section change, coordinate misalignment, or another technical factor"
            if language == "en"
            else "切面生物学变化、配准坐标差异或技术因素"
        ),
    }
    if keys & {"cell_density", "n_locations"}:
        causes.append(translations["density"])
    if keys & {"fragmentation", "largest_component_fraction"}:
        causes.append(translations["fragmentation"])
    if keys & {"hole_fraction", "knn_tail_ratio"}:
        causes.append(translations["holes"])
    if keys & {
        "median_total_counts",
        "median_n_genes",
        "detected_gene_fraction",
        "library_complexity",
        "zero_fraction",
        "regional_low_depth_cluster_fraction",
        "local_low_depth_fraction",
    }:
        causes.append(translations["expression"])
    if "median_pct_mito" in keys:
        causes.append(translations["damage"])
    if not causes:
        causes.append(translations["other"])
    return ("; " if language == "en" else "；").join(dict.fromkeys(causes))


def build_directional_slice_explanations(
    metrics: pd.DataFrame, *, language: str = "zh"
) -> pd.DataFrame:
    """Return one plain-language binary explanation and evidence list per slice."""
    if language not in {"zh", "en"}:
        raise ValueError("language must be 'zh' or 'en'")
    rows: list[dict[str, Any]] = []
    for record in metrics.to_dict("records"):
        call = str(record.get("final_call", record.get("recommendation", "")))
        score = _number(record.get("quality_anomaly_score"))
        evidence = _evidence_for_row(record, language)
        adverse = [item for item in evidence if item["adverse"]]
        adaptive_resolution = str(record.get("decision_basis")) in {
            "adaptive_multidomain_resolution",
            "tiered_multidomain_resolution",
        }
        tier = str(record.get("review_resolution_tier", "")).strip()
        stability = _number(
            record.get(
                "review_tier_window_stability", record.get("adaptive_window_stability")
            )
        )
        windows = str(record.get("adaptive_windows_tested", "3|5|7")).replace("|", "/")
        if language == "en":
            if call == "exclude":
                detail = "; ".join(item["sentence"] for item in adverse[:4])
                if not detail:
                    detail = "multi-domain anomaly scores increased under two-sided slice comparison"
                if adaptive_resolution:
                    reason = (
                        f"Exclude: {detail}. Its score ({_fmt(score)}) entered the calibrated "
                        f"{tier or 'review'} tier, which required two-sided context, no partial-structure "
                        f"protection, severe corroborated multi-domain evidence, and {_fmt(stability * 100 if stability is not None else None)}% "
                        f"agreement across {windows}-slice windows. Plausible causes include "
                        f"{_likely_causes(adverse[:4], language)}; the H5AD alone cannot identify the specific "
                        "experimental cause."
                    )
                else:
                    reason = (
                        f"Exclude: {detail}. The slice had two-sided local-window evidence, did not trigger "
                        f"partial-structure protection, and its total anomaly score ({_fmt(score)}) passed the "
                        "independently validated high-specificity exclusion threshold. Plausible causes include "
                        f"{_likely_causes(adverse[:4], language)}; the H5AD alone cannot identify the specific "
                        "experimental cause."
                    )
            elif str(record.get("decision_basis")) == "independently_certified":
                reason = (
                    f"Keep: the anomaly score ({_fmt(score)}) and detector evidence were within the "
                    "validated retention range, with no corroborated multi-domain evidence for exclusion."
                )
            elif 'calibrated keep-only score band' in str(record.get('review_resolution_reason', '')):
                reason = (f"Keep after stage-2 tier and permission review ({tier}): this locked score band is keep-only; "
                          "the band's multi-domain and multi-window exclusion tests are not executed. "
                          "This operational retention is not evidence that the slice is free of defects.")
            elif str(record.get("decision_basis")) == "review_resolved_keep":
                reason = (
                    f"Keep after fine screen: stage-1 threshold triage placed this slice in the internal "
                    f"review queue, but the second-stage multi-domain, anatomical, confidence, and 3/5/7-window "
                    f"checks did not all support exclusion. No additional user decision is required."
                )
            else:
                reason = (
                    f"Keep: the anomaly score ({_fmt(score)}) did not jointly pass the high-specificity "
                    "exclusion threshold and detector guardrails. The final binary policy retained the slice "
                    "to avoid destructive false exclusion; no additional user decision is required."
                )
        else:
            if call == "exclude":
                detail = "；".join(item["sentence"] for item in adverse[:4])
                if not detail:
                    detail = "多维异常得分在双侧相邻切片比较中显著升高"
                if adaptive_resolution:
                    reason = (
                        f"排除：{detail}。其总分 {_fmt(score)} 进入校准后的 {tier or 'review'} 分层，因此解析器要求"
                        "双侧上下文、未触发部分结构保护、严重且相互印证的多域证据，以及"
                        f"在 {windows} 张窗口下 {_fmt(stability * 100 if stability is not None else None)}% 的一致结论。"
                        f"可能原因包括：{_likely_causes(adverse[:4], language)}；仅凭 H5AD 无法区分具体实验原因。"
                    )
                else:
                    reason = (
                        f"排除：{detail}。该切片有双侧局部窗口证据、未触发部分结构保护，"
                        f"且总异常分数 {_fmt(score)} 通过经独立验证的高特异度排除门槛。"
                        f"可能原因包括：{_likely_causes(adverse[:4], language)}；"
                        "仅凭 H5AD 无法区分具体实验原因。"
                    )
            elif str(record.get("decision_basis")) == "independently_certified":
                reason = (
                    f"保留：异常分数 {_fmt(score)} 与检测器证据均处于经验证的保留范围，"
                    "未见足以支持排除的多域局部异常。"
                )
            elif 'calibrated keep-only score band' in str(record.get('review_resolution_reason', '')):
                reason = (f"第二阶段已完成分层与排除权限审查（{tier}）：该锁定区间为仅保留，判为keep；"
                          "不会执行该层的多域与跨窗口排除检查；保留不证明切片没有缺陷。")
            elif str(record.get("decision_basis")) == "review_resolved_keep":
                reason = (
                    "细筛后保留：第一阶段阈值初筛将该切片放入内部 review 队列，但第二阶段的多域证据、"
                    "解剖保护、置信度和 3/5/7 窗口稳定性条件未能同时支持排除。无需用户再次判断。"
                )
            else:
                reason = (
                    f"保留：异常分数 {_fmt(score)} 未同时通过高特异度排除门槛与检测器保护条件；"
                    "最终二元策略为避免误删而保留该切片，不要求用户再次作质量判定。"
                )
        rows.append(
            {
                "slice_id": str(record.get("slice_id")),
                "user_reason": reason,
                "primary_reason": (
                    ("; " if language == "en" else "；").join(
                        item["sentence"] for item in adverse[:2]
                    )
                    if adverse
                    else (
                        "No reproducible major anomaly"
                        if language == "en"
                        else "未见可重复的主要异常"
                    )
                ),
                "likely_causes": (
                    _likely_causes(adverse[:4], language)
                    if adverse
                    else (
                        "No specific technical anomaly indicated"
                        if language == "en"
                        else "无明确技术异常指征"
                    )
                ),
                "directional_evidence": evidence,
            }
        )
    return pd.DataFrame(rows)


def _fixed_bandwidth_kde(
    coords: np.ndarray, bandwidth: float, evaluate: Optional[np.ndarray] = None
) -> np.ndarray:
    if evaluate is None:
        evaluate = coords
    if not len(coords) or not len(evaluate):
        return np.zeros(len(evaluate), dtype=float)
    bandwidth = max(float(bandwidth), math.sqrt(_EPS))
    result = np.zeros(len(evaluate), dtype=float)
    chunk = 512
    normalizer = 2.0 * math.pi * bandwidth * bandwidth
    for start in range(0, len(evaluate), chunk):
        points = evaluate[start : start + chunk]
        delta = points[:, None, :] - coords[None, :, :]
        squared = np.einsum("ijk,ijk->ij", delta, delta)
        result[start : start + chunk] = np.exp(
            -0.5 * squared / (bandwidth * bandwidth)
        ).sum(axis=1)
    return result / normalizer


def _sample_indices(n: int, maximum: int) -> np.ndarray:
    if n <= maximum:
        return np.arange(n, dtype=int)
    return np.linspace(0, n - 1, maximum, dtype=int)


def _stable_sample_indices(
    slice_id: str, n: int, maximum: int, seed: int = 13
) -> np.ndarray:
    if n <= maximum:
        return np.arange(n, dtype=int)
    stable = int(hashlib.sha256(slice_id.encode("utf-8")).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed + stable)
    return np.sort(rng.choice(n, maximum, replace=False))


def _component_ranks(
    coords: np.ndarray, median_nn: float
) -> tuple[np.ndarray, int, float]:
    if len(coords) < 2 or not np.isfinite(median_nn) or median_nn <= 0:
        return np.ones(len(coords), dtype=int), 1, 1.0
    tree = cKDTree(coords)
    graph = tree.sparse_distance_matrix(
        tree,
        max_distance=max(float(median_nn) * 3.25, math.sqrt(_EPS)),
        output_type="coo_matrix",
    )
    graph.data[:] = 1
    count, labels = csgraph.connected_components(graph.tocsr(), directed=False)
    sizes = np.bincount(labels, minlength=count)
    label_order = np.argsort(-sizes, kind="stable")
    label_to_rank = np.empty(count, dtype=int)
    label_to_rank[label_order] = np.arange(1, count + 1)
    ranks = label_to_rank[labels]
    largest = float(sizes.max() / max(len(coords), 1))
    return ranks, int(count), largest


def _manifest_source_contract(
    manifest: Mapping[str, Any],
    source_h5ad: Optional[Union[str, Path]],
    slice_key: Optional[str],
    spatial_key: Optional[str],
) -> tuple[Optional[Path], Optional[str], Optional[str]]:
    sources = manifest.get("sources") or []
    record = sources[0] if sources else {}
    source_value = source_h5ad or record.get("path")
    source_path = Path(source_value).expanduser() if source_value else None
    if not slice_key:
        match = re.match(r"obs:(.+)", str(record.get("slice_source", "")))
        if match:
            slice_key = match.group(1)
        else:
            match = re.match(
                r"obsm:([^\[]+)\[:,\s*2\]", str(record.get("slice_source", ""))
            )
            slice_key = f"obsm:{match.group(1)}[:,2]" if match else None
    if not spatial_key:
        match = re.match(r"obsm:([^\[]+)", str(record.get("coordinate_source", "")))
        spatial_key = match.group(1) if match else None
    return source_path, slice_key, spatial_key


def _first_numeric_obs(
    adata: Any, requested: Optional[str], candidates: Sequence[str]
) -> tuple[Optional[str], Optional[np.ndarray]]:
    keys = [requested] if requested else list(candidates)
    for key in keys:
        if key and key in adata.obs:
            values = pd.to_numeric(adata.obs[key], errors="coerce").to_numpy(
                dtype=float
            )
            if np.isfinite(values).any():
                return str(key), values
    return None, None


def _prepare_source_evidence_points(
    source_h5ad: Path,
    order: Sequence[str],
    bandwidth: float,
    max_points_per_slice: int,
    *,
    slice_key: str,
    spatial_key: str,
    total_counts_key: Optional[str] = None,
    n_genes_key: Optional[str] = None,
    mito_key: Optional[str] = None,
) -> tuple[dict[str, dict[str, Any]], dict[str, np.ndarray], dict[str, Any]]:
    import anndata as ad

    adata = ad.read_h5ad(source_h5ad, backed="r")
    try:
        if spatial_key not in adata.obsm:
            raise KeyError(
                f"spatial key {spatial_key!r} is absent from source H5AD obsm"
            )
        coordinates_all = np.asarray(adata.obsm[spatial_key], dtype=float)
        if coordinates_all.ndim != 2 or coordinates_all.shape[1] < 2:
            raise ValueError(
                f"obsm[{spatial_key!r}] must contain at least two coordinate columns"
            )
        if slice_key in adata.obs:
            labels = adata.obs[slice_key].astype(str).to_numpy()
            resolved_slice_source = f"obs:{slice_key}"
        else:
            z_match = re.fullmatch(r"obsm:([^\[]+)\[:,\s*2\]", slice_key)
            z_key = z_match.group(1) if z_match else spatial_key
            if z_key not in adata.obsm:
                raise KeyError(
                    f"slice key {slice_key!r} is absent from source H5AD obs and "
                    f"obsm[{z_key!r}] is unavailable"
                )
            z_coordinates = np.asarray(adata.obsm[z_key], dtype=float)
            if z_coordinates.ndim != 2 or z_coordinates.shape[1] < 3:
                raise ValueError(
                    "A source without an obs slice key needs a 3D obsm array "
                    "whose third column contains discrete z values"
                )
            z_values = z_coordinates[:, 2]
            labels = np.asarray(
                [f"z{value:g}" if np.isfinite(value) else "zNA" for value in z_values],
                dtype=str,
            )
            resolved_slice_source = f"obsm:{z_key}[:,2]"
        coordinates_all = coordinates_all[:, :2]
        counts_name, counts_all = _first_numeric_obs(
            adata,
            total_counts_key,
            ("total_counts", "n_counts", "total_counts_X", "UMI_count", "nCount_RNA"),
        )
        genes_name, genes_all = _first_numeric_obs(
            adata,
            n_genes_key,
            ("n_genes_by_counts", "n_genes", "nFeature_RNA", "gene_count"),
        )
        mito_name, mito_all = _first_numeric_obs(
            adata,
            mito_key,
            ("pct_counts_mt", "pct_counts_mito", "percent_mito", "mito_percent"),
        )
        display: dict[str, dict[str, Any]] = {}
        full: dict[str, np.ndarray] = {}
        for slice_id in order:
            row_indices = np.flatnonzero(labels == slice_id)
            coords = coordinates_all[row_indices]
            valid = np.isfinite(coords).all(axis=1)
            row_indices = row_indices[valid]
            coords = coords[valid]
            full[slice_id] = coords
            chosen = _stable_sample_indices(slice_id, len(coords), max_points_per_slice)
            if len(coords) >= 2:
                nearest = cKDTree(coords).query(coords, k=2)[0][:, 1]
                positive = nearest[np.isfinite(nearest) & (nearest > 0)]
                median_nn = (
                    float(np.median(positive)) if len(positive) else float("nan")
                )
            else:
                median_nn = float("nan")
            component_rank, component_count, largest_fraction = _component_ranks(
                coords, median_nn
            )
            sampled_counts = (
                counts_all[row_indices[chosen]]
                if counts_all is not None
                else np.asarray([])
            )
            sampled_genes = (
                genes_all[row_indices[chosen]]
                if genes_all is not None
                else np.asarray([])
            )
            if len(chosen) and (not sampled_counts.size or not sampled_genes.size):
                sampled_matrix = adata.X[row_indices[chosen]]
                if sparse.issparse(sampled_matrix):
                    if not sampled_counts.size:
                        sampled_counts = np.asarray(sampled_matrix.sum(axis=1)).ravel()
                    if not sampled_genes.size:
                        sampled_genes = np.asarray(
                            sampled_matrix.getnnz(axis=1)
                        ).ravel()
                else:
                    sampled_matrix = np.asarray(sampled_matrix)
                    if not sampled_counts.size:
                        sampled_counts = np.asarray(sampled_matrix.sum(axis=1)).ravel()
                    if not sampled_genes.size:
                        sampled_genes = np.count_nonzero(sampled_matrix, axis=1)
                counts_name = counts_name or "X:row_sum"
                genes_name = genes_name or "X:nonzero_genes"
            sampled_mito = (
                mito_all[row_indices[chosen]]
                if mito_all is not None
                else np.asarray([])
            )
            if sampled_counts.size:
                log_counts = np.log1p(np.maximum(sampled_counts, 0))
                lo, hi = np.nanquantile(log_counts, [0.02, 0.98])
                depth = np.clip((log_counts - lo) / max(float(hi - lo), _EPS), 0, 1)
            else:
                log_counts = np.asarray([])
                depth = np.zeros(len(chosen), dtype=float)
            display[slice_id] = {
                "x": np.round(coords[chosen, 0], 3).tolist(),
                "y": np.round(coords[chosen, 1], 3).tolist(),
                "kde": np.round(
                    _fixed_bandwidth_kde(coords, bandwidth, coords[chosen]), 8
                ).tolist(),
                "depth": np.round(depth, 4).tolist(),
                "captured_counts": np.round(sampled_counts, 3).tolist(),
                "log_captured_counts": np.round(log_counts, 5).tolist(),
                "n_genes": np.round(sampled_genes, 3).tolist(),
                "pct_mito": np.round(sampled_mito, 4).tolist(),
                "component_rank": component_rank[chosen].astype(int).tolist(),
                "component_count": component_count,
                "largest_component_fraction": round(largest_fraction, 5),
                "displayed": int(len(chosen)),
                "available": int(len(coords)),
            }
        metadata = {
            "available": True,
            "source_h5ad": str(source_h5ad),
            "slice_key": resolved_slice_source,
            "spatial_key": spatial_key,
            "total_counts_key": counts_name,
            "n_genes_key": genes_name,
            "mito_key": mito_name,
            "absolute_expression_available": bool(counts_name or genes_name),
            "component_evidence_available": True,
        }
        return display, full, metadata
    finally:
        adata.file.close()


def _prepare_kde_points(
    point_samples: Mapping[str, Mapping[str, Sequence[Any]]],
    order: Sequence[str],
    bandwidth: float,
    max_points_per_slice: int,
) -> tuple[dict[str, dict[str, Any]], dict[str, np.ndarray]]:
    display: dict[str, dict[str, Any]] = {}
    full: dict[str, np.ndarray] = {}
    for slice_id in order:
        source = point_samples.get(slice_id, {})
        x = np.asarray(source.get("x", []), dtype=float)
        y = np.asarray(source.get("y", []), dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)
        coords = np.column_stack([x[valid], y[valid]])
        depth = np.asarray(source.get("depth", np.zeros(len(x))), dtype=float)
        depth = (
            depth[valid] if len(depth) == len(x) else np.zeros(len(coords), dtype=float)
        )
        kde = _fixed_bandwidth_kde(coords, bandwidth)
        full[slice_id] = coords
        chosen = _sample_indices(len(coords), max_points_per_slice)
        display[slice_id] = {
            "x": np.round(coords[chosen, 0], 3).tolist(),
            "y": np.round(coords[chosen, 1], 3).tolist(),
            "kde": np.round(kde[chosen], 8).tolist(),
            "depth": np.round(np.clip(depth[chosen], 0, 1), 3).tolist(),
            "displayed": int(len(chosen)),
            "available": int(len(coords)),
        }
    return display, full


def _prepare_preregistered_evidence_points(source, payload, order, bandwidth, maximum):
    """Reuse full raw coordinate caches; fit/score/source data remain untouched."""
    import hashlib
    registration = payload['preregistration']
    points, coordinates = {}, {}
    full_cache_slices = []
    for sid in order:
        sample = payload['point_samples'][sid]
        info = registration['slices'][sid]
        raw = np.column_stack([sample['x'], sample['y']]).astype(float)
        chosen = _sample_indices(len(raw), maximum)
        raw = raw[chosen]
        cache_path = source / 'display_preregistration' / info['cache_file']
        full = raw
        full_available = False
        if cache_path.exists():
            actual = hashlib.sha256(cache_path.read_bytes()).hexdigest()
            if actual != registration['files'][info['cache_file']]:
                raise ValueError('Coordinate cache checksum mismatch: ' + sid)
            with np.load(cache_path, allow_pickle=False) as cache:
                all_raw = cache['raw_xy']
                source_rows = cache['source_row']
                wanted = np.asarray(sample['source_row'])[chosen]
                lookup = {int(row): i for i, row in enumerate(source_rows)}
                offsets = np.asarray([lookup[int(row)] for row in wanted])
                if not np.allclose(all_raw[offsets], raw, atol=1e-10):
                    raise ValueError('Display sample/source row mismatch: ' + sid)
                full = all_raw[np.isfinite(all_raw).all(axis=1)]
                nearest = cKDTree(full).query(full, k=2)[0][:, 1] if len(full) > 1 else np.array([])
                positive = nearest[np.isfinite(nearest) & (nearest > 0)]
                ranks, component_count, largest = _component_ranks(full, float(np.median(positive)) if len(positive) else float('nan'))
                valid_offsets = np.full(len(all_raw), -1, dtype=int)
                valid_offsets[np.flatnonzero(np.isfinite(all_raw).all(axis=1))] = np.arange(len(full))
                sampled_ranks = ranks[valid_offsets[offsets]].tolist()
                full_available = True
                full_cache_slices.append(sid)
        if not full_available:
            sampled_ranks, component_count, largest = [], None, None
        matrix = np.asarray(info['matrix'], dtype=float)
        display = (matrix @ np.column_stack([raw, np.ones(len(raw))]).T).T[:, :2]
        counts = np.asarray(sample.get('counts', [None]*len(sample['x'])), dtype=float)[chosen]
        # Missing capture evidence stays missing, never a zero measurement.
        safe_counts = [float(v) if np.isfinite(v) else None for v in counts]
        points[sid] = dict(x=raw[:,0].tolist(), y=raw[:,1].tolist(),
            display_x=display[:,0].tolist(), display_y=display[:,1].tolist(),
            kde=_fixed_bandwidth_kde(full, bandwidth, raw).tolist(),
            captured_counts=safe_counts,
            log_captured_counts=[float(np.log1p(max(v,0))) if v is not None else None for v in safe_counts],
            n_genes=[], component_rank=sampled_ranks, component_count=component_count,
            largest_component_fraction=largest, displayed=len(raw), available=sample['n_total'],
            source_row=np.asarray(sample['source_row'])[chosen].tolist(),
            obs_names=np.asarray(sample['obs_names'])[chosen].tolist(),
            source_file=info['source_file'], coordinate_source=info['coordinate_source'],
            evidence_geometry_scope='full raw points' if full_available else 'display sample only')
        # Density-deficit boxes in raw unregistered frames would imply false correspondence.
        coordinates[sid] = (matrix @ np.column_stack([full, np.ones(len(full))]).T).T[:, :2]
    return points, coordinates, dict(available=True, absolute_expression_available=any(any(v is not None for v in p['captured_counts']) for p in points.values()),
        component_evidence_available=len(full_cache_slices)==len(order),
        full_coordinate_cache_slices=full_cache_slices,
        source='saved preregistration cache and detector-selected capture sample',
        gene_maps_available=False, coordinate_frame=registration['frame_id'],
        note='KDE and components use full raw geometry when its verified cache is present. '
             'Without that cache, KDE is sample-only and component maps are unavailable. '
             'Captured counts use the QC-selected expression source; per-point gene maps are unavailable.')


def _largest_density_deficit_region(
    focal_id: str,
    order: Sequence[str],
    coordinates: Mapping[str, np.ndarray],
    bandwidth: float,
    half_window: int = 2,
    language: str = "zh",
) -> Optional[dict[str, Any]]:
    index = order.index(focal_id)
    lo = max(0, index - half_window)
    hi = min(len(order), index + half_window + 1)
    neighbor_ids = [
        order[i]
        for i in range(lo, hi)
        if i != index and len(coordinates.get(order[i], []))
    ]
    focal = coordinates.get(focal_id)
    if focal is None or len(focal) < 3 or len(neighbor_ids) < 2:
        return None
    arrays = [focal] + [coordinates[item] for item in neighbor_ids]
    pooled = np.vstack(arrays)
    min_xy = np.nanmin(pooled, axis=0)
    max_xy = np.nanmax(pooled, axis=0)
    span = np.maximum(max_xy - min_xy, bandwidth)
    min_xy -= span * 0.035
    max_xy += span * 0.035
    gx = np.linspace(min_xy[0], max_xy[0], 42)
    gy = np.linspace(min_xy[1], max_xy[1], 32)
    xx, yy = np.meshgrid(gx, gy)
    grid = np.column_stack([xx.ravel(), yy.ravel()])
    focal_density = _fixed_bandwidth_kde(focal, bandwidth, grid).reshape(
        len(gy), len(gx)
    )
    neighbor_density = np.stack(
        [
            _fixed_bandwidth_kde(coordinates[item], bandwidth, grid)
            for item in neighbor_ids
        ],
        axis=0,
    )
    expected = np.mean(neighbor_density, axis=0).reshape(len(gy), len(gx))
    positive = expected[expected > 0]
    if not positive.size:
        return None
    support = expected >= np.quantile(positive, 0.58)
    ratio = focal_density / np.maximum(expected, _EPS)
    mask = support & (ratio <= 0.48)
    labels, count = ndimage.label(mask)
    if not count:
        return None
    sizes = ndimage.sum(mask, labels, index=np.arange(1, count + 1))
    label = int(np.argmax(sizes) + 1)
    region = labels == label
    if int(region.sum()) < 3:
        return None
    iy, ix = np.where(region)
    x0, x1 = float(gx[ix.min()]), float(gx[ix.max()])
    y0, y1 = float(gy[iy.min()]), float(gy[iy.max()])
    support_count = max(int(support.sum()), 1)
    return {
        "x0": round(x0, 3),
        "x1": round(x1, 3),
        "y0": round(y0, 3),
        "y1": round(y1, 3),
        "deficit_fraction": round(float(region.sum()) / support_count, 3),
        "median_density_ratio": round(float(np.median(ratio[region])), 3),
        "neighbor_ids": neighbor_ids,
        "note": (
            (
                "The yellow box marks the largest contiguous candidate region in which neighboring slices "
                "have high expected KDE but the focal-slice density is below 48% of that expectation; it "
                f"covers approximately {100 * float(region.sum()) / support_count:.1f}% of comparable "
                "high-density grid cells."
            )
            if language == "en"
            else (
                f"黄色框为局部窗口 KDE 期望较高、但本切片密度不足其 48% 的最大连续候选区；"
                f"覆盖约 {100 * float(region.sum()) / support_count:.1f}% 的可比高密度网格。"
            )
        ),
    }


def _summary_stat(values: pd.Series) -> dict[str, Optional[float]]:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if not len(numeric):
        return {"median": None, "q25": None, "q75": None}
    return {
        "median": float(np.median(numeric)),
        "q25": float(np.quantile(numeric, 0.25)),
        "q75": float(np.quantile(numeric, 0.75)),
    }


def _metric_supports_group(
    record: Mapping[str, Any], metric: str, threshold: float = 0.5
) -> bool:
    anomaly = _number(record.get(f"{metric}_anomaly"))
    if anomaly is None or anomaly < threshold:
        return False
    if metric in {"expression_profile", "celltype_composition"}:
        return True
    definition = next(
        (item for item in _EVIDENCE_DEFINITIONS if item["key"] == metric), None
    )
    if definition is None:
        return False
    actual = _number(record.get(metric))
    expected = _number(record.get(f"{metric}_expected"))
    if actual is None or expected is None:
        return False
    return actual < expected if definition["bad"] == "low" else actual > expected


def _exclusion_evidence_summary(
    metrics: pd.DataFrame, language: str
) -> list[dict[str, Any]]:
    excluded = metrics.loc[metrics["final_call"].astype(str).eq("exclude")]
    total = int(len(excluded))
    results: list[dict[str, Any]] = []
    for group in _EXCLUSION_EVIDENCE_GROUPS:
        slice_ids = [
            str(record["slice_id"])
            for record in excluded.to_dict("records")
            if any(
                _metric_supports_group(record, metric) for metric in group["metrics"]
            )
        ]
        if not slice_ids:
            continue
        results.append(
            {
                "key": group["key"],
                "label": group["label_en"] if language == "en" else group["label"],
                "count": len(slice_ids),
                "percent": round(100 * len(slice_ids) / max(total, 1), 1),
                "slice_ids": slice_ids,
                "panel": group["panel"],
            }
        )
    return sorted(results, key=lambda item: (-item["count"], item["label"]))


def _visual_panel_findings(
    record: Mapping[str, Any], language: str
) -> dict[str, dict[str, str]]:
    evidence = [item for item in _evidence_for_row(record, language) if item["adverse"]]
    definitions = {
        "kde": {"cell_density", "n_locations", "knn_tail_ratio", "hole_fraction"},
        "expression": {
            "median_total_counts",
            "median_n_genes",
            "detected_gene_fraction",
            "library_complexity",
            "zero_fraction",
            "regional_low_depth_cluster_fraction",
            "local_low_depth_fraction",
        },
        "components": {"fragmentation", "largest_component_fraction"},
    }
    output: dict[str, dict[str, str]] = {}
    excluded = str(record.get("final_call")) == "exclude"
    for panel, keys in definitions.items():
        items = [item for item in evidence if item["key"] in keys]
        strong = [item for item in items if item["anomaly"] >= 0.5]
        supports = bool(excluded and strong)
        selected = (strong or items)[:2]
        if language == "en":
            prefix = "Supports exclusion" if supports else "Context only"
            detail = "; ".join(item["sentence"] for item in selected)
            if not detail:
                detail = "no adverse directional evidence reached the display threshold"
        else:
            prefix = "支持排除" if supports else "背景参照"
            detail = "；".join(item["sentence"] for item in selected)
            if not detail:
                detail = "未见达到展示阈值的不利方向证据"
        output[panel] = {
            "status": "supports" if supports else "context",
            "summary": (
                f"{prefix}：{detail}" if language == "zh" else f"{prefix}: {detail}"
            ),
        }
    return output


def _write_paper_tables(
    metrics: pd.DataFrame,
    output_dir: Path,
    selected_window: Any,
    *,
    dataset_name: str,
    title: str,
    language: str = "zh",
) -> dict[str, str]:
    calls = metrics["final_call"].astype(str)
    summary: dict[str, Any] = {
        "dataset": dataset_name,
        "n_slices": int(len(metrics)),
        "keep": int(calls.eq("keep").sum()),
        "exclude": int(calls.eq("exclude").sum()),
        "exclude_rate_percent": round(float(calls.eq("exclude").mean() * 100), 2),
        "selected_detection_window": selected_window,
        "final_states": "keep / exclude",
    }
    for definition in _REPORT_METRICS:
        key = definition["key"]
        if key not in metrics:
            continue
        statistic = _summary_stat(metrics[key])
        summary[f"{key}_median"] = statistic["median"]
        summary[f"{key}_q25"] = statistic["q25"]
        summary[f"{key}_q75"] = statistic["q75"]
    summary_path = output_dir / "paper_dataset_summary.csv"
    pd.DataFrame([summary]).to_csv(summary_path, index=False)

    explanation = build_directional_slice_explanations(
        metrics, language=language
    ).set_index("slice_id")
    excluded = metrics.loc[calls.eq("exclude")].copy()
    columns = [
        "slice_index",
        "slice_id",
        "final_call",
        "quality_anomaly_score",
        "n_locations",
        "cell_density",
        "cell_density_expected",
        "median_n_genes",
        "median_n_genes_expected",
        "median_total_counts",
        "median_total_counts_expected",
        "median_nn_distance",
        "detected_genes",
    ]
    columns = [column for column in columns if column in excluded]
    excluded = excluded[columns]
    excluded["primary_exclusion_reason"] = excluded["slice_id"].map(
        explanation["primary_reason"]
    )
    excluded["possible_causes_not_diagnosis"] = excluded["slice_id"].map(
        explanation["likely_causes"]
    )
    exclude_path = output_dir / "paper_excluded_slices.csv"
    excluded.to_csv(exclude_path, index=False)

    paper_html = output_dir / "paper_tables.html"
    summary_view = pd.DataFrame(
        [
            {
                "Dataset": summary["dataset"],
                "Slices": summary["n_slices"],
                "Keep": summary["keep"],
                "Exclude": summary["exclude"],
                "Exclude rate": f"{summary['exclude_rate_percent']:.1f}%",
                "Window": selected_window,
                **{
                    definition["label"]: (
                        f"{_fmt(summary.get(definition['key'] + '_median'))} "
                        f"[{_fmt(summary.get(definition['key'] + '_q25'))}, "
                        f"{_fmt(summary.get(definition['key'] + '_q75'))}]"
                    )
                    for definition in _REPORT_METRICS
                    if f"{definition['key']}_median" in summary
                },
            }
        ]
    )
    paper_html.write_text(
        "<!doctype html><meta charset='utf-8'><title>"
        + html.escape(title)
        + " · QC tables</title>"
        "<style>body{font:13px/1.45 Arial,sans-serif;margin:32px;color:#172033}"
        "h1{font-size:20px}h2{font-size:15px;margin-top:28px}table{border-collapse:collapse;width:100%}"
        "th,td{padding:7px 8px;border:1px solid #cfd7e3;text-align:left;vertical-align:top}"
        "th{background:#eef3f8}caption{caption-side:bottom;text-align:left;margin-top:8px;color:#5c687a}</style>"
        + "<h1>"
        + html.escape(title)
        + "</h1><h2>Dataset summary</h2>"
        + summary_view.to_html(index=False, border=0, escape=True)
        + "<p>Values are median [Q1, Q3] across slices. Density and spacing use native coordinate units.</p>"
        + "<h2>Automatically excluded slices</h2>"
        + excluded.to_html(index=False, border=0, escape=True)
        + "<p>Possible causes are interpretations, not experimental diagnoses. H5AD cannot identify raw-read depth, duplication or saturation.</p>",
        encoding="utf-8",
    )
    return {
        "paper_dataset_summary": str(summary_path),
        "paper_excluded_slices": str(exclude_path),
        "paper_tables_html": str(paper_html),
    }


def write_collection_paper_tables(
    dataset_summary_csv: Union[str, Path],
    output_dir: Union[str, Path],
    *,
    title: str = "Spateo referee dataset summary",
) -> dict[str, Any]:
    """Condense a grouped collection report into main- and supplement-table CSVs."""
    source = Path(dataset_summary_csv).expanduser().resolve()
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(source)
    required = {"name", "group", "n_slices", "keep", "exclude"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(
            f"Collection summary is missing required columns: {sorted(missing)}"
        )
    if (
        "review" in data
        and pd.to_numeric(data["review"], errors="coerce").fillna(0).ne(0).any()
    ):
        raise ValueError(
            "Paper tables require the public keep/exclude-only collection summary"
        )
    rows: list[dict[str, Any]] = []
    for record in data.to_dict("records"):
        output: dict[str, Any] = {
            "Biological type": record.get("group"),
            "Subgroup": record.get("subgroup"),
            "Dataset": record.get("name"),
            "Slices": int(record.get("n_slices", 0)),
            "Keep": int(record.get("keep", 0)),
            "Exclude": int(record.get("exclude", 0)),
            "Exclude rate (%)": round(
                100
                * int(record.get("exclude", 0))
                / max(int(record.get("n_slices", 0)), 1),
                2,
            ),
            "Selected window": record.get("selected_window"),
        }
        for definition in _REPORT_METRICS:
            key = definition["key"]
            median = _number(record.get(f"{key}_median"))
            q25 = _number(record.get(f"{key}_q25"))
            q75 = _number(record.get(f"{key}_q75"))
            output[f"{definition['label']}, median [Q1, Q3]"] = (
                f"{_fmt(median)} [{_fmt(q25)}, {_fmt(q75)}]"
                if median is not None
                else ""
            )
        rows.append(output)
    supplementary = pd.DataFrame(rows).sort_values(
        ["Biological type", "Subgroup", "Dataset"], kind="stable", na_position="last"
    )
    supplement_path = destination / "paper_all_datasets_summary.csv"
    supplementary.to_csv(supplement_path, index=False)

    group_rows: list[dict[str, Any]] = []
    for group, frame in data.groupby("group", dropna=False, sort=True):
        slices = int(pd.to_numeric(frame["n_slices"], errors="coerce").fillna(0).sum())
        keep = int(pd.to_numeric(frame["keep"], errors="coerce").fillna(0).sum())
        exclude = int(pd.to_numeric(frame["exclude"], errors="coerce").fillna(0).sum())
        group_rows.append(
            {
                "Biological type": group,
                "Datasets": int(len(frame)),
                "Slices": slices,
                "Keep": keep,
                "Exclude": exclude,
                "Exclude rate (%)": round(100 * exclude / max(slices, 1), 2),
            }
        )
    groups = pd.DataFrame(group_rows)
    group_path = destination / "paper_biological_type_summary.csv"
    groups.to_csv(group_path, index=False)

    html_path = destination / "paper_collection_tables.html"
    html_path.write_text(
        "<!doctype html><meta charset='utf-8'><title>"
        + html.escape(title)
        + "</title><style>body{font:13px/1.45 Arial,sans-serif;margin:30px;color:#172033}"
        "h1{font-size:20px}h2{font-size:15px;margin-top:28px}.table{overflow:auto}"
        "table{border-collapse:collapse;min-width:100%;white-space:nowrap}th,td{padding:7px 8px;"
        "border:1px solid #cfd7e3;text-align:left;vertical-align:top}th{background:#eef3f8}"
        "p{color:#5c687a}</style><h1>"
        + html.escape(title)
        + "</h1><h2>Main table — biological-type summary</h2><div class='table'>"
        + groups.to_html(index=False, border=0, escape=True)
        + "</div><h2>Supplementary table — all datasets</h2><div class='table'>"
        + supplementary.to_html(index=False, border=0, escape=True)
        + "</div><p>Metric cells show the across-slice median [Q1, Q3]. Captured counts are an H5AD "
        "molecule/library-size proxy. Density and observed spacing use native coordinate units and should "
        "not be ranked across incompatible technologies.</p>",
        encoding="utf-8",
    )
    return {
        "datasets": int(len(data)),
        "biological_types": int(data["group"].nunique(dropna=False)),
        "slices": int(pd.to_numeric(data["n_slices"], errors="coerce").fillna(0).sum()),
        "keep": int(pd.to_numeric(data["keep"], errors="coerce").fillna(0).sum()),
        "exclude": int(pd.to_numeric(data["exclude"], errors="coerce").fillna(0).sum()),
        "all_datasets_csv": str(supplement_path),
        "biological_type_csv": str(group_path),
        "html": str(html_path),
    }


def _write_detail_page(page_payload, report_path):
    """Render already prepared scientific evidence without recomputing it."""
    records = page_payload['records']
    if any(row['final_call'] not in {'keep','exclude'} for row in records):
        raise ValueError('Detailed viewer requires complete binary calls')
    if [str(row['slice_id']) for row in records] != page_payload['order']:
        raise ValueError('Detailed viewer slice order mismatch')
    encoded = json.dumps(page_payload, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    text = _report_template(page_payload['language']).replace('__TITLE__', html.escape(page_payload['title'])).replace('__PAYLOAD__', encoded)
    if page_payload.get('binary_application', {}).get('application_scope') == 'new_input_unvalidated':
        note = ('Frozen-policy application to new input — not independently validated on this input.'
                if page_payload['language']=='en' else '已锁定策略在新输入上的应用结果；本输入未经独立验证。')
        text = text.replace('</header>', '<p style="font-weight:700">'+note+'</p></header>')
    if page_payload.get('binary_application', {}).get('application_scope') == 'experimental_policy':
        note = ('Experimental joint-review policy: evaluated with metric perturbations; independent raw-matrix and biological validation is incomplete. These input datasets are not certified.'
                if page_payload['language']=='en' else '实验性联合细筛策略：已进行指标扰动评估，尚未完成独立原始矩阵及生物学验证；当前输入未经认证。')
        text = text.replace('</header>', '<p style="font-weight:700">'+note+'</p></header>')
    policy = page_payload.get('binary_application', {}).get('policy', {})
    tiers = policy.get('review_exclusion_tiers', [])
    if tiers:
        descriptions = []
        for tier in tiers:
            span = f"[{tier['min_score']:.3f}, {tier['max_score']:.3f})" if tier.get('max_score') is not None else f"score ≥ {tier['min_score']:.3f}"
            label = ('Evidence fine-screen band' if tier.get('enable_exclude') else 'Conservative retention band') if page_payload['language']=='en' else ('证据细筛区间' if tier.get('enable_exclude') else '保守保留区间')
            descriptions.append(html.escape(label + ': ' + span))
        heading = 'Fixed policy across datasets' if page_payload['language']=='en' else '所有数据集共用的锁定策略'
        note = ('These score boundaries were selected in historical calibration and then frozen; they are not re-estimated for each dataset. Every internal review receives a stage-2 band assessment and recorded outcome. Exclusion evidence gates run only in the evidence fine-screen band. Local expectations and detection windows can adapt; these boundaries do not.' if page_payload['language']=='en' else '分数边界经历史校准后锁定，不按当前数据集重新估计。每条内部review均进入第二阶段并记录分层审查结论；仅证据细筛区间执行排除证据门控。局部参照与窗口可以自适应，这些分数边界不会自适应变化。')
        if page_payload.get('binary_application', {}).get('application_scope') == 'experimental_policy':
            note = ('Stage-1 thresholds and the 0.540 boundary retain historical values. Both review bands execute evidence checks; their domain thresholds were jointly selected and frozen before held-out metric-stress evaluation. They are not re-estimated per dataset. This is not biological certification.' if page_payload['language']=='en' else '第一阶段与0.540分界沿用历史值；两个review区间均执行证据检查，域门槛经联合选择后冻结，再进行留出指标扰动评估。不按数据集重新估计，亦不代表生物学认证。')
        block='<section class="section"><h2>'+heading+'</h2><p>'+'; '.join(descriptions)+'</p><p>'+note+'</p></section>'
        en = page_payload['language']=='en'
        k, e = policy['keep_max_score'], policy['exclude_min_score']
        operation_rows = [(f'score ≤ {k:.3f}', 'Stage-1 keep candidate; no review evidence pass' if en else '第一阶段保留候选；不进入review证据细筛', 'keep')]
        for tier in tiers:
            lo, hi = tier['min_score'], tier.get('max_score')
            lower = f'{lo:.3f} < score' if lo == k else f'{lo:.3f} ≤ score'
            span = lower + (f' < {hi:.3f}' if hi is not None else '')
            enabled = tier.get('enable_exclude')
            action = (('Evaluate all exclusion-evidence gates; retain if any fail' if enabled else 'Band/rule assessment and audit only; no further exclusion-evidence checks') if en else ('逐项判断全部排除证据条件；任一失败则保留' if enabled else '只做分层规则核对和审计；不再执行排除证据检查'))
            if enabled:
                action += (f". Domains ≥{tier['min_corroborating_domains']}; strongest domain ≥{tier['severe_domain_threshold']:.2f}; confidence field ≥{tier['min_score_confidence']:.2f}; required window support {tier['min_window_stability']:.0%}" if en else f"。支持维度≥{tier['min_corroborating_domains']}；最强维度≥{tier['severe_domain_threshold']:.2f}；邻域支持字段≥{tier['min_score_confidence']:.2f}；窗口一致性≥{tier['min_window_stability']:.0%}")
            operation_rows.append((('Internal review: ' if en else '内部review：')+span, action, 'exclude / keep' if enabled else 'keep'))
        operation_rows.extend([(f'score ≥ {e:.3f}; '+('direct gate passes' if en else '直接排除门控通过'), 'Direct exclusion; no stage-2 pass' if en else '第一阶段直接排除；不进第二阶段', 'exclude'),(f'score ≥ {e:.3f}; '+('direct gate fails' if en else '直接排除门控未通过'), 'Downgrade to review and test all evidence-band gates' if en else '降级review，执行证据细筛区间的全部检查', 'exclude / keep')])
        table='<div class="table-wrap"><table class="data-table"><thead><tr><th>'+('Score and route' if en else '分数与路径')+'</th><th>'+('Actual operation' if en else '实际操作')+'</th><th>'+('Final action' if en else '最终动作')+'</th></tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in row)+'</tr>' for row in operation_rows)+'</tbody></table></div>'
        block=block.replace('</section>', table+'</section>')
        text=text.replace('<section id="methods" class="panel">','<section id="methods" class="panel">'+block)
    if page_payload.get('preregistration'):
        extension = Path(__file__).with_name('slice_quality_detail_roi.js').read_text()
        text = text.replace('</body>', '<script>'+extension+'</script></body>')
    Path(report_path).write_text(text, encoding='utf-8')


def render_binary_slice_quality_appendix(
    input_dir: Union[str, Path],
    output_dir: Union[str, Path],
    *,
    title: str = "Spateo referee · serial-slice QC",
    display_window: int = 5,
    max_points_per_slice: int = 900,
    language: str = "zh",
    source_h5ad: Optional[Union[str, Path]] = None,
    slice_key: Optional[str] = None,
    spatial_key: Optional[str] = None,
    total_counts_key: Optional[str] = None,
    n_genes_key: Optional[str] = None,
    mito_key: Optional[str] = None,
) -> dict[str, Any]:
    """Build a self-contained user report and publication-oriented tables."""
    if display_window not in {3, 5}:
        raise ValueError("display_window must be 3 or 5")
    if language not in {"zh", "en"}:
        raise ValueError("language must be 'zh' or 'en'")
    source = Path(input_dir).expanduser().resolve()
    dataset_name = source.name
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    audit_path = source / "slice_quality_binary_audit.csv"
    payload_path = source / "slice_quality_display_payload.json"
    manifest_path = source / "slice_quality_manifest.json"
    if (
        not audit_path.exists()
        or not payload_path.exists()
        or not manifest_path.exists()
    ):
        raise FileNotFoundError(
            "Detailed report requires binary audit, display payload and manifest"
        )
    metrics = (
        pd.read_csv(audit_path, dtype={"slice_id": str})
        .sort_values("slice_index", kind="stable")
        .reset_index(drop=True)
    )
    calls = metrics["final_call"].astype("string")
    if calls.isna().any() or not calls.isin(["keep", "exclude"]).all():
        raise ValueError(
            "Detailed user report requires a complete keep/exclude final_call for every slice"
        )
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    order = metrics["slice_id"].astype(str).tolist()
    spacing = pd.to_numeric(metrics.get("median_nn_distance"), errors="coerce")
    valid_spacing = spacing[np.isfinite(spacing) & (spacing > 0)]
    bandwidth = float(np.median(valid_spacing) * 2.0) if len(valid_spacing) else 1.0
    points, coordinates = _prepare_kde_points(
        payload.get("point_samples", {}), order, bandwidth, max_points_per_slice
    )
    evidence_source: dict[str, Any] = {
        "available": False,
        "absolute_expression_available": False,
        "component_evidence_available": False,
    }
    explicit_source = source_h5ad is not None
    resolved_h5ad, resolved_slice_key, resolved_spatial_key = _manifest_source_contract(
        manifest, source_h5ad, slice_key, spatial_key
    )
    if payload.get("preregistration"):
        points, coordinates, evidence_source = _prepare_preregistered_evidence_points(
            source, payload, order, bandwidth, max_points_per_slice
        )
    elif resolved_h5ad is not None and resolved_h5ad.exists():
        if not resolved_slice_key or not resolved_spatial_key:
            raise ValueError(
                "Source-H5AD evidence requires resolvable slice and spatial keys"
            )
        points, coordinates, evidence_source = _prepare_source_evidence_points(
            resolved_h5ad,
            order,
            bandwidth,
            max_points_per_slice,
            slice_key=resolved_slice_key,
            spatial_key=resolved_spatial_key,
            total_counts_key=total_counts_key,
            n_genes_key=n_genes_key,
            mito_key=mito_key,
        )
    elif explicit_source:
        raise FileNotFoundError(f"Source H5AD not found: {resolved_h5ad}")

    explanations = build_directional_slice_explanations(
        metrics, language=language
    ).set_index("slice_id")
    regions: dict[str, Any] = {}
    for slice_id in metrics.loc[calls.eq("exclude"), "slice_id"].astype(str):
        if payload.get("preregistration"):
            info = payload["preregistration"]["slices"][slice_id]
            # Avoid presenting uncertain registration as a local dropout finding.
            neighbors = order[max(0,order.index(slice_id)-2):order.index(slice_id)+3]
            if any(payload['preregistration']['slices'][key]['status'] not in {'reference','estimated'}
                   or payload['preregistration']['slices'][key]['warnings'] for key in neighbors):
                continue
        region = _largest_density_deficit_region(
            slice_id, order, coordinates, bandwidth, language=language
        )
        if region is not None:
            regions[slice_id] = region

    allowed_columns = {
        "slice_id",
        "slice_index",
        "final_call",
        "decision_basis",
        "review_resolution_tier",
        "review_tier_window_stability",
        "quality_anomaly_score",
        "density_domain_score",
        "expression_domain_score",
        "damage_domain_score",
        "continuity_domain_score",
        "n_locations",
        "n_locations_expected",
        "cell_density",
        "cell_density_expected",
        "median_n_genes",
        "median_n_genes_expected",
        "median_total_counts",
        "median_total_counts_expected",
        "median_nn_distance",
        "detected_genes",
        "largest_component_fraction",
        "largest_component_fraction_expected",
        "fragmentation",
        "fragmentation_expected",
        "hole_fraction",
        "hole_fraction_expected",
        "knn_tail_ratio",
        "knn_tail_ratio_expected",
        "detected_gene_fraction",
        "detected_gene_fraction_expected",
        "library_complexity",
        "library_complexity_expected",
        "zero_fraction",
        "zero_fraction_expected",
        "local_low_depth_fraction",
        "local_low_depth_fraction_expected",
        "regional_low_depth_cluster_fraction",
        "regional_low_depth_cluster_fraction_expected",
        "median_pct_mito",
        "median_pct_mito_expected",
        "partial_structure_protection",
        "window_context",
        "adaptive_windows_tested",
        "adaptive_window_stability",
        "adaptive_exclude_call_fraction",
        "adaptive_min_score_across_windows",
        "adaptive_median_score_across_windows",
        "adaptive_exclusion_gate",
        "standard_exclusion_gate",
    }
    records: list[dict[str, Any]] = []
    for record in metrics.to_dict("records"):
        slice_id = str(record["slice_id"])
        output = {
            key: (None if pd.isna(value) else value)
            for key, value in record.items()
            if key in allowed_columns
        }
        explanation = explanations.loc[slice_id]
        output["user_reason"] = explanation["user_reason"]
        if record.get('application_scope') == 'new_input_unvalidated':
            output['user_reason'] = output['user_reason'].replace('independently validated high-specificity', 'historically validated, frozen high-specificity').replace('经独立验证的高特异度', '历史验证并锁定的高特异度')
            output['user_reason'] += (' New-input application; no independent validation of this input.' if language == 'en' else ' 本结论为新输入应用结果，未对本输入独立验证。')
        if record.get('application_scope') == 'experimental_policy':
            output['user_reason'] = output['user_reason'].replace('independently validated high-specificity', 'frozen experimental-policy').replace('entered the calibrated', 'entered the experimental').replace('经独立验证的高特异度', '实验性锁定策略的').replace('校准分层', '实验分层')
            output['user_reason'] += (' Experimental policy application: metric-stress evidence only; this input is not independently validated.' if language == 'en' else ' 实验策略应用：仅有指标扰动验证依据，本输入未经独立验证。')
        output["primary_reason"] = explanation["primary_reason"]
        output["likely_causes"] = explanation["likely_causes"]
        output["directional_evidence"] = explanation["directional_evidence"]
        output["visual_findings"] = _visual_panel_findings(record, language)
        records.append(output)

    statistics = []
    for definition in _REPORT_METRICS:
        key = definition["key"]
        if key not in metrics:
            continue
        summary = _summary_stat(metrics[key])
        if summary["median"] is not None:
            statistics.append({**definition, **summary})

    page_payload = {
        "title": title,
        "preregistration": payload.get("preregistration"),
        "binary_application": json.loads((source / "binary_policy_application.json").read_text()) if (source / "binary_policy_application.json").exists() else {},
        "records": records,
        "order": order,
        "points": points,
        "regions": regions,
        "statistics": statistics,
        "exclusion_evidence": _exclusion_evidence_summary(metrics, language),
        "evidence_source": evidence_source,
        "kde_bandwidth": bandwidth,
        "display_window": display_window,
        "language": language,
        "selected_detection_window": manifest.get("selected_window"),
        "measurement_note": (
            (
                "Captured counts are an H5AD molecule-capture/library-size proxy and should be described "
                "as UMIs only when UMI semantics are explicit. Observed spot spacing is a center-to-center "
                "spacing proxy in native coordinates, not the platform's nominal resolution. Point-level "
                "expression panels use absolute observed values and one shared scale per slice window."
            )
            if language == "en"
            else (
                "Captured counts 是 H5AD 中的分子捕获/文库量代理；只有数据明确为 UMI 时才可解释为 UMI。"
                "Observed spot spacing 是原生坐标中的中心间距代理，不是平台标称分辨率。"
                "逐点表达图使用绝对观测值，并在同一切片窗口内共用色标。"
            )
        ),
        "source_note": (
            (
                "KDE uses one fixed bandwidth across the dataset and a shared color scale within each "
                "three- or five-slice window. The yellow box marks the largest contiguous region with high "
                "neighboring-slice KDE expectation but a marked focal-slice deficit."
            )
            if language == "en"
            else (
                "KDE 使用全数据集统一固定带宽，并在每个 3–5 张窗口内使用共同色标。"
                "黄色框表示邻片 KDE 期望较高而焦点切片明显不足的最大连续候选区。"
            )
        ),
        "source": manifest.get("sources", []),
    }
    report_path = destination / "index.html"
    _write_detail_page(page_payload, report_path)
    table_outputs = _write_paper_tables(
        metrics,
        destination,
        manifest.get("selected_window"),
        dataset_name=dataset_name,
        title=title,
        language=language,
    )
    summary_path = destination / "report_summary.json"
    summary = {
        "dataset": dataset_name,
        "slices": int(len(metrics)),
        "keep": int(calls.eq("keep").sum()),
        "exclude": int(calls.eq("exclude").sum()),
        "review": 0,
        "adaptive_resolved": int(
            metrics.get("decision_basis", pd.Series(dtype=str))
            .astype(str)
            .isin(["adaptive_multidomain_resolution", "tiered_multidomain_resolution"])
            .sum()
        ),
        "language": language,
        "kde_bandwidth_native_units": bandwidth,
        "display_points": int(sum(item["displayed"] for item in points.values())),
        "available_points": int(sum(item["available"] for item in points.values())),
        "density_deficit_regions": int(len(regions)),
        "renderer": "slice_quality_visualization.render_binary_slice_quality_appendix",
        "application_scope": page_payload["binary_application"].get("application_scope", "existing_binary_audit"),
        "preregistration_frame_id": payload.get("preregistration", {}).get("frame_id"),
        "evidence_source": evidence_source,
        "absolute_expression_maps": bool(
            evidence_source.get("absolute_expression_available")
        ),
        "connected_component_maps": bool(
            evidence_source.get("component_evidence_available")
        ),
        "report": str(report_path),
        **table_outputs,
    }
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    summary["summary_json"] = str(summary_path)
    return summary


_REPORT_TEMPLATE = r"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
:root{--bg:#f6f8fb;--surface:#fff;--ink:#172033;--muted:#627086;--line:#dce3ed;--navy:#173a5e;--blue:#2f67d8;--cyan:#42b8d4;--yellow:#ffd54a;--orange:#f37a49;--red:#cf3f4f;--green:#14805e;--plot:#07111f}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans SC",sans-serif}header{background:linear-gradient(125deg,#132b47,#285e86);color:#fff;padding:24px clamp(18px,4vw,52px)}header h1{font-size:26px;margin:0 0 5px;font-weight:650}header p{margin:0;color:#d8e6f2}.wrap{width:min(1500px,100%);margin:auto;padding:18px clamp(12px,3vw,34px) 42px}.tabs{display:flex;gap:4px;overflow:auto;border-bottom:1px solid var(--line);margin-bottom:18px}.tab{border:0;background:transparent;color:#536278;padding:10px 15px;font:inherit;white-space:nowrap;cursor:pointer;border-bottom:3px solid transparent}.tab[aria-selected="true"]{color:var(--blue);border-bottom-color:var(--blue);font-weight:650}.panel{display:none}.panel.active{display:block}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:10px;margin-bottom:18px}.card{background:var(--surface);border:1px solid var(--line);border-radius:10px;padding:13px 15px}.card b{font-size:23px;display:block}.muted{color:var(--muted)}h2{font-size:18px;margin:4px 0 12px}h3{font-size:15px;margin:0 0 8px}.section{background:var(--surface);border:1px solid var(--line);border-radius:11px;padding:15px;margin-bottom:14px}.call{display:inline-block;border-radius:999px;color:#fff;font-weight:700;padding:2px 8px;text-transform:uppercase;font-size:12px}.keep{background:var(--green)}.exclude{background:var(--red)}.toolbar{display:flex;align-items:flex-end;gap:10px;flex-wrap:wrap;margin-bottom:12px}.field{display:grid;gap:4px;color:var(--muted)}select,input{font:inherit;border:1px solid #c9d3e0;border-radius:7px;padding:7px 9px;background:#fff;color:var(--ink)}button{font:inherit}.plots{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:12px}.plot{min-width:0}.plot svg{width:100%;height:205px;display:block}.plot-head,.evidence-head{display:flex;justify-content:space-between;gap:14px;align-items:flex-start;flex-wrap:wrap}.metric{color:var(--blue);font-weight:700}.window{display:grid;grid-template-columns:repeat(auto-fit,minmax(185px,1fr));gap:9px}.slice-panel{min-width:0}.slice-panel canvas{width:100%;aspect-ratio:1.5/1;display:block;background:var(--plot);border-radius:7px}.slice-cap{text-align:center;color:var(--muted);margin-top:4px;overflow-wrap:anywhere;font-size:12px}.slice-cap.focus{color:var(--ink);font-weight:700}.legend{display:flex;align-items:center;gap:8px;margin:9px 0;flex-wrap:wrap}.ramp{width:min(260px,55vw);height:10px;border-radius:999px;background:linear-gradient(90deg,#243b82,#18a7c4,#f1dd3f,#f06432)}.component-key{display:flex;gap:12px;align-items:center;flex-wrap:wrap}.component-key i{width:10px;height:10px;border-radius:50%;display:inline-block;margin-right:4px}.finding{max-width:760px;padding:7px 10px;border-radius:7px;font-size:12px}.finding.supports{background:#fff0f1;color:#8b2633;border-left:3px solid var(--red)}.finding.context{background:#eff4f8;color:#536278;border-left:3px solid #9aaabd}.reason-box{border-left:4px solid var(--red);padding:10px 13px;background:#fff6f6}.reason-box.keep-reason{border-left-color:var(--green);background:#f1fbf7}.visual-links{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:9px;margin-top:11px}.visual-link{padding:9px 10px;border:1px solid var(--line);border-radius:8px}.visual-link.supports{border-color:#e9a5ae;background:#fff7f7}.visual-link b{display:block;margin-bottom:3px}.domain-bars{display:grid;gap:7px;margin-top:12px}.domain-row{display:grid;grid-template-columns:110px minmax(120px,1fr) 44px;gap:8px;align-items:center;font-size:12px}.domain-track,.frequency-track{height:10px;background:#e8edf3;border-radius:999px;overflow:hidden;position:relative}.domain-track:after{content:"";position:absolute;left:45%;top:0;bottom:0;border-left:1px dashed #6f7f93}.domain-track i,.frequency-track i{display:block;height:100%;background:var(--blue);border-radius:999px}.domain-track i.strong{background:var(--red)}.frequency-list{display:grid;gap:12px;max-width:1050px}.frequency-row{display:grid;grid-template-columns:minmax(190px,260px) minmax(180px,1fr) 70px minmax(170px,1fr);gap:11px;align-items:center;border:0;background:transparent;text-align:left;color:inherit;font:inherit;padding:3px 0;cursor:pointer}.frequency-row:hover .frequency-track{outline:2px solid #b9cae1}.frequency-track{height:17px}.frequency-track i{background:var(--red)}.frequency-count{font-weight:700;color:var(--red)}.evidence-layout{display:grid;grid-template-columns:minmax(0,1.15fr) minmax(280px,.85fr);gap:13px}.evidence-table,.data-table{width:100%;border-collapse:collapse}.evidence-table th,.evidence-table td,.data-table th,.data-table td{padding:7px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}.evidence-table th,.data-table th{background:#eef3f8}.evidence-table td.num,.data-table td.num{text-align:right;font-variant-numeric:tabular-nums}.table-wrap{overflow:auto;max-height:660px}.data-table{font-size:12px}.data-table th{position:sticky;top:0;z-index:2}.data-table tr{cursor:pointer}.data-table tr:hover{background:#f0f5fb}.reason-cell{min-width:360px;white-space:normal}.summary-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:12px}.exclude-item{border-top:1px solid var(--line);padding:10px 0}.exclude-item:first-child{border-top:0}.note{font-size:12px;color:var(--muted)}.downloads{display:flex;gap:8px;flex-wrap:wrap}.downloads a{display:inline-block;text-decoration:none;border:1px solid #b9c7d8;border-radius:7px;padding:7px 10px;color:#245792;background:#fff}.methods{columns:2 350px;column-gap:30px}.methods section{break-inside:avoid;margin-bottom:18px}
@media(max-width:900px){header h1{font-size:22px}.evidence-layout{grid-template-columns:1fr}.plots{grid-template-columns:1fr}.window{grid-template-columns:repeat(2,minmax(0,1fr))}.visual-links{grid-template-columns:1fr}.frequency-row{grid-template-columns:minmax(150px,220px) 1fr 60px}.frequency-slices{grid-column:1/-1}.reason-cell{min-width:260px}.methods{columns:1}}
@media(max-width:520px){.window{grid-template-columns:1fr}.section{padding:12px}.plot svg{height:185px}.frequency-row{grid-template-columns:1fr 60px}.frequency-row>:first-child,.frequency-slices{grid-column:1/-1}.domain-row{grid-template-columns:90px 1fr 40px}}
</style></head><body>
<header><h1 id="pageTitle"></h1><p>单数据集开发附件 · 全部切片仅保留 keep / exclude · 不需要用户再次判定</p></header>
<main class="wrap">
<nav class="tabs" aria-label="报告页面">
<button class="tab" data-tab="summary" aria-selected="true">结果总览</button>
<button class="tab" data-tab="statistics" aria-selected="false">统计趋势</button>
<button class="tab" data-tab="density" aria-selected="false">切片证据图</button>
<button class="tab" data-tab="slices" aria-selected="false">全部切片</button>
<button class="tab" data-tab="methods" aria-selected="false">方法与口径</button>
</nav>
<section id="summary" class="panel active"><div id="cards" class="cards"></div><div class="summary-grid"><div class="section"><h2>自动排除切片</h2><div id="excludeSummary"></div></div><div class="section"><h2>读图要点</h2><p>红色代表最终自动排除；绿色代表最终保留。每张切片的低质量原因会直接对应到 KDE 密度、表达捕获和连通组织成分三组图。</p><p>证据图在同一页面并排展示连续切片。每组图上方明确标记该证据是否支持排除，避免仅凭颜色主观判断。</p><p class="note" id="measurementNote"></p><div class="downloads"><a href="paper_dataset_summary.csv">文章数据集汇总表 CSV</a><a href="paper_excluded_slices.csv">文章排除切片表 CSV</a><a href="paper_tables.html">文章表格预览</a></div></div></div></section>
<section id="statistics" class="panel"><div class="section"><h2>逐切片统计趋势</h2><p class="muted">五项指标同时展示，虚线为本数据集切片中位数；红色方形为 exclude，绿色圆点为 keep。点击数据点可进入对应证据页面。</p><div id="plots" class="plots"></div></div><div class="section"><h2>排除切片的关键证据</h2><p class="muted">每个横条表示最终 exclude 切片中有多少张出现该类显著不利证据。证据可以重叠，因此横条数量之和不等于排除切片总数。点击横条可查看其中第一张切片。</p><div id="exclusionEvidenceChart" class="frequency-list"></div></div></section>
<section id="density" class="panel">
<div class="section"><h2>切片低质量结论与图像证据</h2><div class="toolbar"><label class="field">焦点切片<select id="sliceSelect"></select></label><label class="field">展示窗口<select id="windowSelect"><option value="5">5 张连续切片</option><option value="3">3 张连续切片</option></select></label><label class="field">表达指标<select id="expressionSelect"><option value="log_captured_counts">Captured counts / spot</option><option value="n_genes">Detected genes / spot</option></select></label><label class="field">空间范围<select id="extentSelect"><option value="shared">窗口共同坐标范围</option><option value="individual">每张切片独立放大</option></select></label></div><div id="reason" class="reason-box"></div><div id="domainBars" class="domain-bars"></div><div id="visualLinks" class="visual-links"></div></div>
<div class="section"><div class="evidence-head"><h2>A · KDE 局部细胞/位置密度</h2><div id="kdeFinding" class="finding"></div></div><div id="kdeWindow" class="window"></div><div class="legend"><span>低</span><span class="ramp"></span><span>高</span><span id="kdeLegend" class="muted"></span></div><p id="regionNote" class="note"></p><p class="note" id="sourceNote"></p></div>
<div class="section"><div class="evidence-head"><h2>B · 表达捕获</h2><div id="expressionFinding" class="finding"></div></div><div id="expressionWindow" class="window"></div><div class="legend"><span>低</span><span class="ramp"></span><span>高</span><span id="expressionLegend" class="muted"></span></div><p class="note">连续切片使用共同的绝对观测值色标；该图显示 H5AD 捕获量，不代表原始 read depth。</p></div>
<div class="section"><div class="evidence-head"><h2>C · 连通组织成分</h2><div id="componentFinding" class="finding"></div></div><div id="componentWindow" class="window"></div><div class="legend component-key"><span><i style="background:#22b7c7"></i>最大连通成分</span><span><i style="background:#f27b4d"></i>第二成分</span><span><i style="background:#b36ac8"></i>其他成分</span></div><p class="note">颜色表示每张切片内部的连通成分排名；小型分离组织块会以不同颜色显示。</p></div>
<div class="section"><h2>焦点切片的方向性指标证据</h2><div class="evidence-layout"><p class="muted">表格给出观测值与局部连续切片期望值。图像负责空间定位，表格负责说明变化方向与幅度。</p><div class="table-wrap"><table class="evidence-table"><thead><tr><th>指标</th><th>方向</th><th>观测</th><th>局部期望</th></tr></thead><tbody id="evidenceRows"></tbody></table></div></div></div>
</section>
<section id="slices" class="panel"><div class="section"><h2>全部切片最终结果</h2><div class="toolbar"><label class="field">状态<select id="callFilter"><option value="all">全部</option><option value="exclude">Exclude</option><option value="keep">Keep</option></select></label><label class="field">切片搜索<input id="sliceSearch" placeholder="例如 slices_4"></label></div><div class="table-wrap"><table class="data-table"><thead><tr><th>#</th><th>Slice</th><th>最终状态</th><th>Score</th><th>位置数</th><th>Cell density</th><th>Genes / spot</th><th>Counts / spot</th><th>Observed spacing</th><th>Total genes</th><th>结论及原因</th></tr></thead><tbody id="sliceRows"></tbody></table></div></div></section>
<section id="methods" class="panel"><div class="section methods"><section><h2>检测流程</h2><p>每张切片与前后连续切片比较，默认检测窗口由 3/5/7 候选中选择；同时使用更宽的稳健趋势防止连续坏片互相掩护。密度/完整性、表达捕获、损伤及跨切片连续性四个证据域共同形成异常分数。</p></section><section><h2>最终二元规则</h2><p>方法分为两阶段：第一阶段按锁定阈值形成 keep / review / exclude 内部初筛，高分但未通过保护条件的切片降级为 review；第二阶段逐条审查全部内部 review，完成锁定分层归属、排除权限核查并记录结论。仅保留层据此给出 keep，不执行该层排除证据门控；只有启用排除的层才检查多域、严重异常、双侧上下文、结构保护、置信度及 3/5/7 窗口稳定性。启用层全部条件通过才转为 exclude，否则为 keep。公开结果只有 keep / exclude。</p></section><section><h2>关键证据条形图</h2><p>横条统计 exclude 切片中某类不利方向指标异常达到 0.50 的切片数；同一切片可同时支持多个证据类别，因此横条数值不能相加。</p></section><section><h2>KDE 图如何生成</h2><p id="kdeMethod"></p></section><section><h2>表达图如何生成</h2><p>Captured counts 和 detected genes 直接来自 H5AD 的逐点观测列；显示时在当前连续切片窗口内使用共同的 4%–96% 色标，不对每张切片单独归一化。</p></section><section><h2>连通成分如何生成</h2><p>在每张切片内以 3.25 倍中位最近邻距离连接空间点，随后计算无向连通成分；颜色表示按点数排序的成分排名，与检测器的 fragmentation 和 largest-component fraction 定义一致。</p></section><section><h2>原因如何解释</h2><p>“比局部窗口期望低/高”描述的是相邻切片趋势中的相对变化，不是跨技术平台阈值。可能原因只用于实验排查，不能从 H5AD 单独确诊。</p></section><section><h2>表达深度限制</h2><p>H5AD 可观察每点捕获计数、检出基因和复杂度，但没有原始 reads、PCR duplication、mapping rate 或 saturation 时，不应把 captured-count proxy 直接称为原始测序深度。</p></section><section><h2>空间分辨率限制</h2><p>本报告使用 observed nearest-neighbour spacing 作为分辨率相关的采样间距代理，单位来自原生坐标；它不是平台标称的 spot diameter 或分辨率。</p></section></div></section>
</main>
<script id="payload" type="application/json">__PAYLOAD__</script>
<script>
const D=JSON.parse(document.getElementById('payload').textContent),R=D.records,byId=Object.fromEntries(R.map(r=>[String(r.slice_id),r])),EN=D.language==='en';
const colors={keep:'#14805e',exclude:'#cf3f4f'},componentColors=['#22b7c7','#f27b4d','#b36ac8','#e2bd36','#6687d9','#94a0ad'];
const T=EN?{all:'All slices',selectedWindow:'Selected detection window',adaptiveResolved:'Review-tier excludes',causes:'Plausible causes: ',noneExcluded:'No slice was automatically excluded in this dataset.',density:'Density / amount',expression:'Expression capture',components:'Connected components',strong:'strong evidence ≥ 0.45',noRegion:'No stable contiguous KDE density-deficit region was identified; the call is supported by the other evidence shown below.',noEvidence:'No directional anomaly reached the display threshold',unavailable:'Point-level evidence unavailable',locations:'locations',largest:'largest component',componentsN:'components',counts:'captured counts / spot',genes:'detected genes / spot'}:{all:'全部切片',selectedWindow:'自动选择检测窗口',adaptiveResolved:'分层重审排除',causes:'可能原因：',noneExcluded:'本数据集没有自动排除切片。',density:'密度 / 组织量',expression:'表达捕获',components:'连通组织成分',strong:'强证据阈值 ≥ 0.45',noRegion:'该切片未形成稳定的最大连续 KDE 密度缺失候选区；结论来自下方其他证据。',noEvidence:'未见达到展示阈值的方向性异常',unavailable:'逐点证据不可用',locations:'位置',largest:'最大成分',componentsN:'成分数',counts:'每点捕获计数',genes:'每点检出基因'};
const fmt=x=>{if(x===null||x===undefined)return 'NA';x=Number(x);if(!Number.isFinite(x))return 'NA';let a=Math.abs(x);if(a&&a<.01)return x.toExponential(2);if(a>=1e6)return x.toExponential(2);return new Intl.NumberFormat(EN?'en-US':'zh-CN',{maximumFractionDigits:a<10?3:a<100?1:0}).format(x)},pct=x=>x!==null&&x!==undefined&&Number.isFinite(Number(x))?`${Math.abs(Number(x)*100).toFixed(1)}%`:'NA';
const q=(values,p)=>{let a=values.filter(Number.isFinite).sort((x,y)=>x-y);if(!a.length)return 0;let i=(a.length-1)*p,l=Math.floor(i),h=Math.ceil(i);return a[l]+(a[h]-a[l])*(i-l)},ramp=t=>{t=Math.max(0,Math.min(1,t));let stops=[[36,59,130],[24,167,196],[241,221,63],[240,100,50]],u=t*(stops.length-1),i=Math.min(stops.length-2,Math.floor(u)),f=u-i,a=stops[i],b=stops[i+1];return `rgb(${a.map((v,k)=>Math.round(v+(b[k]-v)*f)).join(',')})`};
document.getElementById('pageTitle').textContent=D.title;document.getElementById('measurementNote').textContent=D.measurement_note;document.getElementById('sourceNote').textContent=D.source_note;document.getElementById('kdeMethod').textContent=EN?`Two-dimensional Gaussian KDE uses a fixed bandwidth of ${fmt(D.kde_bandwidth)} native coordinate units. A shared scale is used within each slice window. The page displays ${Object.values(D.points).reduce((a,p)=>a+p.displayed,0).toLocaleString()} / ${Object.values(D.points).reduce((a,p)=>a+p.available,0).toLocaleString()} available points.`:`二维 Gaussian KDE 使用统一固定带宽 ${fmt(D.kde_bandwidth)} 个原生坐标单位；窗口内共同色标便于直接比较。页面展示 ${Object.values(D.points).reduce((a,p)=>a+p.displayed,0).toLocaleString()} / ${Object.values(D.points).reduce((a,p)=>a+p.available,0).toLocaleString()} 个可用点。`;
document.querySelectorAll('.tab').forEach(button=>button.onclick=()=>{document.querySelectorAll('.tab').forEach(x=>x.setAttribute('aria-selected',String(x===button)));document.querySelectorAll('.panel').forEach(x=>x.classList.toggle('active',x.id===button.dataset.tab));if(button.dataset.tab==='density')requestAnimationFrame(drawEvidence)});
const keeps=R.filter(r=>r.final_call==='keep').length,excludes=R.length-keeps,adaptiveResolved=R.filter(r=>['adaptive_multidomain_resolution','tiered_multidomain_resolution'].includes(r.decision_basis)).length;document.getElementById('cards').innerHTML=`<div class="card"><b>${R.length}</b><span class="muted">${T.all}</span></div><div class="card"><b style="color:${colors.keep}">${keeps}</b><span class="muted">Keep</span></div><div class="card"><b style="color:${colors.exclude}">${excludes}</b><span class="muted">Exclude</span></div><div class="card"><b>${adaptiveResolved}</b><span class="muted">${T.adaptiveResolved}</span></div><div class="card"><b>${D.selected_detection_window}</b><span class="muted">${T.selectedWindow}</span></div>`;
document.getElementById('excludeSummary').innerHTML=R.filter(r=>r.final_call==='exclude').map(r=>`<div class="exclude-item"><span class="call exclude">exclude</span> <b>${r.slice_id}</b> · score ${fmt(r.quality_anomaly_score)}<br><span>${r.primary_reason}</span><br><span class="note">${T.causes}${r.likely_causes}</span></div>`).join('')||`<p>${T.noneExcluded}</p>`;
function showEvidence(id){document.getElementById('sliceSelect').value=id;document.querySelector('[data-tab="density"]').click();drawEvidence()}
function makePlots(){let container=document.getElementById('plots');container.innerHTML=D.statistics.map(s=>{let vals=R.map((r,i)=>({r,i,v:r[s.key]==null?NaN:Number(r[s.key])})).filter(o=>Number.isFinite(o.v)),lo=Math.min(...vals.map(o=>o.v)),hi=Math.max(...vals.map(o=>o.v));if(!(hi>lo)){lo-=1;hi+=1}let w=520,h=205,l=52,rr=12,t=18,b=34,iw=w-l-rr,ih=h-t-b,sx=i=>l+i*iw/Math.max(R.length-1,1),sy=v=>t+(hi-v)*ih/(hi-lo),grid=[lo,(lo+hi)/2,hi].map(v=>`<line x1="${l}" x2="${w-rr}" y1="${sy(v)}" y2="${sy(v)}" stroke="#dfe6ef"/><text x="2" y="${sy(v)+4}" fill="#66758a" font-size="11">${fmt(v)}</text>`).join(''),median=`<line x1="${l}" x2="${w-rr}" y1="${sy(s.median)}" y2="${sy(s.median)}" stroke="#2f67d8" stroke-dasharray="5 4"/>`,line=`<polyline points="${vals.map(o=>`${sx(o.i)},${sy(o.v)}`).join(' ')}" fill="none" stroke="#8294aa" stroke-width="1.2"/>`,dots=vals.map(o=>o.r.final_call==='exclude'?`<rect data-id="${o.r.slice_id}" x="${sx(o.i)-4}" y="${sy(o.v)-4}" width="8" height="8" fill="${colors.exclude}"><title>${o.r.slice_id}: ${fmt(o.v)} ${s.unit}</title></rect>`:`<circle data-id="${o.r.slice_id}" cx="${sx(o.i)}" cy="${sy(o.v)}" r="3.2" fill="${colors.keep}"><title>${o.r.slice_id}: ${fmt(o.v)} ${s.unit}</title></circle>`).join('');return `<div class="section plot"><div class="plot-head"><h3>${s.label}</h3><span class="metric">median ${fmt(s.median)}</span></div><svg role="img" aria-label="${s.label} across slices" viewBox="0 0 ${w} ${h}" preserveAspectRatio="none">${grid}${median}${line}${dots}<text x="${l}" y="${h-7}" fill="#66758a" font-size="11">slice order →</text></svg><div class="note">${s.unit} · IQR ${fmt(s.q25)}–${fmt(s.q75)}</div></div>`}).join('');container.querySelectorAll('[data-id]').forEach(n=>n.onclick=()=>showEvidence(n.dataset.id))}
function makeEvidenceFrequency(){let max=Math.max(1,...D.exclusion_evidence.map(x=>x.count));document.getElementById('exclusionEvidenceChart').innerHTML=D.exclusion_evidence.map(x=>`<button type="button" class="frequency-row" data-id="${x.slice_ids[0]}" aria-label="${x.label}: ${x.count} of ${excludes} excluded slices"><b>${x.label}</b><span class="frequency-track"><i style="width:${100*x.count/max}%"></i></span><span class="frequency-count">${x.count} / ${excludes}</span><span class="frequency-slices note">${x.slice_ids.join(', ')}</span></button>`).join('');document.querySelectorAll('.frequency-row').forEach(n=>n.onclick=()=>showEvidence(n.dataset.id))}
makePlots();makeEvidenceFrequency();
const select=document.getElementById('sliceSelect');select.innerHTML=R.map(r=>`<option value="${r.slice_id}">${r.slice_id} · ${r.final_call}</option>`).join('');let firstExclude=R.find(r=>r.final_call==='exclude');select.value=firstExclude?firstExclude.slice_id:R[0].slice_id;document.getElementById('windowSelect').value=String(D.display_window);['sliceSelect','windowSelect','expressionSelect','extentSelect'].forEach(id=>document.getElementById(id).onchange=drawEvidence);
function currentIds(){let id=select.value,index=D.order.indexOf(id),width=Number(document.getElementById('windowSelect').value),half=Math.floor(width/2),start=Math.max(0,Math.min(index-half,D.order.length-width));return D.order.slice(start,Math.min(D.order.length,start+width))}
function sharedExtent(ids){if(document.getElementById('extentSelect').value==='individual')return null;let xs=ids.flatMap(id=>D.points[id]?.x||[]),ys=ids.flatMap(id=>D.points[id]?.y||[]);return xs.length?[Math.min(...xs),Math.max(...xs),Math.min(...ys),Math.max(...ys)]:null}
function drawCanvas(canvas,id,extent,mode,vlo,vhi,focus,region){let p=D.points[id],ratio=Math.min(window.devicePixelRatio||1,1.35),w=canvas.clientWidth,h=canvas.clientHeight;canvas.width=Math.round(w*ratio);canvas.height=Math.round(h*ratio);let c=canvas.getContext('2d');c.setTransform(ratio,0,0,ratio,0,0);c.fillStyle='#07111f';c.fillRect(0,0,w,h);if(!p||!p.x.length)return;let e=extent||[Math.min(...p.x),Math.max(...p.x),Math.min(...p.y),Math.max(...p.y)],sx=(w-18)/Math.max(e[1]-e[0],1e-9),sy=(h-18)/Math.max(e[3]-e[2],1e-9),s=Math.min(sx,sy),ox=(w-(e[1]-e[0])*s)/2,oy=(h-(e[3]-e[2])*s)/2,values=mode==='kde'?p.kde:p[mode]||[];for(let i=0;i<p.x.length;i++){let v=values[i]==null?NaN:Number(values[i]),fill='#8793a3';if(mode==='component_rank'&&Number.isFinite(v)){fill=componentColors[Math.min(Math.max(Math.round(v)-1,0),componentColors.length-1)]}else if(Number.isFinite(v)){fill=ramp((v-vlo)/Math.max(vhi-vlo,1e-12))}c.fillStyle=fill;c.globalAlpha=.9;c.beginPath();c.arc(ox+(p.x[i]-e[0])*s,h-oy-(p.y[i]-e[2])*s,Math.max(1.4,Math.min(2.45,390/Math.max(p.x.length,160))),0,Math.PI*2);c.fill()}c.globalAlpha=1;if(mode==='kde'&&region){let x=ox+(region.x0-e[0])*s,y=h-oy-(region.y1-e[2])*s,rw=(region.x1-region.x0)*s,rh=(region.y1-region.y0)*s;c.strokeStyle='#ffd54a';c.lineWidth=focus?2.3:1.2;c.setLineDash(focus?[]:[5,4]);c.strokeRect(x,y,rw,rh);c.setLineDash([])}}
function panelCaption(id,mode,focus){let r=byId[id],p=D.points[id],prefix=id===focus?'● ':'';if(mode==='kde')return `${prefix}${id}<br>n=${fmt(r.n_locations)} · density=${fmt(r.cell_density)}`;if(mode==='component_rank')return `${prefix}${id}<br>${T.componentsN}=${p?.component_count??'NA'} · ${T.largest}=${p?.largest_component_fraction==null?'NA':pct(p.largest_component_fraction)}`;if(mode==='n_genes')return `${prefix}${id}<br>${T.genes}: ${fmt(r.median_n_genes)} · expected ${fmt(r.median_n_genes_expected)}`;return `${prefix}${id}<br>${T.counts}: ${fmt(r.median_total_counts)} · expected ${fmt(r.median_total_counts_expected)}`}
function drawRow(containerId,ids,mode,extent,vlo,vhi,focus,region){let root=document.getElementById(containerId);root.innerHTML=ids.map((id,i)=>`<div class="slice-panel"><canvas id="${containerId}${i}" role="img" aria-label="${id} ${mode}"></canvas><div class="slice-cap ${id===focus?'focus':''}">${panelCaption(id,mode,focus)}</div></div>`).join('');ids.forEach((id,i)=>drawCanvas(document.getElementById(`${containerId}${i}`),id,extent,mode,vlo,vhi,id===focus,region))}
function setFinding(id,finding){let node=document.getElementById(id);node.className=`finding ${finding.status}`;node.textContent=finding.summary}
function drawEvidence(){let ids=currentIds(),focus=select.value,r=byId[focus],extent=sharedExtent(ids),region=D.regions[focus]||null,kdeValues=ids.flatMap(id=>D.points[id]?.kde||[]).filter(v=>v!==null&&v!==undefined).map(Number).filter(Number.isFinite),klo=q(kdeValues,.04),khi=q(kdeValues,.96),expressionMode=document.getElementById('expressionSelect').value,expressionValues=ids.flatMap(id=>D.points[id]?.[expressionMode]||[]).filter(v=>v!==null&&v!==undefined).map(Number).filter(Number.isFinite),elo=q(expressionValues,.04),ehi=q(expressionValues,.96);document.getElementById('reason').className=`reason-box ${r.final_call==='keep'?'keep-reason':''}`;document.getElementById('reason').innerHTML=`<span class="call ${r.final_call}">${r.final_call}</span> <b>${r.slice_id}</b><p>${r.user_reason}</p><span class="note">${T.causes}${r.likely_causes}</span>`;let domains=[['density',r.density_domain_score],['expression',r.expression_domain_score],['damage',r.damage_domain_score],['continuity',r.continuity_domain_score]];document.getElementById('domainBars').innerHTML=domains.map(([label,value])=>`<div class="domain-row"><span>${label}</span><span class="domain-track"><i class="${Number(value)>=.45?'strong':''}" style="width:${100*Math.max(0,Math.min(1,Number(value)||0))}%"></i></span><b>${fmt(value)}</b></div>`).join('')+`<div class="note">${T.strong}</div>`;let panels=[['kde',T.density],['expression',T.expression],['components',T.components]];document.getElementById('visualLinks').innerHTML=panels.map(([key,label])=>{let f=r.visual_findings[key];return `<div class="visual-link ${f.status}"><b>${label}</b><span>${f.summary}</span></div>`}).join('');setFinding('kdeFinding',r.visual_findings.kde);setFinding('expressionFinding',r.visual_findings.expression);setFinding('componentFinding',r.visual_findings.components);drawRow('kdeWindow',ids,'kde',extent,klo,khi,focus,region);drawRow('expressionWindow',ids,expressionMode,extent,elo,ehi,focus,null);drawRow('componentWindow',ids,'component_rank',extent,0,1,focus,null);document.getElementById('kdeLegend').textContent=`fixed-bandwidth KDE ${fmt(klo)}–${fmt(khi)}`;document.getElementById('expressionLegend').textContent=expressionValues.length?(expressionMode==='log_captured_counts'?`${fmt(Math.expm1(elo))}–${fmt(Math.expm1(ehi))} ${T.counts}`:`${fmt(elo)}–${fmt(ehi)} ${T.genes}`):T.unavailable;document.getElementById('regionNote').textContent=region?region.note:T.noRegion;document.getElementById('evidenceRows').innerHTML=(r.directional_evidence||[]).slice(0,9).map(e=>`<tr><td>${e.label}<br><span class="note">${e.domain}</span></td><td>${e.direction} ${e.delta==null?'':pct(e.delta)}</td><td class="num">${fmt(e.actual)}</td><td class="num">${fmt(e.expected)}</td></tr>`).join('')||`<tr><td colspan="4">${T.noEvidence}</td></tr>`}
window.addEventListener('resize',()=>{if(document.getElementById('density').classList.contains('active'))drawEvidence()});
function table(){let f=document.getElementById('callFilter').value,q=document.getElementById('sliceSearch').value.trim().toLowerCase(),rows=R.filter(r=>(f==='all'||r.final_call===f)&&String(r.slice_id).toLowerCase().includes(q));document.getElementById('sliceRows').innerHTML=rows.map(r=>`<tr data-id="${r.slice_id}"><td>${Number(r.slice_index)+1}</td><td><b>${r.slice_id}</b></td><td><span class="call ${r.final_call}">${r.final_call}</span></td><td class="num">${fmt(r.quality_anomaly_score)}</td><td class="num">${fmt(r.n_locations)}</td><td class="num">${fmt(r.cell_density)}</td><td class="num">${fmt(r.median_n_genes)}</td><td class="num">${fmt(r.median_total_counts)}</td><td class="num">${fmt(r.median_nn_distance)}</td><td class="num">${fmt(r.detected_genes)}</td><td class="reason-cell">${r.user_reason}</td></tr>`).join('');document.querySelectorAll('#sliceRows tr').forEach(row=>row.onclick=()=>showEvidence(row.dataset.id))}
document.getElementById('callFilter').onchange=table;document.getElementById('sliceSearch').oninput=table;table();drawEvidence();
</script></body></html>"""


_REPORT_TRANSLATIONS_EN: tuple[tuple[str, str], ...] = (
    (
        "单数据集开发附件 · 全部切片仅保留 keep / exclude · 不需要用户再次判定",
        "Single-dataset supplementary report · final binary keep/exclude calls · no additional user adjudication required",
    ),
    (
        "红色代表最终自动排除；绿色代表最终保留。每张切片的低质量原因会直接对应到 KDE 密度、表达捕获和连通组织成分三组图。",
        "Red denotes final automatic exclusion and green denotes retention. Each slice-level rationale is linked directly to the KDE density, expression-capture, and connected-component panels.",
    ),
    (
        "证据图在同一页面并排展示连续切片。每组图上方明确标记该证据是否支持排除，避免仅凭颜色主观判断。",
        "All evidence panels show the same consecutive slices on one page. Each panel explicitly states whether its evidence supports exclusion, avoiding subjective interpretation from color alone.",
    ),
    (
        "五项指标同时展示，虚线为本数据集切片中位数；红色方形为 exclude，绿色圆点为 keep。点击数据点可进入对应证据页面。",
        "Five metrics are shown across slice order. The dashed line is the dataset median; red squares are exclude calls and green circles are keep calls. Select a point to open its evidence page.",
    ),
    (
        "每个横条表示最终 exclude 切片中有多少张出现该类显著不利证据。证据可以重叠，因此横条数量之和不等于排除切片总数。点击横条可查看其中第一张切片。",
        "Each bar reports how many excluded slices contain that type of substantial adverse evidence. Evidence categories may overlap, so bar counts do not sum to the number of excluded slices. Select a bar to inspect its first listed slice.",
    ),
    (
        "连续切片使用共同的绝对观测值色标；该图显示 H5AD 捕获量，不代表原始 read depth。",
        "Consecutive slices use one shared scale of absolute observed values. The panel shows captured abundance in the H5AD, not raw read depth.",
    ),
    (
        "颜色表示每张切片内部的连通成分排名；小型分离组织块会以不同颜色显示。",
        "Colors indicate connected-component rank within each slice; small detached tissue components are shown in distinct colors.",
    ),
    (
        "表格给出观测值与局部连续切片期望值。图像负责空间定位，表格负责说明变化方向与幅度。",
        "The table compares observed values with local serial-slice expectations. Images localize the evidence, while the table states its direction and magnitude.",
    ),
    (
        "Captured counts 和 detected genes 直接来自 H5AD 的逐点观测列；显示时在当前连续切片窗口内使用共同的 4%–96% 色标，不对每张切片单独归一化。",
        "Captured counts and detected genes come directly from point-level H5AD observation columns. The current serial-slice window shares a 4th–96th percentile color scale; slices are not normalized independently.",
    ),
    (
        "在每张切片内以 3.25 倍中位最近邻距离连接空间点，随后计算无向连通成分；颜色表示按点数排序的成分排名，与检测器的 fragmentation 和 largest-component fraction 定义一致。",
        "Within each slice, spatial points are linked at 3.25 times the median nearest-neighbour distance before undirected connected components are computed. Colors show components ranked by point count, matching the detector definitions of fragmentation and largest-component fraction.",
    ),
    (
        "横条统计 exclude 切片中某类不利方向指标异常达到 0.50 的切片数；同一切片可同时支持多个证据类别，因此横条数值不能相加。",
        "Bars count excluded slices with an adverse directional metric anomaly of at least 0.50 in each evidence category. One slice may support multiple categories, so bar values must not be summed.",
    ),
    (
        "直方图汇总每项指标在全部切片中的分布，并将 keep（绿色）和 exclude（红色）堆叠显示；蓝色虚线为中位数。分箱使用 Freedman–Diaconis 规则并限制为 6–14 箱，跨度超过两个数量级时使用 log10 横轴。",
        "Each histogram summarizes the across-slice distribution, with keep (green) and exclude (red) counts stacked within each bin; the dashed blue line marks the dataset median. Bins follow the Freedman–Diaconis rule with 6–14-bin guardrails, and a log10 x-axis is used when the range spans at least two orders of magnitude.",
    ),
    (
        "五项指标同时展示，虚线为本数据集切片中位数；红色方形为 exclude，绿色圆点为 keep。点击数据点可进入对应 KDE 窗口。",
        "Five metrics are shown across slice order. The dashed line is the dataset median; red squares are exclude calls and green circles are keep calls. Select a point to open its KDE window.",
    ),
    (
        "每张切片与前后连续切片比较，默认检测窗口由 3/5/7 候选中选择；同时使用更宽的稳健趋势防止连续坏片互相掩护。密度/完整性、表达捕获、损伤及跨切片连续性四个证据域共同形成异常分数。",
        "Each slice is compared with adjacent serial sections. The detection window is selected from 3-, 5-, and 7-slice candidates, while a wider robust trend limits mutual masking by consecutive poor-quality sections. Density/integrity, expression capture, damage, and cross-slice continuity jointly define the anomaly score.",
    ),
    (
        "方法分为两阶段：第一阶段按锁定阈值形成 keep / review / exclude 内部初筛，高分但未通过保护条件的切片降级为 review；第二阶段用连续分数段评估全部 review。只有在校准实验中兼具安全性和额外检出增益的分数段才启用排除；缺少证据支持的较低分段明确设为 keep-only。启用的排除路径仍必须满足双侧上下文、未触发部分结构保护、多域证据、置信度及 3/5/7 窗口稳定性。公开结果只有 keep / exclude。",
        "The method has two stages. Stage 1 uses locked score thresholds to create an internal keep/review/exclude triage and downgrades unsafe high-score exclusions to review. Stage 2 evaluates every review row with contiguous calibrated score bands spanning the full uncertainty interval. A band is exclude-enabled only when calibration demonstrates positive safe detection gain; an unsupported lower-score band is explicitly keep-only. Enabled lower-score rules require stronger multi-domain evidence, and every enabled route also requires two-sided context, no partial-structure protection, adequate confidence, and stable 3/5/7-window evidence. A row becomes exclude only when every enabled condition passes; otherwise it becomes keep. Public output contains only keep/exclude.",
    ),
    (
        "“比局部窗口期望低/高”描述的是相邻切片趋势中的相对变化，不是跨技术平台阈值。可能原因只用于实验排查，不能从 H5AD 单独确诊。",
        "Lower or higher than the local-window expectation describes a relative change along adjacent sections, not a threshold transferable across technologies. Plausible causes are provided for experimental troubleshooting and cannot be diagnosed from the H5AD alone.",
    ),
    (
        "H5AD 可观察每点捕获计数、检出基因和复杂度，但没有原始 reads、PCR duplication、mapping rate 或 saturation 时，不应把 captured-count proxy 直接称为原始测序深度。",
        "The H5AD provides captured counts per spot, detected genes, and complexity. Without raw reads, PCR duplication, mapping rate, or saturation, the captured-count proxy must not be presented as raw sequencing depth.",
    ),
    (
        "本报告使用 observed nearest-neighbour spacing 作为分辨率相关的采样间距代理，单位来自原生坐标；它不是平台标称的 spot diameter 或分辨率。",
        "This report uses observed nearest-neighbour spacing as a resolution-related sampling-distance proxy in native coordinate units; it is not the platform's nominal spot diameter or stated resolution.",
    ),
    (
        "二维 Gaussian KDE 使用统一固定带宽 ${fmt(D.kde_bandwidth)} 个原生坐标单位；窗口内共同色标便于直接比较。页面展示 ${Object.values(D.points).reduce((a,p)=>a+p.displayed,0).toLocaleString()} / ${Object.values(D.points).reduce((a,p)=>a+p.available,0).toLocaleString()} 个可用点。",
        "Two-dimensional Gaussian KDE uses a fixed bandwidth of ${fmt(D.kde_bandwidth)} native coordinate units. A shared color scale supports direct comparison within each window. The report displays ${Object.values(D.points).reduce((a,p)=>a+p.displayed,0).toLocaleString()} / ${Object.values(D.points).reduce((a,p)=>a+p.available,0).toLocaleString()} available points.",
    ),
    (
        "红色代表最终自动排除；绿色代表最终保留。异常原因写明了相对局部连续切片期望值的方向、幅度和观测值。",
        "Red denotes a final automatic exclusion and green denotes retention. Each anomaly explanation states its direction, magnitude, and observed value relative to the local serial-section expectation.",
    ),
    (
        "KDE 图使用同一检测窗口中的共同色标；黄色框把“邻片通常存在、焦点切片明显稀疏”的区域直接标出来。",
        "KDE panels share a color scale within each detection window. The yellow box directly marks an area that is supported by neighboring slices but markedly sparse in the focal slice.",
    ),
    (
        "该切片未形成稳定的最大连续 KDE 密度缺失候选区；结论来自表中其他多域证据。",
        "No stable, largest contiguous KDE density-deficit region was identified in this slice; the call is supported by the other multi-domain evidence shown in the table.",
    ),
    (
        "本数据集没有自动排除切片。",
        "No slice was automatically excluded in this dataset.",
    ),
    (
        "该切片未形成稳定的最大连续 KDE 密度缺失候选区；结论来自下方其他证据。",
        "No stable contiguous KDE density-deficit region was identified; the call is supported by the other evidence shown below.",
    ),
    (
        "未见达到展示阈值的方向性异常",
        "No directional anomaly reached the display threshold",
    ),
    (
        "未见达到展示阈值的不利方向证据",
        "No adverse directional evidence reached the display threshold",
    ),
    ("逐点证据不可用", "Point-level evidence unavailable"),
    ("局部密度缺失候选区", "Candidate local density deficit"),
    ("切片低质量结论与图像证据", "Slice-level quality call and visual evidence"),
    ("焦点切片的方向性指标证据", "Directional metric evidence for the focal slice"),
    ("排除切片的关键证据", "Key evidence among excluded slices"),
    ("切片证据图", "Slice evidence"),
    ("A · KDE 局部细胞/位置密度", "A · Local cell/spot density by KDE"),
    ("B · 表达捕获", "B · Expression capture"),
    ("C · 连通组织成分", "C · Connected tissue components"),
    ("最大连通成分", "Largest component"),
    ("第二成分", "Second component"),
    ("其他成分", "Other components"),
    ("表达图如何生成", "Expression-capture visualization"),
    ("连通成分如何生成", "Connected-component visualization"),
    ("关键证据条形图", "Key-evidence frequency chart"),
    ("表达指标", "Expression metric"),
    ("密度 / 组织量", "Density / tissue amount"),
    ("分层重审排除", "Review-tier excludes"),
    ("表达捕获", "Expression capture"),
    ("连通组织成分", "Connected components"),
    ("强证据阈值", "Strong-evidence threshold"),
    ("最大成分", "largest component"),
    ("成分数", "components"),
    ("每点捕获计数", "captured counts per spot"),
    ("每点检出基因", "detected genes per spot"),
    ("切片指标分布直方图", "Slice-metric distribution histograms"),
    ("KDE 局部细胞/位置密度", "Local cell/spot density by KDE"),
    ("焦点切片的可读证据", "Interpretable evidence for the focal slice"),
    ("全部切片最终结果", "Final calls for all slices"),
    ("自动排除切片", "Automatically excluded slices"),
    ("逐切片统计趋势", "Per-slice metric trends"),
    ("窗口共同坐标范围", "Shared coordinate extent within the window"),
    ("每张切片独立放大", "Fit each slice independently"),
    ("文章数据集汇总表 CSV", "Dataset summary table (CSV)"),
    ("文章排除切片表 CSV", "Excluded-slice table (CSV)"),
    ("文章表格预览", "Publication table preview"),
    ("KDE 密度窗口", "KDE density window"),
    ("方法与口径", "Methods and interpretation"),
    ("最终二元规则", "Final binary decision rule"),
    ("KDE 图如何生成", "KDE visualization"),
    ("原因如何解释", "Interpretation of causes"),
    ("表达深度限制", "Expression-depth limitation"),
    ("空间分辨率限制", "Spatial-resolution limitation"),
    ("检测流程", "Detection workflow"),
    ("结果总览", "Overview"),
    ("统计趋势", "Statistics"),
    ("全部切片", "All slices"),
    ("读图要点", "How to read this report"),
    ("焦点切片", "Focal slice"),
    ("展示窗口", "Display window"),
    ("着色指标", "Color metric"),
    ("空间范围", "Spatial extent"),
    ("5 张连续切片", "5 consecutive slices"),
    ("3 张连续切片", "3 consecutive slices"),
    ("KDE 细胞密度", "KDE cell/spot density"),
    ("最终状态", "Final call"),
    ("结论及原因", "Call and rationale"),
    ("切片搜索", "Slice search"),
    ("例如 slices_4", "e.g., slices_4"),
    ("自动选择检测窗口", "Selected detection window"),
    ("可能原因：", "Plausible causes: "),
    ("位置数", "Locations"),
    ("位置", "locations"),
    ("局部期望", "Local expectation"),
    ("指标", "Metric"),
    ("方向", "Direction"),
    ("观测", "Observed"),
    ("状态", "Call"),
    ("全部", "All"),
    ("报告页面", "Report sections"),
    ("低", "Low"),
    ("高", "High"),
)


def _report_template(language: str) -> str:
    if language == "zh":
        return _REPORT_TEMPLATE
    rendered = _REPORT_TEMPLATE.replace('lang="zh-CN"', 'lang="en"')
    rendered = rendered.replace('"Noto Sans SC"', "Arial")
    rendered = rendered.replace("'zh-CN'", "'en-US'")
    for source, target in _REPORT_TRANSLATIONS_EN:
        rendered = rendered.replace(source, target)
    return rendered


__all__ = [
    "build_directional_slice_explanations",
    "render_binary_slice_quality_appendix",
    "write_collection_paper_tables",
]
