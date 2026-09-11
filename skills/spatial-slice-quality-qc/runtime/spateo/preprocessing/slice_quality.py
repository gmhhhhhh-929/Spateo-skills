"""Optional post-I/O quality control for serial spatial transcriptomics slices.

The routines in this module are deliberately non-destructive.  They calculate
slice-level evidence, compare each slice with a configurable local window, and
return ``keep``, ``review`` or ``exclude`` recommendations.  They never remove
observations from the caller's :class:`~anndata.AnnData` object.

Two public entry points cover the common in-memory use case:

``calculate_slice_quality``
    Calculate per-slice metrics and recommendations for one AnnData object.

``simulate_slice_quality_artifacts``
    Create a copy containing a separate simulated count layer and optional
    cell-density defects for controlled ground-truth evaluation.

The file-oriented helper ``scan_h5ad_series`` accepts either one multi-slice
H5AD or several H5AD files belonging to one ordered series.  Use
``scan_h5ad_collection`` for multiple independent datasets.  HTML and other
viewer concerns intentionally live in the companion skill rather than this
scientific-computation module.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad
from scipy import ndimage, sparse
from scipy.sparse import csgraph
from scipy.spatial import ConvexHull, cKDTree, distance

EPS = np.finfo(float).eps
SCHEMA_VERSION = "1.0"


@dataclass
class SliceQCConfig:
    """Configuration for serial-slice quality assessment.

    ``window`` is the total odd window width, including the focal slice.  With
    the default value 3, the expected value for an internal slice is derived
    from its immediate predecessor and successor.  Set it to ``"auto"`` to
    choose among ``window_candidates`` using controlled metric perturbations
    and baseline stability.
    """

    window: Union[int, str] = 3
    window_candidates: tuple[int, ...] = (3, 5, 7)
    k_neighbors: int = 12
    min_points_geometry: int = 20
    max_points_per_slice_report: int = 1200
    random_seed: int = 13
    review_threshold: float = 0.38
    exclude_threshold: float = 0.64
    severe_domain_threshold: float = 0.78
    minimum_corrob_domains: int = 2
    profile_genes: int = 3000
    mito_prefixes: tuple[str, ...] = (
        "MT-",
        "mt-",
        "Mt-",
        "mt:",
        "mitochondrion_genome",
    )
    report_title: str = "Pre-alignment spatial slice quality control"


@dataclass
class SliceSeriesResult:
    """Calculated metrics, calls, display samples and provenance."""

    metrics: pd.DataFrame
    point_samples: dict[str, dict[str, list[Any]]]
    provenance: dict[str, Any]
    profiles: Optional[np.ndarray] = None
    profile_genes: tuple[str, ...] = ()


@dataclass(frozen=True)
class ReviewEvidenceTier:
    """One calibrated evidence route for resolving a stage-1 review slice.

    Tiers partition the review score interval. Exclude-enabled lower-score
    tiers require progressively stronger evidence. A calibrated interval with
    no safe incremental exclusion rule can be declared keep-only by setting
    ``enable_exclude=False``. Tiers are evaluated in declared order and never
    bypass two-sided-context or anatomical-protection guardrails.
    """

    name: str
    min_score: float
    max_score: Optional[float] = None
    min_corroborating_domains: int = 2
    severe_domain_threshold: float = 0.85
    min_window_stability: float = 1.0
    min_score_confidence: float = 0.90
    window_detector_mode: str = "candidate"
    enable_exclude: bool = True

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "ReviewEvidenceTier":
        return cls(
            name=str(value["name"]),
            min_score=float(value["min_score"]),
            max_score=(None if value.get("max_score") is None else float(value["max_score"])),
            min_corroborating_domains=int(value.get("min_corroborating_domains", 2)),
            severe_domain_threshold=float(value.get("severe_domain_threshold", 0.85)),
            min_window_stability=float(value.get("min_window_stability", 1.0)),
            min_score_confidence=float(value.get("min_score_confidence", 0.90)),
            window_detector_mode=str(value.get("window_detector_mode", "candidate")),
            enable_exclude=bool(value.get("enable_exclude", True)),
        )


@dataclass
class HighConfidencePolicy:
    """Two-stage publication policy for binary keep/exclude outputs.

    Stage 1 assigns an internal threshold band (keep/review/exclude) and
    downgrades unsafe high-score exclusions to review.  Stage 2 resolves only
    the review queue with multiscale, multi-domain evidence and anatomical
    guardrails.  Public output contains only keep/exclude; the intermediate
    state remains in the audit table.
    """

    keep_max_score: float
    exclude_min_score: float
    calibration_id: str = "unversioned"
    require_two_sided: bool = True
    require_detector_call: bool = True
    enable_keep: bool = True
    enable_exclude: bool = True
    unresolved_action: str = "withhold"
    adaptive_exclude_min_score: Optional[float] = None
    adaptive_min_corroborating_domains: int = 2
    adaptive_severe_domain_threshold: float = 0.85
    adaptive_min_window_stability: float = 0.90
    adaptive_min_score_confidence: float = 0.90
    review_exclusion_tiers: tuple[ReviewEvidenceTier, ...] = ()
    benchmark_summary: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "HighConfidencePolicy":
        fields = {
            "keep_max_score": float(value["keep_max_score"]),
            "exclude_min_score": float(value["exclude_min_score"]),
            "calibration_id": str(value.get("calibration_id", "unversioned")),
            "require_two_sided": bool(value.get("require_two_sided", True)),
            "require_detector_call": bool(value.get("require_detector_call", True)),
            "enable_keep": bool(value.get("enable_keep", True)),
            "enable_exclude": bool(value.get("enable_exclude", True)),
            "unresolved_action": str(value.get("unresolved_action", "withhold")),
            "adaptive_exclude_min_score": (
                None if value.get("adaptive_exclude_min_score") is None else float(value["adaptive_exclude_min_score"])
            ),
            "adaptive_min_corroborating_domains": int(value.get("adaptive_min_corroborating_domains", 2)),
            "adaptive_severe_domain_threshold": float(value.get("adaptive_severe_domain_threshold", 0.85)),
            "adaptive_min_window_stability": float(value.get("adaptive_min_window_stability", 0.90)),
            "adaptive_min_score_confidence": float(value.get("adaptive_min_score_confidence", 0.90)),
            "review_exclusion_tiers": tuple(
                ReviewEvidenceTier.from_mapping(item) for item in value.get("review_exclusion_tiers", ())
            ),
            "benchmark_summary": dict(value.get("benchmark_summary", {})),
        }
        return cls(**fields)


@dataclass
class _FileSliceData:
    metrics: list[dict[str, Any]] = field(default_factory=list)
    point_samples: dict[str, dict[str, list[Any]]] = field(default_factory=dict)
    pseudobulk: dict[str, tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    celltype_profiles: dict[str, tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    source_info: dict[str, Any] = field(default_factory=dict)


_LOW_BAD_RULES: dict[str, tuple[float, float]] = {
    "n_locations": (0.35, 0.68),
    "cell_density": (0.25, 0.58),
    "largest_component_fraction": (0.08, 0.32),
    "median_total_counts": (0.25, 0.62),
    "median_n_genes": (0.22, 0.55),
    "library_complexity": (0.15, 0.42),
    "detected_gene_fraction": (0.15, 0.45),
}

_HIGH_BAD_RULES: dict[str, tuple[float, float]] = {
    "hole_fraction": (0.04, 0.24),
    "fragmentation": (0.03, 0.24),
    "knn_tail_ratio": (0.35, 1.50),
    "local_low_depth_fraction": (0.05, 0.32),
    "regional_low_depth_cluster_fraction": (0.03, 0.20),
    "zero_fraction": (0.08, 0.32),
    "median_pct_mito": (0.05, 0.25),
}


def _natural_key(value: Any) -> tuple[Any, ...]:
    return tuple(int(token) if token.isdigit() else token.lower() for token in re.split(r"(\d+)", str(value)))


def _jsonable(value: Any) -> Any:
    if isinstance(value, float):
        return None if not math.isfinite(value) else value
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.ndarray):
        return [_jsonable(x) for x in value.tolist()]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(x) for x in value]
    return value


_REPORT_STATISTICS: tuple[dict[str, str], ...] = (
    {
        "key": "median_n_genes",
        "label": "Median genes / spot",
        "short_label": "Genes / spot",
        "unit": "genes",
        "help": "Median number of detected genes per observed location in each slice.",
    },
    {
        "key": "median_total_counts",
        "label": "Median captured counts / spot",
        "short_label": "Counts / spot",
        "unit": "captured counts",
        "help": (
            "Median captured-molecule total per observed location. This is an H5AD count/UMI proxy, "
            "not raw read depth."
        ),
    },
    {
        "key": "cell_density",
        "label": "Cell / spot density",
        "short_label": "Density",
        "unit": "locations / coordinate-area",
        "help": "Observed locations divided by the two-dimensional convex-hull area in native coordinate units.",
    },
    {
        "key": "median_nn_distance",
        "label": "Observed spot spacing",
        "short_label": "Spot spacing",
        "unit": "coordinate units",
        "help": (
            "Median nearest-neighbour center-to-center distance in native coordinates. It is a sampling-spacing "
            "proxy, not the platform's nominal spatial resolution."
        ),
    },
    {
        "key": "detected_genes",
        "label": "Total detected genes / slice",
        "short_label": "Detected genes",
        "unit": "genes",
        "help": "Number of genes with non-zero captured expression anywhere in the slice.",
    },
)


def _summarize_report_statistics(metrics: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """Summarize available descriptive slice statistics without inventing missing values."""
    summaries: dict[str, dict[str, Any]] = {}
    for definition in _REPORT_STATISTICS:
        key = definition["key"]
        if key not in metrics:
            continue
        values = pd.to_numeric(metrics[key], errors="coerce").to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if not values.size:
            continue
        summaries[key] = {
            **definition,
            "median": float(np.median(values)),
            "q25": float(np.quantile(values, 0.25)),
            "q75": float(np.quantile(values, 0.75)),
            "minimum": float(np.min(values)),
            "maximum": float(np.max(values)),
            "n_slices": int(values.size),
        }
    return summaries


def _report_measurement_context(provenance: Mapping[str, Any]) -> dict[str, Any]:
    """Describe count and coordinate semantics conservatively for report readers."""
    sources = [item for item in provenance.get("sources", []) if isinstance(item, Mapping)]
    count_flags = [bool(item["count_like"]) for item in sources if item.get("count_like") is not None]
    if count_flags and all(count_flags):
        count_semantics = "count-like H5AD values; captured-count / UMI proxy"
    elif count_flags and not any(count_flags):
        count_semantics = "non-count-like H5AD values; expression-value scale"
    elif count_flags:
        count_semantics = "mixed count-like and non-count-like H5AD values"
    else:
        count_semantics = "H5AD value semantics not recorded"
    coordinate_sources = sorted(
        {str(item.get("coordinate_source")) for item in sources if item.get("coordinate_source")}
    )
    return {
        "count_semantics": count_semantics,
        "coordinate_unit": "native coordinate units",
        "coordinate_sources": coordinate_sources,
        "cross_dataset_note": (
            "Expression summaries are most comparable when count semantics match. Density and spot-spacing "
            "summaries are most comparable within the same platform and coordinate system."
        ),
    }


def _sha256(path: Union[str, Path], block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def _as_matrix(matrix: Any) -> Union[np.ndarray, sparse.csr_matrix]:
    """Materialize an AnnData-backed matrix without densifying sparse data."""
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()
    elif not sparse.issparse(matrix) and not isinstance(matrix, np.ndarray):
        try:
            matrix = matrix[:]
        except Exception:
            matrix = np.asarray(matrix)
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return np.asarray(matrix)


def _matrix_row_sum(matrix: Union[np.ndarray, sparse.spmatrix]) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel().astype(float)


def _matrix_row_nnz(matrix: Union[np.ndarray, sparse.spmatrix]) -> np.ndarray:
    if sparse.issparse(matrix):
        clean = matrix.tocsr(copy=True)
        clean.sum_duplicates()
        clean.eliminate_zeros()
        return np.asarray(clean.getnnz(axis=1)).ravel().astype(float)
    return np.count_nonzero(matrix, axis=1).astype(float)


def _matrix_data_sample(matrix: Union[np.ndarray, sparse.spmatrix], limit: int = 10000) -> np.ndarray:
    if sparse.issparse(matrix):
        data = matrix.data
    else:
        data = np.asarray(matrix).ravel()
    if data.size <= limit:
        return np.asarray(data, dtype=float)
    indices = np.linspace(0, data.size - 1, limit, dtype=int)
    return np.asarray(data[indices], dtype=float)


def _is_count_like(matrix: Union[np.ndarray, sparse.spmatrix]) -> bool:
    sample = _matrix_data_sample(matrix)
    if sample.size == 0:
        return True
    sample = sample[np.isfinite(sample)]
    if sample.size == 0 or np.nanmin(sample) < 0:
        return False
    return bool(np.mean(np.abs(sample - np.rint(sample)) < 1e-6) >= 0.98)


def _choose_layer(adata: AnnData, requested: Optional[str]) -> tuple[str, Union[np.ndarray, sparse.csr_matrix], bool]:
    if requested and requested.lower() != "auto":
        name = requested
    else:
        priorities = ("counts", "counts_X", "raw_counts", "count", "raw")
        name = next((key for key in priorities if key in adata.layers), "X")
    if name == "X":
        matrix = _as_matrix(adata.X)
    else:
        if name not in adata.layers:
            raise KeyError(f"Count layer {name!r} is missing. Available layers: {list(adata.layers.keys())}")
        matrix = _as_matrix(adata.layers[name])
    return name, matrix, _is_count_like(matrix)


def _first_obs_column(adata: AnnData, candidates: Sequence[str]) -> Optional[str]:
    return next((key for key in candidates if key in adata.obs.columns), None)


def _resolve_slice_labels(
    adata: AnnData,
    slice_key: Optional[str],
    spatial_key: Optional[str],
    file_stem: str,
) -> tuple[np.ndarray, str, Optional[np.ndarray]]:
    z_values: Optional[np.ndarray] = None
    key = slice_key
    if key is None or key == "auto":
        key = _first_obs_column(
            adata,
            (
                "slice_id",
                "slices",
                "slice",
                "sample_order",
                "section_id",
                "section",
                "library_id",
            ),
        )
    if key:
        labels = np.asarray(adata.obs[key].astype(str))
        z_key = _first_obs_column(adata, ("align_z", "z", "z_coord", "spatial_z"))
        if z_key:
            z_values = pd.to_numeric(adata.obs[z_key], errors="coerce").to_numpy(dtype=float)
        return labels, f"obs:{key}", z_values

    candidate_keys: list[str] = []
    if spatial_key and spatial_key != "auto":
        candidate_keys.append(spatial_key)
    candidate_keys.extend(["spatial_3d", "spatial", "align_spatial", "tdr_spatial"])
    for key_candidate in dict.fromkeys(candidate_keys):
        if key_candidate not in adata.obsm:
            continue
        coords = np.asarray(adata.obsm[key_candidate])
        if coords.ndim == 2 and coords.shape[1] >= 3:
            z_values = np.asarray(coords[:, 2], dtype=float)
            unique = np.unique(z_values[np.isfinite(z_values)])
            if 1 < unique.size <= max(500, int(math.sqrt(max(adata.n_obs, 1))) * 4):
                mapping = {value: f"z{value:g}" for value in sorted(unique)}
                labels = np.asarray([mapping.get(value, "zNA") for value in z_values], dtype=object)
                return labels.astype(str), f"obsm:{key_candidate}[:,2]", z_values
    return np.repeat(file_stem, adata.n_obs).astype(str), "file_stem", None


def _resolve_xy(
    adata: AnnData,
    spatial_key: Optional[str],
    x_key: Optional[str],
    y_key: Optional[str],
) -> tuple[np.ndarray, np.ndarray, str]:
    if x_key and y_key:
        if x_key not in adata.obs or y_key not in adata.obs:
            raise KeyError(f"Requested coordinate columns {x_key!r}/{y_key!r} are missing from adata.obs")
        return (
            pd.to_numeric(adata.obs[x_key], errors="coerce").to_numpy(dtype=float),
            pd.to_numeric(adata.obs[y_key], errors="coerce").to_numpy(dtype=float),
            f"obs:{x_key},{y_key}",
        )

    if spatial_key and spatial_key != "auto":
        if spatial_key not in adata.obsm:
            raise KeyError(f"Requested spatial_key={spatial_key!r} is missing from adata.obsm")
        coords = np.asarray(adata.obsm[spatial_key])
        if coords.ndim != 2 or coords.shape[1] < 2:
            raise ValueError(f"adata.obsm[{spatial_key!r}] must have at least two coordinate columns")
        return (
            coords[:, 0].astype(float),
            coords[:, 1].astype(float),
            f"obsm:{spatial_key}[:,0:2]",
        )

    if "raw_spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["raw_spatial"])
        if coords.ndim == 2 and coords.shape[1] >= 2:
            return (
                coords[:, 0].astype(float),
                coords[:, 1].astype(float),
                "obsm:raw_spatial[:,0:2]",
            )

    for pair in (
        ("align_x", "align_y"),
        ("spatial_X", "spatial_Y"),
        ("x", "y"),
        ("x_coord", "y_coord"),
    ):
        if pair[0] in adata.obs and pair[1] in adata.obs:
            return (
                pd.to_numeric(adata.obs[pair[0]], errors="coerce").to_numpy(dtype=float),
                pd.to_numeric(adata.obs[pair[1]], errors="coerce").to_numpy(dtype=float),
                f"obs:{pair[0]},{pair[1]}",
            )

    candidate_keys: list[str] = []
    if spatial_key and spatial_key != "auto":
        candidate_keys.append(spatial_key)
    # Raw per-slice coordinates are preferable for a pre-alignment check.
    candidate_keys.extend(
        [
            "raw_spatial",
            "aligned_spatial_3D",
            "3d_align_spatial",
            "spatial",
            "align_spatial",
            "spatial_3d",
            "tdr_spatial",
        ]
    )
    for key in dict.fromkeys(candidate_keys):
        if key in adata.obsm:
            coords = np.asarray(adata.obsm[key])
            if coords.ndim == 2 and coords.shape[1] >= 2:
                return (
                    coords[:, 0].astype(float),
                    coords[:, 1].astype(float),
                    f"obsm:{key}[:,0:2]",
                )
    raise KeyError("Could not infer spatial x/y coordinates. Pass spatial_key or x_key/y_key explicitly.")


def _resolve_order(
    labels: np.ndarray,
    z_values: Optional[np.ndarray],
    order_values: Optional[np.ndarray],
) -> tuple[list[str], dict[str, float], str]:
    unique = list(pd.unique(labels.astype(str)))
    if order_values is not None:
        table = pd.DataFrame({"slice": labels.astype(str), "order": order_values})
        med = table.groupby("slice", observed=True)["order"].median()
        ordered = sorted(
            unique,
            key=lambda value: (float(med.get(value, np.inf)), _natural_key(value)),
        )
        return (
            ordered,
            {value: float(med.get(value, np.nan)) for value in unique},
            "explicit_order",
        )
    if z_values is not None and np.isfinite(z_values).any():
        table = pd.DataFrame({"slice": labels.astype(str), "z": z_values})
        med = table.groupby("slice", observed=True)["z"].median()
        if med.nunique(dropna=True) > 1:
            ordered = sorted(
                unique,
                key=lambda value: (float(med.get(value, np.inf)), _natural_key(value)),
            )
            return (
                ordered,
                {value: float(med.get(value, np.nan)) for value in unique},
                "z_median",
            )
    ordered = sorted(unique, key=_natural_key)
    return (
        ordered,
        {value: float(i) for i, value in enumerate(ordered)},
        "natural_label",
    )


def _safe_hull_area(xy: np.ndarray) -> float:
    if xy.shape[0] < 3:
        return float("nan")
    try:
        return float(ConvexHull(xy).volume)
    except Exception:
        spans = np.nanmax(xy, axis=0) - np.nanmin(xy, axis=0)
        return float(np.prod(spans))


def _component_metrics(xy: np.ndarray, median_nn: float) -> tuple[int, float]:
    if xy.shape[0] < 2 or not np.isfinite(median_nn) or median_nn <= 0:
        return 1, 1.0
    radius = max(median_nn * 3.25, EPS)
    tree = cKDTree(xy)
    graph = tree.sparse_distance_matrix(tree, max_distance=radius, output_type="coo_matrix")
    graph.data[:] = 1
    n_comp, labels = csgraph.connected_components(graph.tocsr(), directed=False)
    counts = np.bincount(labels, minlength=n_comp)
    return int(n_comp), float(counts.max() / max(xy.shape[0], 1))


def _raster_hole_fraction(xy: np.ndarray, median_nn: float) -> float:
    if xy.shape[0] < 20 or not np.isfinite(median_nn) or median_nn <= 0:
        return 0.0
    mins = np.nanmin(xy, axis=0)
    spans = np.nanmax(xy, axis=0) - mins
    bin_size = max(median_nn * 1.8, float(np.nanmax(spans)) / 220.0, EPS)
    shape = np.maximum(np.ceil(spans / bin_size).astype(int) + 3, 4)
    if int(np.prod(shape)) > 300_000:
        bin_size *= math.sqrt(float(np.prod(shape)) / 300_000.0)
        shape = np.maximum(np.ceil(spans / bin_size).astype(int) + 3, 4)
    indices = np.floor((xy - mins) / bin_size).astype(int) + 1
    indices[:, 0] = np.clip(indices[:, 0], 0, shape[0] - 1)
    indices[:, 1] = np.clip(indices[:, 1], 0, shape[1] - 1)
    occupied = np.zeros(tuple(shape), dtype=bool)
    occupied[indices[:, 0], indices[:, 1]] = True
    support = ndimage.binary_dilation(occupied, iterations=1)
    support = ndimage.binary_closing(support, iterations=2)
    filled = ndimage.binary_fill_holes(support)
    denom = int(filled.sum())
    if denom == 0:
        return 0.0
    holes = filled & ~support
    labels, n_labels = ndimage.label(holes)
    if n_labels:
        sizes = np.bincount(labels.ravel())
        # Ignore one-bin raster speckles; retain true enclosed gaps.
        keep = np.where(sizes >= 3)[0]
        keep = keep[keep != 0]
        hole_count = int(sum(int(sizes[i]) for i in keep))
    else:
        hole_count = 0
    return float(hole_count / denom)


def _largest_low_depth_cluster(xy: np.ndarray, low: np.ndarray, median_nn: float) -> float:
    indices = np.flatnonzero(low)
    if indices.size == 0:
        return 0.0
    if indices.size == 1 or not np.isfinite(median_nn) or median_nn <= 0:
        return float(1 / max(xy.shape[0], 1))
    low_xy = xy[indices]
    tree = cKDTree(low_xy)
    graph = tree.sparse_distance_matrix(tree, max_distance=median_nn * 2.75, output_type="coo_matrix")
    graph.data[:] = 1
    n_comp, comp = csgraph.connected_components(graph.tocsr(), directed=False)
    size = np.bincount(comp, minlength=n_comp).max()
    return float(size / max(xy.shape[0], 1))


def _local_depth_metrics(xy: np.ndarray, total_counts: np.ndarray, k: int, median_nn: float) -> tuple[float, float]:
    n = xy.shape[0]
    if n < max(8, k + 2):
        return 0.0, 0.0
    kk = min(k + 1, n)
    tree = cKDTree(xy)
    _, neighbors = tree.query(xy, k=kk)
    neighbors = np.asarray(neighbors)
    if neighbors.ndim == 1:
        return 0.0, 0.0
    log_counts = np.log1p(np.maximum(total_counts, 0))
    local_med = np.median(log_counts[neighbors[:, 1:]], axis=1)
    delta = log_counts - local_med
    center = np.median(delta)
    mad = np.median(np.abs(delta - center)) * 1.4826
    if not np.isfinite(mad) or mad <= EPS:
        local_low = delta < -math.log(2.5)
    else:
        local_low = (delta - center) / mad < -2.75
    # A broad contiguous dryspot can contain most of each point's neighbours,
    # so a purely local comparison masks it.  Add a robust within-slice global
    # floor; persistent biological regions cancel in the later cross-slice
    # comparison, while a newly damaged region raises both fraction and cluster
    # evidence in the focal section.
    global_center = float(np.median(log_counts))
    global_low = log_counts < global_center - math.log(2.5)
    low = local_low | global_low
    return float(np.mean(low)), _largest_low_depth_cluster(xy, low, median_nn)


def _geometry_metrics(
    x: np.ndarray, y: np.ndarray, total_counts: np.ndarray, config: SliceQCConfig
) -> dict[str, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    xy = np.column_stack([x[valid], y[valid]])
    totals = np.asarray(total_counts)[valid]
    n = xy.shape[0]
    if n < 2:
        return {
            "hull_area": float("nan"),
            "cell_density": float("nan"),
            "median_nn_distance": float("nan"),
            "knn_tail_ratio": float("nan"),
            "component_count": 1,
            "largest_component_fraction": 1.0,
            "fragmentation": 0.0,
            "hole_fraction": 0.0,
            "local_low_depth_fraction": 0.0,
            "regional_low_depth_cluster_fraction": 0.0,
            "shape_aspect_ratio": float("nan"),
            "radial_q50": float("nan"),
            "radial_q90": float("nan"),
        }
    tree = cKDTree(xy)
    distances, _ = tree.query(xy, k=min(2, n))
    nn = distances[:, -1]
    positive = nn[np.isfinite(nn) & (nn > 0)]
    median_nn = float(np.median(positive)) if positive.size else float("nan")
    q95 = float(np.quantile(positive, 0.95)) if positive.size else float("nan")
    knn_tail = q95 / max(median_nn, EPS) - 1 if np.isfinite(median_nn) else float("nan")
    area = _safe_hull_area(xy)
    components, largest = _component_metrics(xy, median_nn)
    centered = xy - np.mean(xy, axis=0)
    covariance = np.cov(centered.T) if n > 2 else np.eye(2)
    eig = np.sort(np.maximum(np.linalg.eigvalsh(covariance), 0))[::-1]
    aspect = float(math.sqrt(eig[0] / max(eig[1], EPS))) if eig.size >= 2 else 1.0
    scale = math.sqrt(max(eig.sum(), EPS))
    radial = np.linalg.norm(centered, axis=1) / scale
    local_low, regional_low = _local_depth_metrics(xy, totals, config.k_neighbors, median_nn)
    return {
        "hull_area": area,
        "cell_density": (float(n / area) if np.isfinite(area) and area > 0 else float("nan")),
        "median_nn_distance": median_nn,
        "knn_tail_ratio": (float(max(knn_tail, 0)) if np.isfinite(knn_tail) else float("nan")),
        "component_count": components,
        "largest_component_fraction": largest,
        "fragmentation": float(1.0 - largest),
        "hole_fraction": _raster_hole_fraction(xy, median_nn),
        "local_low_depth_fraction": local_low,
        "regional_low_depth_cluster_fraction": regional_low,
        "shape_aspect_ratio": aspect,
        "radial_q50": float(np.quantile(radial, 0.50)),
        "radial_q90": float(np.quantile(radial, 0.90)),
    }


def _entropy(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    total = values.sum()
    if total <= 0:
        return float("nan")
    p = values[values > 0] / total
    return float(-np.sum(p * np.log(p)) / max(np.log(max(p.size, 2)), EPS))


def _sample_points(
    slice_id: str,
    x: np.ndarray,
    y: np.ndarray,
    total_counts: np.ndarray,
    max_points: int,
    seed: int,
) -> dict[str, list[Any]]:
    n = len(x)
    if n > max_points > 0:
        stable = int(hashlib.sha256(slice_id.encode("utf-8")).hexdigest()[:8], 16)
        rng = np.random.default_rng(seed + stable)
        idx = np.sort(rng.choice(n, max_points, replace=False))
    else:
        idx = np.arange(n)
    counts = np.log1p(np.maximum(np.asarray(total_counts)[idx], 0))
    lo, hi = np.quantile(counts, [0.02, 0.98]) if counts.size else (0.0, 1.0)
    scaled = np.clip((counts - lo) / max(hi - lo, EPS), 0, 1)
    return {
        "x": np.asarray(x)[idx].astype(float).round(5).tolist(),
        "y": np.asarray(y)[idx].astype(float).round(5).tolist(),
        "depth": scaled.astype(float).round(5).tolist(),
        "counts": _jsonable(np.asarray(total_counts)[idx].astype(float)),
    }


def _obs_metric(adata: AnnData, candidates: Sequence[str]) -> Optional[np.ndarray]:
    key = _first_obs_column(adata, candidates)
    if not key:
        return None
    values = pd.to_numeric(adata.obs[key], errors="coerce").to_numpy(dtype=float)
    return values if np.isfinite(values).any() else None


def _extract_one_file(
    path: Union[str, Path],
    config: SliceQCConfig,
    slice_key: Optional[str],
    spatial_key: Optional[str],
    x_key: Optional[str],
    y_key: Optional[str],
    order_key: Optional[str],
    layer: Optional[str],
    celltype_key: Optional[str],
    file_prefix: bool,
) -> _FileSliceData:
    path = Path(path).expanduser().resolve()
    adata = read_h5ad(path, backed="r")
    try:
        labels, label_source, z_values = _resolve_slice_labels(adata, slice_key, spatial_key, path.stem)
        x, y, coordinate_source = _resolve_xy(adata, spatial_key, x_key, y_key)
        order_values = None
        if order_key:
            if order_key not in adata.obs:
                raise KeyError(f"order_key={order_key!r} is missing from {path.name}")
            order_values = pd.to_numeric(adata.obs[order_key], errors="coerce").to_numpy(dtype=float)
        ordered, order_map, order_source = _resolve_order(labels, z_values, order_values)
        if file_prefix and len(ordered) > 1:
            remap = {value: f"{path.stem}::{value}" for value in ordered}
            labels = np.asarray([remap[value] for value in labels], dtype=str)
            ordered = [remap[value] for value in ordered]
            order_map = {remap[key]: value for key, value in order_map.items()}

        resolved_layer, matrix, count_like = _choose_layer(adata, layer)
        representation = str(adata.uns.get("blind_input_contract", {}).get("expression_representation", ""))
        annotation_matrix = resolved_layer == "X" and "one-hot" in representation.lower()
        matrix_totals = _matrix_row_sum(matrix)
        matrix_genes = _matrix_row_nnz(matrix)
        observed_totals = _obs_metric(adata, ("total_counts", "nCounts", "nCount_RNA", "total_umi", "UMI_count"))
        observed_genes = _obs_metric(adata, ("n_genes_by_counts", "nGenes", "nFeature_RNA", "gene_count"))
        total_counts = (
            observed_totals
            if observed_totals is not None and not count_like and resolved_layer != "slice_qc_simulated_counts"
            else matrix_totals
        )
        n_genes = (
            observed_genes
            if observed_genes is not None and not count_like and resolved_layer != "slice_qc_simulated_counts"
            else matrix_genes
        )
        if annotation_matrix:
            total_counts = np.full(adata.n_obs, np.nan)
            n_genes = np.full(adata.n_obs, np.nan)
        mito = _obs_metric(adata, ("pct_counts_mt", "pMito", "percent_mito", "mito_ratio"))
        if mito is None and count_like:
            names = np.asarray(adata.var_names.astype(str))
            mito_mask = np.zeros(adata.n_vars, dtype=bool)
            for prefix in config.mito_prefixes:
                mito_mask |= np.char.startswith(names.astype(str), prefix)
            if mito_mask.any():
                mito_counts = _matrix_row_sum(matrix[:, mito_mask])
                mito = np.divide(mito_counts, np.maximum(matrix_totals, EPS))
        if mito is not None and np.nanmedian(mito) > 1.5:
            mito = mito / 100.0

        requested_celltype = celltype_key
        if not requested_celltype or requested_celltype == "auto":
            requested_celltype = _first_obs_column(
                adata,
                (
                    "celltype",
                    "annotation",
                    "anno",
                    "Annotation_2_tissue",
                    "lineage",
                    "cluster",
                ),
            )
        celltypes = np.asarray(adata.obs[requested_celltype].astype(str)) if requested_celltype else None

        codes = pd.Categorical(labels, categories=ordered, ordered=True).codes
        if np.any(codes < 0):
            raise ValueError("Failed to assign one or more observations to a slice")
        indicator = sparse.csr_matrix(
            (np.ones(adata.n_obs, dtype=np.float32), (codes, np.arange(adata.n_obs))),
            shape=(len(ordered), adata.n_obs),
        )
        bulk = indicator @ matrix
        bulk = bulk.tocsr() if sparse.issparse(bulk) else np.asarray(bulk)

        output = _FileSliceData()
        output.source_info = {
            "path": str(path),
            "sha256": _sha256(path),
            "shape": [int(adata.n_obs), int(adata.n_vars)],
            "slice_source": label_source,
            "coordinate_source": coordinate_source,
            "order_source": order_source,
            "count_layer": resolved_layer,
            "count_like": bool(count_like),
            "depth_source": (
                "obs"
                if observed_totals is not None and not count_like and resolved_layer != "slice_qc_simulated_counts"
                else resolved_layer
            ),
            "celltype_key": requested_celltype,
            "expression_capture_available": not annotation_matrix,
            "expression_representation": representation or "selected expression matrix",
            "obs_totals_match_selected_matrix": (bool(np.allclose(observed_totals, matrix_totals, equal_nan=True)) if observed_totals is not None else None),
            "obs_genes_match_selected_matrix": (bool(np.allclose(observed_genes, matrix_genes, equal_nan=True)) if observed_genes is not None else None),
        }
        var_names = np.asarray(adata.var_names.astype(str))
        for position, slice_id in enumerate(ordered):
            idx = np.flatnonzero(codes == position)
            sx = np.asarray(x[idx], dtype=float)
            sy = np.asarray(y[idx], dtype=float)
            stotal = np.asarray(total_counts[idx], dtype=float)
            sgenes = np.asarray(n_genes[idx], dtype=float)
            geometry = _geometry_metrics(sx, sy, stotal, config)
            if sparse.issparse(bulk):
                profile = np.asarray(bulk.getrow(position).toarray()).ravel().astype(float)
            else:
                profile = np.asarray(bulk[position]).ravel().astype(float)
            detected = int(np.count_nonzero(profile))
            median_total = float(np.nanmedian(stotal)) if np.isfinite(stotal).any() else float("nan")
            median_genes = float(np.nanmedian(sgenes)) if np.isfinite(sgenes).any() else float("nan")
            complexity_values = np.divide(sgenes, np.maximum(np.log1p(stotal), EPS))
            row = {
                "slice_id": slice_id,
                "source_file": str(path),
                "source_slice_id": slice_id.split("::", 1)[-1],
                "slice_order_value": order_map.get(slice_id, float(position)),
                "n_locations": int(idx.size),
                "median_total_counts": median_total,
                "mean_total_counts": (float(np.nanmean(stotal)) if np.isfinite(stotal).any() else float("nan")),
                "median_n_genes": median_genes,
                "zero_fraction": float(1.0 - np.nansum(sgenes) / max(idx.size * adata.n_vars, 1)),
                "library_complexity": (
                    float(np.nanmedian(complexity_values)) if np.isfinite(complexity_values).any() else float("nan")
                ),
                "detected_genes": detected,
                "detected_gene_fraction": float(detected / max(adata.n_vars, 1)),
                "median_pct_mito": (float(np.nanmedian(mito[idx])) if mito is not None and idx.size else float("nan")),
                "celltype_entropy": float("nan"),
                **geometry,
            }
            if celltypes is not None:
                categories, counts = np.unique(celltypes[idx], return_counts=True)
                row["celltype_entropy"] = _entropy(counts)
                output.celltype_profiles[slice_id] = (
                    categories.astype(str),
                    counts.astype(float),
                )
            row["expression_capture_available"] = not annotation_matrix
            row["coordinate_valid_fraction"] = float(np.mean(np.isfinite(sx) & np.isfinite(sy)))
            row["input_warning"] = "annotation_one_hot_no_capture_evidence" if annotation_matrix else ""
            if annotation_matrix:
                for key in ("median_total_counts", "mean_total_counts", "median_n_genes", "zero_fraction", "library_complexity", "detected_genes", "detected_gene_fraction", "median_pct_mito", "local_low_depth_fraction", "regional_low_depth_cluster_fraction"):
                    row[key] = float("nan")
            output.metrics.append(row)
            output.point_samples[slice_id] = _sample_points(
                slice_id,
                sx,
                sy,
                stotal,
                config.max_points_per_slice_report,
                config.random_seed,
            )
            if not annotation_matrix:
                output.pseudobulk[slice_id] = (var_names, profile)
        return output
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def _align_profiles(
    pseudobulk: Mapping[str, tuple[np.ndarray, np.ndarray]],
    slice_order: Sequence[str],
    max_genes: int,
    fixed_genes: Optional[Sequence[str]] = None,
) -> tuple[Optional[np.ndarray], tuple[str, ...]]:
    if not pseudobulk or any(s not in pseudobulk for s in slice_order):
        return None, ()
    gene_sets = [set(map(str, pseudobulk[slice_id][0])) for slice_id in slice_order]
    common = set.intersection(*gene_sets) if gene_sets else set()
    if not common:
        return None, ()
    if fixed_genes is not None:
        common_sorted = [str(gene) for gene in fixed_genes if str(gene) in common]
        if not common_sorted:
            return None, ()
    else:
        common_sorted = sorted(common)
    raw = np.zeros((len(slice_order), len(common_sorted)), dtype=float)
    for i, slice_id in enumerate(slice_order):
        genes, values = pseudobulk[slice_id]
        lookup = {str(gene): j for j, gene in enumerate(genes)}
        raw[i] = np.asarray([values[lookup[gene]] for gene in common_sorted], dtype=float)
    if fixed_genes is None and raw.shape[1] > max_genes > 0:
        normalized = np.log1p(raw / np.maximum(raw.sum(axis=1, keepdims=True), EPS) * 10_000)
        variance = np.var(normalized, axis=0)
        keep = np.argsort(variance)[-max_genes:]
        raw = raw[:, keep]
        common_sorted = [common_sorted[i] for i in keep]
    normalized = np.log1p(raw / np.maximum(raw.sum(axis=1, keepdims=True), EPS) * 10_000)
    return normalized, tuple(common_sorted)


def _profile_continuity(profiles: Optional[np.ndarray], window: int) -> tuple[np.ndarray, np.ndarray]:
    if profiles is None:
        return np.asarray([]), np.asarray([])
    n = profiles.shape[0]
    divergence = np.full(n, np.nan)
    neighbor_agreement = np.full(n, np.nan)
    half = max(window // 2, 1)
    for i in range(n):
        left = list(range(max(0, i - half), i))
        right = list(range(i + 1, min(n, i + half + 1)))
        neighbors = left + right
        if not neighbors:
            continue
        expected = np.mean(profiles[neighbors], axis=0)
        denom = np.linalg.norm(profiles[i]) * np.linalg.norm(expected)
        divergence[i] = 1.0 - float(np.dot(profiles[i], expected) / max(denom, EPS))
        if left and right:
            left_mean = np.mean(profiles[left], axis=0)
            right_mean = np.mean(profiles[right], axis=0)
            denom_lr = np.linalg.norm(left_mean) * np.linalg.norm(right_mean)
            neighbor_agreement[i] = float(np.dot(left_mean, right_mean) / max(denom_lr, EPS))
    return divergence, neighbor_agreement


def _celltype_js_profiles(
    celltype_profiles: Mapping[str, tuple[np.ndarray, np.ndarray]],
    slice_order: Sequence[str],
    window: int,
) -> np.ndarray:
    if not celltype_profiles:
        return np.full(len(slice_order), np.nan)
    categories = sorted(set().union(*(set(values[0]) for values in celltype_profiles.values())))
    lookup = {category: i for i, category in enumerate(categories)}
    matrix = np.zeros((len(slice_order), len(categories)), dtype=float)
    for i, slice_id in enumerate(slice_order):
        if slice_id not in celltype_profiles:
            continue
        names, counts = celltype_profiles[slice_id]
        for name, count in zip(names, counts):
            matrix[i, lookup[str(name)]] = count
    matrix /= np.maximum(matrix.sum(axis=1, keepdims=True), EPS)
    out = np.full(len(slice_order), np.nan)
    half = max(window // 2, 1)
    for i in range(len(slice_order)):
        neighbors = list(range(max(0, i - half), i)) + list(range(i + 1, min(len(slice_order), i + half + 1)))
        if neighbors:
            expected = np.mean(matrix[neighbors], axis=0)
            out[i] = float(distance.jensenshannon(matrix[i] + EPS, expected + EPS, base=2.0) ** 2)
    return out


def _robust_scale(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 1.0
    mad = np.median(np.abs(values - np.median(values))) * 1.4826
    if np.isfinite(mad) and mad > EPS:
        return float(mad)
    std = np.std(values)
    return float(std) if np.isfinite(std) and std > EPS else 1.0


def _local_expected(values: np.ndarray, window: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(values)
    half = max(window // 2, 1)
    expected = np.full(n, np.nan)
    two_sided = np.zeros(n, dtype=bool)
    neighbor_disagreement = np.full(n, np.nan)
    positions = np.arange(n, dtype=float)
    for i in range(n):
        left = list(range(max(0, i - half), i))
        right = list(range(i + 1, min(n, i + half + 1)))
        valid_left = [j for j in left if np.isfinite(values[j])]
        valid_right = [j for j in right if np.isfinite(values[j])]
        neighbors = valid_left + valid_right
        if not neighbors:
            continue
        if valid_left and valid_right:
            two_sided[i] = True
            x = positions[neighbors]
            y = values[neighbors]
            if len(neighbors) >= 2 and np.ptp(x) > 0:
                # Median pairwise slopes are resistant to a short consecutive
                # run of damaged neighbours; ordinary least squares lets those
                # slices mask one another.
                dx = x[np.newaxis, :] - x[:, np.newaxis]
                dy = y[np.newaxis, :] - y[:, np.newaxis]
                upper = np.triu(np.ones_like(dx, dtype=bool), k=1) & (np.abs(dx) > EPS)
                slopes = dy[upper] / dx[upper]
                slope = float(np.median(slopes)) if slopes.size else 0.0
                intercept = float(np.median(y - slope * x))
                expected[i] = intercept + slope * positions[i]
            else:
                expected[i] = float(np.mean(y))
            left_value = float(np.median(values[valid_left]))
            right_value = float(np.median(values[valid_right]))
            neighbor_disagreement[i] = abs(left_value - right_value)
        else:
            expected[i] = float(np.median(values[neighbors]))
    return expected, two_sided, neighbor_disagreement


def _metric_anomaly(
    values: np.ndarray,
    window: int,
    direction: str,
    mild: float,
    severe: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    expected, two_sided, disagreement = _local_expected(values, window)
    residual = expected - values if direction == "low" else values - expected
    finite_residual = residual[np.isfinite(residual)]
    center = float(np.median(finite_residual)) if finite_residual.size else 0.0
    scale = _robust_scale(residual)
    z = np.maximum((residual - center) / max(scale, EPS), 0)
    global_scale = _robust_scale(values)
    agreement_penalty = np.exp(-np.nan_to_num(disagreement, nan=0.0) / max(global_scale * 2.5, EPS))

    magnitude = np.zeros_like(values, dtype=float)
    valid = np.isfinite(values) & np.isfinite(expected)
    if direction == "low":
        magnitude[valid] = np.maximum(
            (expected[valid] - values[valid]) / np.maximum(np.abs(expected[valid]), EPS),
            0,
        )
    else:
        # For bounded fractions, absolute increases are more interpretable than ratios near zero.
        magnitude[valid] = np.maximum(values[valid] - expected[valid], 0)
    ratio_score = np.clip((magnitude - mild) / max(severe - mild, EPS), 0, 1)
    z_score = np.clip((z - 1.5) / 3.0, 0, 1)
    if direction == "high":
        # A standardized residual can become arbitrarily large when a bounded
        # metric and its neighbours are all close to zero.  Require a minimum
        # absolute effect before that residual contributes fully.  This keeps,
        # for example, a 0.4% low-depth region from looking as severe as a 40%
        # low-depth region merely because the local expectation was exactly 0.
        z_score *= np.clip(magnitude / max(mild, EPS), 0, 1)
    score = np.maximum(ratio_score, z_score)
    score *= np.where(two_sided, 1.0, 0.72)
    score *= np.clip(agreement_penalty, 0.35, 1.0)
    score[~valid] = 0.0
    return score, expected, two_sided


def _weighted_available(frame: pd.DataFrame, weights: Mapping[str, float]) -> np.ndarray:
    numerator = np.zeros(len(frame), dtype=float)
    denominator = np.zeros(len(frame), dtype=float)
    for key, weight in weights.items():
        if key not in frame:
            continue
        values = pd.to_numeric(frame[key], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(values)
        numerator[valid] += values[valid] * weight
        denominator[valid] += weight
    return np.divide(numerator, np.maximum(denominator, EPS))


def _top_two_evidence(frame: pd.DataFrame, columns: Sequence[str]) -> np.ndarray:
    available = [key for key in columns if key in frame]
    if not available:
        return np.zeros(len(frame), dtype=float)
    values = frame[available].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(dtype=float)
    values.sort(axis=1)
    strongest = values[:, -1]
    second = values[:, -2] if values.shape[1] >= 2 else np.zeros(len(frame), dtype=float)
    return 0.65 * strongest + 0.35 * second


def _score_metrics(
    metrics: pd.DataFrame,
    profiles: Optional[np.ndarray],
    celltype_profiles: Mapping[str, tuple[np.ndarray, np.ndarray]],
    config: SliceQCConfig,
    window: int,
) -> pd.DataFrame:
    if window < 3 or window % 2 == 0:
        raise ValueError("window must be an odd integer >= 3")
    out = metrics.copy().reset_index(drop=True)
    support_window = min(9, len(out) if len(out) % 2 == 1 else len(out) - 1)
    support_window = max(support_window, 3)
    support_window = max(window, support_window)
    out["support_window_size"] = int(support_window)
    two_sided_matrix = []
    for metric, (mild, severe) in _LOW_BAD_RULES.items():
        if metric not in out:
            continue
        score, expected, two_sided = _metric_anomaly(
            pd.to_numeric(out[metric], errors="coerce").to_numpy(dtype=float),
            window,
            "low",
            mild,
            severe,
        )
        out[f"{metric}_expected"] = expected
        support_expected = expected
        if support_window > window:
            support_score, support_expected, _ = _metric_anomaly(
                pd.to_numeric(out[metric], errors="coerce").to_numpy(dtype=float),
                support_window,
                "low",
                mild,
                severe,
            )
            out[f"{metric}_support_expected"] = support_expected
            out[f"{metric}_support_anomaly"] = support_score
            score = np.maximum(score, support_score)
        if metric in {"median_total_counts", "median_n_genes"}:
            values = pd.to_numeric(out[metric], errors="coerce").to_numpy(dtype=float)
            valid = np.isfinite(values) & np.isfinite(expected)
            relative_drop = np.zeros(len(out), dtype=float)
            relative_drop[valid] = np.maximum(
                (expected[valid] - values[valid]) / np.maximum(np.abs(expected[valid]), EPS),
                0.0,
            )
            support_valid = np.isfinite(values) & np.isfinite(support_expected)
            relative_drop[support_valid] = np.maximum(
                relative_drop[support_valid],
                np.maximum(
                    (support_expected[support_valid] - values[support_valid])
                    / np.maximum(np.abs(support_expected[support_valid]), EPS),
                    0.0,
                ),
            )
            direct_effect = np.clip((relative_drop - mild) / max(severe - mild, EPS), 0, 1)
            direct_effect *= np.where(two_sided, 1.0, 0.72)
            score = np.maximum(score, direct_effect)
        out[f"{metric}_anomaly"] = score
        two_sided_matrix.append(two_sided)
    for metric, (mild, severe) in _HIGH_BAD_RULES.items():
        if metric not in out:
            continue
        score, expected, two_sided = _metric_anomaly(
            pd.to_numeric(out[metric], errors="coerce").to_numpy(dtype=float),
            window,
            "high",
            mild,
            severe,
        )
        out[f"{metric}_expected"] = expected
        if support_window > window:
            support_score, support_expected, _ = _metric_anomaly(
                pd.to_numeric(out[metric], errors="coerce").to_numpy(dtype=float),
                support_window,
                "high",
                mild,
                severe,
            )
            out[f"{metric}_support_expected"] = support_expected
            out[f"{metric}_support_anomaly"] = support_score
            score = np.maximum(score, support_score)
        out[f"{metric}_anomaly"] = score
        two_sided_matrix.append(two_sided)

    expression_divergence, expression_neighbor_similarity = _profile_continuity(profiles, window)
    if expression_divergence.size:
        out["expression_profile_divergence"] = expression_divergence
        out["neighbor_expression_similarity"] = expression_neighbor_similarity
        div_score, _, _ = _metric_anomaly(expression_divergence, window, "high", 0.03, 0.18)
        # Absolute divergence matters even if several adjacent slices are affected.
        absolute = np.clip((np.nan_to_num(expression_divergence) - 0.06) / 0.24, 0, 1)
        agreement = np.clip(
            (np.nan_to_num(expression_neighbor_similarity, nan=0.5) - 0.25) / 0.6,
            0.25,
            1,
        )
        out["expression_profile_anomaly"] = np.maximum(div_score, absolute) * agreement
    elif "expression_profile_anomaly" in out:
        # Re-scoring a saved metric table (used by the repeated calibration
        # benchmark) must retain the real profile-continuity evidence.
        out["expression_profile_anomaly"] = pd.to_numeric(out["expression_profile_anomaly"], errors="coerce").fillna(
            0.0
        )
    else:
        out["expression_profile_divergence"] = np.nan
        out["neighbor_expression_similarity"] = np.nan
        out["expression_profile_anomaly"] = 0.0

    slice_order = out["slice_id"].astype(str).tolist()
    celltype_js = _celltype_js_profiles(celltype_profiles, slice_order, window)
    if np.isfinite(celltype_js).any():
        out["celltype_js_divergence"] = celltype_js
        out["celltype_composition_anomaly"] = np.clip((np.nan_to_num(celltype_js) - 0.08) / 0.32, 0, 1)
    elif "celltype_composition_anomaly" in out:
        out["celltype_composition_anomaly"] = pd.to_numeric(
            out["celltype_composition_anomaly"], errors="coerce"
        ).fillna(0.0)
    else:
        out["celltype_js_divergence"] = np.nan
        out["celltype_composition_anomaly"] = 0.0

    out["density_domain_score"] = _weighted_available(
        out,
        {
            "cell_density_anomaly": 0.28,
            "n_locations_anomaly": 0.10,
            "hole_fraction_anomaly": 0.27,
            "knn_tail_ratio_anomaly": 0.13,
            "fragmentation_anomaly": 0.12,
            "largest_component_fraction_anomaly": 0.10,
        },
    )
    out["density_domain_score"] = np.maximum(
        out["density_domain_score"].to_numpy(),
        _top_two_evidence(
            out,
            (
                "cell_density_anomaly",
                "n_locations_anomaly",
                "hole_fraction_anomaly",
                "knn_tail_ratio_anomaly",
                "fragmentation_anomaly",
                "largest_component_fraction_anomaly",
            ),
        ),
    )
    out["expression_domain_score"] = _weighted_available(
        out,
        {
            "median_total_counts_anomaly": 0.27,
            "median_n_genes_anomaly": 0.24,
            "zero_fraction_anomaly": 0.14,
            "library_complexity_anomaly": 0.11,
            "detected_gene_fraction_anomaly": 0.08,
            "local_low_depth_fraction_anomaly": 0.08,
            "regional_low_depth_cluster_fraction_anomaly": 0.08,
        },
    )
    out["expression_domain_score"] = np.maximum(
        out["expression_domain_score"].to_numpy(),
        _top_two_evidence(
            out,
            (
                "median_total_counts_anomaly",
                "median_n_genes_anomaly",
                "zero_fraction_anomaly",
                "library_complexity_anomaly",
                "detected_gene_fraction_anomaly",
                "local_low_depth_fraction_anomaly",
                "regional_low_depth_cluster_fraction_anomaly",
            ),
        ),
    )
    out["damage_domain_score"] = _weighted_available(
        out,
        {"median_pct_mito_anomaly": 0.65, "fragmentation_anomaly": 0.35},
    )
    out["continuity_domain_score"] = _weighted_available(
        out,
        {
            "expression_profile_anomaly": 0.58,
            "celltype_composition_anomaly": 0.22,
            "hole_fraction_anomaly": 0.20,
        },
    )

    domains = out[
        [
            "density_domain_score",
            "expression_domain_score",
            "damage_domain_score",
            "continuity_domain_score",
        ]
    ].to_numpy()
    out["corroborating_domains"] = np.sum(domains >= 0.45, axis=1)
    raw_score = (
        0.31 * out["density_domain_score"].to_numpy()
        + 0.39 * out["expression_domain_score"].to_numpy()
        + 0.12 * out["damage_domain_score"].to_numpy()
        + 0.18 * out["continuity_domain_score"].to_numpy()
    )
    severe_single = np.max(domains, axis=1)
    raw_score = np.maximum(raw_score, np.clip((severe_single - 0.68) * 0.78 + 0.48, 0, 1))
    if two_sided_matrix:
        internal = np.any(np.vstack(two_sided_matrix), axis=0)
    else:
        internal = np.ones(len(out), dtype=bool)
    out["window_context"] = np.where(internal, "two_sided", "one_sided_endpoint")
    out["score_confidence"] = np.where(internal, 1.0, 0.72)

    density = out["density_domain_score"].to_numpy()
    expression = out["expression_domain_score"].to_numpy()
    damage = out["damage_domain_score"].to_numpy()
    continuity = out["continuity_domain_score"].to_numpy()
    holes = out.get("hole_fraction_anomaly", pd.Series(np.zeros(len(out)))).to_numpy(dtype=float)
    density_only = (density >= 0.48) & (expression < 0.32) & (damage < 0.38) & (continuity < 0.45)
    taper_like = (
        (out.get("n_locations_anomaly", pd.Series(np.zeros(len(out)))).to_numpy(dtype=float) >= 0.45)
        & (out.get("cell_density_anomaly", pd.Series(np.zeros(len(out)))).to_numpy(dtype=float) < 0.30)
        & (holes < 0.30)
        & (expression < 0.35)
    )
    out["partial_structure_protection"] = density_only | taper_like
    final_score = raw_score.copy()
    final_score[out["partial_structure_protection"].to_numpy()] = np.minimum(
        final_score[out["partial_structure_protection"].to_numpy()],
        config.exclude_threshold - 0.02,
    )
    final_score[~internal] *= 0.82
    out["quality_anomaly_score"] = np.clip(final_score, 0, 1)

    recommendations: list[str] = []
    reasons: list[str] = []
    domain_names = ["density", "expression", "damage", "continuity"]
    metric_columns = [key for key in out.columns if key.endswith("_anomaly") and key != "quality_anomaly"]
    for i in range(len(out)):
        domain_values = dict(zip(domain_names, domains[i]))
        strong = [name for name, value in domain_values.items() if value >= 0.45]
        severe = [name for name, value in domain_values.items() if value >= config.severe_domain_threshold]
        expression_severe = expression[i] >= config.severe_domain_threshold
        severe_metric = max((float(out.loc[i, key]) for key in metric_columns), default=0.0) >= 0.75
        can_exclude = internal[i] and not bool(out.loc[i, "partial_structure_protection"])
        if can_exclude and (
            (final_score[i] >= config.exclude_threshold and len(strong) >= config.minimum_corrob_domains)
            or expression_severe
            or (density[i] >= 0.82 and holes[i] >= 0.55 and continuity[i] >= 0.42)
        ):
            recommendation = "exclude"
        elif (
            final_score[i] >= config.review_threshold
            or bool(strong)
            or severe
            or severe_metric
            or bool(out.loc[i, "partial_structure_protection"])
        ):
            recommendation = "review"
        else:
            recommendation = "keep"
        detail = []
        if strong:
            detail.append("strong domains: " + ", ".join(strong))
        if severe:
            detail.append("severe: " + ", ".join(severe))
        if bool(out.loc[i, "partial_structure_protection"]):
            detail.append("geometry-only/partial-structure protection applied")
        if not internal[i]:
            detail.append("endpoint: one-sided evidence only")
        if not detail:
            detail.append("no corroborated local anomaly")
        recommendations.append(recommendation)
        reasons.append("; ".join(detail))
    out["recommendation"] = recommendations
    out["low_quality_flag"] = out["recommendation"].eq("exclude")
    out["reason"] = reasons
    out["window_size"] = int(window)
    return out


def add_multiscale_exclusion_evidence(
    metrics: pd.DataFrame,
    *,
    config: Optional[SliceQCConfig] = None,
    windows: Sequence[int] = (3, 5, 7),
    minimum_corroborating_domains: int = 2,
    severe_domain_threshold: float = 0.85,
) -> pd.DataFrame:
    """Add auditable multi-window evidence for resolving near-boundary slices.

    The saved slice metric table is rescored independently at every requested
    odd window width. It stores the score, detector call, four domain scores,
    context, and anatomical-protection state for every window. The legacy
    ``supports_adaptive_exclusion`` flag still requires an ``exclude`` detector
    call, but a calibrated tier may independently require candidate-level or
    exclude-level detector support. The original score and detector call are
    not replaced.

    This helper is deliberately separate from policy application so a future
    dataset- or technology-level calibration can lock the acceptable score and
    stability thresholds without recomputing the primary scan.
    """
    if "slice_id" not in metrics:
        raise KeyError("metrics are missing required column: slice_id")
    if metrics["slice_id"].astype(str).duplicated().any():
        raise ValueError("slice_id must be unique for multi-window stability scoring")
    if minimum_corroborating_domains < 2:
        raise ValueError("minimum_corroborating_domains must be at least 2")
    if not 0 <= severe_domain_threshold <= 1:
        raise ValueError("severe_domain_threshold must be between 0 and 1")

    available_windows = sorted(
        {int(window) for window in windows if int(window) >= 3 and int(window) % 2 == 1 and int(window) <= len(metrics)}
    )
    if not available_windows:
        raise ValueError("windows must contain at least one odd width between 3 and the series length")

    active_config = config or SliceQCConfig(window=available_windows[0])
    domain_columns = [
        "density_domain_score",
        "expression_domain_score",
        "damage_domain_score",
        "continuity_domain_score",
    ]
    slice_ids = metrics["slice_id"].astype(str).tolist()
    score_matrix: list[np.ndarray] = []
    support_matrix: list[np.ndarray] = []
    exclude_matrix: list[np.ndarray] = []
    detail_by_slice: dict[str, list[dict[str, Any]]] = {slice_id: [] for slice_id in slice_ids}

    for window in available_windows:
        scored = _score_metrics(metrics, None, {}, active_config, window)
        scored = scored.set_index(scored["slice_id"].astype(str)).loc[slice_ids]
        domains = scored[domain_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        max_domain = domains.max(axis=1).to_numpy(dtype=float)
        corroborating = pd.to_numeric(scored["corroborating_domains"], errors="coerce").fillna(0).to_numpy(dtype=float)
        detector_exclude = scored["recommendation"].astype(str).eq("exclude").to_numpy()
        supported = (
            detector_exclude
            & scored["window_context"].astype(str).eq("two_sided").to_numpy()
            & ~scored["partial_structure_protection"].fillna(False).astype(bool).to_numpy()
            & (corroborating >= minimum_corroborating_domains)
            & (max_domain >= severe_domain_threshold)
        )
        scores = pd.to_numeric(scored["quality_anomaly_score"], errors="coerce").to_numpy(dtype=float)
        score_matrix.append(scores)
        support_matrix.append(supported)
        exclude_matrix.append(detector_exclude)
        for index, slice_id in enumerate(slice_ids):
            detail_by_slice[slice_id].append(
                {
                    "window": int(window),
                    "score": (float(scores[index]) if np.isfinite(scores[index]) else None),
                    "detector_call": str(scored.iloc[index]["recommendation"]),
                    "corroborating_domains": int(corroborating[index]),
                    "maximum_domain_score": float(max_domain[index]),
                    "density_domain_score": float(domains.iloc[index]["density_domain_score"]),
                    "expression_domain_score": float(domains.iloc[index]["expression_domain_score"]),
                    "damage_domain_score": float(domains.iloc[index]["damage_domain_score"]),
                    "continuity_domain_score": float(domains.iloc[index]["continuity_domain_score"]),
                    "window_context": str(scored.iloc[index]["window_context"]),
                    "partial_structure_protection": bool(scored.iloc[index]["partial_structure_protection"]),
                    "supports_adaptive_exclusion": bool(supported[index]),
                }
            )

    stacked_scores = np.vstack(score_matrix)
    out = metrics.copy()
    out["adaptive_windows_tested"] = "|".join(str(window) for window in available_windows)
    out["adaptive_window_stability"] = np.mean(np.vstack(support_matrix), axis=0)
    out["adaptive_exclude_call_fraction"] = np.mean(np.vstack(exclude_matrix), axis=0)
    out["adaptive_min_score_across_windows"] = np.nanmin(stacked_scores, axis=0)
    out["adaptive_median_score_across_windows"] = np.nanmedian(stacked_scores, axis=0)
    out["adaptive_window_details"] = [
        json.dumps(detail_by_slice[slice_id], ensure_ascii=False, separators=(",", ":")) for slice_id in slice_ids
    ]
    return out


def _choose_window(
    metrics: pd.DataFrame,
    profiles: Optional[np.ndarray],
    celltype_profiles: Mapping[str, tuple[np.ndarray, np.ndarray]],
    config: SliceQCConfig,
) -> tuple[int, list[dict[str, Any]]]:
    n = len(metrics)
    candidates = sorted(
        {
            int(value)
            for value in config.window_candidates
            if int(value) >= 3 and int(value) % 2 == 1 and int(value) <= n
        }
    )
    if not candidates:
        return 3, [{"window": 3, "utility": None, "note": "fewer slices than candidate widths"}]
    trials: list[dict[str, Any]] = []
    for window in candidates:
        baseline = _score_metrics(metrics, profiles, celltype_profiles, config, window)
        base_flags = float(np.mean(baseline["recommendation"].ne("keep")))
        eligible = np.arange(max(window // 2, 1), max(n - window // 2, 1))
        if eligible.size == 0:
            eligible = np.arange(n)
        if eligible.size > 8:
            eligible = eligible[np.linspace(0, eligible.size - 1, 8, dtype=int)]
        hits = []
        deltas = []
        for target in eligible:
            perturbed = metrics.copy()
            for key in (
                "median_total_counts",
                "median_n_genes",
                "library_complexity",
                "detected_gene_fraction",
            ):
                if key in perturbed:
                    perturbed.loc[target, key] = float(perturbed.loc[target, key]) * 0.24
            if "zero_fraction" in perturbed:
                perturbed.loc[target, "zero_fraction"] = min(
                    0.999, float(perturbed.loc[target, "zero_fraction"]) + 0.35
                )
            scored = _score_metrics(perturbed, profiles, celltype_profiles, config, window)
            hits.append(scored.loc[target, "recommendation"] in {"review", "exclude"})
            deltas.append(
                float(scored.loc[target, "quality_anomaly_score"] - baseline.loc[target, "quality_anomaly_score"])
            )
        recall = float(np.mean(hits)) if hits else 0.0
        mean_delta = float(np.mean(deltas)) if deltas else 0.0
        utility = recall + 0.35 * mean_delta - 0.35 * base_flags - 0.012 * (window - 3)
        trials.append(
            {
                "window": int(window),
                "synthetic_recall": recall,
                "mean_target_score_delta": mean_delta,
                "baseline_nonkeep_rate": base_flags,
                "utility": utility,
            }
        )
    best = max(trials, key=lambda item: (item["utility"], -item["window"]))
    return int(best["window"]), trials


def scan_h5ad_series(
    inputs: Sequence[Union[str, Path]],
    *,
    config: Optional[SliceQCConfig] = None,
    sort_inputs: bool = True,
    slice_key: Optional[str] = "auto",
    spatial_key: Optional[str] = "auto",
    x_key: Optional[str] = None,
    y_key: Optional[str] = None,
    order_key: Optional[str] = None,
    layer: Optional[str] = "auto",
    celltype_key: Optional[str] = "auto",
    profile_gene_names: Optional[Sequence[str]] = None,
) -> SliceSeriesResult:
    """Scan one multi-slice H5AD or a series of H5AD files.

    Multiple file paths are naturally sorted by filename unless
    ``sort_inputs=False`` is used for a caller-defined anatomical order. When
    several files contain more than one internal slice, slice identifiers
    are prefixed with the filename to prevent collisions.  For one-file-per-
    slice input, the filename stem is used when no explicit slice annotation is
    present.
    """
    config = config or SliceQCConfig()
    paths = [Path(value).expanduser().resolve() for value in inputs]
    if not paths:
        raise ValueError("At least one H5AD input is required")
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing input(s): " + ", ".join(missing))
    if sort_inputs and len(paths) > 1:
        paths = sorted(paths, key=lambda path: _natural_key(path.name))

    all_metrics: list[dict[str, Any]] = []
    point_samples: dict[str, dict[str, list[Any]]] = {}
    pseudobulk: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    celltype_profiles: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    sources: list[dict[str, Any]] = []
    for file_index, path in enumerate(paths):
        extracted = _extract_one_file(
            path,
            config,
            slice_key,
            spatial_key,
            x_key,
            y_key,
            order_key,
            layer,
            celltype_key,
            file_prefix=len(paths) > 1,
        )
        for row in extracted.metrics:
            row["_input_file_index"] = int(file_index)
        all_metrics.extend(extracted.metrics)
        overlap = set(point_samples).intersection(extracted.point_samples)
        if overlap:
            raise ValueError(f"Duplicate slice identifiers after input resolution: {sorted(overlap)}")
        point_samples.update(extracted.point_samples)
        pseudobulk.update(extracted.pseudobulk)
        celltype_profiles.update(extracted.celltype_profiles)
        sources.append(extracted.source_info)

    metrics = pd.DataFrame(all_metrics)
    metrics["_natural"] = metrics["slice_id"].map(_natural_key)
    metrics = metrics.sort_values(["_input_file_index", "slice_order_value", "_natural"], kind="stable").drop(
        columns=["_input_file_index", "_natural"]
    )
    metrics = metrics.reset_index(drop=True)
    metrics["slice_index"] = np.arange(len(metrics), dtype=int)
    slice_order = metrics["slice_id"].astype(str).tolist()
    profiles, profile_genes = _align_profiles(
        pseudobulk, slice_order, config.profile_genes, fixed_genes=profile_gene_names
    )

    if isinstance(config.window, str):
        if config.window.lower() != "auto":
            raise ValueError("config.window must be an odd integer or 'auto'")
        selected_window, window_trials = _choose_window(metrics, profiles, celltype_profiles, config)
    else:
        selected_window = int(config.window)
        window_trials = [{"window": selected_window, "selected": True, "mode": "fixed"}]
    scored = _score_metrics(metrics, profiles, celltype_profiles, config, selected_window)
    provenance = {
        "schema_version": SCHEMA_VERSION,
        "method": "multi-domain local-window pre-alignment slice QC",
        "config": _jsonable(asdict(config)),
        "selected_window": int(selected_window),
        "window_trials": _jsonable(window_trials),
        "sources": sources,
        "expression_profile_gene_policy": (
            "fixed_from_reference" if profile_gene_names is not None else "selected_by_variance"
        ),
        "expression_profile_gene_count": int(len(profile_genes)),
        "observable_depth_note": (
            "H5AD count matrices provide library-size/captured-molecule and detected-gene proxies. "
            "Read depth, duplication and sequencing saturation are unavailable unless separately supplied."
        ),
        "recommendation_policy": {
            "keep": "no corroborated local anomaly",
            "review": "ambiguous, single-domain, partial-structure, or endpoint evidence",
            "exclude": "strong expression evidence or corroborated multi-domain evidence",
            "automatic_deletion": False,
        },
    }
    return SliceSeriesResult(scored, point_samples, provenance, profiles, profile_genes)


def calculate_slice_quality(
    adata: AnnData,
    *,
    slice_key: Optional[str] = "auto",
    spatial_key: Optional[str] = "auto",
    x_key: Optional[str] = None,
    y_key: Optional[str] = None,
    order_key: Optional[str] = None,
    layer: Optional[str] = "auto",
    celltype_key: Optional[str] = "auto",
    profile_gene_names: Optional[Sequence[str]] = None,
    config: Optional[SliceQCConfig] = None,
    inplace: bool = False,
) -> pd.DataFrame:
    """Calculate pre-alignment slice quality for an in-memory AnnData object.

    The implementation uses a temporary H5AD only when the caller passes an
    in-memory object, keeping file and in-memory behavior identical.  No input
    coordinates or observations are removed.  With ``inplace=True``, only a
    JSON-compatible summary is added to ``adata.uns['slice_quality_qc']``.
    """
    import tempfile

    config = config or SliceQCConfig()
    with tempfile.TemporaryDirectory(prefix="spateo_slice_qc_") as tmp:
        path = Path(tmp) / "input.h5ad"
        adata.write_h5ad(path)
        result = scan_h5ad_series(
            [path],
            config=config,
            slice_key=slice_key,
            spatial_key=spatial_key,
            x_key=x_key,
            y_key=y_key,
            order_key=order_key,
            layer=layer,
            celltype_key=celltype_key,
            profile_gene_names=profile_gene_names,
        )
    if inplace:
        adata.uns["slice_quality_qc"] = {
            "provenance": _jsonable(result.provenance),
            "slice_calls": _jsonable(
                result.metrics[["slice_id", "quality_anomaly_score", "recommendation", "reason"]].to_dict("records")
            ),
        }
    return result.metrics


def scan_h5ad_collection(
    datasets: Mapping[str, Union[str, Path, Sequence[Union[str, Path]]]],
    *,
    config: Optional[SliceQCConfig] = None,
    dataset_options: Optional[Mapping[str, Mapping[str, Any]]] = None,
) -> dict[str, SliceSeriesResult]:
    """Scan multiple independent serial-slice datasets without merging them.

    Args:
        datasets: Ordered mapping from a dataset id to one multi-slice H5AD or
            a sequence of per-slice H5AD files belonging to that dataset.
        config: Common detector configuration. A dataset may override it with
            ``dataset_options[dataset_id]["config"]``.
        dataset_options: Optional per-dataset keyword arguments accepted by
            :func:`scan_h5ad_series`, such as ``slice_key``, ``spatial_key``,
            ``order_key``, ``layer``, ``celltype_key`` and ``sort_inputs``.

    Returns:
        One independent :class:`SliceSeriesResult` per dataset id, preserving
        input order.

    Notes:
        Dataset boundaries are never crossed when computing local expectations.
        Concatenating unrelated datasets would create invalid neighbor pairs.
    """
    if not isinstance(datasets, Mapping) or not datasets:
        raise ValueError("datasets must be a non-empty mapping of dataset_id to H5AD path(s)")
    options_by_dataset = {str(key).strip(): dict(value) for key, value in (dataset_options or {}).items()}
    dataset_ids = {str(key).strip() for key in datasets}
    unknown_dataset_options = set(options_by_dataset).difference(dataset_ids)
    if unknown_dataset_options:
        raise KeyError(f"dataset_options contains unknown dataset ids: {sorted(unknown_dataset_options)}")
    allowed_options = {
        "slice_key",
        "spatial_key",
        "x_key",
        "y_key",
        "order_key",
        "layer",
        "celltype_key",
        "profile_gene_names",
        "sort_inputs",
        "config",
    }
    results: dict[str, SliceSeriesResult] = {}
    for raw_dataset_id, raw_paths in datasets.items():
        dataset_id = str(raw_dataset_id).strip()
        if not dataset_id:
            raise ValueError("dataset ids must be non-empty strings")
        if dataset_id in results:
            raise ValueError(f"duplicate dataset id after string normalization: {dataset_id!r}")
        paths = [raw_paths] if isinstance(raw_paths, (str, Path)) else list(raw_paths)
        if not paths:
            raise ValueError(f"dataset {dataset_id!r} has no H5AD inputs")
        options = dict(options_by_dataset.get(dataset_id, {}))
        unexpected = set(options).difference(allowed_options)
        if unexpected:
            raise TypeError(f"unsupported options for dataset {dataset_id!r}: {sorted(unexpected)}")
        dataset_config = options.pop("config", config)
        results[dataset_id] = scan_h5ad_series(
            paths,
            config=dataset_config or SliceQCConfig(),
            **options,
        )
    return results


def _binomial_thin_matrix(
    matrix: Union[np.ndarray, sparse.spmatrix],
    row_mask: np.ndarray,
    rate: float,
    rng: np.random.Generator,
) -> Union[np.ndarray, sparse.csr_matrix]:
    if not 0 <= rate <= 1:
        raise ValueError("thinning rate must be between 0 and 1")
    if sparse.issparse(matrix):
        out = matrix.tocsr(copy=True)
        rows = np.flatnonzero(row_mask)
        for row in rows:
            start, end = out.indptr[row], out.indptr[row + 1]
            if start < end:
                counts = np.rint(np.maximum(out.data[start:end], 0)).astype(np.int64)
                out.data[start:end] = rng.binomial(counts, rate)
        out.eliminate_zeros()
        return out
    out = np.array(matrix, copy=True)
    counts = np.rint(np.maximum(out[row_mask], 0)).astype(np.int64)
    out[row_mask] = rng.binomial(counts, rate)
    return out


def simulate_slice_quality_artifacts(
    adata: AnnData,
    artifacts: Sequence[Mapping[str, Any]],
    *,
    slice_key: str,
    spatial_key: str = "spatial",
    layer: Optional[str] = "auto",
    output_layer: str = "slice_qc_simulated_counts",
    random_seed: int = 13,
) -> AnnData:
    """Return a copy with controlled slice defects and provenance.

    Supported artifact records are:

    ``{"slice_id": ..., "kind": "depth", "rate": 0.2}``
        Binomially thin every count in the target slice.

    ``{"slice_id": ..., "kind": "regional_depth", "rate": 0.1,
       "radius_quantile": 0.35}``
        Thin counts only inside a deterministic central region.

    ``{"slice_id": ..., "kind": "cell_dropout", "keep_rate": 0.35}``
        Randomly retain a fraction of locations.

    ``{"slice_id": ..., "kind": "spatial_hole", "radius_quantile": 0.35}``
        Remove locations in a central elliptical region.

    Source layers are not overwritten.  The simulated counts are stored in
    ``output_layer`` and any removed observations are applied to the returned
    copy only.
    """
    if slice_key not in adata.obs:
        raise KeyError(f"slice_key={slice_key!r} is missing from adata.obs")
    if spatial_key not in adata.obsm:
        raise KeyError(f"spatial_key={spatial_key!r} is missing from adata.obsm")
    resolved_layer, matrix, count_like = _choose_layer(adata, layer)
    if not count_like:
        raise ValueError("Simulation requires a non-negative integer-like count matrix/layer")
    labels = np.asarray(adata.obs[slice_key].astype(str))
    coords = np.asarray(adata.obsm[spatial_key])[:, :2].astype(float)
    rng = np.random.default_rng(random_seed)
    keep = np.ones(adata.n_obs, dtype=bool)
    expression_actions: list[tuple[np.ndarray, float]] = []
    normalized_plan: list[dict[str, Any]] = []
    for artifact in artifacts:
        item = dict(artifact)
        slice_id = str(item["slice_id"])
        kind = str(item["kind"])
        target = labels == slice_id
        if not target.any():
            raise KeyError(f"Simulation target slice {slice_id!r} was not found")
        if kind == "depth":
            rate = float(item.get("rate", 0.2))
            expression_actions.append((target.copy(), rate))
        elif kind in {"regional_depth", "spatial_hole"}:
            target_idx = np.flatnonzero(target)
            xy = coords[target_idx]
            center = np.median(xy, axis=0)
            scale = np.maximum(np.quantile(np.abs(xy - center), 0.75, axis=0), EPS)
            radius = np.sqrt(np.sum(((xy - center) / scale) ** 2, axis=1))
            quantile = float(item.get("radius_quantile", 0.35))
            region_local = radius <= np.quantile(radius, np.clip(quantile, 0.05, 0.9))
            region = np.zeros(adata.n_obs, dtype=bool)
            region[target_idx[region_local]] = True
            item["affected_locations"] = int(region.sum())
            if kind == "regional_depth":
                expression_actions.append((region, float(item.get("rate", 0.1))))
            else:
                keep[region] = False
        elif kind == "cell_dropout":
            keep_rate = float(item.get("keep_rate", 0.35))
            if not 0 < keep_rate <= 1:
                raise ValueError("keep_rate must be in (0, 1]")
            target_idx = np.flatnonzero(target)
            retained = rng.choice(
                target_idx,
                max(1, int(round(target_idx.size * keep_rate))),
                replace=False,
            )
            keep[target_idx] = False
            keep[retained] = True
            item["retained_locations"] = int(retained.size)
        else:
            raise ValueError(f"Unsupported simulation artifact kind: {kind!r}")
        normalized_plan.append(_jsonable(item))

    simulated = matrix
    for row_mask, rate in expression_actions:
        simulated = _binomial_thin_matrix(simulated, row_mask, rate, rng)
    out = adata[keep].copy()
    # AnnData treats keys ending in ``_colors`` specially while slicing.  Some
    # valid source files store a color mapping dict there, which older AnnData
    # releases wrap into a one-element object array that can no longer be
    # written.  Restore only mappings affected by that conversion in the
    # synthetic copy; the source object remains untouched.
    import copy

    for uns_key, original_value in adata.uns.items():
        copied_value = out.uns.get(uns_key)
        if isinstance(original_value, Mapping) and isinstance(copied_value, np.ndarray):
            if copied_value.dtype == object and copied_value.size == 1 and isinstance(copied_value.flat[0], Mapping):
                out.uns[uns_key] = copy.deepcopy(original_value)
    raw_var_index_repaired = False
    if out.raw is not None and "_index" in out.raw.var.columns:
        raw_adata = out.raw.to_adata()
        recovered_names = raw_adata.var["_index"].astype(str).to_numpy()
        raw_adata.var = raw_adata.var.drop(columns="_index")
        raw_adata.var_names = pd.Index(recovered_names)
        raw_adata.var_names_make_unique()
        out.raw = raw_adata
        raw_var_index_repaired = True
    out.layers[output_layer] = simulated[keep]
    sim_totals = _matrix_row_sum(out.layers[output_layer])
    sim_genes = _matrix_row_nnz(out.layers[output_layer])
    out.obs["slice_qc_sim_total_counts"] = sim_totals
    out.obs["slice_qc_sim_n_genes"] = sim_genes.astype(int)
    out.uns["slice_quality_simulation"] = {
        "schema_version": SCHEMA_VERSION,
        "random_seed": int(random_seed),
        "source_layer": resolved_layer,
        "output_layer": output_layer,
        "artifacts_json": json.dumps(normalized_plan, ensure_ascii=False, sort_keys=True),
        "raw_var_index_repaired_for_h5ad": bool(raw_var_index_repaired),
        "ground_truth_slices": np.asarray(sorted({str(item["slice_id"]) for item in normalized_plan}), dtype=str),
    }
    return out


def evaluate_slice_calls(metrics: pd.DataFrame, ground_truth: Mapping[str, bool]) -> dict[str, Any]:
    """Evaluate suspicious (review or exclude) calls against explicit ground truth."""
    truth = np.asarray([bool(ground_truth.get(str(value), False)) for value in metrics["slice_id"]])
    calls = metrics["recommendation"].astype(str)
    predicted = calls.ne("keep").to_numpy()
    tp = int(np.sum(truth & predicted))
    fp = int(np.sum(~truth & predicted))
    fn = int(np.sum(truth & ~predicted))
    tn = int(np.sum(~truth & ~predicted))
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, EPS)
    exclude = calls.eq("exclude").to_numpy()
    return {
        "n_slices": int(len(metrics)),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "positive_definition": "recommendation in {review, exclude}",
        "exclude_tp": int(np.sum(truth & exclude)),
        "exclude_fp": int(np.sum(~truth & exclude)),
        "exclude_candidate_count": int(np.sum(exclude)),
    }


def evaluate_paired_simulation(
    baseline_metrics: pd.DataFrame,
    simulated_metrics: pd.DataFrame,
    artifacts: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, Any], pd.DataFrame]:
    """Compare synthetic calls with the same slices before perturbation.

    Real input data can already contain native QC candidates.  A paired
    comparison therefore avoids treating every non-injected slice as a known
    high-quality negative.  It reports whether injected slices are detected,
    whether their anomaly scores increase, and whether any non-injected slice
    becomes newly suspicious after the simulation.
    """
    required = {"slice_id", "quality_anomaly_score", "recommendation"}
    for name, frame in (
        ("baseline", baseline_metrics),
        ("simulated", simulated_metrics),
    ):
        missing = required.difference(frame.columns)
        if missing:
            raise KeyError(f"{name} metrics are missing required columns: {sorted(missing)}")
        if frame["slice_id"].astype(str).duplicated().any():
            raise ValueError(f"{name} metrics contain duplicate slice_id values")

    baseline = baseline_metrics[list(required)].copy()
    simulated = simulated_metrics[list(required)].copy()
    baseline["slice_id"] = baseline["slice_id"].astype(str)
    simulated["slice_id"] = simulated["slice_id"].astype(str)
    baseline = baseline.rename(
        columns={
            "quality_anomaly_score": "baseline_score",
            "recommendation": "baseline_recommendation",
        }
    )
    simulated = simulated.rename(
        columns={
            "quality_anomaly_score": "simulated_score",
            "recommendation": "simulated_recommendation",
        }
    )
    paired = baseline.merge(simulated, on="slice_id", how="outer", validate="one_to_one", indicator=True)
    if not paired["_merge"].eq("both").all():
        missing = paired.loc[paired["_merge"].ne("both"), ["slice_id", "_merge"]].to_dict("records")
        raise ValueError(f"Baseline and simulated slice sets differ: {missing}")
    paired = paired.drop(columns="_merge")

    artifact_map: dict[str, list[str]] = {}
    for artifact in artifacts:
        artifact_map.setdefault(str(artifact["slice_id"]), []).append(str(artifact["kind"]))
    paired["injected"] = paired["slice_id"].isin(artifact_map)
    paired["artifact"] = paired["slice_id"].map(lambda value: "+".join(artifact_map.get(str(value), [])) or "none")
    paired["score_delta"] = paired["simulated_score"] - paired["baseline_score"]
    paired["baseline_suspicious"] = paired["baseline_recommendation"].ne("keep")
    paired["detected_after"] = paired["simulated_recommendation"].ne("keep")
    paired["newly_flagged"] = ~paired["baseline_suspicious"] & paired["detected_after"]
    rank = {"keep": 0, "review": 1, "exclude": 2}
    paired["call_worsened"] = paired["simulated_recommendation"].map(rank).fillna(-1) > paired[
        "baseline_recommendation"
    ].map(rank).fillna(-1)

    injected = paired[paired["injected"]]
    eligible = injected[~injected["baseline_suspicious"]]
    noninjected = paired[~paired["injected"]]
    median_delta = float(injected["score_delta"].median()) if len(injected) else float("nan")
    summary = {
        "comparison": "paired baseline versus simulated counts using the same slice window",
        "injected_total": int(len(injected)),
        "injected_detected_after": int(injected["detected_after"].sum()),
        "injected_detection_rate_after": (float(injected["detected_after"].mean()) if len(injected) else 0.0),
        "injected_preexisting_nonkeep": int(injected["baseline_suspicious"].sum()),
        "injected_baseline_keep": int(len(eligible)),
        "injected_newly_flagged": int(eligible["newly_flagged"].sum()),
        "injected_new_flag_rate": (float(eligible["newly_flagged"].mean()) if len(eligible) else None),
        "injected_score_increased": int((injected["score_delta"] > 0).sum()),
        "injected_median_score_delta": median_delta,
        "injected_call_worsened": int(injected["call_worsened"].sum()),
        "noninjected_newly_flagged": int(noninjected["newly_flagged"].sum()),
        "baseline_nonkeep": int(paired["baseline_suspicious"].sum()),
        "simulated_nonkeep": int(paired["detected_after"].sum()),
        "interpretation": (
            "Native baseline candidates are not labelled false positives. "
            "Review injected detection, score deltas, and non-injected newly flagged slices separately."
        ),
    }
    return _jsonable(summary), paired.sort_values("slice_id", key=lambda series: series.map(_natural_key))


def _validate_review_tiers(
    tiers: Sequence[ReviewEvidenceTier],
    keep_max_score: float,
) -> None:
    if tiers and abs(tiers[0].min_score - keep_max_score) > EPS:
        raise ValueError("the first review evidence tier must start at keep_max_score")
    names: set[str] = set()
    previous_upper: Optional[float] = None
    previous_tier: Optional[ReviewEvidenceTier] = None
    detector_strength = {"ignore": 0, "candidate": 1, "exclude": 2}
    for index, tier in enumerate(tiers):
        if not tier.name or tier.name in names:
            raise ValueError("review evidence tier names must be non-empty and unique")
        names.add(tier.name)
        if not keep_max_score <= tier.min_score <= 1:
            raise ValueError("review evidence tier min_score must be between keep_max_score and 1")
        if tier.max_score is not None and not tier.min_score < tier.max_score <= 1:
            raise ValueError("review evidence tier max_score must be greater than min_score and at most 1")
        if tier.min_corroborating_domains < 2 or tier.min_corroborating_domains > 4:
            raise ValueError("review evidence tiers require between two and four corroborating domains")
        for name, value in {
            "severe_domain_threshold": tier.severe_domain_threshold,
            "min_window_stability": tier.min_window_stability,
            "min_score_confidence": tier.min_score_confidence,
        }.items():
            if not 0 <= value <= 1:
                raise ValueError(f"review evidence tier {name} must be between 0 and 1")
        if tier.window_detector_mode not in {"candidate", "exclude", "ignore"}:
            raise ValueError("window_detector_mode must be candidate, exclude, or ignore")
        if index > 0 and previous_upper is None:
            raise ValueError("only the final review evidence tier may omit max_score")
        if previous_upper is not None and abs(tier.min_score - previous_upper) > EPS:
            raise ValueError("review evidence tiers must form contiguous score bands")
        if previous_tier is not None:
            if previous_tier.enable_exclude and not tier.enable_exclude:
                raise ValueError("an exclude-enabled tier cannot precede a higher-score keep-only tier")
            if not previous_tier.enable_exclude and not tier.enable_exclude:
                raise ValueError("adjacent keep-only review tiers must be merged")
            if previous_tier.enable_exclude:
                lower_not_weaker = (
                    previous_tier.min_corroborating_domains >= tier.min_corroborating_domains
                    and previous_tier.severe_domain_threshold + EPS >= tier.severe_domain_threshold
                    and previous_tier.min_window_stability + EPS >= tier.min_window_stability
                    and previous_tier.min_score_confidence + EPS >= tier.min_score_confidence
                    and detector_strength[previous_tier.window_detector_mode]
                    >= detector_strength[tier.window_detector_mode]
                )
                lower_strictly_stronger = (
                    previous_tier.min_corroborating_domains > tier.min_corroborating_domains
                    or previous_tier.severe_domain_threshold > tier.severe_domain_threshold + EPS
                    or previous_tier.min_window_stability > tier.min_window_stability + EPS
                    or previous_tier.min_score_confidence > tier.min_score_confidence + EPS
                    or detector_strength[previous_tier.window_detector_mode]
                    > detector_strength[tier.window_detector_mode]
                )
                if not lower_not_weaker or not lower_strictly_stronger:
                    raise ValueError(
                        "each lower-score review tier must require strictly stronger evidence "
                        "than the adjacent higher-score tier"
                    )
        previous_upper = tier.max_score
        previous_tier = tier
    if tiers and tiers[-1].max_score is not None:
        raise ValueError("the final review evidence tier must have max_score=None")


def _detector_mode_passes(call: str, mode: str, required: bool) -> bool:
    if not required or mode == "ignore":
        return True
    if mode == "exclude":
        return call == "exclude"
    return call in {"review", "exclude"}


def _evaluate_review_tier_row(
    row: pd.Series,
    tier: ReviewEvidenceTier,
    *,
    require_two_sided: bool,
    require_detector_call: bool,
) -> tuple[bool, float, list[str]]:
    if not tier.enable_exclude:
        return False, 0.0, ["calibrated keep-only score band"]
    failed: list[str] = []
    score = pd.to_numeric(pd.Series([row.get("quality_anomaly_score")]), errors="coerce").iloc[0]
    if not np.isfinite(score) or float(score) < tier.min_score:
        failed.append("score below tier floor")
    if tier.max_score is not None and np.isfinite(score) and float(score) >= tier.max_score:
        failed.append("score outside tier band")
    call = str(row.get("recommendation", ""))
    if not _detector_mode_passes(call, tier.window_detector_mode, require_detector_call):
        failed.append("primary detector call does not support tier")
    if require_two_sided and str(row.get("window_context", "")) != "two_sided":
        failed.append("two-sided context unavailable")
    if bool(row.get("partial_structure_protection", False)):
        failed.append("partial/anatomical structure protection")
    corroborating = pd.to_numeric(pd.Series([row.get("corroborating_domains")]), errors="coerce").iloc[0]
    if not np.isfinite(corroborating) or int(corroborating) < tier.min_corroborating_domains:
        failed.append("fewer than tier corroborating domains")
    domain_values = pd.to_numeric(
        pd.Series(
            [
                row.get("density_domain_score"),
                row.get("expression_domain_score"),
                row.get("damage_domain_score"),
                row.get("continuity_domain_score"),
            ]
        ),
        errors="coerce",
    ).fillna(0.0)
    if float(domain_values.max()) < tier.severe_domain_threshold:
        failed.append("no domain reaches tier severity")
    confidence = pd.to_numeric(pd.Series([row.get("score_confidence")]), errors="coerce").iloc[0]
    if not np.isfinite(confidence) or float(confidence) < tier.min_score_confidence:
        failed.append("insufficient score confidence")

    try:
        raw_details = row.get("adaptive_window_details", "[]")
        details = json.loads(raw_details) if isinstance(raw_details, str) else list(raw_details)
    except (TypeError, ValueError, json.JSONDecodeError):
        details = []
    supported = 0
    for detail in details:
        window_score = pd.to_numeric(pd.Series([detail.get("score")]), errors="coerce").iloc[0]
        window_ok = bool(np.isfinite(window_score) and float(window_score) >= tier.min_score)
        window_ok &= _detector_mode_passes(
            str(detail.get("detector_call", "")),
            tier.window_detector_mode,
            require_detector_call,
        )
        if require_two_sided:
            window_ok &= str(detail.get("window_context", "")) == "two_sided"
        window_ok &= not bool(detail.get("partial_structure_protection", False))
        window_ok &= int(detail.get("corroborating_domains", 0)) >= tier.min_corroborating_domains
        maximum_domain = pd.to_numeric(pd.Series([detail.get("maximum_domain_score")]), errors="coerce").iloc[0]
        window_ok &= bool(np.isfinite(maximum_domain) and float(maximum_domain) >= tier.severe_domain_threshold)
        supported += int(window_ok)
    stability = supported / len(details) if details else 0.0
    if stability + EPS < tier.min_window_stability:
        failed.append("multi-window support below tier minimum")
    return not failed, float(stability), failed


def audit_review_evidence(
    metrics: pd.DataFrame, tier: ReviewEvidenceTier
) -> pd.DataFrame:
    """Inspect a review band without granting it publication permission.

    Reuses the complete production evidence checks, including every available
    confirming window. A candidate may trade severity for corroboration and is
    deliberately separate from the legacy monotonic publication validator.
    Existing final calls and scores are preserved verbatim. Passing this audit
    is experimental evidence, never certification or a final exclude action.
    """
    required = {"threshold_triage_call", "quality_anomaly_score", "recommendation",
                "window_context", "partial_structure_protection", "corroborating_domains",
                "score_confidence", "adaptive_window_details", "density_domain_score",
                "expression_domain_score", "damage_domain_score", "continuity_domain_score"}
    missing = required.difference(metrics.columns)
    if missing:
        raise ValueError(f"Review evidence audit requires columns: {sorted(missing)}")
    if not tier.enable_exclude:
        raise ValueError("Supply an enabled evidence candidate; keep-only short circuits are not evidence audits")
    score = pd.to_numeric(metrics["quality_anomaly_score"], errors="coerce")
    selected = metrics["threshold_triage_call"].eq("review") & score.ge(tier.min_score)
    if tier.max_score is not None:
        selected &= score.lt(tier.max_score)
    out = metrics.loc[selected].copy(deep=True)
    evaluated = [_evaluate_review_tier_row(row, tier, require_two_sided=True,
                                         require_detector_call=True) for _, row in out.iterrows()]
    out["evidence_audit_candidate"] = tier.name
    out["evidence_audit_status"] = "experimental_not_certified"
    out["evidence_audit_performed"] = True
    out["evidence_candidate_passes"] = [item[0] for item in evaluated]
    out["evidence_window_support"] = [item[1] for item in evaluated]
    out["evidence_failed_conditions"] = [json.dumps(item[2]) for item in evaluated]
    out["evidence_policy_permission"] = False
    return out


def apply_high_confidence_policy(
    metrics: pd.DataFrame,
    policy: Union[HighConfidencePolicy, Mapping[str, Any]],
) -> pd.DataFrame:
    """Apply threshold triage, resolve review rows, and publish binary actions.

    The audit table records four distinct concepts:

    - ``threshold_band``: pure score-based keep/review/exclude partition;
    - ``threshold_triage_call``: the stage-1 partition after exclusion
      guardrails can downgrade an unsafe exclude to review;
    - ``review_resolution``: the stage-2 keep/exclude/withhold result applied
      only to rows in the review queue;
    - ``final_call``: the public operational action.

    With ``unresolved_action='keep'``, every row receives a public keep/exclude
    action.  Intermediate review evidence remains auditable but is never a
    third user-facing action.
    """
    if not isinstance(policy, HighConfidencePolicy):
        policy = HighConfidencePolicy.from_mapping(policy)
    required = {"slice_id", "quality_anomaly_score", "recommendation"}
    missing = required.difference(metrics.columns)
    if missing:
        raise KeyError(f"metrics are missing required columns: {sorted(missing)}")
    if not 0 <= policy.keep_max_score < policy.exclude_min_score <= 1:
        raise ValueError("policy thresholds must satisfy 0 <= keep < exclude <= 1")
    if policy.unresolved_action not in {"withhold", "keep"}:
        raise ValueError("unresolved_action must be either 'withhold' or 'keep'")
    tiered_enabled = bool(policy.review_exclusion_tiers)
    legacy_adaptive_enabled = policy.adaptive_exclude_min_score is not None
    if tiered_enabled:
        _validate_review_tiers(policy.review_exclusion_tiers, policy.keep_max_score)
        tiered_required = {
            "corroborating_domains",
            "score_confidence",
            "adaptive_window_details",
            "density_domain_score",
            "expression_domain_score",
            "damage_domain_score",
            "continuity_domain_score",
        }
        tiered_missing = tiered_required.difference(metrics.columns)
        if tiered_missing:
            raise KeyError("tiered review policy requires multi-window evidence columns: " f"{sorted(tiered_missing)}")
    elif legacy_adaptive_enabled:
        assert policy.adaptive_exclude_min_score is not None
        if not 0 <= policy.adaptive_exclude_min_score < policy.exclude_min_score:
            raise ValueError("adaptive exclude threshold must be lower than the standard exclude threshold")
        if policy.adaptive_min_corroborating_domains < 2:
            raise ValueError("adaptive exclusion requires at least two corroborating domains")
        for name, value in {
            "adaptive_severe_domain_threshold": policy.adaptive_severe_domain_threshold,
            "adaptive_min_window_stability": policy.adaptive_min_window_stability,
            "adaptive_min_score_confidence": policy.adaptive_min_score_confidence,
        }.items():
            if not 0 <= value <= 1:
                raise ValueError(f"{name} must be between 0 and 1")
        adaptive_required = {
            "corroborating_domains",
            "score_confidence",
            "adaptive_window_stability",
            "adaptive_min_score_across_windows",
            "density_domain_score",
            "expression_domain_score",
            "damage_domain_score",
            "continuity_domain_score",
        }
        adaptive_missing = adaptive_required.difference(metrics.columns)
        if adaptive_missing:
            raise KeyError("adaptive policy requires multi-window evidence columns: " f"{sorted(adaptive_missing)}")

    out = metrics.copy()
    score = pd.to_numeric(out["quality_anomaly_score"], errors="coerce").to_numpy(dtype=float)
    recommendation = out["recommendation"].astype(str).to_numpy()
    if "window_context" in out:
        two_sided = out["window_context"].astype(str).eq("two_sided").to_numpy()
    else:
        two_sided = np.ones(len(out), dtype=bool)
    if "partial_structure_protection" in out:
        protected = out["partial_structure_protection"].fillna(False).astype(bool).to_numpy()
    else:
        protected = np.zeros(len(out), dtype=bool)
    context_ok = two_sided if policy.require_two_sided else np.ones(len(out), dtype=bool)
    detector_keep = recommendation == "keep" if policy.require_detector_call else np.ones(len(out), dtype=bool)
    detector_exclude = recommendation == "exclude" if policy.require_detector_call else np.ones(len(out), dtype=bool)
    finite_score = np.isfinite(score)
    threshold_band = np.full(len(out), "review", dtype=object)
    threshold_band[finite_score & (score <= policy.keep_max_score)] = "keep"
    threshold_band[finite_score & (score >= policy.exclude_min_score)] = "exclude"
    publish_keep = (
        policy.enable_keep & finite_score & (score <= policy.keep_max_score) & detector_keep & context_ok & ~protected
    )
    standard_exclude = (
        policy.enable_exclude
        & finite_score
        & (score >= policy.exclude_min_score)
        & detector_exclude
        & context_ok
        & ~protected
    )
    threshold_triage_call = threshold_band.copy()
    threshold_triage_call[(threshold_band == "exclude") & ~standard_exclude] = "review"
    review_queue = threshold_triage_call == "review"
    adaptive_exclude = np.zeros(len(out), dtype=bool)
    matched_tier = np.full(len(out), "", dtype=object)
    matched_tier_stability = np.zeros(len(out), dtype=float)
    tier_diagnostics: list[list[dict[str, Any]]] = [[] for _ in range(len(out))]
    if tiered_enabled:
        for index in np.flatnonzero(review_queue):
            row = out.iloc[index]
            for tier in policy.review_exclusion_tiers:
                passed, stability, failed = _evaluate_review_tier_row(
                    row,
                    tier,
                    require_two_sided=policy.require_two_sided,
                    require_detector_call=policy.require_detector_call,
                )
                tier_diagnostics[index].append(
                    {
                        "tier": tier.name,
                        "passed": bool(passed),
                        "window_stability": float(stability),
                        "failed": failed,
                    }
                )
                score_in_band = bool(
                    finite_score[index]
                    and score[index] >= tier.min_score
                    and (tier.max_score is None or score[index] < tier.max_score)
                )
                if score_in_band:
                    matched_tier[index] = tier.name
                    matched_tier_stability[index] = stability
                    adaptive_exclude[index] = bool(policy.enable_exclude and tier.enable_exclude and passed)
                    break
        adaptive_exclude &= ~standard_exclude
    elif legacy_adaptive_enabled:
        assert policy.adaptive_exclude_min_score is not None
        domain_values = (
            out[
                [
                    "density_domain_score",
                    "expression_domain_score",
                    "damage_domain_score",
                    "continuity_domain_score",
                ]
            ]
            .apply(pd.to_numeric, errors="coerce")
            .fillna(0.0)
            .to_numpy(dtype=float)
        )
        corroborating = pd.to_numeric(out["corroborating_domains"], errors="coerce").fillna(0).to_numpy()
        score_confidence = pd.to_numeric(out["score_confidence"], errors="coerce").fillna(0).to_numpy()
        window_stability = pd.to_numeric(out["adaptive_window_stability"], errors="coerce").fillna(0).to_numpy()
        multiscale_min_score = (
            pd.to_numeric(out["adaptive_min_score_across_windows"], errors="coerce").fillna(-np.inf).to_numpy()
        )
        detector_candidate = (
            np.isin(recommendation, ["review", "exclude"])
            if policy.require_detector_call
            else np.ones(len(out), dtype=bool)
        )
        adaptive_exclude = (
            policy.enable_exclude
            & finite_score
            & (score >= policy.adaptive_exclude_min_score)
            & (multiscale_min_score >= policy.adaptive_exclude_min_score)
            & review_queue
            & detector_candidate
            & context_ok
            & ~protected
            & (corroborating >= policy.adaptive_min_corroborating_domains)
            & (np.max(domain_values, axis=1) >= policy.adaptive_severe_domain_threshold)
            & (score_confidence >= policy.adaptive_min_score_confidence)
            & (window_stability >= policy.adaptive_min_window_stability)
        )
        adaptive_exclude &= ~standard_exclude
    publish_exclude = standard_exclude | adaptive_exclude
    if np.any(publish_keep & publish_exclude):
        raise RuntimeError("High-confidence keep and exclude masks overlap")

    certified_call = np.full(len(out), None, dtype=object)
    certified_call[publish_keep] = "keep"
    certified_call[publish_exclude] = "exclude"
    if policy.unresolved_action == "keep":
        final_call = np.full(len(out), "keep", dtype=object)
        final_call[publish_exclude] = "exclude"
        published = np.ones(len(out), dtype=bool)
    else:
        final_call = certified_call.copy()
        published = publish_keep | publish_exclude
    review_resolution = np.full(len(out), "not_applicable", dtype=object)
    review_resolution[review_queue & adaptive_exclude] = "exclude"
    review_resolution[review_queue & ~adaptive_exclude] = "keep" if policy.unresolved_action == "keep" else "withhold"
    confidence = np.zeros(len(out), dtype=float)
    confidence[publish_keep] = np.clip(
        (policy.keep_max_score - score[publish_keep]) / max(policy.keep_max_score, EPS),
        0,
        1,
    )
    confidence[publish_exclude] = np.clip(
        (score[publish_exclude] - policy.exclude_min_score) / max(1 - policy.exclude_min_score, EPS),
        0,
        1,
    )
    if tiered_enabled and np.any(adaptive_exclude):
        primary_confidence = (
            pd.to_numeric(out.loc[adaptive_exclude, "score_confidence"], errors="coerce")
            .fillna(0.0)
            .to_numpy(dtype=float)
        )
        confidence[adaptive_exclude] = np.minimum(
            matched_tier_stability[adaptive_exclude],
            primary_confidence,
        )
    elif legacy_adaptive_enabled and np.any(adaptive_exclude):
        assert policy.adaptive_exclude_min_score is not None
        confidence[adaptive_exclude] = np.minimum(
            pd.to_numeric(out.loc[adaptive_exclude, "adaptive_window_stability"], errors="coerce")
            .fillna(0.0)
            .to_numpy(dtype=float),
            np.clip(
                (score[adaptive_exclude] - policy.adaptive_exclude_min_score)
                / max(policy.exclude_min_score - policy.adaptive_exclude_min_score, EPS),
                0,
                1,
            ),
        )
    reasons = np.full(
        len(out),
        "withheld inside calibrated uncertainty interval or detector guardrail",
        dtype=object,
    )
    if policy.unresolved_action == "keep":
        reasons[:] = (
            "operational keep: exclusion evidence did not pass the independently validated "
            "high-specificity exclusion gate"
        )
    reasons[publish_keep] = "published keep: calibrated low anomaly and detector agreement"
    reasons[standard_exclude] = "published exclude: calibrated high anomaly and detector agreement"
    if policy.unresolved_action == "keep":
        reasons[review_queue & ~adaptive_exclude] = (
            "published keep: internal review was resolved to keep because the fine-screen "
            "exclusion criteria were not all satisfied"
        )
    reasons[adaptive_exclude] = (
        "published exclude: near-boundary score resolved by severe corroborated multi-domain evidence "
        "stable across slice-window scales"
    )
    if tiered_enabled:
        for index in np.flatnonzero(adaptive_exclude):
            reasons[index] = (
                "published exclude: calibrated review tier "
                f"{matched_tier[index]} passed its multi-domain and multi-window evidence gate"
            )

    threshold_triage_reason = np.full(len(out), "score is inside the review interval", dtype=object)
    threshold_triage_reason[~finite_score] = "review: anomaly score is unavailable"
    threshold_triage_reason[threshold_band == "keep"] = f"threshold keep: score <= {policy.keep_max_score:.3f}"
    threshold_triage_reason[standard_exclude] = (
        f"threshold exclude: score >= {policy.exclude_min_score:.3f} and detector/guardrails agree"
    )
    downgraded_exclude = (threshold_band == "exclude") & ~standard_exclude
    threshold_triage_reason[downgraded_exclude] = (
        "review: high-score threshold band was downgraded because detector agreement or "
        "anatomical exclusion guardrails were not satisfied"
    )

    review_resolution_reason = np.full(len(out), "not applicable: stage-1 triage did not assign review", dtype=object)
    if tiered_enabled:
        for index in np.flatnonzero(review_queue):
            tier_name = str(matched_tier[index]) or "unmatched"
            selected = next(
                (item for item in tier_diagnostics[index] if item["tier"] == tier_name),
                None,
            )
            if adaptive_exclude[index]:
                review_resolution_reason[index] = (
                    f"exclude after review fine screen: calibrated tier {tier_name} passed "
                    "its score-band, multi-domain, multi-window, context, and protection checks"
                )
            else:
                action = "keep" if policy.unresolved_action == "keep" else "withhold"
                failures = selected["failed"] if selected is not None else ["no calibrated tier matched"]
                review_resolution_reason[index] = f"{action} after review fine screen ({tier_name}): " + "; ".join(
                    failures or ["exclusion gate not satisfied"]
                )
    elif legacy_adaptive_enabled:
        review_resolution_reason[review_queue & adaptive_exclude] = (
            "exclude after review fine screen: two-sided unprotected context, severe corroborated "
            "multi-domain evidence, and stable 3/5/7-window support"
        )
        assert policy.adaptive_exclude_min_score is not None
        for index in np.flatnonzero(review_queue & ~adaptive_exclude):
            failed: list[str] = []
            if not policy.enable_exclude:
                failed.append("exclude direction disabled")
            if not finite_score[index] or score[index] < policy.adaptive_exclude_min_score:
                failed.append("score below fine-screen floor")
            if multiscale_min_score[index] < policy.adaptive_exclude_min_score:
                failed.append("minimum multi-window score below floor")
            if not detector_candidate[index]:
                failed.append("detector did not flag a candidate")
            if not context_ok[index]:
                failed.append("two-sided context unavailable")
            if protected[index]:
                failed.append("partial/anatomical structure protection")
            if corroborating[index] < policy.adaptive_min_corroborating_domains:
                failed.append("fewer than required corroborating domains")
            if np.max(domain_values[index]) < policy.adaptive_severe_domain_threshold:
                failed.append("no severe evidence domain")
            if score_confidence[index] < policy.adaptive_min_score_confidence:
                failed.append("insufficient score confidence")
            if window_stability[index] < policy.adaptive_min_window_stability:
                failed.append("3/5/7-window support not stable")
            action = "keep" if policy.unresolved_action == "keep" else "withhold"
            review_resolution_reason[index] = f"{action} after review fine screen: " + "; ".join(
                failed or ["exclusion gate not satisfied"]
            )
    else:
        action = "keep" if policy.unresolved_action == "keep" else "withhold"
        review_resolution_reason[review_queue] = (
            f"{action} after review: no second-stage fine-screen resolver is configured"
        )

    out["internal_recommendation"] = out["recommendation"].astype(str)
    out["threshold_band"] = pd.Series(threshold_band, dtype="string")
    out["threshold_triage_call"] = pd.Series(threshold_triage_call, dtype="string")
    out["threshold_triage_reason"] = threshold_triage_reason
    out["review_resolution"] = pd.Series(review_resolution, dtype="string")
    out["review_resolution_reason"] = review_resolution_reason
    out["review_resolution_tier"] = pd.Series(matched_tier, dtype="string")
    out["review_tier_window_stability"] = matched_tier_stability
    out["review_tier_diagnostics"] = [
        json.dumps(items, ensure_ascii=False, separators=(",", ":")) for items in tier_diagnostics
    ]
    out["certified_call"] = pd.Series(certified_call, dtype="string")
    out["final_call"] = pd.Series(final_call, dtype="string")
    decision_basis = np.full(
        len(out),
        ("conservative_keep_default" if policy.unresolved_action == "keep" else "withheld"),
        dtype=object,
    )
    decision_basis[publish_keep | standard_exclude] = "independently_certified"
    if policy.unresolved_action == "keep":
        decision_basis[review_queue & ~adaptive_exclude] = "review_resolved_keep"
    decision_basis[adaptive_exclude] = (
        "tiered_multidomain_resolution" if tiered_enabled else "adaptive_multidomain_resolution"
    )
    out["decision_basis"] = decision_basis
    out["standard_exclusion_gate"] = standard_exclude
    out["adaptive_exclusion_gate"] = adaptive_exclude
    out["publication_status"] = np.where(published, "published", "withheld")
    out["binary_confidence_margin"] = confidence
    out["binary_reason"] = reasons
    out["binary_policy_id"] = policy.calibration_id
    return out


def write_high_confidence_outputs(
    metrics: pd.DataFrame,
    policy: Union[HighConfidencePolicy, Mapping[str, Any]],
    output_dir: Union[str, Path],
    *,
    application_scope: str = "certified",
) -> dict[str, Any]:
    """Write auditable binary outputs, optionally with complete operational coverage."""
    if not isinstance(policy, HighConfidencePolicy):
        policy = HighConfidencePolicy.from_mapping(policy)
    output_path = Path(output_dir).expanduser().resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    applied = apply_high_confidence_policy(metrics, policy)
    if application_scope not in {"certified", "new_input_unvalidated", "experimental_policy"}:
        raise ValueError("Unknown binary application scope")
    if application_scope in {"new_input_unvalidated", "experimental_policy"}:
        applied["policy_gate_call"] = applied["certified_call"]
        applied["policy_gate_basis"] = applied["decision_basis"]
        applied["certified_call"] = pd.Series(pd.NA, index=applied.index, dtype="string")
        applied["decision_basis"] = applied["decision_basis"].replace({"independently_certified": "frozen_policy_direct"})
        applied["binary_reason"] = applied["binary_reason"].astype(str).str.replace("certified", "frozen_policy_gate", regex=False)
    applied["application_scope"] = application_scope
    audit_path = output_path / "slice_quality_binary_audit.csv"
    calls_path = output_path / "slice_quality_binary_calls.csv"
    withheld_path = output_path / "slice_quality_withheld_audit.csv"
    summary_path = output_path / "binary_policy_application.json"
    applied.to_csv(audit_path, index=False)
    public_calls = applied.loc[applied["publication_status"].eq("published")].copy()
    public_calls["recommendation"] = public_calls["final_call"].astype(str)
    public_calls.drop(
        columns=[
            "internal_recommendation",
            "threshold_band",
            "threshold_triage_call",
            "threshold_triage_reason",
            "review_resolution",
            "review_resolution_reason",
            "review_resolution_tier",
            "review_tier_window_stability",
            "review_tier_diagnostics",
        ],
        errors="ignore",
    ).to_csv(calls_path, index=False)
    applied.loc[applied["publication_status"].eq("withheld")].to_csv(withheld_path, index=False)
    counts = applied["final_call"].value_counts(dropna=True).to_dict()
    certified = applied["decision_basis"].eq("independently_certified")
    adaptive_resolved = applied["decision_basis"].isin(
        ["adaptive_multidomain_resolution", "tiered_multidomain_resolution"]
    )
    review_resolved_keep = applied["decision_basis"].eq("review_resolved_keep")
    operational_default = applied["decision_basis"].isin(["conservative_keep_default", "review_resolved_keep"])
    summary = {
        "schema_version": SCHEMA_VERSION,
        "application_scope": application_scope,
        "new_input_independently_validated": application_scope == "certified",
        "policy_validation_scope": ("metric_stress_tests_only" if application_scope == "experimental_policy"
                                    else "historical_policy_validation"),
        "policy": _jsonable(asdict(policy)),
        "n_slices": int(len(applied)),
        "published": int(applied["publication_status"].eq("published").sum()),
        "withheld": int(applied["publication_status"].eq("withheld").sum()),
        "coverage": (float(applied["publication_status"].eq("published").mean()) if len(applied) else 0.0),
        "keep": int(counts.get("keep", 0)),
        "exclude": int(counts.get("exclude", 0)),
        "certified": int(certified.sum()),
        "standard_exclude": int(applied["standard_exclusion_gate"].sum()),
        "adaptive_exclude": int(applied["adaptive_exclusion_gate"].sum()),
        "adaptive_resolved": int(adaptive_resolved.sum()),
        "threshold_triage": {
            key: int(value) for key, value in applied["threshold_triage_call"].value_counts().to_dict().items()
        },
        "review_resolved_keep": int(review_resolved_keep.sum()),
        "review_resolved_exclude": int(adaptive_resolved.sum()),
        "operational_default_keep": int(operational_default.sum()),
        "complete_binary": bool(policy.unresolved_action == "keep"),
        "enabled_directions": {
            "keep": bool(policy.enable_keep),
            "exclude": bool(policy.enable_exclude),
        },
        "interpretation": (
            "Stage 1 assigns threshold-based keep/review/exclude triage. Stage 2 resolves only review rows "
            "with multiscale, multi-domain evidence and anatomical guardrails. Every public slice has one "
            "final keep/exclude action; decision_basis and the audit-only stage columns preserve the path."
            if policy.unresolved_action == "keep"
            else "Only independently certified keep/exclude calls are published; unresolved rows are withheld."
        ),
        "outputs": {
            "audit": audit_path.name,
            "calls": calls_path.name,
            "withheld_audit": withheld_path.name,
        },
    }
    summary_path.write_text(json.dumps(_jsonable(summary), indent=2, ensure_ascii=False), encoding="utf-8")
    return {
        "audit": str(audit_path),
        "calls": str(calls_path),
        "withheld_audit": str(withheld_path),
        "summary": str(summary_path),
        **{
            key: summary[key]
            for key in (
                "n_slices",
                "published",
                "withheld",
                "coverage",
                "keep",
                "exclude",
            )
        },
    }


def write_slice_quality_outputs(
    result: SliceSeriesResult,
    output_dir: Union[str, Path],
    *,
    title: Optional[str] = None,
    ground_truth: Optional[Mapping[str, bool]] = None,
    write_display_payload: bool = False,
) -> dict[str, str]:
    """Write core QC tables and provenance without embedding a viewer.

    The optional bounded display payload is consumed by the companion
    ``spatial-slice-quality-qc`` skill. Keeping HTML generation outside the
    package prevents UI templates from becoming a scientific runtime
    dependency.
    """
    output_path = Path(output_dir).expanduser().resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    metrics_path = output_path / "slice_quality_metrics.csv"
    calls_path = output_path / "slice_quality_calls.csv"
    manifest_path = output_path / "slice_quality_manifest.json"
    policy_path = output_path / "alignment_slice_policy.json"
    template_path = output_path / "manual_ground_truth_template.csv"
    display_payload_path = output_path / "slice_quality_display_payload.json"

    result.metrics.to_csv(metrics_path, index=False)
    call_columns = [
        "slice_index",
        "slice_id",
        "quality_anomaly_score",
        "recommendation",
        "low_quality_flag",
        "partial_structure_protection",
        "window_context",
        "reason",
    ]
    result.metrics[[column for column in call_columns if column in result.metrics]].to_csv(calls_path, index=False)
    template = result.metrics[["slice_index", "slice_id", "recommendation"]].copy()
    template["manual_label"] = "unreviewed"
    template["manual_reason"] = ""
    template.to_csv(template_path, index=False)

    detector_policy = {
        "schema_version": SCHEMA_VERSION,
        "automatic_apply": False,
        "keep_slices": result.metrics.loc[result.metrics.recommendation == "keep", "slice_id"].astype(str).tolist(),
        "review_slices": result.metrics.loc[result.metrics.recommendation == "review", "slice_id"].astype(str).tolist(),
        "exclude_candidate_slices": result.metrics.loc[result.metrics.recommendation == "exclude", "slice_id"]
        .astype(str)
        .tolist(),
        "note": (
            "Internal detector triage only. Apply an independently validated two-stage binary policy "
            "before using exclusions for alignment."
        ),
    }
    policy_path.write_text(
        json.dumps(_jsonable(detector_policy), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    resolved_title = title or result.provenance.get("config", {}).get("report_title") or "Slice QC"
    output_names = {
        "metrics": metrics_path.name,
        "calls": calls_path.name,
        "alignment_policy": policy_path.name,
        "manual_ground_truth_template": template_path.name,
    }
    if write_display_payload:
        output_names["display_payload"] = display_payload_path.name
    manifest = dict(result.provenance)
    manifest["outputs"] = output_names
    if ground_truth is not None:
        manifest["ground_truth_evaluation"] = evaluate_slice_calls(result.metrics, ground_truth)
    manifest_path.write_text(json.dumps(_jsonable(manifest), indent=2, ensure_ascii=False), encoding="utf-8")

    if write_display_payload:
        display_payload_path.write_text(
            json.dumps(
                _jsonable(
                    {
                        "schema_version": SCHEMA_VERSION,
                        "title": resolved_title,
                        "point_samples": result.point_samples,
                        "provenance": result.provenance,
                        "ground_truth": ground_truth or {},
                    }
                ),
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
    outputs = {
        "metrics": str(metrics_path),
        "calls": str(calls_path),
        "manifest": str(manifest_path),
        "alignment_policy": str(policy_path),
        "manual_ground_truth_template": str(template_path),
    }
    if write_display_payload:
        outputs["display_payload"] = str(display_payload_path)
    return outputs


def write_slice_quality_collection_outputs(
    results: Mapping[str, SliceSeriesResult],
    output_root: Union[str, Path],
    *,
    titles: Optional[Mapping[str, str]] = None,
    ground_truth: Optional[Mapping[str, Mapping[str, bool]]] = None,
    write_display_payload: bool = False,
) -> dict[str, dict[str, str]]:
    """Write one isolated core-output directory per dataset result."""
    if not isinstance(results, Mapping) or not results:
        raise ValueError("results must be a non-empty mapping")
    root = Path(output_root).expanduser().resolve()
    outputs: dict[str, dict[str, str]] = {}
    for raw_dataset_id, result in results.items():
        dataset_id = str(raw_dataset_id).strip()
        relative = Path(dataset_id)
        if not dataset_id or relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"dataset id is not a safe relative output path: {dataset_id!r}")
        outputs[dataset_id] = write_slice_quality_outputs(
            result,
            root / relative,
            title=(titles or {}).get(dataset_id),
            ground_truth=(ground_truth or {}).get(dataset_id),
            write_display_payload=write_display_payload,
        )
    return outputs
