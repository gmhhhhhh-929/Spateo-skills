#!/usr/bin/env python3
"""CLI for pre-alignment spatial slice quality control."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Any


def natural_key(value: Any) -> tuple[Any, ...]:
    return tuple(
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", str(value))
    )


def bootstrap_spateo(source: str | None) -> None:
    bundled = Path(__file__).resolve().parents[1] / 'runtime/spateo/preprocessing/slice_quality.py'
    if bundled.is_file():
        if source:
            os.environ['SPATEO_REFEREE_SOURCE'] = str(Path(source).expanduser().resolve())
        from activate_candidate import activate
        activate()
        return
    candidates: list[Path] = []
    if source:
        candidates.append(Path(source).expanduser())
    env_source = os.environ.get("SPATEO_SOURCE")
    if env_source:
        candidates.append(Path(env_source).expanduser())
    candidates.extend([Path.cwd() / "spateo-code", Path.cwd()])
    for candidate in candidates:
        if (candidate / "spateo" / "preprocessing" / "slice_quality.py").exists():
            sys.path.insert(0, str(candidate.resolve()))
            break
    try:
        from spateo.preprocessing.slice_quality import SliceQCConfig  # noqa: F401
    except Exception as exc:
        raise SystemExit(
            "Could not import spateo.preprocessing.slice_quality. "
            "Pass --spateo-source /path/to/the/directory/containing/spateo. "
            f"Original error: {exc}"
        ) from exc


def resolve_inputs(values: list[str], manifest: str | None) -> list[Path]:
    paths: list[Path] = []
    if manifest:
        import pandas as pd

        table = pd.read_csv(manifest)
        if "path" not in table:
            raise SystemExit("Manifest CSV must contain a 'path' column")
        if "order" in table:
            table = table.sort_values("order", kind="stable")
        paths.extend(Path(value).expanduser() for value in table["path"].astype(str))
    for value in values:
        path = Path(value).expanduser()
        if path.is_dir():
            paths.extend(
                sorted(path.glob("*.h5ad"), key=lambda item: natural_key(item.name))
            )
        else:
            paths.append(path)
    deduplicated: list[Path] = []
    seen: set[str] = set()
    for path in paths:
        resolved = path.resolve()
        if str(resolved) not in seen:
            deduplicated.append(resolved)
            seen.add(str(resolved))
    if not deduplicated:
        raise SystemExit("No H5AD inputs were resolved")
    return deduplicated


def parse_window(value: str) -> int | str:
    if value.lower() == "auto":
        return "auto"
    parsed = int(value)
    if parsed < 3 or parsed % 2 == 0:
        raise argparse.ArgumentTypeError("window must be an odd integer >= 3 or 'auto'")
    return parsed


def make_config(args: argparse.Namespace):
    from spateo.preprocessing.slice_quality import SliceQCConfig

    return SliceQCConfig(
        window=args.window,
        window_candidates=tuple(args.window_candidates),
        k_neighbors=args.k_neighbors,
        max_points_per_slice_report=args.report_points_per_slice,
        random_seed=args.seed,
        review_threshold=args.review_threshold,
        exclude_threshold=args.exclude_threshold,
        severe_domain_threshold=args.severe_domain_threshold,
        minimum_corrob_domains=args.minimum_corroborating_domains,
        profile_genes=args.profile_genes,
        report_title=args.title,
    )


def add_detection_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--slice-key", default="auto")
    parser.add_argument("--spatial-key", default="auto")
    parser.add_argument("--x-key")
    parser.add_argument("--y-key")
    parser.add_argument("--order-key")
    parser.add_argument("--layer", default="auto")
    parser.add_argument("--celltype-key", default="auto")
    parser.add_argument("--window", type=parse_window, default=3)
    parser.add_argument("--window-candidates", type=int, nargs="+", default=[3, 5, 7])
    parser.add_argument("--k-neighbors", type=int, default=12)
    parser.add_argument("--report-points-per-slice", type=int, default=1200)
    parser.add_argument("--profile-genes", type=int, default=3000)
    parser.add_argument("--review-threshold", type=float, default=0.38)
    parser.add_argument("--exclude-threshold", type=float, default=0.64)
    parser.add_argument("--severe-domain-threshold", type=float, default=0.78)
    parser.add_argument("--minimum-corroborating-domains", type=int, default=2)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument(
        "--title", default="Pre-alignment spatial slice quality control"
    )


def load_ground_truth(path: str | None) -> dict[str, bool] | None:
    if not path:
        return None
    import pandas as pd

    table = pd.read_csv(path)
    if "slice_id" not in table:
        raise SystemExit("Ground-truth CSV must contain slice_id")
    label_key = next(
        (
            key
            for key in ("low_quality", "is_low_quality", "label", "manual_label")
            if key in table
        ),
        None,
    )
    if label_key is None:
        raise SystemExit(
            "Ground-truth CSV must contain low_quality/is_low_quality/label/manual_label"
        )
    truthy = {"1", "true", "yes", "low_quality", "exclude", "bad"}
    return {
        str(row.slice_id): str(getattr(row, label_key)).strip().lower() in truthy
        for row in table.itertuples()
    }


def run_scan(args: argparse.Namespace) -> None:
    from spateo.preprocessing.slice_quality import (
        scan_h5ad_series,
        write_slice_quality_outputs,
    )

    paths = resolve_inputs(args.input, args.manifest)
    result = scan_h5ad_series(
        paths,
        config=make_config(args),
        sort_inputs=args.manifest is None,
        slice_key=args.slice_key,
        spatial_key=args.spatial_key,
        x_key=args.x_key,
        y_key=args.y_key,
        order_key=args.order_key,
        layer=args.layer,
        celltype_key=args.celltype_key,
    )
    ground_truth = load_ground_truth(args.ground_truth)
    outputs = write_slice_quality_outputs(
        result,
        args.output_dir,
        title=args.title,
        ground_truth=ground_truth,
        write_display_payload=args.write_display_payload or not args.no_report,
    )
    if not args.no_report:
        from slice_quality_report import render_slice_quality_report_from_outputs

        outputs["report"] = render_slice_quality_report_from_outputs(
            args.output_dir,
            title=args.title,
        )
    summary = result.metrics["recommendation"].value_counts().to_dict()
    print(
        json.dumps(
            {
                "inputs": [str(path) for path in paths],
                "calls": summary,
                "outputs": outputs,
            },
            indent=2,
        )
    )


def _load_dataset_manifest(
    path: str,
) -> tuple[dict[str, list[Path]], dict[str, dict[str, Any]]]:
    import pandas as pd

    table = pd.read_csv(Path(path).expanduser().resolve())
    required = {"dataset_id", "path"}
    missing = required.difference(table.columns)
    if missing:
        raise SystemExit(f"Dataset manifest is missing columns: {sorted(missing)}")
    table = table.copy()
    table["dataset_id"] = table["dataset_id"].astype(str).str.strip()
    if table["dataset_id"].eq("").any():
        raise SystemExit("Dataset manifest contains an empty dataset_id")
    table["__row_order"] = range(len(table))
    datasets: dict[str, list[Path]] = {}
    options: dict[str, dict[str, Any]] = {}
    option_columns = (
        "slice_key",
        "spatial_key",
        "x_key",
        "y_key",
        "order_key",
        "layer",
        "celltype_key",
    )
    for dataset_id, group in table.groupby("dataset_id", sort=False):
        ordered = (
            group.sort_values("order", kind="stable")
            if "order" in group
            else group.sort_values("__row_order")
        )
        datasets[dataset_id] = [
            Path(value).expanduser().resolve() for value in ordered["path"].astype(str)
        ]
        dataset_options: dict[str, Any] = {}
        for column in option_columns:
            if column not in group:
                continue
            values = [
                str(value).strip()
                for value in group[column].dropna().unique()
                if str(value).strip()
            ]
            if len(values) > 1:
                raise SystemExit(
                    f"Dataset {dataset_id!r} has multiple values for {column}: {values}. "
                    "Use one setting per dataset."
                )
            if values:
                dataset_options[column] = values[0]
        options[dataset_id] = dataset_options
    return datasets, options


def run_scan_batch(args: argparse.Namespace) -> None:
    from spateo.preprocessing.slice_quality import (
        scan_h5ad_collection,
        write_slice_quality_collection_outputs,
    )

    datasets, manifest_options = _load_dataset_manifest(args.dataset_manifest)
    defaults = {
        "slice_key": args.slice_key,
        "spatial_key": args.spatial_key,
        "x_key": args.x_key,
        "y_key": args.y_key,
        "order_key": args.order_key,
        "layer": args.layer,
        "celltype_key": args.celltype_key,
        "sort_inputs": False,
    }
    dataset_options = {
        dataset_id: {**defaults, **manifest_options.get(dataset_id, {})}
        for dataset_id in datasets
    }
    results = scan_h5ad_collection(
        datasets,
        config=make_config(args),
        dataset_options=dataset_options,
    )
    titles = {dataset_id: f"{args.title} · {dataset_id}" for dataset_id in results}
    outputs = write_slice_quality_collection_outputs(
        results,
        args.output_dir,
        titles=titles,
        write_display_payload=args.write_display_payload or not args.no_report,
    )
    if not args.no_report:
        from slice_quality_report import render_slice_quality_report_from_outputs

        for dataset_id, dataset_outputs in outputs.items():
            dataset_outputs["report"] = render_slice_quality_report_from_outputs(
                Path(args.output_dir).expanduser().resolve() / dataset_id,
                title=titles[dataset_id],
            )
    calls = {
        dataset_id: result.metrics["recommendation"].value_counts().to_dict()
        for dataset_id, result in results.items()
    }
    print(
        json.dumps(
            {"datasets": list(results), "calls": calls, "outputs": outputs}, indent=2
        )
    )


def load_plan(value: str) -> list[dict[str, Any]]:
    path = Path(value).expanduser()
    obj = (
        json.loads(path.read_text(encoding="utf-8"))
        if path.exists()
        else json.loads(value)
    )
    if isinstance(obj, dict) and "artifacts" in obj:
        obj = obj["artifacts"]
    if not isinstance(obj, list):
        raise SystemExit(
            "Simulation plan must be a JSON list or an object with an 'artifacts' list"
        )
    return [dict(item) for item in obj]


def run_simulate(args: argparse.Namespace) -> None:
    import anndata as ad
    from spateo.preprocessing.slice_quality import simulate_slice_quality_artifacts

    source = Path(args.input).expanduser().resolve()
    adata = ad.read_h5ad(source)
    simulated = simulate_slice_quality_artifacts(
        adata,
        load_plan(args.plan),
        slice_key=args.slice_key,
        spatial_key=args.spatial_key,
        layer=args.layer,
        output_layer=args.output_layer,
        random_seed=args.seed,
    )
    output = Path(args.output_h5ad).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    simulated.write_h5ad(output, compression="gzip")
    print(
        json.dumps(
            {
                "output_h5ad": str(output),
                "simulation": simulated.uns["slice_quality_simulation"],
            },
            indent=2,
        )
    )


def _auto_benchmark_plan(labels: list[str]) -> list[dict[str, Any]]:
    if len(labels) < 9:
        raise SystemExit(
            "Automatic benchmark planning requires at least 9 slices; provide --plan for a shorter series"
        )
    positions = [
        round((len(labels) - 1) * fraction) for fraction in (0.22, 0.42, 0.62, 0.80)
    ]
    positions = [min(max(value, 2), len(labels) - 3) for value in positions]
    # Preserve order while avoiding collisions from short series.
    unique_positions: list[int] = []
    for value in positions:
        while value in unique_positions and value < len(labels) - 3:
            value += 1
        unique_positions.append(value)
    return [
        {"slice_id": labels[unique_positions[0]], "kind": "depth", "rate": 0.18},
        {
            "slice_id": labels[unique_positions[1]],
            "kind": "regional_depth",
            "rate": 0.08,
            "radius_quantile": 0.42,
        },
        {
            "slice_id": labels[unique_positions[2]],
            "kind": "cell_dropout",
            "keep_rate": 0.28,
        },
        {
            "slice_id": labels[unique_positions[3]],
            "kind": "spatial_hole",
            "radius_quantile": 0.38,
        },
    ]


def _prepare_benchmark_slices(
    adata,
    slice_key: str,
    spatial_key: str,
    x_key: str | None,
    y_key: str | None,
) -> tuple[str, str, list[str]]:
    """Resolve an explicit simulation slice key, including discrete-z inputs."""
    import numpy as np
    import pandas as pd

    def resolve_simulation_spatial() -> str:
        if spatial_key != "auto" and spatial_key in adata.obsm:
            return spatial_key
        inferred = next(
            (
                key
                for key in ("raw_spatial", "spatial", "spatial_3d", "align_spatial")
                if key in adata.obsm
            ),
            None,
        )
        if inferred is not None:
            return inferred
        if x_key and y_key and x_key in adata.obs and y_key in adata.obs:
            derived_spatial = "__slice_qc_benchmark_spatial"
            adata.obsm[derived_spatial] = np.column_stack(
                [
                    pd.to_numeric(adata.obs[x_key], errors="coerce").to_numpy(
                        dtype=float
                    ),
                    pd.to_numeric(adata.obs[y_key], errors="coerce").to_numpy(
                        dtype=float
                    ),
                ]
            )
            return derived_spatial
        raise SystemExit(
            "Benchmark simulation requires an obsm spatial key or explicit --x-key/--y-key coordinates"
        )

    if slice_key != "auto":
        if slice_key not in adata.obs:
            raise SystemExit(f"slice key {slice_key!r} is missing")
        resolved_spatial = resolve_simulation_spatial()
        labels = sorted(pd.unique(adata.obs[slice_key].astype(str)), key=natural_key)
        return slice_key, resolved_spatial, labels

    obs_key = next(
        (
            key
            for key in (
                "slice_id",
                "slices",
                "slice",
                "sample_order",
                "section_id",
                "section",
                "library_id",
            )
            if key in adata.obs
        ),
        None,
    )
    if obs_key is not None:
        resolved_spatial = resolve_simulation_spatial()
        labels = sorted(pd.unique(adata.obs[obs_key].astype(str)), key=natural_key)
        return obs_key, resolved_spatial, labels

    resolved_spatial = spatial_key
    if resolved_spatial == "auto":
        resolved_spatial = next(
            (
                key
                for key in ("spatial_3d", "spatial", "align_spatial", "tdr_spatial")
                if key in adata.obsm
            ),
            "auto",
        )
    if resolved_spatial == "auto" or resolved_spatial not in adata.obsm:
        raise SystemExit(
            "benchmark --slice-key auto requires an inferable 3D spatial coordinate key"
        )
    coords = np.asarray(adata.obsm[resolved_spatial])
    if coords.ndim != 2 or coords.shape[1] < 3:
        raise SystemExit(
            "benchmark --slice-key auto requires discrete z values in the third spatial column"
        )
    z_values = np.asarray(coords[:, 2], dtype=float)
    unique = np.unique(z_values[np.isfinite(z_values)])
    if unique.size < 2:
        raise SystemExit(
            "Could not resolve multiple discrete z slices for benchmark simulation"
        )
    mapping = {value: f"z{value:g}" for value in sorted(unique)}
    derived_key = "__slice_qc_benchmark_slice_id"
    adata.obs[derived_key] = [mapping.get(value, "zNA") for value in z_values]
    labels = [mapping[value] for value in sorted(unique)]
    return derived_key, resolved_spatial, labels


def run_benchmark(args: argparse.Namespace) -> None:
    from dataclasses import replace

    import anndata as ad
    import pandas as pd
    from spateo.preprocessing.slice_quality import (
        evaluate_paired_simulation,
        evaluate_slice_calls,
        scan_h5ad_series,
        simulate_slice_quality_artifacts,
        write_slice_quality_outputs,
    )

    source = Path(args.input).expanduser().resolve()
    adata = ad.read_h5ad(source)
    simulation_slice_key, simulation_spatial_key, labels = _prepare_benchmark_slices(
        adata, args.slice_key, args.spatial_key, args.x_key, args.y_key
    )
    plan = load_plan(args.plan) if args.plan else _auto_benchmark_plan(labels)
    benchmark_config = make_config(args)
    baseline_result = scan_h5ad_series(
        [source],
        config=benchmark_config,
        slice_key=args.slice_key,
        spatial_key=args.spatial_key,
        x_key=args.x_key,
        y_key=args.y_key,
        order_key=args.order_key,
        layer=args.layer,
        celltype_key=args.celltype_key,
    )
    selected_window = int(baseline_result.provenance["selected_window"])
    simulated = simulate_slice_quality_artifacts(
        adata,
        plan,
        slice_key=simulation_slice_key,
        spatial_key=simulation_spatial_key,
        layer=args.layer,
        output_layer=args.output_layer,
        random_seed=args.seed,
    )
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    synthetic_path = output_dir / "synthetic_low_quality.h5ad"
    simulated.write_h5ad(synthetic_path, compression="gzip")
    ground_truth = {str(value): False for value in labels}
    for item in plan:
        ground_truth[str(item["slice_id"])] = True
    gt_path = output_dir / "synthetic_ground_truth.csv"
    pd.DataFrame(
        [
            {
                "slice_id": key,
                "low_quality": value,
                "artifact": next(
                    (x["kind"] for x in plan if str(x["slice_id"]) == key), "none"
                ),
            }
            for key, value in ground_truth.items()
        ]
    ).to_csv(gt_path, index=False)
    (output_dir / "simulation_plan.json").write_text(
        json.dumps({"artifacts": plan}, indent=2), encoding="utf-8"
    )
    baseline_outputs = write_slice_quality_outputs(
        baseline_result,
        output_dir / "baseline_reference",
        title=f"{args.title} — unmodified baseline reference",
        write_display_payload=args.write_display_payload or not args.no_report,
    )
    if not args.no_report:
        from slice_quality_report import render_slice_quality_report_from_outputs

        baseline_outputs["report"] = render_slice_quality_report_from_outputs(
            output_dir / "baseline_reference",
            title=f"{args.title} — unmodified baseline reference",
        )
    result = scan_h5ad_series(
        [synthetic_path],
        config=replace(benchmark_config, window=selected_window),
        slice_key=simulation_slice_key,
        spatial_key=simulation_spatial_key,
        x_key=args.x_key,
        y_key=args.y_key,
        order_key=args.order_key,
        layer=args.output_layer,
        celltype_key=args.celltype_key,
        profile_gene_names=baseline_result.profile_genes,
    )
    outputs = write_slice_quality_outputs(
        result,
        output_dir,
        title=args.title,
        ground_truth=ground_truth,
        write_display_payload=args.write_display_payload or not args.no_report,
    )
    if not args.no_report:
        from slice_quality_report import render_slice_quality_report_from_outputs

        outputs["report"] = render_slice_quality_report_from_outputs(
            output_dir, title=args.title
        )
    evaluation = evaluate_slice_calls(result.metrics, ground_truth)
    paired_evaluation, paired_table = evaluate_paired_simulation(
        baseline_result.metrics, result.metrics, plan
    )
    paired_path = output_dir / "paired_detection_comparison.csv"
    paired_table.to_csv(paired_path, index=False)
    evaluation["simulation_plan"] = plan
    evaluation["synthetic_h5ad"] = str(synthetic_path)
    evaluation["paired_evaluation"] = paired_evaluation
    evaluation["paired_comparison_csv"] = str(paired_path)
    evaluation["baseline_outputs"] = baseline_outputs
    evaluation["selected_window_locked_from_baseline"] = selected_window
    evaluation["naive_confusion_matrix_warning"] = (
        "The tp/fp/fn/tn block assumes every non-injected slice is high quality. "
        "Use paired_evaluation when the real baseline may contain native QC candidates."
    )
    evaluation_path = output_dir / "synthetic_benchmark_evaluation.json"
    evaluation_path.write_text(json.dumps(evaluation, indent=2), encoding="utf-8")
    print(json.dumps({"evaluation": evaluation, "outputs": outputs}, indent=2))


def run_summarize(args: argparse.Namespace) -> None:
    from slice_quality_report import write_slice_quality_collection_report

    catalog_records = None
    if args.catalog:
        catalog_payload = json.loads(
            Path(args.catalog).expanduser().resolve().read_text(encoding="utf-8")
        )
        catalog_records = catalog_payload.get("records", catalog_payload)
    summary = write_slice_quality_collection_report(
        [Path(value).expanduser().resolve() for value in args.run_dir],
        args.output_html,
        title=args.title,
        binary_only=args.binary_only,
        catalog_records=catalog_records,
    )
    print(json.dumps(summary, indent=2))


def run_render(args: argparse.Namespace) -> None:
    from slice_quality_report import (
        render_binary_slice_quality_report_from_outputs,
        render_slice_quality_report_from_outputs,
    )

    renderer = (
        render_binary_slice_quality_report_from_outputs
        if args.binary_only
        else render_slice_quality_report_from_outputs
    )
    report = renderer(
        args.input_dir,
        args.output_html,
        title=args.title,
    )
    print(json.dumps({"report": report}, indent=2))


def run_render_detail(args: argparse.Namespace) -> None:
    from slice_quality_visualization import render_binary_slice_quality_appendix

    outputs = render_binary_slice_quality_appendix(
        args.input_dir,
        args.output_dir,
        title=args.title,
        display_window=args.display_window,
        max_points_per_slice=args.max_points_per_slice,
        language=args.language,
        source_h5ad=args.source_h5ad,
        slice_key=args.slice_key,
        spatial_key=args.spatial_key,
        total_counts_key=args.total_counts_key,
        n_genes_key=args.n_genes_key,
        mito_key=args.mito_key,
    )
    print(json.dumps(outputs, ensure_ascii=False, indent=2))


def run_summarize_paper(args: argparse.Namespace) -> None:
    from slice_quality_visualization import write_collection_paper_tables

    outputs = write_collection_paper_tables(
        args.dataset_summary_csv,
        args.output_dir,
        title=args.title,
    )
    print(json.dumps(outputs, ensure_ascii=False, indent=2))


def run_publish(args: argparse.Namespace) -> None:
    from dataclasses import fields

    import pandas as pd
    from spateo.preprocessing.slice_quality import (
        HighConfidencePolicy,
        SliceQCConfig,
        add_multiscale_exclusion_evidence,
        write_high_confidence_outputs,
    )

    input_dir = Path(args.input_dir).expanduser().resolve()
    metrics_path = input_dir / "slice_quality_metrics.csv"
    if not metrics_path.exists():
        raise SystemExit(f"Missing baseline metrics: {metrics_path}")
    policy_path = Path(args.policy).expanduser().resolve()
    policy_payload = json.loads(policy_path.read_text(encoding="utf-8"))
    if "policy" in policy_payload:
        policy_payload = policy_payload["policy"]
    benchmark = policy_payload.get("benchmark_summary", {})
    base_validation_passed = bool(
        benchmark.get(
            "base_policy_validation_passed", benchmark.get("validation_passed", False)
        )
    )
    resolver_validation_passed = bool(
        benchmark.get("adaptive_extension_validation_passed", False)
    )
    application_scope = getattr(args, "application_scope", "certified")
    experimental_metric_application = application_scope == "experimental_policy"
    if experimental_metric_application and not benchmark.get("metric_stress_validation_passed", False):
        raise SystemExit("Experimental-policy application requires a passed, frozen metric-stress evaluation; it does not establish biological certification.")
    if not resolver_validation_passed and application_scope == "new_input_unvalidated":
        raise SystemExit("new_input_unvalidated transfers a historically validated policy. Use experimental_policy for a new metric-tested resolver.")
    if not base_validation_passed and not args.allow_unvalidated_policy:
        raise SystemExit(
            "Refusing to publish an unvalidated binary policy. "
            "The policy benchmark_summary.base_policy_validation_passed (or legacy validation_passed) "
            "field must be true; "
            "use --allow-unvalidated-policy only for an explicitly labelled development audit."
        )
    certifications = policy_payload.get("benchmark_summary", {}).get(
        "dataset_direction_certification", {}
    )
    application_scope = getattr(args, "application_scope", "certified")
    if certifications and application_scope == "certified":
        dataset_id = f"{input_dir.parent.name}/{input_dir.name}"
        certification = certifications.get(dataset_id)
        if certification is None:
            raise SystemExit(
                f"Refusing to publish dataset absent from the direction-certification map: {dataset_id}"
            )
        policy_payload["enable_keep"] = bool(certification.get("keep_certified", False))
        policy_payload["enable_exclude"] = bool(
            certification.get("exclude_certified", False)
        )
    if args.complete_binary:
        resolver_configured = bool(policy_payload.get("review_exclusion_tiers")) or (
            policy_payload.get("adaptive_exclude_min_score") is not None
        )
        if not resolver_configured:
            raise SystemExit(
                "Complete two-stage publication requires calibrated review_exclusion_tiers "
                "(or the legacy adaptive_exclude_min_score) and their review-resolution settings. "
                "Refusing to relabel review rows without a configured fine screen."
            )
        if not resolver_validation_passed and not args.allow_unvalidated_policy and not experimental_metric_application:
            raise SystemExit(
                "Refusing complete binary publication because the frozen review resolver has not "
                "passed independent validation. Expected "
                "benchmark_summary.adaptive_extension_validation_passed=true."
            )
        policy_payload["unresolved_action"] = "keep"
    policy = HighConfidencePolicy.from_mapping(policy_payload)
    metrics = pd.read_csv(metrics_path, dtype={"slice_id": str})
    tiered = bool(policy.review_exclusion_tiers)
    adaptive_required = (
        {"adaptive_window_details"}
        if tiered
        else {"adaptive_window_stability", "adaptive_min_score_across_windows"}
    )
    if (
        tiered or policy.adaptive_exclude_min_score is not None
    ) and not adaptive_required.issubset(metrics.columns):
        manifest_path = input_dir / "slice_quality_manifest.json"
        if not manifest_path.exists():
            raise SystemExit(
                "Adaptive review resolution needs slice_quality_manifest.json to reconstruct "
                "the detector configuration and compute 3/5/7-window evidence."
            )
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        saved_config = dict(manifest.get("config", {}))
        allowed = {item.name for item in fields(SliceQCConfig)}
        config_values = {
            key: value for key, value in saved_config.items() if key in allowed
        }
        for key in ("window_candidates", "mito_prefixes"):
            if key in config_values:
                config_values[key] = tuple(config_values[key])
        config = SliceQCConfig(**config_values)
        adaptive_summary = policy.benchmark_summary.get("adaptive_extension", {})
        windows = tuple(
            int(value) for value in adaptive_summary.get("windows", (3, 5, 7))
        )
        minimum_domains = (
            min(
                (
                    tier.min_corroborating_domains
                    for tier in policy.review_exclusion_tiers
                    if tier.enable_exclude
                ),
                default=policy.adaptive_min_corroborating_domains,
            )
            if tiered
            else policy.adaptive_min_corroborating_domains
        )
        severe_threshold = (
            min(
                (
                    tier.severe_domain_threshold
                    for tier in policy.review_exclusion_tiers
                    if tier.enable_exclude
                ),
                default=policy.adaptive_severe_domain_threshold,
            )
            if tiered
            else policy.adaptive_severe_domain_threshold
        )
        metrics = add_multiscale_exclusion_evidence(
            metrics,
            config=config,
            windows=windows,
            minimum_corroborating_domains=minimum_domains,
            severe_domain_threshold=severe_threshold,
        )
    outputs = write_high_confidence_outputs(
        metrics,
        policy,
        Path(args.output_dir).expanduser().resolve() if args.output_dir else input_dir,
        application_scope=application_scope,
    )
    import hashlib
    summary_path = Path(outputs["summary"])
    summary = json.loads(summary_path.read_text())
    summary.update(policy_source=str(policy_path), policy_sha256=hashlib.sha256(policy_path.read_bytes()).hexdigest(), input_dataset=str(input_dir))
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(outputs, indent=2))


def run_preregister(args):
    from slice_quality_preregistration import preregister_outputs
    from spateo.preprocessing.slice_preregistration import PreregistrationConfig
    result = preregister_outputs(args.input_dir, args.dataset_id,
        PreregistrationConfig(sample_size=args.fit_points, seed=args.seed,annotation_weight=args.annotation_weight), args.display_points)
    print(json.dumps({"frame_id":result["frame_id"],"reference":result["reference_slice"],"seconds":result["total_elapsed_seconds"]}))


def run_export_roi(args):
    from slice_quality_preregistration import export_roi
    print(json.dumps(export_roi(args.input_dir,args.roi,args.output_csv)))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Detect low-quality serial spatial transcriptomics slices before alignment"
    )
    parser.add_argument(
        "--spateo-source", help="Directory containing the spateo Python package"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    prereg = subparsers.add_parser("preregister", help="Display-only rigid preregistration of an existing QC run")
    prereg.add_argument("--input-dir", required=True)
    prereg.add_argument("--dataset-id", required=True)
    prereg.add_argument("--fit-points", type=int, default=1000)
    prereg.add_argument("--display-points", type=int, default=1800)
    prereg.add_argument("--seed", type=int, default=13)
    prereg.add_argument("--annotation-weight", type=float, default=0.0, help="Optional observed-cell-type candidate ranking; 0 preserves geometry-only fitting")
    prereg.set_defaults(func=run_preregister)
    roi = subparsers.add_parser("export-roi", help="Exact full-data replay of a display-frame ROI")
    roi.add_argument("--input-dir", required=True)
    roi.add_argument("--roi", required=True)
    roi.add_argument("--output-csv", required=True)
    roi.set_defaults(func=run_export_roi)


    scan = subparsers.add_parser(
        "scan", help="Scan one multi-slice H5AD or multiple per-slice H5AD files"
    )
    scan.add_argument("--input", nargs="+", default=[])
    scan.add_argument("--manifest", help="CSV with path and optional order columns")
    scan.add_argument("--output-dir", required=True)
    scan.add_argument("--ground-truth")
    scan.add_argument(
        "--no-report",
        action="store_true",
        help="Do not generate HTML on the compute host",
    )
    scan.add_argument(
        "--write-display-payload",
        action="store_true",
        help="Write a compact point-sample payload for later local HTML rendering",
    )
    add_detection_options(scan)
    scan.set_defaults(func=run_scan)

    scan_batch = subparsers.add_parser(
        "scan-batch",
        help="Scan multiple independent datasets from one manifest without crossing dataset boundaries",
    )
    scan_batch.add_argument(
        "--dataset-manifest",
        required=True,
        help=(
            "CSV with dataset_id and path; optional order and per-dataset slice_key, spatial_key, "
            "x_key, y_key, order_key, layer, and celltype_key columns"
        ),
    )
    scan_batch.add_argument("--output-dir", required=True)
    scan_batch.add_argument(
        "--no-report", action="store_true", help="Do not generate per-dataset HTML"
    )
    scan_batch.add_argument(
        "--write-display-payload",
        action="store_true",
        help="Write compact point samples for later local HTML rendering",
    )
    add_detection_options(scan_batch)
    scan_batch.set_defaults(func=run_scan_batch)

    simulate = subparsers.add_parser(
        "simulate", help="Create a non-destructive synthetic low-quality H5AD"
    )
    simulate.add_argument("--input", required=True)
    simulate.add_argument("--output-h5ad", required=True)
    simulate.add_argument("--plan", required=True, help="JSON file or inline JSON list")
    simulate.add_argument("--slice-key", required=True)
    simulate.add_argument("--spatial-key", default="spatial")
    simulate.add_argument("--layer", default="auto")
    simulate.add_argument("--output-layer", default="slice_qc_simulated_counts")
    simulate.add_argument("--seed", type=int, default=13)
    simulate.set_defaults(func=run_simulate)

    benchmark = subparsers.add_parser(
        "benchmark", help="Inject known defects and calculate slice-level accuracy"
    )
    benchmark.add_argument("--input", required=True)
    benchmark.add_argument("--output-dir", required=True)
    benchmark.add_argument(
        "--plan", help="Optional JSON plan; otherwise four internal slices are selected"
    )
    benchmark.add_argument("--output-layer", default="slice_qc_simulated_counts")
    benchmark.add_argument(
        "--no-report",
        action="store_true",
        help="Do not generate HTML on the compute host",
    )
    benchmark.add_argument(
        "--write-display-payload",
        action="store_true",
        help="Write compact point samples for later local HTML rendering",
    )
    add_detection_options(benchmark)
    benchmark.set_defaults(func=run_benchmark)

    summarize = subparsers.add_parser(
        "summarize", help="Bundle completed QC runs into one interactive HTML"
    )
    summarize.add_argument("--run-dir", nargs="+", required=True)
    summarize.add_argument("--output-html", required=True)
    summarize.add_argument("--title", default="Spateo referee")
    summarize.add_argument(
        "--catalog", help="Optional grouped official-file catalogue JSON"
    )
    summarize.add_argument(
        "--binary-only",
        action="store_true",
        help="Build the public collection from validated keep/exclude publication tables only",
    )
    summarize.set_defaults(func=run_summarize)

    render = subparsers.add_parser(
        "render", help="Render HTML from synced metrics and display payload"
    )
    render.add_argument("--input-dir", required=True)
    render.add_argument("--output-html")
    render.add_argument("--title")
    render.add_argument(
        "--binary-only",
        action="store_true",
        help="Render only validated published keep/exclude rows from the binary audit",
    )
    render.set_defaults(func=run_render)

    render_detail = subparsers.add_parser(
        "render-detail",
        help="Build a detailed keep/exclude appendix with KDE density windows and paper tables",
    )
    render_detail.add_argument("--input-dir", required=True)
    render_detail.add_argument("--output-dir", required=True)
    render_detail.add_argument("--title", default="Spateo referee · serial-slice QC")
    render_detail.add_argument("--display-window", type=int, choices=[3, 5], default=5)
    render_detail.add_argument("--max-points-per-slice", type=int, default=900)
    render_detail.add_argument(
        "--language",
        choices=["zh", "en"],
        default="zh",
        help="Report language: zh for Chinese or en for an English submission appendix",
    )
    render_detail.add_argument(
        "--source-h5ad",
        help="Optional source H5AD for absolute expression and connected-component evidence maps",
    )
    render_detail.add_argument(
        "--slice-key", help="Source-H5AD obs slice key; inferred from the run manifest"
    )
    render_detail.add_argument(
        "--spatial-key", help="Source-H5AD obsm key; inferred from the run manifest"
    )
    render_detail.add_argument(
        "--total-counts-key", help="Optional obs column containing captured counts"
    )
    render_detail.add_argument(
        "--n-genes-key", help="Optional obs column containing detected genes"
    )
    render_detail.add_argument(
        "--mito-key", help="Optional obs column containing mitochondrial percentage"
    )
    render_detail.set_defaults(func=run_render_detail)

    summarize_paper = subparsers.add_parser(
        "summarize-paper",
        help="Condense a keep/exclude-only collection into main and supplementary paper tables",
    )
    summarize_paper.add_argument("--dataset-summary-csv", required=True)
    summarize_paper.add_argument("--output-dir", required=True)
    summarize_paper.add_argument("--title", default="Spateo referee dataset summary")
    summarize_paper.set_defaults(func=run_summarize_paper)

    publish = subparsers.add_parser(
        "publish",
        help="Publish calibrated keep/exclude calls, optionally with a complete operational action for every slice",
    )
    publish.add_argument(
        "--input-dir", required=True, help="Completed baseline QC run directory"
    )
    publish.add_argument(
        "--policy", required=True, help="Calibrated binary policy JSON"
    )
    publish.add_argument("--output-dir", help="Defaults to the baseline run directory")
    publish.add_argument(
        "--complete-binary",
        action="store_true",
        help=(
            "Run the required review-only fine screen, resolve every slice to keep/exclude, and retain "
            "the internal decision path in the audit"
        ),
    )
    publish.add_argument(
        "--allow-unvalidated-policy",
        action="store_true",
        help="Development-only override; never use for a public keep/exclude result",
    )
    publish.add_argument("--application-scope", choices=["certified", "new_input_unvalidated", "experimental_policy"], default="certified", help="Choose certified publication, transfer of a historical policy, or an explicitly experimental metric-tested policy")
    publish.set_defaults(func=run_publish)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    bootstrap_spateo(args.spateo_source)
    args.func(args)


if __name__ == "__main__":
    main()
