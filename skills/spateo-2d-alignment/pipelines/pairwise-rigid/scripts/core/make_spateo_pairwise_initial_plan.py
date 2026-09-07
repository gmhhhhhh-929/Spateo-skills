#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from pairwise_core import (
    build_adjacent_edges,
    json_dump,
    now_iso,
    sanitize_token,
    select_order_rows,
    shell_quote,
    write_csv_rows,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate adjacent pairwise Spateo SN-S jobs.")
    parser.add_argument("--dataset-root", required=True)
    parser.add_argument("--runner-script", required=True)
    parser.add_argument("--slice-order-csv")
    parser.add_argument("--stage")
    parser.add_argument("--window-start")
    parser.add_argument("--window-end")
    parser.add_argument("--coordinate-key", default="spatial")
    parser.add_argument("--project-account", required=True)
    parser.add_argument("--conda-env", required=True)
    parser.add_argument("--stage1-mode", default="SN-S")
    parser.add_argument("--stage2-mode", default="none")
    parser.add_argument("--stage1-max-iter", type=int, default=300)
    parser.add_argument("--stage2-max-iter", type=int, default=300)
    parser.add_argument("--expression-mode", choices=["normal", "spatial_only"], default="normal")
    parser.add_argument("--dummy-rep-dim", type=int, default=30)
    parser.add_argument("--sigma2-init-scale", default="")
    parser.add_argument("--sigma2-end", default="")
    parser.add_argument("--n-sampling-ref", type=int, default=40000)
    parser.add_argument("--max-cells-per-slice", type=int, default=0)
    parser.add_argument("--batch-size", type=int, default=800)
    parser.add_argument("--save-full-assignment", action="store_true")
    parser.add_argument("--no-require-cuda", dest="require_cuda", action="store_false")
    parser.set_defaults(require_cuda=True)
    parser.add_argument("--rng-seed", type=int, default=20260612)
    parser.add_argument("--exclude-celltypes-regex", default="")
    parser.add_argument("--dsub-resource", default="cpu=64;mem=240000;gpu=1")
    parser.add_argument("--walltime", type=int, default=345600)
    parser.add_argument("--job-name-prefix", default="spateo_pair")
    parser.add_argument("--runtime-output-dir")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    outdir = Path(args.output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    runtime_outdir = Path(args.runtime_output_dir).expanduser() if args.runtime_output_dir else outdir

    def runtime_path(local_path: Path) -> Path:
        if not args.runtime_output_dir:
            return local_path
        return runtime_outdir / local_path.relative_to(outdir)

    dataset_root = Path(args.dataset_root).expanduser()
    order_rows = select_order_rows(
        dataset_root,
        Path(args.slice_order_csv).expanduser().resolve() if args.slice_order_csv else None,
        args.stage,
        args.window_start,
        args.window_end,
    )
    if len(order_rows) < 2:
        raise SystemExit("Need at least two ordered slices to generate pairwise jobs.")
    edges = build_adjacent_edges(order_rows)
    runner_script = Path(args.runner_script).expanduser()
    rows = []
    for edge in edges:
        edge_local = outdir / "edges" / edge["edge_id"]
        edge_runtime = runtime_path(edge_local)
        payload_local = outdir / "payloads" / f"{sanitize_token(edge['edge_id'])}.sh"
        payload_runtime = runtime_path(payload_local)
        edge_local.mkdir(parents=True, exist_ok=True)
        (edge_local / "logs").mkdir(parents=True, exist_ok=True)
        job_name = (args.job_name_prefix + "_" + sanitize_token(edge["edge_id"]))[:120]
        provenance = {
            "version": "spateo_pairwise_initial_edge_v1",
            "created_at": now_iso(),
            "status": "planned",
            **edge,
            "dataset_root": str(dataset_root),
            "coordinate_key": args.coordinate_key,
            "runner_script": str(runner_script),
            "runner_parameters": {
                "stage1_mode": args.stage1_mode,
                "stage2_mode": args.stage2_mode,
                "stage1_max_iter": args.stage1_max_iter,
                "stage2_max_iter": args.stage2_max_iter,
                "expression_mode": args.expression_mode,
                "dummy_rep_dim": args.dummy_rep_dim,
                "sigma2_init_scale": args.sigma2_init_scale,
                "sigma2_end": args.sigma2_end,
                "n_sampling_ref": args.n_sampling_ref,
                "max_cells_per_slice": args.max_cells_per_slice,
                "batch_size": args.batch_size,
                "save_full_assignment": args.save_full_assignment,
                "require_cuda": args.require_cuda,
                "rng_seed": args.rng_seed,
                "exclude_celltypes_regex": args.exclude_celltypes_regex,
            },
            "dsub": {
                "project_account": args.project_account,
                "resource": args.dsub_resource,
                "walltime": args.walltime,
                "job_name": job_name,
                "job_id": "",
            },
            "local_edge_dir": str(edge_local),
            "runtime_edge_dir": str(edge_runtime),
        }
        json_dump(edge_local / "provenance.json", provenance)
        env_lines = [
            f"export DATASET_ROOT={shell_quote(dataset_root)}",
            f"export OUTDIR={shell_quote(edge_runtime)}",
            f"export PAIR_START_ORDER={int(edge['pair_start_order'])}",
            "export N_SLICES=2",
            f"export STAGE1_MODE={shell_quote(args.stage1_mode)}",
            f"export STAGE2_MODE={shell_quote(args.stage2_mode)}",
            f"export STAGE1_MAX_ITER={int(args.stage1_max_iter)}",
            f"export STAGE2_MAX_ITER={int(args.stage2_max_iter)}",
            f"export EXPRESSION_MODE={shell_quote(args.expression_mode)}",
            f"export DUMMY_REP_DIM={int(args.dummy_rep_dim)}",
            f"export N_SAMPLING_REF={int(args.n_sampling_ref)}",
            f"export MAX_CELLS_PER_SLICE={int(args.max_cells_per_slice)}",
            f"export BATCH_SIZE={int(args.batch_size)}",
            f"export SAVE_FULL_ASSIGNMENT={1 if args.save_full_assignment else 0}",
            f"export REQUIRE_CUDA={1 if args.require_cuda else 0}",
            f"export RNG_SEED={int(args.rng_seed)}",
            f"export EXCLUDE_CELLTYPES_REGEX={shell_quote(args.exclude_celltypes_regex)}",
            f"export COORDINATE_KEY={shell_quote(args.coordinate_key)}",
            "export INITIAL_PRETRANSFORM_RECIPE=''",
            "export INITIAL_COORDINATES_CSV=''",
        ]
        if args.sigma2_init_scale != "":
            env_lines.append(f"export SIGMA2_INIT_SCALE={shell_quote(args.sigma2_init_scale)}")
        if args.sigma2_end != "":
            env_lines.append(f"export SIGMA2_END={shell_quote(args.sigma2_end)}")
        payload = [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "date",
            "hostname",
            "nvidia-smi || true",
            f"export RUNNER_SCRIPT={shell_quote(runner_script)}",
            *env_lines,
            "mkdir -p \"$OUTDIR\"",
            "set +u",
            "source ~/.bashrc || true",
            "if [ -f /home/HPCBase/tools/Anaconda3/etc/profile.d/conda.sh ]; then source /home/HPCBase/tools/Anaconda3/etc/profile.d/conda.sh; fi",
            "set -u",
            f"conda activate {shell_quote(args.conda_env)}",
            "python \"$RUNNER_SCRIPT\"",
            "date",
        ]
        payload_local.parent.mkdir(parents=True, exist_ok=True)
        payload_local.write_text("\n".join(payload) + "\n", encoding="utf-8")
        payload_local.chmod(0o755)
        rows.append(
            {
                "edge_id": edge["edge_id"],
                "stage": edge["stage"],
                "fixed_slice_id": edge["fixed_slice_id"],
                "moving_slice_id": edge["moving_slice_id"],
                "fixed_sl_number": edge["fixed_sl_number"],
                "moving_sl_number": edge["moving_sl_number"],
                "pair_start_order": edge["pair_start_order"],
                "edge_dir": str(edge_runtime),
                "local_edge_dir": str(edge_local),
                "payload_script": str(payload_runtime),
                "local_payload_script": str(payload_local),
                "job_name": job_name,
                "stage1_mode": args.stage1_mode,
                "stage2_mode": args.stage2_mode,
                "expression_mode": args.expression_mode,
                "dummy_rep_dim": args.dummy_rep_dim,
                "sigma2_init_scale": args.sigma2_init_scale,
                "sigma2_end": args.sigma2_end,
                "dsub_resource": args.dsub_resource,
            }
        )

    plan = {
        "version": "spateo_pairwise_initial_plan_v2_clean",
        "created_at": now_iso(),
        "dataset_root": str(dataset_root),
        "slice_order_csv": args.slice_order_csv or "",
        "stage": args.stage or "",
        "window_start": args.window_start or "",
        "window_end": args.window_end or "",
        "runner_script": str(runner_script),
        "project_account": args.project_account,
        "conda_env": args.conda_env,
        "dsub_resource": args.dsub_resource,
        "walltime": args.walltime,
        "stage1_mode": args.stage1_mode,
        "stage2_mode": args.stage2_mode,
        "expression_mode": args.expression_mode,
        "dummy_rep_dim": args.dummy_rep_dim,
        "sigma2_init_scale": args.sigma2_init_scale,
        "sigma2_end": args.sigma2_end,
        "local_output_dir": str(outdir),
        "runtime_output_dir": str(runtime_outdir),
        "n_selected_slices": len(order_rows),
        "n_edges": len(edges),
        "selected_slices": order_rows,
        "edges": edges,
    }
    json_dump(outdir / "pairwise_initial_plan.json", plan)
    write_csv_rows(outdir / "pairwise_jobs.csv", rows)

    submit = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        f"PLAN_DIR={shell_quote(runtime_outdir)}",
        f"PROJECT_ACCOUNT={shell_quote(args.project_account)}",
        f"DSUB_RESOURCE={shell_quote(args.dsub_resource)}",
        f"WALLTIME={int(args.walltime)}",
        "",
        "mkdir -p \"$PLAN_DIR/logs\" \"$PLAN_DIR/payloads\" \"$PLAN_DIR/edges\"",
    ]
    for row in rows:
        edge_runtime = Path(row["edge_dir"])
        submit.extend(
            [
                "",
                f"# {row['edge_id']}: {row['fixed_slice_id']} -> {row['moving_slice_id']}",
                f"mkdir -p {shell_quote(edge_runtime / 'logs')}",
                "dsub \\",
                f"  -A {shell_quote(args.project_account)} \\",
                f"  -n {shell_quote(row['job_name'])} \\",
                "  -T \"$WALLTIME\" \\",
                f"  -oo {shell_quote(edge_runtime / 'logs' / (row['job_name'] + '.%J.out'))} \\",
                f"  -eo {shell_quote(edge_runtime / 'logs' / (row['job_name'] + '.%J.err'))} \\",
                "  -R \"$DSUB_RESOURCE\" \\",
                f"  -s {shell_quote(row['payload_script'])}",
            ]
        )
    submit_path = outdir / "submit_pairwise_spateo.sh"
    submit_path.write_text("\n".join(submit) + "\n", encoding="utf-8")
    submit_path.chmod(0o755)
    print(outdir / "pairwise_initial_plan.json")


if __name__ == "__main__":
    main()
