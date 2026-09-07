#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from pairwise_core import json_dump, now_iso, read_csv_rows, sanitize_token, shell_quote, write_csv_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Create dsub scripts for pairwise Spateo refine reruns from candidate recipes.")
    parser.add_argument("--candidate-ranking", required=True)
    parser.add_argument("--dataset-root", required=True)
    parser.add_argument("--runner-script", required=True)
    parser.add_argument("--project-account", required=True)
    parser.add_argument("--conda-env", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--runtime-output-dir")
    parser.add_argument("--top-k", type=int, default=3)
    parser.add_argument("--stage1-mode", default="SN-S")
    parser.add_argument("--stage2-mode", default="none")
    parser.add_argument("--stage1-max-iter", type=int, default=300)
    parser.add_argument("--stage2-max-iter", type=int, default=300)
    parser.add_argument("--n-sampling-ref", type=int, default=40000)
    parser.add_argument("--max-cells-per-slice", type=int, default=0)
    parser.add_argument("--batch-size", type=int, default=800)
    parser.add_argument("--no-require-cuda", dest="require_cuda", action="store_false")
    parser.set_defaults(require_cuda=True)
    parser.add_argument("--rng-seed", type=int, default=20260612)
    parser.add_argument("--dsub-resource", default="cpu=64;mem=240000;gpu=1")
    parser.add_argument("--walltime", type=int, default=345600)
    parser.add_argument("--job-name-prefix", default="spateo_pair_refine")
    args = parser.parse_args()

    outdir = Path(args.output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    runtime_outdir = Path(args.runtime_output_dir).expanduser() if args.runtime_output_dir else outdir

    def runtime_path(local_path: Path) -> Path:
        if not args.runtime_output_dir:
            return local_path
        return runtime_outdir / local_path.relative_to(outdir)

    candidates = read_csv_rows(Path(args.candidate_ranking).expanduser().resolve())[: args.top_k]
    if not candidates:
        raise SystemExit("No candidates found.")

    jobs = []
    for row in candidates:
        candidate_id = row["candidate_id"]
        edge_id = row["edge_id"]
        edge_dir = outdir / "reruns" / candidate_id
        edge_runtime = runtime_path(edge_dir)
        payload_local = outdir / "payloads" / f"{sanitize_token(candidate_id)}.sh"
        payload_runtime = runtime_path(payload_local)
        recipe_path = Path(row["transform_recipe"]).expanduser()
        recipe_runtime = recipe_path
        job_name = (args.job_name_prefix + "_" + sanitize_token(edge_id) + "_" + sanitize_token(candidate_id))[:120]
        edge_dir.mkdir(parents=True, exist_ok=True)
        (edge_dir / "logs").mkdir(parents=True, exist_ok=True)
        provenance = {
            "version": "spateo_pairwise_refine_rerun_edge_v1",
            "created_at": now_iso(),
            "status": "planned",
            "edge_id": edge_id,
            "candidate_id": candidate_id,
            "target_celltype": row.get("target_celltype", ""),
            "fixed_slice_id": row.get("fixed_slice_id", ""),
            "moving_slice_id": row.get("moving_slice_id", ""),
            "fixed_sl_number": row.get("fixed_sl_number", ""),
            "moving_sl_number": row.get("moving_sl_number", ""),
            "pair_start_order": row.get("pair_start_order", ""),
            "dataset_root": args.dataset_root,
            "runner_script": args.runner_script,
            "initial_pretransform_recipe": str(recipe_runtime),
            "candidate": row,
            "runner_parameters": {
                "stage1_mode": args.stage1_mode,
                "stage2_mode": args.stage2_mode,
                "stage1_max_iter": args.stage1_max_iter,
                "stage2_max_iter": args.stage2_max_iter,
                "n_sampling_ref": args.n_sampling_ref,
                "max_cells_per_slice": args.max_cells_per_slice,
                "batch_size": args.batch_size,
                "require_cuda": args.require_cuda,
                "rng_seed": args.rng_seed,
            },
            "dsub": {"project_account": args.project_account, "resource": args.dsub_resource, "walltime": args.walltime, "job_name": job_name, "job_id": ""},
        }
        json_dump(edge_dir / "provenance.json", provenance)
        start_order = row.get("pair_start_order") or row.get("fixed_global_order_by_SL") or row.get("fixed_order") or ""
        env_lines = [
            f"export DATASET_ROOT={shell_quote(args.dataset_root)}",
            f"export OUTDIR={shell_quote(edge_runtime)}",
            f"export PAIR_START_ORDER={shell_quote(start_order)}",
            "export N_SLICES=2",
            f"export STAGE1_MODE={shell_quote(args.stage1_mode)}",
            f"export STAGE2_MODE={shell_quote(args.stage2_mode)}",
            f"export STAGE1_MAX_ITER={int(args.stage1_max_iter)}",
            f"export STAGE2_MAX_ITER={int(args.stage2_max_iter)}",
            f"export N_SAMPLING_REF={int(args.n_sampling_ref)}",
            f"export MAX_CELLS_PER_SLICE={int(args.max_cells_per_slice)}",
            f"export BATCH_SIZE={int(args.batch_size)}",
            f"export REQUIRE_CUDA={1 if args.require_cuda else 0}",
            f"export RNG_SEED={int(args.rng_seed)}",
            f"export INITIAL_PRETRANSFORM_RECIPE={shell_quote(recipe_runtime)}",
            "export INITIAL_COORDINATES_CSV=''",
        ]
        payload = [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "date",
            "hostname",
            "nvidia-smi || true",
            f"export RUNNER_SCRIPT={shell_quote(args.runner_script)}",
            *env_lines,
            "if [ -z \"${PAIR_START_ORDER}\" ]; then echo 'PAIR_START_ORDER is empty; regenerate plan with pair_start_order in candidate_ranking.csv' >&2; exit 2; fi",
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
        jobs.append(
            {
                "candidate_id": candidate_id,
                "edge_id": edge_id,
                "target_celltype": row.get("target_celltype", ""),
                "rank_score": row.get("rank_score", ""),
                "pair_start_order": start_order,
                "rerun_edge_dir": str(edge_runtime),
                "local_rerun_edge_dir": str(edge_dir),
                "payload_script": str(payload_runtime),
                "local_payload_script": str(payload_local),
                "job_name": job_name,
                "initial_pretransform_recipe": str(recipe_runtime),
            }
        )

    write_csv_rows(outdir / "pairwise_refine_rerun_jobs.csv", jobs)
    plan = {
        "version": "spateo_pairwise_refine_rerun_plan_v1",
        "created_at": now_iso(),
        "candidate_ranking": str(Path(args.candidate_ranking).expanduser().resolve()),
        "dataset_root": args.dataset_root,
        "runner_script": args.runner_script,
        "project_account": args.project_account,
        "conda_env": args.conda_env,
        "dsub_resource": args.dsub_resource,
        "stage1_mode": args.stage1_mode,
        "stage2_mode": args.stage2_mode,
        "n_jobs": len(jobs),
        "jobs": jobs,
    }
    json_dump(outdir / "pairwise_refine_rerun_plan.json", plan)

    submit = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        f"DSUB_RESOURCE={shell_quote(args.dsub_resource)}",
        f"WALLTIME={int(args.walltime)}",
    ]
    for job in jobs:
        edge_runtime = Path(job["rerun_edge_dir"])
        submit.extend(
            [
                "",
                f"# {job['candidate_id']} {job['edge_id']}",
                f"mkdir -p {shell_quote(edge_runtime / 'logs')}",
                "dsub \\",
                f"  -A {shell_quote(args.project_account)} \\",
                f"  -n {shell_quote(job['job_name'])} \\",
                "  -T \"$WALLTIME\" \\",
                f"  -oo {shell_quote(edge_runtime / 'logs' / (job['job_name'] + '.%J.out'))} \\",
                f"  -eo {shell_quote(edge_runtime / 'logs' / (job['job_name'] + '.%J.err'))} \\",
                "  -R \"$DSUB_RESOURCE\" \\",
                f"  -s {shell_quote(job['payload_script'])}",
            ]
        )
    submit_path = outdir / "submit_pairwise_refine_spateo.sh"
    submit_path.write_text("\n".join(submit) + "\n", encoding="utf-8")
    submit_path.chmod(0o755)
    print(outdir / "pairwise_refine_rerun_plan.json")


if __name__ == "__main__":
    main()
