---
name: remote-workflow-intake
description: Collect the minimal remote execution context before spatial alignment work starts. Use when a task will run on a remote cluster or needs SSH host, remote workdir, source h5ad/raw-data path, and conda environment before using alignment, Spateo, or visualization skills.
---

# Remote Workflow Intake

## Package paths

The configuration template is bundled in `references/`.
This intake skill can be installed independently.


Use this skill at the start of any remote spatial alignment task. Its job is to
collect the minimum context needed to begin safely without hardcoded paths.

Always apply the repository workflow principles:

- source h5ad/raw data stays remote;
- source h5ad obs identity is the authority for cell identities. Before
  alignment work, confirm whether `obs_names`, `obs["cell_id"]`,
  `obs["CellID"]`, and slice fields are present; downstream clean outputs must
  derive cell ids from these fields, not from newly invented row numbers;
- Spateo jobs, ROI/drop refinement, and recipe extraction run on the remote
  server because they depend on remote raw data and cluster resources;
- local review uses an immutable local coordinate cache plus remote-produced
  `recipe.json`/`recipe_chain.json` and small QC metadata, not expression
  matrices or raw h5ad;
- deterministic recipe application may run both remotely and locally, but both
  sides must use the same recipe chain, same transform code version, same input
  coordinate columns, and the same output format;
- avoid transferring full coordinate outputs during routine iteration; transfer
  recipes, provenance, audit JSON, logs, and small QC tables instead;
- remote workflows default to `GENERATE_HTML=0`.
- local coordinate outputs are allowed only as replay outputs in a new directory
  from a preserved original local coordinate cache. Never overwrite the local
  original cache.
- verify local replay against the remote reference result with `sha256` when
  byte-identical output is expected, or with a coordinate hash plus numeric
  tolerance when CSV formatting can differ.

## Workflow

1. Identify the intended task: diagnosis, candidate search, recipe replay,
   Spateo rerun, local viewer generation, or data/cache setup.
2. Inspect existing local files or provided paths first. If SSH access is
   available and the user has asked to run remotely, inspect remote paths with
   read-only commands before asking for details.
3. Build a Minimal Remote Context summary with the four required fields below.
4. Mark each field as `confirmed`, `missing`, or `assumed`.
5. Ask only for missing required fields that cannot be discovered safely.
6. Defer task-specific fields until the downstream skill actually needs them.

## Minimal Required Context

Collect only these fields before starting:

- remote host or SSH alias;
- remote repo/workdir path;
- source h5ad/raw-data root on remote;
- conda environment path or name;

For spatial h5ad inputs, also note the expected identity fields if they are
known: `obs_names`, `cell_id`, `CellID`, slice column, and celltype column.
These fields should be verified by read-only inspection before final export.

Everything else is task-specific. Ask for it later only when needed.

## Task-Specific Fields

Collect these only when a downstream task requires them:

- full coordinate table path on remote;
- local original coordinate cache path, treated as read-only;
- coordinate columns or keys;
- recipe path or candidate root;
- local and remote transform script paths or code hashes;
- output directory;
- scheduler/account/resource settings;
- slice/stage window;
- transfer policy;
- explicit sampling caps;
- `GENERATE_HTML`, which defaults to `0`.

## Config Output

When preparing only the initial context, produce a concise block:

```bash
SSH_HOST=HOST_OR_ALIAS
WORKDIR=/remote/path/to/spateo_interactive
DATASET_ROOT=/remote/path/to/input_h5ad_directory
CONDA_ENV=/remote/path/to/conda/env
```

If the next task needs a full `project_paths.env`, extend this block using
`references/project_paths.env.example`. Use placeholders only in
drafts. Do not submit jobs with placeholders.

## Read-Only Checks

Use read-only checks before running or writing:

```bash
ssh HOST 'test -d DATASET_ROOT && find DATASET_ROOT -name "*.h5ad" | head'
ssh HOST 'test -d WORKDIR && test -x "$(command -v bash)"'
ssh HOST 'source ~/.bashrc >/dev/null 2>&1 || true; conda env list | head'
```

Adapt commands to the actual host and shell. Avoid commands that copy, delete,
or rewrite data during intake.

## Do Not Do

- Do not infer private project paths from older conversations when a current
  path is needed.
- Do not transfer `.h5ad`, expression matrices, raw imaging data, generated
  HTML, or large corrected coordinate copies during routine iteration unless
  the user explicitly asks to share or inspect a final remote-produced cache.
- Do not run raw-data or cluster-dependent alignment operations locally.
  Pairwise Spateo jobs, ROI/drop refinement, and recipe extraction must happen
  on the remote server.
- Do not overwrite the local original coordinate cache. Local recipe replay must
  write to a new output directory and preserve input hashes.
- Do not claim a local corrected CSV/h5ad is consistent with remote unless a
  recorded checksum or numeric replay comparison passes.
- Do not generate HTML remotely unless the user explicitly sets
  `GENERATE_HTML=1` for an interactive/local run.
- Do not silently downsample visualization. Use full coordinates by default.
- Do not submit cluster jobs until required path/env/account fields are
  confirmed for that specific task.

## Handoff

After intake, hand off to the task-specific skill with the Remote Context
summary. Example handoffs:

- `spateo-pairwise-run`: raw h5ad root, runner script, conda env,
  scheduler/account, stage/window, run mode, and output dir;
- `spateo-pairwise-qc`: remote pairwise run dir, remote full-coordinate CSV,
  coordinate columns, slice column, and cell-type column;
- `spateo-roi-refine`: suspect edge id, ROI/drop condition, pairwise run dir,
  remote full-coordinate CSV, raw h5ad root, runner script, conda env, and
  scheduler/account;
- `spatial-alignment-compose`: confirmed remote recipes, remote and/or local
  original coordinate cache, output directories, transform code path/hash, and
  consistency-check target;
- `spatial-pointcloud-viewer`: local corrected replay output or copied remote
  coordinate cache and viewer output path.
