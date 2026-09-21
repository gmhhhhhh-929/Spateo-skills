# Configuration v3

The canonical schema is `spateo-4d/v3`; the template contains all accepted fields. Sections merge with explicit defaults. Unknown keys fail. Old notebook/v2 JSON requires deliberate migration rather than silent compatibility.

| Section | Meaning / invalidation |
| --- | --- |
| inputs | Absolute or config-relative stage1/stage2 H5AD and declared coordinate_unit; changed contents invalidate alignment even at the same filename. |
| labels | Source/target labels recorded in the config; changes conservatively invalidate alignment. |
| alignment | spatial_key, counts_layer, explicit x_is_counts assertion, log_layer, aligned_key, target_sum, mode SN-S/SN-N, n_sampling, sampling_method, max_iter, nonrigid_start_iter. Changes invalidate all stages. |
| subset | annotation_key and group (null means all cells). Changes invalidate mapping and descendants. |
| mapping | key, alpha, numItermax, numItermaxEmd, normalization target_sum, max_pairs. Changes invalidate mapping and descendants. |
| morphofield | key, M, lambda_, restart_num, restart_seed, MaxIter. Fit at source points; these points also define NX evaluation, avoiding mandatory mesh reconstruction. |
| trajectory | enabled, key, positive t_end, interpolation_num and direction. Changes rerun trajectory and dashboard only. Native implementation currently ignores cores/layer options; the runner does not promise multiprocessing. |
| metrics | enabled, selected geometric quantities, glm_metrics, explicit glm_genes, qval_threshold, llf_threshold. GLMs use normalized (not log) expression and a real interpolated formula. |
| gp | enabled, verified genes, training_iter, method, inducing_num. Space→expression interpolation evaluated at source points and saved independently. |
| dashboard | enabled, max_points, max_target_points, max_vectors, default_feature, explicit cdn. Display changes rebuild only dashboard. |
| runtime | device and seed. Changes conservatively invalidate all computations. |

`alignment.aligned_key` is used consistently throughout. `mode=SN-S` selects the rigid output; `SN-N` selects nonrigid. The source also produces `_rigid` and `_nonrigid` keys, but they are not silently substituted. Normalization uses target_sum=10000 by default; the protocol's target_sum=None (library median) can be explicitly chosen if reproducing that preprocessing, subject to validation.

Optional trajectory/metrics/gp/dashboard stages can be disabled. GP defaults off with no example organism genes. GLMs default to an empty selection; enabling them requires verified gene IDs and sufficiently variable features. The runner does not hard-code CNS, GPU 0, arbitrary z offsets, mesh smoothing or notebook cameras.

Checkpoints are stage-specific. Mapping depends on alignment; morphofield on mapping; trajectory and metrics independently on morphofield; GP on mapping. Dashboard depends on their current states. Changes to source files, Spateo Python files, skill Python files or recorded package versions invalidate reuse. Reuse additionally checks every output's SHA256; missing/modified artifacts force that stage and descendants to run.

`--stop-after alignment` produces an explicitly partial run. A child using `--parent-manifest` resumes from valid completed checkpoints. All runs use fresh directories and no overwrite option. Old v2 manifests cannot be reused as v3 checkpoints.

`max_iter` must exceed `nonrigid_start_iter + 1`: the pinned reference alignment source initializes deformation coefficients only after that boundary. Tiny smoke tests explicitly lower both settings; production defaults retain 200/80.

AnnData 0.10 cannot write native preprocessing history lists of dictionaries directly. The runner preserves lists/tuples/None as reversible tagged mappings in `.uns` (`__spateo_skill_type__`); `load_pair()` restores them before native calls. Integer mapping keys are persisted as strings, including trajectory indices. Numeric scientific arrays are unchanged. Use the helper when resuming manually; do not discard preprocessing history to make H5AD writing succeed.
