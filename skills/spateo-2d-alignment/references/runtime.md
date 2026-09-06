# Runtime compatibility

Direct external imports used by the bundled execution paths are `anndata`, `h5py`, `numpy`, `pandas`, `scipy`, `matplotlib`, `torch`, `POT` (`ot`) and `spateo-release` (`spateo`). There is no direct external Dynamo requirement in the retained pipeline closure. Use the dependency versions supported by the installed Spateo distribution rather than installing arbitrary latest packages into a working environment.

`scripts/check_runtime.py --device cpu` verifies imports, distribution metadata and explicit alignment parameters. For a CUDA run use the intended device index, for example `--device 0`. It checks both public `morpho_align`/`morpho_align_ref` and the `Morpho_pairwise` constructor because many public keywords are forwarded to that constructor. A passing signature check is not proof of numerical or GPU compatibility.

The intended environment skill's selected source revision was inspected statically: its public alignment functions and constructor accept the preserved kwargs, and alignment downsampling uses native sampling. The exact revision, file hashes and inspection scope are in `provenance/api_compatibility.json`. No code is rewritten to a different alignment API, and no new-main GPU accuracy claim is made.

Useful checks, run with the same Python interpreter intended for the job:

```bash
python /path/to/spateo-2d-alignment/scripts/verify_bundle.py
python /path/to/spateo-2d-alignment/scripts/check_runtime.py --device 0
python /path/to/spateo-2d-alignment/scripts/smoke_test.py
```

The smoke test uses temporary synthetic cells to exercise loaders, feature validation, ID conventions, dependency imports, a proper rigid transform and continuity scores. It calls neither Spateo alignment nor GPU kernels. If a new Spateo installation fails at runtime, preserve the failure and investigate the specific API/dependency problem without silently replacing the algorithm or claiming a repeated benchmark.

The bounded synthetic CPU execution report is in [cpu_smoke.json](../provenance/cpu_smoke.json). It covers two slices and 160 cells with three iterations for expression PCA, spatial-only and continuity modes in the recorded existing environment. It does not establish GPU accuracy or runtime compatibility with the newer pinned main checkout.
