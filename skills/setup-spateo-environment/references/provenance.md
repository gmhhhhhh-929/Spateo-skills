# Source and compatibility profile

This skill is adapted from [`skills/setup-spateo-environment`](https://github.com/gmhhhhhh-929/spateo-release/tree/d6aa68addc475dd0b56f69cebe7823b1f79933a9/skills/setup-spateo-environment) in `gmhhhhhh-929/spateo-release`, commit `d6aa68addc475dd0b56f69cebe7823b1f79933a9`.

The source `environment.yml`, `requirements.txt`, `setup.py`, and package APIs at that revision define the compatibility profile. This publication does not claim that an arbitrary PyPI Spateo version or a future `main` revision satisfies it. Dependencies are constrained ranges, not a fully resolved lockfile; record installed versions after setup.

The original skill's preprocessing/native-vector-field and mesh smoke computations are retained. The standalone adaptation replaces its automatic ancestor-directory import override with explicit installed-package and optional editable-checkout identity checks, adds a metadata-only diagnostic mode, and distinguishes runtime checks from complete pipeline validation.

Static review at this revision finds the two packaged legacy alignment runners' used keywords in `morpho_align`, `morpho_align_ref`, and `Morpho_pairwise.__init__`. Alignment downsampling uses `spateo._native.sample`, with no direct external Dynamo import. The runtime checker supplied with the alignment skill verifies imports and signatures on the actual environment. Static compatibility and small smoke tests do not establish end-to-end alignment or GPU performance.

The source [LICENSE](https://github.com/gmhhhhhh-929/spateo-release/blob/d6aa68addc475dd0b56f69cebe7823b1f79933a9/LICENSE) contains the BSD 2-Clause license, copyright © 2022 Aristotle. The full license is retained inside this skill as [LICENSE](../LICENSE), so it accompanies a standalone copy. The repository NOTICE records the adaptation and other source material. The source package metadata says BSD-3-Clause, which differs from the actual LICENSE; this skill follows the LICENSE text. Preserve the copyright notice, both conditions, and disclaimer when redistributing source or binary forms.
