# Source and licensing notices

The root [MIT license](LICENSE) applies to original repository material except where the notices below identify imported material with a different or unspecified source license.

## Spateo environment skill

`skills/setup-spateo-environment/` is adapted from the identically named directory in `gmhhhhhh-929/spateo-release` at commit `d6aa68addc475dd0b56f69cebe7823b1f79933a9`.

That source repository's `LICENSE` contains the BSD 2-Clause license, copyright 2022 Aristotle. Its full text is preserved in [licenses/spateo-source-BSD-2-Clause.txt](licenses/spateo-source-BSD-2-Clause.txt). The source package metadata labels its license BSD-3-Clause; this repository preserves the actual license file rather than silently replacing its terms.

## Data IO skill

The Data IO instructions and helper code were prepared from a review of that same Spateo commit. References identify the maintained source APIs; the complete Spateo implementation and its datasets are not vendored into this skill.

## Alignment snapshots

The imported alignment implementation files are enumerated in [source_migration.json](skills/spateo-2d-alignment/provenance/source_migration.json). They originate from the user-provided pairwise and continuity snapshots. No explicit license file was found at either source repository root, and the selected Python files contained no license header.

Their inclusion preserves that source licensing status; this notice does not assert that the imported snapshot files have been relicensed under the root MIT license. Source identifiers, hashes, preserved scientific symbols, and packaging changes remain available in the migration manifest. Newly authored wrappers and instructions are covered by the root license unless separately noted.

## Restored workflow companions

The 14 restored entrypoints and fixed support tools originate from the user-provided server skill collection. Imported files and their source/package hashes are recorded in [companion_migration.json](skills/spateo-2d-alignment/provenance/companion_migration.json). Their unspecified source licensing status is preserved under the same alignment-snapshot notice above. The continuity-first entry is a new compatibility adapter to the existing published pipeline, not a copy of an additional historical engine.

## Spateo Referee QC and viewer

The two bundled scientific Referee modules derive from the user-provided Spateo checkout at baseline commit `94fb2f4b71b6e809d3c2341d0e209c30b937573a`, with subsequent QC, display-preregistration and experimental-policy changes. The source BSD-2-Clause license is preserved in [runtime/LICENSE](skills/spatial-slice-quality-qc/runtime/LICENSE). The current skill wrappers and renderers originate from the user-provided Referee skills; packaging hashes and provenance are recorded in [release.json](skills/spatial-slice-quality-qc/provenance/release.json). This snapshot is not claimed to be identical to the baseline commit. No biological point clouds or source matrices are distributed.

## Native IO and 4D migration (2026-09-21)

IO/4D API contracts are source-reviewed against gmhhhhhh-929/spateo-release at `615644f88613bea8ceb2e2df1e2391d16de55ec1`. Scientific algorithms are invoked from that separately installed library; they are not vendored by the new runner. The 4D workflow and portable dashboard were migrated from the user's local `spateo-skills` collection, with the user's Spateo-protocol-files notebooks (`b11ae99fbdc4ae46d41880e9306ab7e5c2751ac5`) used as scientific workflow references. No notebook dataset or image output is redistributed. The existing source notices and licenses remain in force.
