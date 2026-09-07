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
