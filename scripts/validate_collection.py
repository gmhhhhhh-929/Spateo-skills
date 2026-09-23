#!/usr/bin/env python3
"""Offline packaging checks. Does not run alignment or read biological data."""
from pathlib import Path
import ast
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / 'skills/spateo-2d-alignment'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    manifest = json.loads((MAIN / 'provenance/companion_migration.json').read_text())
    skills = sorted((ROOT / 'skills').rglob('SKILL.md'))
    assert len(list((ROOT / 'skills').glob('*/SKILL.md'))) == 6
    assert len(list((MAIN / 'subskills').glob('*/SKILL.md'))) == 14
    assert manifest['total_skill_count'] == 17
    assert len(skills) == 28
    assert len(list((ROOT / 'skills/spateo-4d-pipeline/subskills').glob('*/SKILL.md'))) == 5
    assert len(list((ROOT / 'skills/spateo-3d-pipeline/subskills').glob('*/SKILL.md'))) == 2
    assert (ROOT / 'skills/spatial-slice-quality-qc/subskills/spatial-slice-quality-viewer/SKILL.md').is_file()
    for skill in skills:
        text = skill.read_text()
        assert text.startswith('---\n')
        assert re.search(r'^name: ' + re.escape(skill.parent.name) + r'$', text, re.M), skill
        assert re.search(r'^description: .+', text, re.M), skill
    for row in manifest['inventory']:
        assert (ROOT / row['path']).is_file(), row
    for row in manifest['files']:
        assert digest(ROOT / row['path']) == row['packaged_sha256'], row['path']
    # The historical migration preserved IO too; IO is now explicitly rewritten.
    preserved = [row for row in manifest['preserved_published_python']
                 if row['path'] not in {'skills/spateo-data-io/scripts/spateo_io.py',
                                        'skills/spateo-data-io/scripts/smoke_source.py'}]
    for row in preserved:
        assert digest(ROOT / row['path']) == row['sha256'], row['path']
    runtime = MAIN / 'pipelines/pairwise-rigid'
    lock = json.loads((runtime / 'skill.lock.yaml').read_text())
    for mode, entry in lock['allowed_entrypoints'].items():
        assert digest(runtime / entry['relative_path']) == entry['sha256'], mode
    for path in lock['validators'].values():
        assert (runtime / path).is_file(), path
    for line in (runtime / 'release_manifest.sha256').read_text().splitlines():
        expected, path = line.split(maxsplit=1)
        assert digest(runtime / path.lstrip('*')) == expected, path
    for name in ('spatial-before-after-viewer', 'spatial-pointcloud-viewer'):
        folder = MAIN / 'subskills' / name
        viewer = json.loads((folder / 'provenance/viewer.json').read_text())
        assert digest(folder / viewer['script']) == viewer['sha256'], name
    referee = json.loads((ROOT / 'skills/spatial-slice-quality-qc/provenance/release.json').read_text())
    for row in referee['files']:
        assert digest(ROOT / row['path']) == row['sha256'], row['path']
    pyfiles = list(ROOT.glob('skills/**/*.py'))
    for path in pyfiles:
        ast.parse(path.read_text(), filename=str(path))
    broken = []
    for path in ROOT.rglob('*.md'):
        if '.git' in path.parts:
            continue
        for link in re.findall(r'\]\(([^)]+)\)', path.read_text()):
            if '://' in link or link.startswith(('#', 'mailto:')):
                continue
            target = link.split('#')[0]
            if target and not (path.parent / target).exists():
                broken.append((str(path.relative_to(ROOT)), target))
    assert not broken, broken
    # Subskill command references must resolve from the documented cwd.
    for skill in (MAIN / 'subskills').glob('*/SKILL.md'):
        for relative in re.findall(r'(?<![\w./])(?:\.\./)+(?:[\w./-]+\.(?:py|yaml))', skill.read_text()):
            assert (skill.parent / relative).is_file(), (skill, relative)
    assert sorted(p.name for p in (MAIN / 'pipelines').iterdir() if p.is_dir()) == ['continuity-guided', 'pairwise-rigid']
    print(json.dumps({'status': 'pass', 'skills': len(skills),
        'top_level_skills': 6, 'three_d_subskills': 2, 'four_d_subskills': 5,
        'qc_subskills': 1, 'alignment_subskills': 14,
        'source_entries_covered': len(manifest['inventory']),
        'locked_entrypoints': len(lock['allowed_entrypoints']),
        'python_syntax_checked': len(pyfiles),
        'preserved_published_python': len(preserved),
        'viewer_hashes': 'pass', 'markdown_links': 'pass'}, indent=2))

if __name__ == '__main__':
    main()
