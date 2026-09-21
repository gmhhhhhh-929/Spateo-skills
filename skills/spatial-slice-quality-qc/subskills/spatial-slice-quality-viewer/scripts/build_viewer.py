"""Use the QC skill's canonical detailed renderer; keep computation and UI separate."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', required=True, nargs='+')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--qc-skill', default=str(Path(__file__).resolve().parents[3]))
    parser.add_argument('--spateo-source')
    parser.add_argument('--policy', help='Optional existing frozen policy; run complete two-stage publication before rendering')
    parser.add_argument('--application-scope', choices=['certified','new_input_unvalidated','experimental_policy'], default='certified')
    parser.add_argument('--title', default='Spateo Referee')
    parser.add_argument('--language', choices=['en', 'zh'], default='en')
    parser.add_argument('--display-window', type=int, choices=[3, 5], default=5)
    parser.add_argument('--max-points-per-slice', type=int, default=900)
    parser.add_argument('--source-h5ad', help='Optional single-dataset source for expression and full component evidence')
    args = parser.parse_args()
    inputs = [Path(p).expanduser().resolve() for p in args.input_dir]
    output = Path(args.output_dir).expanduser().resolve()
    if output.exists():
        parser.error('Output must be a new directory; preserve previous reports.')
    if len({p.name for p in inputs}) != len(inputs):
        parser.error('Batch input directory names must be distinct.')
    if args.source_h5ad and len(inputs) != 1:
        parser.error('--source-h5ad applies to one dataset; batch uses recorded source contracts/caches.')
    scripts = Path(args.qc_skill).expanduser().resolve() / 'scripts'
    if not (scripts/'slice_quality_visualization.py').is_file():
        parser.error('The canonical QC skill renderer was not found: '+str(scripts))
    sys.path.insert(0, str(scripts))
    if args.spateo_source:
        os.environ['SPATEO_REFEREE_SOURCE'] = str(Path(args.spateo_source).expanduser().resolve())
    source = os.environ.get('SPATEO_REFEREE_SOURCE')
    if source or (scripts.parent/'runtime/spateo/preprocessing/slice_quality.py').is_file():
        from activate_candidate import activate
        activate()
    from slice_quality_visualization import render_binary_slice_quality_appendix
    from slice_quality_report import write_slice_quality_collection_report
    import pandas as pd
    for directory in inputs:
        for name in ['slice_quality_metrics.csv', 'slice_quality_manifest.json', 'slice_quality_display_payload.json']:
            if not (directory/name).is_file():
                parser.error('Missing QC prerequisite: '+str(directory/name))
        if not args.policy and not (directory/'slice_quality_binary_audit.csv').is_file():
            parser.error('Missing complete binary audit: '+str(directory)+'. Use the QC skill publish --complete-binary with an applicable frozen policy, or provide --policy. The viewer never relabels review.')
    output.mkdir(parents=True)
    runs=[];ledger=[]
    for directory in inputs:
        dest=output/directory.name if len(inputs)>1 else output
        dest.mkdir(exist_ok=True)
        # Independent report package; coordinate caches remain immutable and are reused.
        for name in ['slice_quality_metrics.csv','slice_quality_manifest.json','slice_quality_display_payload.json','slice_quality_binary_audit.csv','slice_quality_binary_calls.csv','binary_policy_application.json']:
            if (directory/name).exists():shutil.copy2(directory/name,dest/name)
        if (directory/'display_preregistration').exists():
            (dest/'display_preregistration').symlink_to(directory/'display_preregistration',target_is_directory=True)
        if args.policy:
            from run_slice_quality_qc import run_publish
            run_publish(SimpleNamespace(input_dir=str(directory),output_dir=str(dest),policy=args.policy,complete_binary=True,allow_unvalidated_policy=False,application_scope=args.application_scope))
        expected=pd.read_csv(dest/'slice_quality_metrics.csv',dtype={'slice_id':str}).sort_values('slice_index')
        audit=pd.read_csv(dest/'slice_quality_binary_audit.csv',dtype={'slice_id':str}).sort_values('slice_index')
        if audit.slice_id.duplicated().any() or audit.slice_id.tolist()!=expected.slice_id.tolist():
            raise ValueError('Binary audit must cover each input slice exactly once in the original order.')
        if not audit.final_call.isin(['keep','exclude']).all():
            raise ValueError('The final viewer accepts only complete keep/exclude calls.')
        result=render_binary_slice_quality_appendix(dest,dest,title=args.title+(' · '+directory.name if len(inputs)>1 else ''),language=args.language,display_window=args.display_window,max_points_per_slice=args.max_points_per_slice,source_h5ad=args.source_h5ad)
        # The existing collection builder resolves this conventional report name.
        (dest/'slice_quality_report.html').symlink_to('index.html')
        runs.append(dest)
        ledger.append({'input_dir':str(directory),'output':result,'binary_audit_sha256':hashlib.sha256((dest/'slice_quality_binary_audit.csv').read_bytes()).hexdigest(),'full_roi_replay_input':str(directory)})
    if len(runs)>1:
        write_slice_quality_collection_report(runs,output/'index.html',title=args.title,binary_only=True)
        if any(item['output'].get('application_scope') == 'new_input_unvalidated' for item in ledger):
            index=output/'index.html'
            index.write_text(index.read_text().replace('</header>', '<p><b>Frozen-policy application to new inputs; these inputs have not been independently validated.</b></p></header>'))
        if any(item['output'].get('application_scope') == 'experimental_policy' for item in ledger):
            index=output/'index.html'
            index.write_text(index.read_text().replace('</header>', '<p><b>Experimental joint-review policy: metric-stress evaluation only; no independent biological validation of this policy or these inputs.</b></p></header>'))
    provenance={'renderer':'slice_quality_visualization.render_binary_slice_quality_appendix','qc_skill':str(scripts.parent),'scripts':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [scripts/'slice_quality_visualization.py',scripts/'slice_quality_detail_roi.js',Path(__file__)]},'runs':ledger,'source_h5ad_modified':False}
    (output/'viewer_workflow.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2))
    print(json.dumps({'viewer':str(output/'index.html'),'datasets':len(runs),'provenance':str(output/'viewer_workflow.json')},ensure_ascii=False))

if __name__=='__main__':main()
