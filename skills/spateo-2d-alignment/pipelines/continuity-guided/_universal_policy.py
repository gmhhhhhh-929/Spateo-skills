"""Anonymous, cross-label-checked terminal proposals for continuity-guided alignment."""
import copy
import hashlib
import math
import numpy as np

VERSION='anonymous-terminal-v1'

def install(native):
    original=native.priority_lineage_terminal_rescue
    def rescue(slices,models,transforms,terminal_indices,args):
        shared=set(slices[terminal_indices[0]].annotations.astype(str))
        for i in terminal_indices[1:]:shared.intersection_update(slices[i].annotations.astype(str))
        eligible=[]
        for label in shared:
            counts=[np.sum(slices[i].annotations.astype(str)==label) for i in terminal_indices]
            if min(counts)<40:continue
            supports=[];spreads=[];member_ids=[]
            for i in terminal_indices:
                mask=slices[i].sample_annotations.astype(str)==label
                group=slices[i].sample_xy[mask]
                if len(group)<8:break
                supports.append(len(group));spreads.append(native.core.robust_diagonal(group)/native.core.robust_diagonal(slices[i].sample_xy))
                member_ids.extend(slices[i].cell_ids[slices[i].annotations.astype(str)==label].astype(str).tolist())
            if len(supports)!=len(terminal_indices):continue
            # ID membership resolves ties without relying on annotation names.
            identity=hashlib.sha256('\n'.join(sorted(member_ids)).encode()).hexdigest()
            reliability=math.sqrt(min(supports))/(0.15+float(np.median(spreads)))
            eligible.append((reliability,identity,label))
        ranked=sorted(eligible,key=lambda x:(-x[0],x[1]))[:3]
        candidates=[];accepted=[]
        if len(eligible)<3:
            return transforms,{'operation':'anonymous_terminal_consensus','accepted':False,'reason':'fewer_than_three_supported_independent_labels','candidate_count':0}
        def other_scores(matrices,priority):
            values={}
            for _,identity,label in eligible:
                if label==priority:continue
                pairs=[]
                for l,r in zip(terminal_indices[:-1],terminal_indices[1:]):
                    v=native.lineage_pair_distance(native.core.aligned_sample(slices[l],matrices[l]),slices[l].sample_annotations,native.core.aligned_sample(slices[r],matrices[r]),slices[r].sample_annotations,label)
                    if np.isfinite(v):pairs.append(v)
                if len(pairs)==len(terminal_indices)-1:values[identity]=float(np.mean(pairs))
            return values
        for rank,(_,identity,label) in enumerate(ranked):
            proposal_args=copy.copy(args);proposal_args.priority_annotation=label
            # Ask the old constructor for a proposal; acceptance is exclusively
            # the multi-label gate below, never its learned_probe bypass.
            proposal_args.repair_policy='learned_probe'
            try:
                trial,info=original(slices,models,[m.copy() for m in transforms],terminal_indices,proposal_args)
                before=other_scores(transforms,label);after=other_scores(trial,label)
                rel=[(before[k]-after[k])/max(before[k],1e-12) for k in sorted(before.keys()&after.keys())]
                full=(info['full_score_before']-info['proposal_full_score'])/max(info['full_score_before'],1e-12)
                bounded=info['max_delta_rotation_deg']<=12 and info['max_delta_translation_ratio']<=.18
                ok=bool(bounded and full>=.02 and len(rel)>=2 and np.median(rel)>=.01 and min(rel)>=-.05 and info['proposal_lineage_score']<=info['lineage_score_before']*.85)
                row={'rank':rank,'support_identity':identity,'selected_annotation':label,'accepted':ok,'full_improvement_fraction':full,'other_label_improvements':rel,'bounded':bounded,'proposal':info}
                candidates.append(row)
                if ok:accepted.append((full,identity,trial,row))
            except Exception as error:candidates.append({'rank':rank,'support_identity':identity,'accepted':False,'error':f'{type(error).__name__}: {error}'})
        if not accepted:
            return transforms,{'operation':'anonymous_terminal_consensus','accepted':False,'reason':'no_candidate_passed_full_and_cross_label_gates','candidate_count':len(candidates),'candidates':candidates}
        _,_,trial,row=sorted(accepted,key=lambda x:(-x[0],x[1]))[0]
        return trial,{'operation':'anonymous_terminal_consensus','accepted':True,'policy':VERSION,'candidate_count':len(candidates),'candidates':candidates,'selected_support_identity':row['support_identity'],'affected_slices':','.join(str(slices[i].slice_id) for i in terminal_indices[:-1]),'selected_annotation':row['selected_annotation']}
    native.priority_lineage_terminal_rescue=rescue
    return original
