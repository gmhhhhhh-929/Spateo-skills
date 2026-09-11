import json

import pandas as pd
import pytest

from spateo.preprocessing.slice_quality import (
    HighConfidencePolicy, ReviewEvidenceTier, apply_high_confidence_policy,
    write_high_confidence_outputs,
)


def policy(low=.65, high=.60):
    return HighConfidencePolicy(.129,.700,unresolved_action='keep',review_exclusion_tiers=(
        ReviewEvidenceTier('low',.129,.54,min_corroborating_domains=2,severe_domain_threshold=low,min_window_stability=1.),
        ReviewEvidenceTier('high',.54,None,min_corroborating_domains=2,severe_domain_threshold=high,min_window_stability=1.)))


def row(sid,score,maximum):
    details=[dict(window=w,score=score,detector_call='review',window_context='two_sided',
                  partial_structure_protection=False,corroborating_domains=2,maximum_domain_score=maximum) for w in [3,5,7]]
    return dict(slice_id=sid,quality_anomaly_score=score,recommendation='review',window_context='two_sided',
        partial_structure_protection=False,corroborating_domains=2,score_confidence=1.,
        density_domain_score=maximum,expression_domain_score=.5,damage_domain_score=.1,continuity_domain_score=.1,
        adaptive_window_details=json.dumps(details))


def test_joint_review_runs_both_bands_and_low_exclusion_requires_more_evidence():
    frame=pd.DataFrame([row('low_pass',.51,.67),row('low_fail',.51,.62),row('high_pass',.55,.62)])
    result=apply_high_confidence_policy(frame,policy()).set_index('slice_id')
    assert result.final_call.to_dict()=={'low_pass':'exclude','low_fail':'keep','high_pass':'exclude'}
    assert result.threshold_triage_call.eq('review').all()
    assert result.review_resolution_tier.to_dict()=={'low_pass':'low','low_fail':'low','high_pass':'high'}
    with pytest.raises(ValueError,match='strictly stronger'):
        apply_high_confidence_policy(frame,policy(.6,.65))


def test_joint_review_keeps_protection_and_each_window_veto():
    protected=row('protected',.51,.67);protected['partial_structure_protection']=True
    unstable=row('unstable',.55,.62);details=json.loads(unstable['adaptive_window_details'])
    details[-1]['maximum_domain_score']=.59;unstable['adaptive_window_details']=json.dumps(details)
    result=apply_high_confidence_policy(pd.DataFrame([protected,unstable]),policy())
    assert result.final_call.eq('keep').all()


def test_experimental_application_cannot_claim_certified_calls(tmp_path):
    result=write_high_confidence_outputs(pd.DataFrame([row('low',.51,.67)]),policy(),tmp_path,
                                          application_scope='experimental_policy')
    audit=pd.read_csv(result['audit']);summary=json.loads(open(result['summary']).read())
    assert audit.final_call.tolist()==['exclude']
    assert audit.certified_call.isna().all()
    assert summary['certified']==0 and not summary['new_input_independently_validated']
    assert summary['policy_validation_scope']=='metric_stress_tests_only'
