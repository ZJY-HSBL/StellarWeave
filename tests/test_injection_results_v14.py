"""Synthetic report fixtures: verify charts and full trial retention, not science."""
import json
import pytest
from stellarweave.injection_results import aggregate,validate
from stellarweave.real_injection import CLUSTERS,NEIGHBOR_CLASSES,FLUX_RMS_MULTIPLIERS


def fake(path,cluster):
    trials=[]
    groups=[]
    for cls in NEIGHBOR_CLASSES:
        for mult in FLUX_RMS_MULTIPLIERS:
            for i in range(3):
                trials.append({'cluster':cluster,'tile_index':0,
                   'class':cls,'flux_over_robust_tile_rms':mult,
                   'detected':i!=0,
                   'detector_peak_position_error_px':None if i==0 else .4,
                   'forced_local_abs_rel_error':.1+i*.02,
                   'forced_center_abs_rel_error':.14+i*.02})
            groups.append({'class':cls,'flux_over_tile_rms':mult})
    report={
        'cluster':cluster,'kind':'REAL_HST_PLANTED_HELDOUT_STAR_RECOVERY',
        'science_calibrated':False,'real_background_images':True,
        'models_fitted_on_injection_sources':False,
        'matched_source_position_known_for_forced_flux':True,
        'fits_sha256':'0'*64,'prf_sha256':'1'*64,
        'total_attempted':len(trials),'total_recovered':sum(t['detected'] for t in trials),
        'trials':trials,'groups':groups,
        'tiles':[{'tile':i,'status':'EVALUATED'} for i in (0,1,2,3,5,6,7,8)]
    }
    path.write_text(json.dumps(report))
    return path


def test_render_frozen_injection_evidence(tmp_path):
    files=[fake(tmp_path/f'{c}.json',c) for c in CLUSTERS]
    out=tmp_path/'plots'
    d=aggregate(files,out,run_id='synthetic-test')
    assert d['total_trials']==54
    assert d['total_recovered']==36
    assert len(d['groups'])==18
    assert d['source_noise_is_real_HST'] is True  # checked as metadata flag only
    assert d['stellar_catalog_completeness_proven'] is False
    for name in ('recovery_by_flux_crowding.png','forced_flux_error.png',
                 'paired_forced_flux_difference.png','per_trial.csv',
                 'groups.csv','summary.json'):
        assert (out/name).stat().st_size>100
    assert len(list((out/'source_reports').glob('*.json')))==3
    assert aggregate(sorted((out/'source_reports').glob('*.json')),
                     tmp_path/'replay')['total_trials']==54


def test_missing_injection_trials_rejected(tmp_path):
    p=fake(tmp_path/'one.json','ngc6397')
    d=json.loads(p.read_text())
    d['trials'].pop()
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='Missing'):
        validate(p)
