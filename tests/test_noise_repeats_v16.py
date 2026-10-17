"""Synthetic v0.16 noise-repeat evidence tests (never HST performance)."""
import json
import numpy as np
import pytest
from stellarweave.noise_repeats import REPLICATES, ADDED_MODES, validate_inputs
from stellarweave.noise_repeat_results import validate_report, aggregate


def synthetic(tmp_path,c):
    trials=[];audit=[]
    for rep in REPLICATES:
        for mode in ADDED_MODES:
            for tile in (0,5):
                audit.append({'tile_index':tile,'replicate_seed_id':rep,'mode':mode,
                    'added_noise_pixel_rms':1.5,'expected_pixel_rms':1.5,
                    'lag1_added_field_correlation':.78 if mode=='correlated_gaussian' else 0.})
                for k in range(3):
                    trials.append({'cluster':c,'tile_index':tile,
                        'class':'near_detected_peak','location_number':k,
                        'flux_over_robust_tile_rms':80.,
                        'x':30+k*12,'y':50+k*12,
                        'source_holdout_index':k,
                        'injected_flux_image_units':800.,
                        'replicate_seed_id':rep,'mode':mode,
                        'preexisting_noise_peak_at_truth':False,
                        'detected':k>0,'local_forced_abs_rel_error':.10+k*.02,
                        'center_forced_abs_rel_error':.14+k*.02})
    d={'kind':'REAL_HST_FROZEN_PLANTED_STARS_REPEATED_NOISE_SENSITIVITY',
        'cluster':c,'replicate_seeds':list(REPLICATES),
        'modes':list(ADDED_MODES),'calibrated_instrument_noise':False,
        'physical_poisson_source_noise_added':False,
        'star_positions_flux_shapes_frozen':True,
        'no_source_model_retrained_on_injected_stars':True,
        'original_sites':6,'repeated_noise_trials':len(trials),
        'repeated_noise_audits':len(audit),
        'fits_sha256':'a'*64,'v015_sha256':'b'*64,
        'trials':trials,'noise_audits':audit,
        'observation_metadata':{'header_raw':{'BUNIT':'ELECTRONS'},
                                'electron_counts_verified':False}}
    path=tmp_path/f'{c}.json'
    path.write_text(json.dumps(d))
    return path


def test_complete_noise_seed_dashboard_and_figures(tmp_path):
    files=[synthetic(tmp_path,c) for c in ('ngc6397','ngc6752','ngc0104')]
    out=tmp_path/'plots'
    d=aggregate(files,out,run_id='synthetic')
    assert d['total_original_sites']==18
    assert d['total_added_noise_trials']==180
    assert len(d['per_seed_groups'])==30
    assert d['HST_noise_physically_calibrated'] is False
    assert d['seed_ranges_are_descriptive_not_confidence_intervals']
    for name in ('recovery_across_seeds.png','forced_flux_across_seeds.png',
                 'paired_model_seed_spread.png','per_trial.csv','summary.json'):
        assert (out/name).stat().st_size>200
    replay=aggregate(sorted((out/'source_reports').glob('*.json')),tmp_path/'again')
    assert replay['total_added_noise_trials']==180


def test_missing_site_and_modified_flux_fail_closed(tmp_path):
    p=synthetic(tmp_path,'ngc6397')
    d=json.loads(p.read_text())
    d['trials'][0]['injected_flux_image_units']=900.
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='identity or flux'):
        validate_report(p)
    d['trials'][0]['injected_flux_image_units']=800.
    d['trials'].pop()
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='Missing'):
        validate_report(p)


def test_noise_audit_and_false_poisson_flag_rejected(tmp_path):
    p=synthetic(tmp_path,'ngc0104')
    d=json.loads(p.read_text())
    d['noise_audits'][0]['added_noise_pixel_rms']=9.
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='variance'):
        validate_report(p)
    d['noise_audits'][0]['added_noise_pixel_rms']=1.5
    d['physical_poisson_source_noise_added']=True
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='science provenance'):
        validate_report(p)


def test_original_v014_v015_history_checks():
    plate={'cluster':'ngc6397','grid':3,'tile_size':512,
           'anchor_index':4,'anchor_status':'HELDOUT_QC_PASS',
           'heldout_tiles_physically_nonoverlapping':True,
           'fits_sha256':'f'}
    prf={'cluster':'ngc6397','kind':'REAL_HST_SPATIAL_NATIVE_PIXEL_PRF_HOLDOUT',
         'seed':42,'no_heldout_star_used_to_fit_a_scored_template':True,
         'fits_sha256':'f','plate_sha256':'p'}
    original={'cluster':'ngc6397','kind':'REAL_HST_PLANTED_HELDOUT_STAR_RECOVERY',
              'seed':42,'models_fitted_on_injection_sources':False,
              'fits_sha256':'f','prf_sha256':'r',
              'trials':[{'tile_index':0,'class':'near_detected_peak',
                         'location_number':0,'flux_over_robust_tile_rms':80.,
                         'source_holdout_index':0,'detected':True}],
              'total_attempted':1,'total_recovered':1}
    prior={'cluster':'ngc6397','kind':'REAL_HST_FROZEN_INJECTION_CONTROLLED_NOISE_STRESS',
           'true_HST_detector_noise_model':False,
           'fits_sha256':'f','v014_baseline_sha256':'o',
           'trials':[{'tile_index':0,'class':'near_detected_peak','location_number':0,
                      'flux_over_robust_tile_rms':80.,'mode':'baseline',
                      'source_holdout_index':0,'detected':True}]}
    assert len(validate_inputs(plate,prf,original,prior,
         cluster='ngc6397',seed=42,fits_sha='f',
         plate_sha='p',prf_sha='r',original_sha='o'))==1
    prior['trials'][0]['source_holdout_index']=8
    with pytest.raises(ValueError,match='heldout shape'):
        validate_inputs(plate,prf,original,prior,
            cluster='ngc6397',seed=42,fits_sha='f',
            plate_sha='p',prf_sha='r',original_sha='o')
