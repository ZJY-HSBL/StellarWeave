"""Synthetic fixtures only; never report them as real-HST performance."""
import json

import numpy as np
import pytest

from stellarweave.noise_stress import (
    added_noise, noise_horizontal_correlation,
    poisson_electron_fluctuation, audit_science_noise_metadata,
    one_augmented_trial, MODES,
)
from stellarweave.spatial_prf import gaussian_reference
from stellarweave.noise_results import aggregate, _validate


def test_equal_pixel_rms_but_distinct_spatial_covariance():
    shape=(256,256)
    mask=np.ones(shape,bool)
    rng=16.5
    a=added_noise(shape,mask,mode='iid_gaussian',tile_rms=rng,seed=77)
    b=added_noise(shape,mask,mode='correlated_gaussian',tile_rms=rng,seed=77)
    np.testing.assert_allclose(np.std(a),rng*.5,atol=1e-11)
    np.testing.assert_allclose(np.std(b),rng*.5,atol=1e-11)
    assert abs(noise_horizontal_correlation(a,mask))<.04
    assert noise_horizontal_correlation(b,mask)>.5
    np.testing.assert_array_equal(
        a,added_noise(shape,mask,mode='iid_gaussian',tile_rms=rng,seed=77))
    assert not np.any(added_noise(shape,mask,mode='baseline',
                                  tile_rms=rng,seed=77))
    with pytest.raises(ValueError,match='Unknown'):
        added_noise(shape,mask,mode='poisson_unverified',tile_rms=rng,seed=77)


def test_no_noise_where_invalid_and_physical_poisson_fail_closed():
    mask=np.ones((100,110),bool)
    mask[:20,:]=False
    a=added_noise(mask.shape,mask,mode='correlated_gaussian',
                  tile_rms=3,seed=5)
    assert not a[~mask].any()
    assert np.std(a[mask])==pytest.approx(1.5,rel=1e-10)
    with pytest.raises(ValueError):
        poisson_electron_fluctuation([-1,1],np.random.default_rng(1))
    rng=np.random.default_rng(10)
    deviations=poisson_electron_fluctuation(np.full(100000,13.),rng)
    assert np.abs(deviations.mean())<.04
    assert np.var(deviations)==pytest.approx(13,rel=.03)
    audit=audit_science_noise_metadata(
        {'BUNIT':'ELECTRONS/S','EXPTIME':2000},err_present=True,dq_present=False)
    assert audit['electron_counts_verified'] is False
    assert audit['poisson_shot_noise_applied_to_real_HST'] is False


def test_noise_trial_never_calls_existing_peak_recovery():
    from stellarweave.detection import scipy_peak_detector
    psf=gaussian_reference()
    img=np.random.default_rng(2).normal(10,1,(160,160))
    img[33-7:33+8,34-7:34+8]+=1000*psf
    before=scipy_peak_detector(img,fwhm=2.6,nsigma=5.)
    assert len(before)>0
    trial={'x':34,'y':33,'injected_flux_image_units':300}
    score=one_augmented_trial(img,np.ones_like(img,bool),trial,psf,
                              psf,psf,before)
    assert score['preexisting_noise_peak_at_truth'] is True
    assert score['detected'] is False
    assert np.isfinite(score['local_forced_abs_rel_error'])


def _fake(path,cluster):
    trials=[]
    for tile in (0,1):
        for cls in ('near_detected_peak','relatively_isolated'):
            for m in (80.,200.):
                for k in range(2):
                    for mode in MODES:
                        trials.append({
                            'cluster':cluster,'tile_index':tile,'class':cls,
                            'location_number':k,'flux_over_robust_tile_rms':m,
                            'x':30+k*10,'y':40+k*15,
                            'source_holdout_index':k,
                            'injected_flux_image_units':1000*m,
                            'mode':mode,'detected':k==1,
                            'preexisting_noise_peak_at_truth':False,
                            'local_forced_abs_rel_error':.12,
                            'center_forced_abs_rel_error':.15})
    report={'cluster':cluster,
            'kind':'REAL_HST_FROZEN_INJECTION_CONTROLLED_NOISE_STRESS',
            'modes':list(MODES),'true_HST_detector_noise_model':False,
            'no_heldout_source_in_training':True,
            'photon_shot_noise_applied_to_real_HST':False,
            'all_original_injection_positions_and_shapes_frozen':True,
            'expected_original_trials':len(trials)//3,'trials':trials,
            'noise_audits':[
                {'tile_index':tile,'mode':mode,
                 'target_added_pixel_rms':0.0 if mode=='baseline' else 1.0,
                 'actual_added_pixel_rms':0.0 if mode=='baseline' else 1.0,
                 'added_field_lag1_correlation':.8 if mode=='correlated_gaussian'
                         else .0 if mode=='iid_gaussian' else None,
                 'baseline_peak_count_after_augmentation':15}
                for tile in (0,1) for mode in MODES],
            'fits_sha256':'f'*64,'v014_baseline_sha256':'b'*64}
    path.write_text(json.dumps(report))
    return path


def test_full_threemode_frozen_figure_reproduction(tmp_path):
    names=('ngc6397','ngc6752','ngc0104')
    paths=[_fake(tmp_path/f'{c}.json',c) for c in names]
    d=aggregate(paths,tmp_path/'figures',run_id='synthetic-fixture')
    assert d['total_mode_trials']==96*3//2
    assert d['total_original_sites']==48
    assert len(d['groups'])==9
    assert not d['calibrated_HST_detector_noise']
    assert len(list((tmp_path/'figures'/'source_reports').glob('*.json')))==3
    for p in ('recovery_noise_sensitivity.png','local_forced_flux_noise.png',
              'paired_prf_noise_robustness.png','per_trial.csv','groups.csv'):
        assert (tmp_path/'figures'/p).stat().st_size>100
    replay=aggregate(sorted((tmp_path/'figures'/'source_reports').glob('*.json')),
                     tmp_path/'replay')
    assert replay['total_mode_trials']==d['total_mode_trials']


def test_missing_control_trial_rejected(tmp_path):
    p=_fake(tmp_path/'bad.json','ngc6397')
    d=json.loads(p.read_text())
    d['trials'].pop()
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='Incomplete'):
        _validate(p)
