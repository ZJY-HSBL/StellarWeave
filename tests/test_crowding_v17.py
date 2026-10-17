"""Synthetic-only v0.17 tests. Never interpret fixtures as HST results."""
import json
import numpy as np
import pytest
from stellarweave.crowding_injections import (
    select_distance_bands, tile_experiment, BANDS, PER_TILE_BAND)
from stellarweave.crowding_results import aggregate, validate
from stellarweave.spatial_prf import gaussian_reference


def test_predeclared_distance_bands_are_reproducible_and_disjoint():
    peaks=np.array([(x,y) for x in range(45,270,25)
                     for y in range(45,270,25)],dtype=float)
    valid=np.ones((320,320),bool)
    sites=select_distance_bands(valid.shape,valid,peaks,seed=53)
    assert sites==select_distance_bands(valid.shape,valid,peaks,seed=53)
    selected=[]
    for band,(low,high) in BANDS.items():
        assert len(sites[band])==PER_TILE_BAND
        for s in sites[band]:
            assert low<=s['nearest_original_detected_peak_px']<high
            selected.append((s['x'],s['y']))
    assert len(set(selected))==len(selected)
    for i,(x,y) in enumerate(selected):
        for x2,y2 in selected[:i]:
            assert np.hypot(x-x2,y-y2)>=12.


def test_missing_band_is_retained():
    valid=np.ones((100,100),bool)
    peaks=np.array([[50.,50.]])
    sites=select_distance_bands(valid.shape,valid,peaks,seed=8,per_band=1)
    assert set(sites)==set(BANDS)
    assert all(isinstance(sites[k],list) for k in BANDS)


def _synthetic_cluster_report(tmp_path,cluster):
    trials=[]
    tiles=[]
    for idx in (0,1,2,3,5,6,7,8):
        tiles.append({'tile':idx,'status':'EVALUATED','bins':{x:1 for x in BANDS}})
        for bi,band in enumerate(BANDS):
            for flux in (80.,200.,500.):
                trials.append({
                    'cluster':cluster,'tile_index':idx,'band':band,
                    'site_index':0,'x':45+25*bi,'y':65+idx,
                    'nearest_original_detected_peak_px':(BANDS[band][0]+BANDS[band][1])/2,
                    'flux_over_tile_robust_pixel_rms':flux,
                    'detected':flux!=80.,
                    'forced_local_abs_rel_error':.08,'forced_center_abs_rel_error':.1,
                    'detector_peak_position_error_px':0.5 if flux!=80. else None})
    report={
        'cluster':cluster,
        'kind':'REAL_HST_EXPANDED_PREDECLARED_CROWDING_INJECTIONS',
        'all_injection_sites_chosen_before_outcomes':True,
        'full_test_set_retained':True,
        'all_models_fitted_on_noninjected_training_stars_only':True,
        'stellar_population_completeness_estimated':False,
        'injected_poisson_noise':False,
        'mosaic_noise_physically_calibrated':False,
        'forced_flux_uses_true_injection_position':True,
        'bands':{x:list(y) for x,y in BANDS.items()},
        'distinct_selected_injection_positions':24,
        'injection_attempts':len(trials),
        'recovered':sum(x['detected'] for x in trials),
        'fits_sha256':'0'*64,'prf_sha256':'1'*64,
        'lineage_audit':{'underlying_exposure_ids_verified':False},
        'tiles':tiles,'trials':trials}
    path=tmp_path/f'{cluster}.json'
    path.write_text(json.dumps(report))
    return path


def test_real_evidence_plot_structure_and_reproduction(tmp_path):
    paths=[_synthetic_cluster_report(tmp_path,c)
           for c in ('ngc6397','ngc6752','ngc0104')]
    out=tmp_path/'plots'
    report=aggregate(paths,out,run_id='synthetic-test-only')
    assert report['total_unique_pixel_sites']==72
    assert report['total_trials']==216
    assert report['total_recovered']==144
    assert len(report['groups'])==27
    for name in ('summary.json','per_trial.csv','groups.csv','tile_sites.csv',
                 'recovery_crowding_distance.png','paired_flux_crowding.png',
                 'injection_site_coverage.png'):
        assert (out/name).stat().st_size>200
    assert aggregate(sorted((out/'source_reports').glob('*.json')),
                     tmp_path/'replay')['total_trials']==216


def test_cherry_picking_or_moved_site_is_rejected(tmp_path):
    p=_synthetic_cluster_report(tmp_path,'ngc0104')
    d=json.loads(p.read_text())
    d['trials'].pop()
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='Incomplete'):
        validate(p)
    d=json.loads(_synthetic_cluster_report(tmp_path,'ngc0104').read_text())
    d['trials'][0]['x']=999
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='Position moved'):
        validate(p)


def test_injection_uses_independent_heldout_shapes_and_no_new_training():
    rng=np.random.default_rng(13)
    im=rng.normal(20,1,(320,320))
    psf=gaussian_reference()
    for x in range(45,270,25):
        for y in range(45,270,25):
            im[y-7:y+8,x-7:x+8]+=900*psf
    stamps=np.repeat(psf[None],48,axis=0).copy()
    xy=np.column_stack((np.arange(48,dtype=float),
                        np.arange(48,dtype=float)*2))
    rows,meta=tile_experiment(im,np.ones_like(im,bool),stamps,xy,psf,
                              cluster='ngc0104',tile_index=0)
    assert meta['status']=='EVALUATED'
    assert meta['heldout_stars']>=3
    assert len(rows)>0
    assert len(rows)==3*sum(meta['bins'].values())
    assert all(0<=row['heldout_stamp_index']<meta['heldout_stars'] for row in rows)
    assert all(row['forcing_uses_known_position'] for row in rows)
