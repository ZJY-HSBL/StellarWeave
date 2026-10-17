"""Synthetic-fixture tests; these are NOT telescope performance measurements."""
import json

import numpy as np
import pytest

from stellarweave.prf_stratification import measure_tile_proxies
from stellarweave.prf_stratification_summary import (
    validate_frozen_reports, summarize_groups, render_diagnostics,
)
from stellarweave.real_hst import sha256_file


def test_real_image_proxy_algorithm_on_synthetic_star_field():
    rng = np.random.default_rng(7)
    image = rng.normal(10, 1, (260, 260))
    y, x = np.mgrid[-8:9, -8:9]
    for cx,cy in [(30,30),(65,35),(100,100),(140,160),(190,210)]:
        image[cy-8:cy+9,cx-8:cx+9] += (
            700*np.exp(-(x*x+y*y)/(2*1.2**2)))
    holdout = [{'tile_x':30.,'tile_y':30.,'paired_gain_pct':-2.},
               {'tile_x':65.,'tile_y':35.,'paired_gain_pct':3.}]
    rows, info = measure_tile_proxies(
        image,np.ones_like(image,bool),holdout)
    assert len(rows)==2
    assert info['detected_peaks']>=5
    assert rows[0]['nearest_detected_peak_px']>30
    assert rows[0]['peak_over_tile_rms']>50
    assert rows[0]['annulus_mad_over_tile_rms'] is not None
    assert rows[0]['paired_gain_pct']==-2.


def test_proxy_fail_closed_on_missing_frozen_star():
    rng=np.random.default_rng(3)
    image=rng.normal(0,1,(100,100))
    image[30,30]+=1000
    image[60,60]+=1000
    with pytest.raises(ValueError,match='not present'):
        measure_tile_proxies(image,np.ones_like(image,bool),
                             [{'tile_x':50,'tile_y':50}])


def _fixture(root,cluster):
    original=[]
    modified=[]
    for tile in (0,1,2,3,5,6,7,8):
        for i in range(4):
            old={'tile_index':tile,'holdout_index':i,
                 'tile_x':float(20+i*30),'tile_y':float(20+tile*20),
                 'local_l1':float(.2+i*.01),
                 'anchor_l1':float(.25+i*.01),
                 'gaussian_l1':float(.3+i*.01),
                 'paired_gain_pct':float(100*.05/(.25+i*.01))}
            original.append(old)
            modified.append(dict(old,
                peak_over_tile_rms=float(10+i*10),
                nearest_detected_peak_px=float(8+i*3),
                neighbors_within_25px=int(i+tile),
                annulus_mad_over_tile_rms=float(1+i*.5)))
    p=root/f'{cluster}_original.json'
    p.write_text(json.dumps({'cluster':cluster,'fits_sha256':'0'*64,
                             'plate_sha256':'1'*64,
                             'kind':'REAL_HST_SPATIAL_NATIVE_PIXEL_PRF_HOLDOUT',
                             'heldout_star_scores':original}))
    strat=root/f'{cluster}_strat.json'
    strat.write_text(json.dumps({
        'cluster':cluster,'kind':'FROZEN_HST_PRF_CROWDING_BRIGHTNESS_EXPLORATORY',
        'models_retrained':False,'stars_removed_for_stratification':False,
        'frozen_prf_scores_unchanged':True,'analysis_is_posthoc_exploratory':True,
        'science_calibrated':False,'tile_count':8,
        'original_prf_report_sha256':sha256_file(p),
        'original_prf_kind':'REAL_HST_SPATIAL_NATIVE_PIXEL_PRF_HOLDOUT',
        'fits_sha256':'0'*64,'plate_sha256':'1'*64,
        'all_original_outer_heldout_stars':len(modified),
        'tiles':[{'tile_index':t} for t in (0,1,2,3,5,6,7,8)],
        'heldout_star_scores_with_proxies':modified}))
    return strat,p


def test_frozen_holdout_provenance_and_all_strata(tmp_path):
    files=[_fixture(tmp_path,c) for c in ('ngc6397','ngc6752','ngc0104')]
    strat,prf=zip(*files)
    reports,_=validate_frozen_reports(strat,prf)
    groups,thresholds=summarize_groups(reports)
    assert len(groups)==24
    assert len(thresholds)==12
    assert sum(g['n'] for g in groups)==3*4*32
    out=tmp_path/'figures'
    result=render_diagnostics(strat,prf,out,run_id='synthetic-fixture')
    assert result['total_original_outer_heldout_stars']==96
    assert result['all_frozen_test_stars_retained'] is True
    assert result['causal_inference_supported'] is False
    for name in ('summary.json','groups.csv','per_star.csv','thresholds.csv',
                 'stratified_prf_gains.png','ngc0104_failure_map.png',
                 'ngc0104_proxy_scatter.png'):
        assert (out/name).stat().st_size>200
    assert len(list((out/'source_reports').glob('*.json')))==3
    assert len(list((out/'original_prf_reports').glob('*.json')))==3
    again=render_diagnostics(
        sorted((out/'source_reports').glob('*.json')),
        sorted((out/'original_prf_reports').glob('*.json')),
        tmp_path/'replay')
    assert again['total_original_outer_heldout_stars']==96


def test_tampering_and_cherry_picking_rejected(tmp_path):
    files=[_fixture(tmp_path,c) for c in ('ngc6397','ngc6752','ngc0104')]
    strat,prf=zip(*files)
    d=json.loads(strat[0].read_text())
    d['heldout_star_scores_with_proxies'][0]['local_l1']+=.01
    strat[0].write_text(json.dumps(d))
    with pytest.raises(ValueError,match='score modified'):
        validate_frozen_reports(strat,prf)
    d['heldout_star_scores_with_proxies'].pop()
    strat[0].write_text(json.dumps(d))
    with pytest.raises(ValueError,match='count mismatch'):
        validate_frozen_reports(strat,prf)
