"""PRF reporting tests with SYNTHETIC fixtures only."""
import csv
import json
import numpy as np
import pytest
from stellarweave.prf_visualization import render_prf_summary,validate,CLUSTERS


def _fixture(root,name):
    root.mkdir(parents=True)
    tiles=[]
    stars=[]
    for idx in range(9):
        tiles.append({'index':idx,'status':'EVALUATED','train_stars':12,
                      'heldout_stars':3,'median_local_l1':.09,
                      'median_anchor_l1':.12 if idx!=4 else .09,
                      'median_gaussian_l1':.15,
                      'median_paired_gain_pct':25. if idx!=4 else 0.})
        for j in range(3):
            stars.append({'tile_index':idx,'holdout_index':j,'tile_x':j*2.,
                          'tile_y':j*3.,'local_l1':.09,
                          'anchor_l1':.12,'gaussian_l1':.15,
                          'paired_gain_pct':25.})
    report={'cluster':name,'kind':'REAL_HST_SPATIAL_NATIVE_PIXEL_PRF_HOLDOUT',
            'science_calibrated':False,'oversampled_instrument_epsf':False,
            'calibrated_photometric_accuracy':False,
            'no_heldout_star_used_to_fit_a_scored_template':True,
            'tiles':tiles,'evaluated_heldout_tiles':8,
            'total_heldout_stars':len(stars),
            'heldout_star_scores':stars,'local_prf_better_tiles':8,
            'median_tile_paired_gain_pct':25.,
            'bootstrap_95ci_median_tile_gain_pct':[20,30],
            'fits_sha256':'0'*64,'plate_sha256':'1'*64}
    p=root/'spatial_prf.json'
    p.write_text(json.dumps(report))
    with (root/'heldout_stars.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(stars[0]))
        w.writeheader();w.writerows(stars)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots()
    ax.imshow(np.ones((15,15)))
    fig.savefig(root/'spatial_prf_templates.png')
    plt.close(fig)
    return p


def test_aggregate_and_provenance(tmp_path):
    paths=[_fixture(tmp_path/c,c) for c in CLUSTERS]
    out=tmp_path/'report'
    d=render_prf_summary(paths,out,run_id='synthetic-fixture',commit='abc')
    assert d['heldout_tile_slots']==24
    assert d['evaluated_heldout_tiles']==24
    assert len(d['clusters'])==3
    assert d['calibrated_photometry'] is False
    for name in ('prf_tile_gains.png','prf_heldout_comparison.png',
                 'summary.json','per_tile.csv','per_star.csv'):
        assert (out/name).stat().st_size>100
    assert len(list((out/'source_reports').glob('*.json')))==3
    # README reproduction command must work from the committed snapshot.
    archived=sorted((out/'source_reports').glob('*.json'))
    replay=render_prf_summary(archived,tmp_path/'reproduced')
    assert replay['evaluated_heldout_tiles']==24
    assert (tmp_path/'reproduced'/'prf_tile_gains.png').stat().st_size>100


def test_missing_and_false_claim_rejected(tmp_path):
    p=_fixture(tmp_path/'a','ngc6397')
    with pytest.raises(ValueError,match='Expected three'):
        render_prf_summary([p],tmp_path/'out')
    d=json.loads(p.read_text())
    d['calibrated_photometric_accuracy']=True
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='scientific claim'):
        validate(p)
