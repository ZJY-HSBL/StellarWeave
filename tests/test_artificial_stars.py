"""Artificial-star regression and reproducibility tests."""
import json
import numpy as np
import pytest
from scipy.spatial.distance import pdist
from stellarweave.artificial_stars import (
    gaussian_prf, normalized_psf, inject_stars, sample_positions,
    recover_injected, wilson_interval, experiment)
from stellarweave.simulation import simulated_field
from stellarweave.detection import scipy_peak_detector


def test_psf_integral():
    k=gaussian_prf(2.6,17)
    assert k.shape==(17,17) and np.isclose(k.sum(),1)
    assert k[8,8]>k[7,8]


def test_subpixel_flux_conservation_and_immutable_original():
    original=np.full((50,50),7.)
    before=original.copy()
    out=inject_stars(original,[[24.28,24.71,1200],[12.15,13.25,500]],gaussian_prf())
    assert np.isclose((out-original).sum(),1700,atol=1e-7)
    np.testing.assert_array_equal(original,before)


def test_bad_footprints():
    a=np.ones((60,60))
    m=np.ones_like(a,bool)
    m[21,21]=False
    with pytest.raises(ValueError,match='masked'):
        inject_stars(a,[[21,21,100]],gaussian_prf(),valid_mask=m)
    with pytest.raises(ValueError,match='border'):
        inject_stars(a,[[1,1,100]],gaussian_prf())


def test_invalid_psfs():
    with pytest.raises(ValueError):
        normalized_psf(-np.ones((7,7)))
    with pytest.raises(ValueError):
        normalized_psf(np.ones((6,6)))


def test_position_sampling():
    m=np.ones((160,160),bool)
    x=sample_positions(m.shape,m,15,8,rng=np.random.default_rng(8),min_separation=15)
    y=sample_positions(m.shape,m,15,8,rng=np.random.default_rng(8),min_separation=15)
    np.testing.assert_allclose(x,y)
    assert pdist(x).min()>=15


def test_recovery_not_credited_for_preexisting_star():
    got,dist,extra=recover_injected([[10,10],[30,30]],[[10,10]],
                                    [[10,10],[30.4,30]],radius=1.5)
    assert got.tolist()==[False,True]
    assert extra==0


def test_wilson_interval():
    assert wilson_interval(0,0)==[None,None]
    assert 0 < wilson_interval(0,10)[1] < .5
    assert .5 < wilson_interval(10,10)[0] < 1


def test_reproducible_experiment_and_added_flux():
    base,_=simulated_field(size=160,n_stars=10,seed=8)
    det=lambda x: scipy_peak_detector(x,nsigma=4)
    kw=dict(detector=det,trials=2,stars_per_flux=4,seed=3,fwhm_px=2.6,
            injection_min_separation=16)
    a=experiment(base,gaussian_prf(),[800,3500],**kw)
    b=experiment(base,gaussian_prf(),[800,3500],**kw)
    assert a['rows']==b['rows']
    assert len(a['rows'])==16
    for row in a['trials']:
        assert np.isclose(row['injected_flux_total'],
                          row['image_flux_difference'],atol=1e-7)
    assert [g['injected'] for g in a['groups']
            if g['crowding']=='all']==[8,8]


def test_cli_smoke(tmp_path):
    from stellarweave.injection_cli import main
    output=tmp_path/'smoke'
    main(['--size','128','--trials','1','--per-flux','2',
          '--flux-levels','800','2500','--seed','7','--out',str(output)])
    result=json.loads((output/'manifest.json').read_text())
    assert result['mode']=='SYNTHETIC_GAUSSIAN_PRF'
    assert (output/'completeness.png').is_file()
    assert len((output/'injections.csv').read_text().splitlines())==5
    assert len(result['groups'])==8


def test_real_fits_requires_empirical_psf_and_unit():
    from stellarweave.injection_cli import main
    with pytest.raises(SystemExit):
        main(['--fits','missing.fits'])
