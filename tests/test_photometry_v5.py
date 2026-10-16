"""Science-integration tests for v0.5 photutils 3.x PSF fitting."""
import csv
import json
import numpy as np
import pytest

pytest.importorskip('astropy')
pytest.importorskip('photutils')

from stellarweave.artificial_stars import gaussian_prf, inject_stars
from stellarweave.photometry import measure_sources
from stellarweave.benchmark import compare_photometry,save_benchmark


def test_pixel_integrated_prf_recovery():
    rng=np.random.default_rng(11)
    base=80+rng.normal(0,1.5,(100,100))
    xy=(48.3,43.7)
    image=inject_stars(base,[[xy[0],xy[1],16000]],gaussian_prf())
    for algorithm in ('psf','iterative'):
        rows=measure_sources(image,algorithm,fwhm=2.6,threshold=6,fit_shape=11)
        assert rows.shape[1]==4
        assert len(rows)>0
        d=np.hypot(rows[:,0]-xy[0],rows[:,1]-xy[1])
        selected=rows[np.argmin(d)]
        assert d.min()<1.0
        assert 0.65 < selected[2]/16000 < 1.35


def test_invalid_configuration():
    image=np.ones((64,64))
    with pytest.raises(ValueError):
        measure_sources(image,'nonexistent')
    with pytest.raises(ValueError):
        measure_sources(image,'psf',fit_shape=8)
    with pytest.raises(ValueError):
        measure_sources(image,'psf',valid_mask=np.ones((10,10),bool))


def test_common_injected_stars_for_each_method(tmp_path):
    rng=np.random.default_rng(91)
    base=70+rng.normal(0,3,(128,128))
    # Strong pre-existing field objects, so baseline exclusion is exercised.
    base=inject_stars(base,[[32.2,40.7,8500],[70.8,83.1,12000]],
                      gaussian_prf())
    result=compare_photometry(
        base,gaussian_prf(),[600.,2400.],
        methods=('dao','psf','iterative'),seed=12,trials=1,
        per_flux=2,fwhm=2.6,threshold=5.,fit_shape=9)
    rows=result['records']
    assert len(rows)==12
    for injection_id in range(4):
        selected=[r for r in rows if r['injection_id']==injection_id]
        assert len(selected)==3
        assert len({(r['x_px'],r['y_px'],r['injected_flux']) for r in selected})==1
    for group in result['groups']:
        if group['injected']:
            assert 0 <= group['completeness'] <= 1
    save_benchmark(result,tmp_path)
    assert (tmp_path/'source_results.csv').exists()
    assert (tmp_path/'method_comparison.png').is_file()
    data=json.loads((tmp_path/'benchmark.json').read_text())
    assert data['metadata']['same_injected_images_all_methods']
    assert len(list(csv.DictReader((tmp_path/'source_results.csv').open())))==12


def test_cli_guard_real_fits():
    from stellarweave.benchmark_cli import main
    with pytest.raises(SystemExit):
        main(['--fits','missing.fits'])
