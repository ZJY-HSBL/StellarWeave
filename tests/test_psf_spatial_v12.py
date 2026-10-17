"""v0.12 spatial PSF diagnostics tests; simulated fixtures are not HST."""
import json
import numpy as np
import pytest

from stellarweave.psf_spatial import psf_shape_metrics, assess_psf_grid
from stellarweave.simulation import simulated_field


def test_nonparametric_psf_shape_is_well_defined():
    y,x=np.mgrid[-7:8,-7:8]
    a=np.exp(-(x*x+y*y)/(2*1.2**2))
    m=psf_shape_metrics(a)
    assert 1.<m['r50_px']<2.5
    assert m['r80_px']>m['r50_px']
    assert 0<m['core_3x3_fraction']<1
    assert m['moment_ellipticity_inner5px']<.03
    assert m['half_light_diameter_px']==2*m['r50_px']


def test_psf_shape_rejects_invalid_inputs():
    with pytest.raises(ValueError,match='square'):
        psf_shape_metrics(np.ones((10,12)))
    with pytest.raises(ValueError,match='finite'):
        psf_shape_metrics(np.full((15,15),np.nan))


def _mosaic():
    n=192
    im=np.zeros((3*n,3*n),float)
    for i in range(9):
        tile,_=simulated_field(size=n,n_stars=50,seed=200+i,min_flux=5000)
        y,x=divmod(i,3)
        im[y*n:(y+1)*n,x*n:(x+1)*n]=tile
    origins=[(i%3*n,i//3*n) for i in range(9)]
    return im,origins


def test_spatial_grid_reports_all_regions_and_failure(tmp_path):
    im,origins=_mosaic()
    valid=np.ones(im.shape,bool)
    valid[:192,:192]=False
    result,psfs=assess_psf_grid(im,valid,origins,tile_size=192,
                               cluster='SYNTHETIC_TEST_ONLY',min_stars=3)
    assert result['tile_count']==9
    assert len(result['tiles'])==9
    assert result['failed_tiles']>=1
    assert result['tiles'][0]['status']=='INSUFFICIENT_ISOLATED_STARS'
    assert result['anchor_usable']
    assert 4 in psfs
    assert result['native_pixel_psf_not_oversampled_ePSF']
    assert result['not_calibrated_photometry']
    assert result['tiles'][4]['distance_to_anchor_total_variation']==0.
    assert result['usable_tiles']>=7
    for p in psfs.values():
        assert p.shape==(15,15)
        assert abs(p.sum()-1)<1e-10


def test_overlapping_tiles_rejected():
    im=np.ones((512,512))
    with pytest.raises(ValueError,match='Repeated'):
        assess_psf_grid(im,np.ones_like(im,bool),[(0,0)]*9,
                        tile_size=128,cluster='test')


def test_fits_integration_uses_frozen_manifest_and_writes_png(tmp_path):
    fits=pytest.importorskip('astropy.io.fits')
    from stellarweave.real_hst import sha256_file
    from stellarweave.psf_spatial import analyze_fits_psf
    im,origins=_mosaic()
    path=tmp_path/'field.fits'
    fits.PrimaryHDU(im.astype('float32')).writeto(path)
    plate={'tile_count':9,'anchor_index':4,'tile_size':192,
           'heldout_tiles_physically_nonoverlapping':True,
           'science_hdu_index':0,'cluster':'SYNTHETIC_TEST_ONLY',
           'fits_sha256':sha256_file(path),
           'tiles':[{'index':i,'tile_xy':list(xy)}
                    for i,xy in enumerate(origins)]}
    manifest=tmp_path/'plate.json'
    manifest.write_text(json.dumps(plate))
    out=tmp_path/'result'
    report=analyze_fits_psf(path,manifest,out,min_stars=3)
    assert report['usable_tiles']>=8
    assert (out/'psf_spatial.json').exists()
    assert (out/'psf_spatial.csv').exists()
    assert (out/'psf_montage.png').stat().st_size>8000
    assert (out/'tile_4_native_psf.npy').exists()
    plate['fits_sha256']='bad'
    manifest.write_text(json.dumps(plate))
    with pytest.raises(ValueError,match='SHA-256'):
        analyze_fits_psf(path,manifest,out,min_stars=3)
