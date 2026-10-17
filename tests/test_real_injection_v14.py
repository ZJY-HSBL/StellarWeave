"""Synthetic fixture tests; NOT real-HST injection performance claims."""
import numpy as np
import pytest

from stellarweave.real_injection import (
    predeclare_positions,forced_amplitude,inject_and_measure,run_injection_tile,
    NEIGHBOR_CLASSES,FLUX_RMS_MULTIPLIERS)
from stellarweave.spatial_prf import gaussian_reference


def synthetic_image(n=256,seed=5):
    rng=np.random.default_rng(seed)
    img=rng.normal(30,1,(n,n))
    psf=gaussian_reference()
    for x,y,f in ((45,45,700),(85,50,800),(120,120,650),(160,180,700),
                  (200,80,750),(55,210,800)):
        img[y-7:y+8,x-7:x+8]+=f*psf
    return img,psf


def test_predeclared_crowding_geometry_and_reproducibility():
    img,_=synthetic_image()
    peaks=np.array([[45.,45.],[85,50],[120,120],[160,180],[200,80],[55,210]])
    good=np.ones_like(img,bool)
    a=predeclare_positions(img,good,peaks,seed=101,per_class=2)
    assert a==predeclare_positions(img,good,peaks,seed=101,per_class=2)
    assert set(a)==set(NEIGHBOR_CLASSES)
    assert len(a['near_detected_peak'])==2
    assert len(a['relatively_isolated'])==2
    for x,y in a['near_detected_peak']:
        d=np.linalg.norm(peaks-np.array([x,y]),axis=1)
        assert 5<=d.min()<11
    for x,y in a['relatively_isolated']:
        d=np.linalg.norm(peaks-np.array([x,y]),axis=1)
        assert d.min()>=16


def test_known_flux_forced_fit_and_detector():
    img,psf=synthetic_image()
    from stellarweave.detection import scipy_peak_detector
    peaks=scipy_peak_detector(img,fwhm=2.6,nsigma=5)
    assert len(peaks)>3
    p=(130,80)
    result=inject_and_measure(img,np.ones_like(img,bool),p,psf,
                             1000.,peaks,psf,psf)
    assert result['detected'] is True
    assert result['detector_peak_position_error_px']<=2.5
    assert result['forced_local_flux']==pytest.approx(result['forced_center_flux'])
    assert result['forcing_uses_known_position'] is True
    with pytest.raises(ValueError,match='already has'):
        inject_and_measure(img,np.ones_like(img,bool),(45,45),psf,
                           1000.,peaks,psf,psf)


def test_tile_trials_reuse_training_templates_and_holdout_shapes():
    img,psf=synthetic_image(seed=12)
    rng=np.random.default_rng(22)
    stamps=np.repeat(psf[None],40,axis=0)
    stamps=np.maximum(stamps+rng.normal(0,.0001,stamps.shape),0)
    stamps/=stamps.sum(axis=(1,2),keepdims=True)
    xy=np.column_stack([np.linspace(20,210,40),np.linspace(20,210,40)])
    rows,meta=run_injection_tile(img,np.ones_like(img,bool),
                                stamps,xy,psf,tile_index=0)
    assert meta['status']=='EVALUATED'
    assert meta['independent_injection_source_stamps']==10
    assert len(rows)<=len(NEIGHBOR_CLASSES)*len(FLUX_RMS_MULTIPLIERS)*2
    assert len(rows)>0
    assert all(row['source_holdout_index']<10 for row in rows)
    assert all(row['injected_flux_image_units']>0 for row in rows)
    assert len({row['flux_over_robust_tile_rms'] for row in rows})==3
    again,_=run_injection_tile(img,np.ones_like(img,bool),
                               stamps,xy,psf,tile_index=0)
    assert rows==again


def test_sparse_tile_retains_failure():
    img,psf=synthetic_image()
    stamps=np.repeat(psf[None],5,axis=0)
    coords=np.arange(10).reshape(5,2)
    rows,meta=run_injection_tile(img,np.ones_like(img,bool),stamps,coords,
                                psf,tile_index=0)
    assert not rows
    assert meta['status']=='INSUFFICIENT_PRF_STARS'
