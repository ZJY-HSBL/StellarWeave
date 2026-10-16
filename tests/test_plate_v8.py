"""v0.8: independent held-out validation for WCS-guided wide-shift fits."""
import numpy as np
import pytest
from stellarweave.astrometry import apply_affine
from stellarweave.robust_plate import (
    unique_links,translation_hypotheses,plate_diagnostics)


def test_global_unique_link_and_no_double_assignment():
    ref=[[10.,10.],[10.2,10.],[40.,40.]]
    detected=[[10.1,10.],[40.,40.]]
    pairs=unique_links(ref,detected,radius=1.)
    assert len(pairs)==2
    assert len(set(pairs[:,1]))==2


def test_wide_translation_votes_with_outliers():
    rng=np.random.default_rng(50)
    refs=rng.uniform(30,460,(110,2))
    shift=np.array([59.,-43.])
    observed=refs+shift+rng.normal(0,.15,(110,2))
    observed=np.vstack((observed,rng.uniform(0,512,(90,2))))
    hypotheses=translation_hypotheses(refs,observed,max_shift_px=100,
                                      bin_width=2,max_hypotheses=10)
    assert hypotheses
    assert any(np.linalg.norm(np.asarray(h['shift'])-shift)<2.5 for h in hypotheses)


def test_wide_shift_affine_with_independent_holdout():
    rng=np.random.default_rng(7)
    source=rng.uniform(20,430,(350,2))
    truth=np.array([[1.0005,-.0007,57.],[.0005,1.0003,-41.]])
    detections=apply_affine(source,truth)+rng.normal(0,.12,(350,2))
    detections=np.vstack((detections,rng.uniform(0,512,(150,2))))
    report=plate_diagnostics(source,detections,shape=(512,512),
                max_shift_px=100,train_limit=190,holdout_limit=150,
                min_train_matches=25,min_holdout_matches=20)
    assert report['status']=='HELDOUT_QC_PASS',report
    assert report['holdout_matches']>80
    assert report['holdout_match_fraction']>.5
    assert report['holdout_vs_null_ratio']>4
    assert report['science_calibrated'] is False
    np.testing.assert_allclose(report['training_affine']['matrix'],truth,atol=.6)


def test_bad_catalog_or_missing_peaks_does_not_pass():
    rng=np.random.default_rng(22)
    ref=rng.uniform(20,490,(350,2))
    det=rng.uniform(700,1000,(170,2))
    report=plate_diagnostics(ref,det,max_shift_px=120)
    assert report['status'] in ('NO_RELIABLE_TRAIN_FIT','HELDOUT_QC_REJECTED')
    assert not report['science_calibrated']


def test_excessive_search_rejected():
    with pytest.raises(ValueError):
        translation_hypotheses([[1,2]],[[1,2]],max_shift_px=-1)


def test_spatial_integration_records_holdout(tmp_path):
    pytest.importorskip('astropy')
    from astropy.io import fits
    from astropy.wcs import WCS
    from stellarweave.simulation import simulated_field
    from stellarweave.spatial_validation import write_spatial_report
    img,truth=simulated_field(size=320,n_stars=95,seed=101,min_flux=4000)
    w=WCS(naxis=2)
    w.wcs.crpix=[1,1]
    w.wcs.crval=[100,-30]
    w.wcs.cdelt=[-0.000014,0.000014]
    w.wcs.ctype=['RA---TAN','DEC--TAN']
    path=tmp_path/'science.fits'
    fits.PrimaryHDU(img.astype('float32'),header=w.to_header()).writeto(path)
    ra,dec=w.all_pix2world(truth[:,0],truth[:,1],0)
    data=[]
    for i,(star,r,d) in enumerate(zip(truth,ra,dec)):
        a=np.zeros(23)
        a[0]=i+1
        a[1:3]=star[:2]
        a[3]=30-2.5*np.log10(star[2])
        a[4]=.02
        a[5]=.9
        a[6]=.02
        a[7]=a[3]-.9
        a[8]=.02
        a[11:13]=[2,2]
        a[13:15]=[1,1]
        a[17:19]=[.01,.01]
        a[19:21]=[.1,.1]
        a[21:23]=[r,d]
        data.append(a)
    catalog=tmp_path/'table.txt'
    np.savetxt(catalog,np.asarray(data),fmt='%.12f')
    result=write_spatial_report(path,catalog,tmp_path/'diag',
                                size=320,min_stars=2,cells=2)
    assert 'wide_shift_heldout_plate_search' in result
    assert result['wide_shift_heldout_plate_search']['science_calibrated'] is False
