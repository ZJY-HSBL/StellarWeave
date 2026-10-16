"""v0.7 spatial diagnostics: synthetic WCS truth and conservative failure tests."""
import json
import numpy as np
import pytest

from stellarweave.simulation import simulated_field
from stellarweave.spatial_validation import (
    match_unique, correspondence_diagnostics, spatial_psf_grid)


def test_match_unique_does_not_reuse_peak():
    ref=np.array([[10.,10.],[10.4,10.],[50.,50.]])
    detections=np.array([[10.1,10.],[50.,50.]])
    pairs=match_unique(ref,detections,radius=1.5)
    assert len(pairs)==2
    assert len(set(pairs[:,1]))==2


def test_match_unique_empty_and_invalid():
    assert match_unique([],[[10.,10.]]).shape==(0,3)
    with pytest.raises(ValueError):
        match_unique([[1,2]],[[1,2]],radius=0)


def test_wcs_catalog_correspondence_recovers_synthetic_stars():
    image,truth=simulated_field(size=256,n_stars=36,seed=140,min_flux=5000)
    report,refs,det,matrix=correspondence_diagnostics(
        image,np.isfinite(image),truth[:,:2],
        np.full(len(truth),18.),radius=3.,min_inliers=10)
    assert report['raw_wcs_match_count']>=12
    assert report['raw_wcs_match_fraction']>.33
    assert report['reference_in_crop']>=20
    assert not report['astrometric_science_validated']
    assert report['status'] in ('CANDIDATE_AFFINE_QC_PASSED','AFFINE_QC_REJECTED')


def test_wrong_catalog_does_not_create_spurious_science_pass():
    image=np.ones((256,256),dtype=float)
    stars=np.column_stack((np.arange(30)*6+10.,np.arange(30)*6+10.))
    report,_,_,_=correspondence_diagnostics(
        image,np.isfinite(image),stars,np.full(30,20.))
    assert report['status']=='INSUFFICIENT_DATA'
    assert report['raw_wcs_match_count']==0


def test_psf_spatial_grid_validates_output():
    image,_=simulated_field(size=320,n_stars=68,seed=21,min_flux=5000)
    cells,stamps=spatial_psf_grid(
        image,np.isfinite(image),cells=2,stamp_size=15,min_stars=2)
    assert len(cells)==4
    assert len(stamps)==4
    successful=[p for p in stamps if p is not None]
    assert len(successful)>=2
    assert all(np.isclose(p.sum(),1) and p.shape==(15,15) for p in successful)


def test_spatial_grid_rejects_invalid_shape():
    with pytest.raises(ValueError):
        spatial_psf_grid(np.zeros((30,30)),np.zeros((10,10),bool))
    with pytest.raises(ValueError):
        spatial_psf_grid(np.zeros((30,30)),np.ones((30,30),bool),cells=5)


def test_real_fits_wcs_spatial_diagnostics(tmp_path):
    pytest.importorskip('astropy')
    from astropy.io import fits
    from astropy.wcs import WCS
    from stellarweave.spatial_validation import write_spatial_report
    image,truth=simulated_field(size=320,n_stars=95,seed=44,min_flux=4000)
    w=WCS(naxis=2)
    w.wcs.crpix=[1,1]
    w.wcs.crval=[100,-30]
    w.wcs.cdelt=[-0.000014,0.000014]
    w.wcs.ctype=['RA---TAN','DEC--TAN']
    header=w.to_header()
    header['BUNIT']='ELECTRONS/S'
    header['FILTER1']='F606W'
    data_path=tmp_path/'image.fits'
    fits.PrimaryHDU(image.astype('float32'),header=header).writeto(data_path)
    ra,dec=w.all_pix2world(truth[:,0],truth[:,1],0)
    rows=[]
    for i,(star,world_ra,world_dec) in enumerate(zip(truth,ra,dec)):
        line=np.zeros(23)
        line[0]=i+1
        line[1:3]=star[:2]
        line[3]=30-2.5*np.log10(star[2])
        line[4]=.02
        line[5]=.9
        line[6]=.03
        line[7]=line[3]-.9
        line[8]=.03
        line[11:13]=[3,3]
        line[13:15]=[1,1]
        line[17:19]=[.01,.01]
        line[19:21]=[.1,.1]
        line[21:23]=[world_ra,world_dec]
        rows.append(line)
    catalog_path=tmp_path/'stars.txt'
    np.savetxt(catalog_path,np.asarray(rows),fmt='%.12f')
    result=write_spatial_report(data_path,catalog_path,tmp_path/'out',
        size=320,cells=2,min_stars=2,matching_radius=3)
    assert result['catalog_rows']==95
    assert result['psf_candidate_cells']>=2
    assert result['catalog_wcs_registration']['raw_wcs_match_count']>=15
    assert result['scientific_astrometry_validated'] is False
    for name in ['spatial_psf_grid.png','astrometry_overlay.png',
                 'spatial_diagnostics.json']:
        assert (tmp_path/'out'/name).exists()
    r=json.loads((tmp_path/'out'/'spatial_diagnostics.json').read_text())
    assert r['kind']=='EXPLORATORY_SPATIAL_HST_DIAGNOSTIC'
