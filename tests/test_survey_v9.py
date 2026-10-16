"""v0.9 disjoint-field transfer protocol and cluster regression tests."""
import json
import numpy as np
import pytest

from stellarweave.simulation import simulated_field
from stellarweave.astrometry import apply_affine
from stellarweave.robust_plate import plate_diagnostics
from stellarweave.survey_validation import (
    tile_origins, global_affine, transfer_score, assess_mosaic, validate_cluster)


def test_disjoint_tiles():
    xy=tile_origins((768,1024),256,3)
    assert len(xy)==9
    assert xy[4]==((1024-256)//2,(768-256)//2)
    assert len(set(xy))==9
    for i,(x,y) in enumerate(xy):
        for u,v in xy[i+1:]:
            assert abs(x-u)>=256 or abs(y-v)>=256
    with pytest.raises(ValueError):
        tile_origins((400,400),256,3)


def test_reference_footprint_drives_disjoint_tile_selection():
    shape=(8000,8000)
    cells=tile_origins(shape,512,grid=3,anchor_center_xy=(4010,3970))
    assert cells[4][0]>3000 and cells[4][0]<4400
    assert cells[4][1]>3000 and cells[4][1]<4400
    assert all(x>=0 and x+512<=8000 and y>=0 and y+512<=8000
               for x,y in cells)
    assert all(abs(x-u)>=512 or abs(y-v)>=512
               for i,(x,y) in enumerate(cells) for u,v in cells[i+1:])
    # A footprint near the mosaic edge must be moved as a whole, not
    # individually clipped into overlapping tiles.
    near_edge=tile_origins(shape,512,grid=3,anchor_center_xy=(25,50))
    assert near_edge[0]==(0,0)
    assert near_edge[4]==(512,512)


def test_affine_lift_to_full_mosaic():
    a=np.array([[1.0002,-.0004,12.1],[.0003,.9998,-7.5]])
    origin=(380,920)
    lifted=global_affine(a,origin)
    p=np.array([[100.,300.],[1024.,1800.]])
    expected=apply_affine(p-np.asarray(origin),a)+origin
    np.testing.assert_allclose(apply_affine(p,lifted),expected,atol=1e-10)


def test_frozen_transfer_rejects_mismatched_regions():
    rng=np.random.default_rng(5)
    ref=rng.uniform(5,251,(95,2))
    mags=np.full(len(ref),20.)
    x=transfer_score(ref,mags,[[1,0,0],[0,1,0]],(0,0),ref,
                     256,max_stars=90)
    assert x['status']=='TRANSFER_QC_PASS'
    assert x['matches']>=85
    bad=transfer_score(ref,mags,[[1,0,2000],[0,1,2000]],(0,0),ref,
                       256,max_stars=90)
    assert bad['status']=='INSUFFICIENT_TEST_REFERENCES'


def test_small_field_dynamic_train_holdout():
    rng=np.random.default_rng(31)
    x=rng.uniform(20,220,(80,2))
    detections=x+np.array([12.,-4.])+rng.normal(0,.1,(80,2))
    d=plate_diagnostics(x,detections,shape=(256,256),max_shift_px=40,
                        train_limit=260,holdout_limit=150,
                        min_train_matches=12,min_holdout_matches=10)
    assert d['status']=='HELDOUT_QC_PASS',d
    assert d['training_reference_count']>=40
    assert d['holdout_reference_count']>=10


def test_synthetic_cross_region_transfer():
    # Dense but high S/N known positions over nine physically disjoint regions.
    image,truth=simulated_field(size=768,n_stars=1200,seed=81,
        min_flux=5500,max_flux=24000)
    report=assess_mosaic(image,np.isfinite(image),truth[:,:2],
                         np.full(len(truth),20.),tile_size=256,grid=3,seed=42)
    assert report['tile_count']==9
    assert report['heldout_tiles_physically_nonoverlapping']
    assert report['transform_frozen_on_heldout_regions']
    assert report['anchor_status']=='HELDOUT_QC_PASS',report['tiles'][4]['local_candidate_fit']
    assert report['transfer_tiles_pass']>=5,report['failed_transfer_tiles']
    assert report['science_calibrated'] is False


def test_validation_cli_on_constructed_fits(tmp_path):
    pytest.importorskip('astropy')
    from astropy.io import fits
    from astropy.wcs import WCS
    image,truth=simulated_field(size=768,n_stars=1200,seed=14,
                                min_flux=5500,max_flux=25000)
    w=WCS(naxis=2)
    w.wcs.crpix=[1,1]
    w.wcs.crval=[100,-30]
    w.wcs.cdelt=[-1.4e-5,1.4e-5]
    w.wcs.ctype=['RA---TAN','DEC--TAN']
    fits_path=tmp_path/'mockscience.fits'
    fits.PrimaryHDU(image.astype('float32'),header=w.to_header()).writeto(fits_path)
    ra,dec=w.all_pix2world(truth[:,0],truth[:,1],0)
    table=[]
    for i,(source,r,d) in enumerate(zip(truth,ra,dec)):
        row=np.zeros(23)
        row[0]=i+1
        row[1:3]=source[:2]
        row[3]=30-2.5*np.log10(source[2])
        row[4]=.02
        row[5]=.9
        row[6]=.03
        row[7]=row[3]-.9
        row[8]=.03
        row[11:13]=[3,3]
        row[13:15]=[1,1]
        row[17:19]=[.02,.02]
        row[19:21]=[.1,.1]
        row[21:23]=[r,d]
        table.append(row)
    catalog_path=tmp_path/'23col.txt'
    np.savetxt(catalog_path,np.array(table),fmt='%.12f')
    data=validate_cluster(fits_path,catalog_path,tmp_path/'out',
                          cluster='fixture',tile_size=256,grid=3)
    assert data['tile_count']==9
    assert data['cluster']=='fixture'
    assert data['anchor_status']=='HELDOUT_QC_PASS'
    assert (tmp_path/'out'/'cross_region.json').exists()
    assert (tmp_path/'out'/'cross_region.csv').exists()
    summary=json.loads((tmp_path/'out'/'cross_region.json').read_text())
    assert not summary['science_calibrated']
