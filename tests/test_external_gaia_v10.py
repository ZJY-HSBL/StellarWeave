"""v0.10 Gaia epoch propagation and independent no-refit integration tests."""
import csv
import json
import numpy as np
import pytest
from stellarweave.external_gaia import (
    read_gaia_dr3_csv,propagate_to_hst_epoch,evaluate_frozen_plate)
from stellarweave.gaia_fetch import query_from_fits
from stellarweave.simulation import simulated_field


def _write_gaia(path,rows):
    fieldnames=('source_id','ra','dec','pmra','pmdec','ref_epoch',
                'phot_g_mean_mag','ruwe','visibility_periods_used')
    with open(path,'w',newline='') as f:
        writer=csv.writer(f);writer.writerow(fieldnames);writer.writerows(rows)


def test_gaia_quality_cuts(tmp_path):
    f=tmp_path/'gaia.csv'
    _write_gaia(f,[(1,100,-30,3,-4,2016,17,1.1,12),
                   (2,100,-30,1,1,2016,18,1.7,12),
                   (3,100,-30,0,0,2016,20,1.1,12)])
    t=read_gaia_dr3_csv(f)
    assert len(t['source_id'])==1
    assert t['source_id'][0]==1
    with pytest.raises(ValueError):
        read_gaia_dr3_csv(f,max_ruwe=0.4)


def test_gaia_motion_epoch_propagation():
    pytest.importorskip('astropy')
    a={'ra':np.array([100.]),'dec':np.array([-30.]),
       'pmra':np.array([100.]),'pmdec':np.array([50.]),
       'ref_epoch':np.array([2016.])}
    from astropy.time import Time
    moved=propagate_to_hst_epoch(a,Time(2006.,format='jyear').mjd)
    # Since pmra includes cos(dec), 10 years with positive pmra moves
    # the past RA left, while positive pmdec moves past DEC down.
    assert moved[0,0]<100
    assert moved[0,1]<-30
    d_ra=(100-moved[0,0])*3600*1000*np.cos(np.deg2rad(-30))
    assert 970<d_ra<1030
    with pytest.raises(ValueError):
        propagate_to_hst_epoch(a,0)


def test_external_validation_on_constructed_fits(tmp_path):
    pytest.importorskip('astropy')
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.time import Time
    image,sources=simulated_field(size=768,n_stars=850,seed=208,
                                  min_flux=5000,max_flux=18000)
    w=WCS(naxis=2)
    w.wcs.crpix=[1,1]
    w.wcs.crval=[100.,-30.]
    w.wcs.cdelt=[-1.4e-5,1.4e-5]
    w.wcs.ctype=['RA---TAN','DEC--TAN']
    science=tmp_path/'science.fits'
    fits.PrimaryHDU(image.astype('float32'),header=w.to_header()).writeto(science)
    subset=((sources[:,0]>12)&(sources[:,0]<242)&
            (sources[:,1]>12)&(sources[:,1]<242))
    star=sources[subset][:85]
    assert len(star)>12
    ra,dec=w.all_pix2world(star[:,0],star[:,1],0)
    gaia=tmp_path/'gaia.csv'
    _write_gaia(gaia,[(i+1,float(r),float(d),0.,0.,2016.,17.,1.,12)
                      for i,(r,d) in enumerate(zip(ra,dec))])
    assert len(read_gaia_dr3_csv(gaia)['source_id'])==len(star)
    tiles=[{'index':0,'tile_xy':[0,0],
            'fixed_anchor_transfer':{'status':'TRANSFER_QC_PASS'}},
           {'index':4,'tile_xy':[256,256],
            'fixed_anchor_transfer':{'status':'ANCHOR_NOT_EVALUATED_AS_TRANSFER'}}]
    plate=tmp_path/'plate.json'
    plate.write_text(json.dumps({'anchor_status':'HELDOUT_QC_PASS',
             'anchor_affine_global':[[1,0,0],[0,1,0]],'anchor_index':4,
             'tiles':tiles,'science_hdu_index':0,'tile_size':256}))
    report=evaluate_frozen_plate(science,gaia,plate,tmp_path/'out.json',
        hst_mjd=Time(2016.,format='jyear').mjd,eval_tile_index=0,
        min_sources=5)
    assert report['status']=='GAIA_DIAGNOSTIC_PASS'
    assert report['gaia_matched']>=10
    assert report['gaia_matched_fraction']>.5
    assert report['fitting_uses_gaia'] is False
    assert not report['absolute_astrometry_certified']
    with pytest.raises(ValueError,match='withheld'):
        evaluate_frozen_plate(science,gaia,plate,tmp_path/'no.json',
             hst_mjd=Time(2016.,format='jyear').mjd,eval_tile_index=4)
    assert (tmp_path/'out.json').is_file()


def test_gaia_adql_is_bounded_and_wcs_derived(tmp_path):
    pytest.importorskip('astropy')
    from astropy.io import fits
    from astropy.wcs import WCS
    w=WCS(naxis=2)
    w.wcs.crpix=[1,1]
    w.wcs.crval=[100.,-30.]
    w.wcs.cdelt=[-0.000014,0.000014]
    w.wcs.ctype=['RA---TAN','DEC--TAN']
    p=tmp_path/'science.fits'
    fits.PrimaryHDU(np.ones((512,512),dtype=np.float32),
                    header=w.to_header()).writeto(p)
    q,meta=query_from_fits(p,max_rows=1000)
    assert 'gaiadr3.gaia_source' in q
    assert 'TOP 1000' in q
    assert 'pmra IS NOT NULL' in q
    assert 0<meta['radius_deg']<.02


def test_bad_gaia_schema(tmp_path):
    f=tmp_path/'wrong.csv'
    f.write_text('ra,dec\n100,-30\n')
    with pytest.raises(ValueError,match='source_id'):
        read_gaia_dr3_csv(f)
