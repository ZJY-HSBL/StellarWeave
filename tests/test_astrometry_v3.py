"""Tests for v0.3 affine registration and paired-band reference export."""
import csv
import numpy as np
import pytest

from stellarweave.astrometry import apply_affine, fit_affine
from stellarweave.dualband import paired_reference, build_dualband
from stellarweave.registration_pipeline import register_catalog


def test_affine_known_geometry_and_outliers():
    rng=np.random.default_rng(102)
    source=rng.uniform(20,900,(220,2))
    truth=np.array([[1.0006,-0.0008,5.4],[0.0007,0.9995,-3.8]])
    observed=apply_affine(source,truth)+rng.normal(0,.12,(220,2))
    false=rng.uniform(0,1000,(100,2))
    report=fit_affine(source,np.vstack((observed,false)))
    np.testing.assert_allclose(report['matrix'],truth,atol=.08)
    assert report['inliers'] >= 190
    assert report['rms_px'] < .4
    assert 0 < report['matched_fraction'] <= 1
    assert len(report['local_residual_grid']) >= 4


def test_affine_rejects_unrelated_fields():
    rng=np.random.default_rng(23)
    source=rng.uniform(0,400,(80,2))
    other=rng.uniform(2000,2500,(90,2))
    with pytest.raises(ValueError):
        fit_affine(source,other,initial_shift=[0,0])


def test_affine_rejects_scale_error():
    rng=np.random.default_rng(19)
    source=rng.uniform(0,600,(120,2))
    observed=source * 1.1 + [3,-2]
    with pytest.raises(ValueError):
        fit_affine(source,observed,initial_shift=[3,-2],match_radius=80)


def test_manual_alignment_is_unverified():
    img=np.zeros((100,100),dtype=float)
    source=np.array([[10.,20.],[30.,40.]])
    adjusted,report,verified=register_catalog(img,source,[19,20],
                                              mode='manual',dx=2,dy=-3)
    np.testing.assert_allclose(adjusted,source+[2,-3])
    assert not verified
    assert report['status']=='UNVERIFIED_USER_SHIFT'


def test_two_band_selection_shared_ids():
    n=30
    cat={'source_id':np.arange(1,n+1,dtype=float),
         'ra_deg':np.linspace(100,100.01,n),
         'dec_deg':np.linspace(-30,-30.01,n)}
    for band in ('f606w','f814w'):
        cat['mag_'+band]=np.full(n,20.0 if band=='f606w' else 19.0)
        cat['magerr_'+band]=np.full(n,.02)
        cat['n_'+band]=np.full(n,2.)
        cat['flag_'+band]=np.full(n,1.)
        cat['crowding_'+band]=np.zeros(n)
        cat['psfq_'+band]=np.full(n,.03)
    m=paired_reference(cat)
    assert m.sum()==n
    cat['flag_f814w'][0]=3
    cat['magerr_f606w'][1]=2
    assert paired_reference(cat).sum()==n-2


def test_two_band_synthetic_fits(tmp_path):
    astropy=pytest.importorskip('astropy')
    from astropy.io import fits
    from astropy.wcs import WCS
    from stellarweave.simulation import simulated_field
    rng=np.random.default_rng(71)
    image,stars=simulated_field(size=512,n_stars=110,seed=71,min_flux=4500)
    second,unused=simulated_field(size=512,n_stars=110,seed=71,min_flux=4500)
    w=WCS(naxis=2)
    w.wcs.crpix=[1,1]
    w.wcs.crval=[100,-30]
    w.wcs.cdelt=[-.000014,.000014]
    w.wcs.ctype=['RA---TAN','DEC--TAN']
    first_path=tmp_path/'first.fits'
    second_path=tmp_path/'second.fits'
    fits.PrimaryHDU(image.astype('float32'),header=w.to_header()).writeto(first_path)
    fits.PrimaryHDU(second.astype('float32'),header=w.to_header()).writeto(second_path)
    ra,dec=w.all_pix2world(stars[:,0],stars[:,1],0)
    values=[]
    for i,(s,ar,d) in enumerate(zip(stars,ra,dec)):
        row=np.zeros(21)
        row[0]=i+1
        row[1:3]=s[:2]
        row[3]=30-2.5*np.log10(s[2])
        row[4]=.02
        row[5]=.9
        row[6]=.03
        row[7]=row[3]-.9
        row[8]=.03
        row[11:13]=[3,3]
        row[13:15]=[1,1]
        row[15:17]=[.1,.1]
        row[17:19]=[.1,.1]
        row[19:21]=[ar,d]
        values.append(row)
    cat=tmp_path/'cat.txt'
    np.savetxt(cat,np.asarray(values),fmt='%.12f')
    report=build_dualband(first_path,second_path,cat,tmp_path/'paired')
    assert report['registered']
    assert report['paired_records'] > 30
    with (tmp_path/'paired'/'paired_reference.csv').open() as f:
        reader=list(csv.DictReader(f))
    assert len(reader)==report['paired_records']
    assert abs(float(reader[0]['color_f606w_minus_f814w'])-.9)<1e-7
