"""v0.6 real FITS audit tests; only constructed fixtures, not HST results."""
import json
import numpy as np
import pytest
from stellarweave.real_hst import empirical_psf
from stellarweave.download_acsggct import filenames
from stellarweave.simulation import simulated_field


def test_archive_inventory_fails_closed_for_unpublished_f814_reference():
    image,cat=filenames('ngc6397','f606w')
    assert image.endswith('f606w_v2_img.fits')
    assert cat.endswith('r.rdviq.cal.adj.zpt')
    with pytest.raises(ValueError,match='F814W'):
        filenames('ngc6397','f814w')


def test_empirical_psf_on_toy_field():
    im,truth=simulated_field(size=320,n_stars=65,seed=11,min_flux=4000)
    psf,xy,report=empirical_psf(im,np.isfinite(im),stamp_size=15,
                                min_stars=3,fwhm=2.6)
    assert len(xy)>=3
    assert psf.shape==(15,15)
    assert psf.min()>=0
    assert abs(psf.sum()-1)<1e-12
    assert psf[7,7] > psf[0,0]
    assert report['no_instrument_calibration']


def test_psf_rejects_too_few_isolated_stars():
    im=np.ones((100,100))
    with pytest.raises(ValueError,match='isolated'):
        empirical_psf(im,np.isfinite(im),min_stars=4)


def test_real_hst_audit_constructed_fits(tmp_path):
    pytest.importorskip('astropy')
    from astropy.io import fits
    from astropy.wcs import WCS
    from stellarweave.real_hst import audit_product,read_science_crop
    im,truth=simulated_field(size=320,n_stars=60,seed=35,min_flux=5000)
    w=WCS(naxis=2)
    w.wcs.crpix=[1,1]
    w.wcs.crval=[100,-30]
    w.wcs.cdelt=[-0.000014,0.000014]
    w.wcs.ctype=['RA---TAN','DEC--TAN']
    hdr=w.to_header()
    hdr['BUNIT']='ELECTRONS/S'
    hdr['FILTER1']='F606W'
    science=fits.ImageHDU(im.astype('float32'),header=hdr,name='SCI')
    errs=fits.ImageHDU(np.full(im.shape,3.,dtype=np.float32),name='ERR')
    dq=fits.ImageHDU(np.zeros(im.shape,np.int16),name='DQ')
    dq.data[1,1]=1
    path=tmp_path/'test.fits'
    fits.HDUList([fits.PrimaryHDU(),science,errs,dq]).writeto(path)
    ra,dec=w.all_pix2world(truth[:,0],truth[:,1],0)
    rows=[]
    for i,((x,y,flux),r,d) in enumerate(zip(truth,ra,dec)):
        line=np.zeros(23)
        line[0]=i+1
        line[1:3]=x,y
        line[3]=30-2.5*np.log10(flux)
        line[4]=.02
        line[5]=.8
        line[6]=.03
        line[7]=line[3]-.8
        line[8]=.02
        line[11:13]=[3,3]
        line[13:15]=[1,1]
        line[15:17]=[.1,.1]
        line[17:19]=[.1,.1]
        line[19:21]=[.1,.1]
        line[21:23]=[r,d]
        rows.append(line)
    catalog=tmp_path/'stars.txt'
    np.savetxt(catalog,np.asarray(rows),fmt='%.12f')
    cropped,mask,err,info=read_science_crop(path,size=320)
    assert mask[1,1]==False
    assert info['has_err_plane'] and info['has_dq_plane']
    assert info['has_celestial_wcs']
    assert 'F606W' in info['filters']
    report=audit_product(path,catalog,tmp_path/'qc',size=320,min_stars=3)
    assert report['kind']=='PRELIMINARY_REAL_DATA_AUDIT'
    assert report['catalog']['row_count']==len(truth)
    assert report['crop_psf']['used_stars']>=3
    assert (tmp_path/'qc'/'empirical_prf_native.npy').is_file()
    assert (tmp_path/'qc'/'psf_preview.png').is_file()
    result=json.loads((tmp_path/'qc'/'audit.json').read_text())
    assert result['catalog']['coordinate_reference_not_verified']


def test_nonexistent_real_file_error(tmp_path):
    from stellarweave.real_hst_cli import main
    with pytest.raises(FileNotFoundError):
        main(['--fits',str(tmp_path/'no.fits'),'--catalog',str(tmp_path/'none.txt')])
