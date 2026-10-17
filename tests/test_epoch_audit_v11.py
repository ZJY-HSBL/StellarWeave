"""HST epoch audit must never convert a header hint to verified Gaia epoch."""
import pytest
from stellarweave.epoch_audit import audit_epoch


def test_epoch_audit_reads_realistic_hst_keywords(tmp_path):
    pytest.importorskip('astropy')
    from astropy.io import fits
    import numpy as np
    hdr=fits.Header()
    hdr['EXPSTART']=53870.125
    hdr['EXPEND']=53870.135
    hdr['DATE-OBS']='2006-05-15'
    hdr['HISTORY']='DRIZZLED MULTIPLE INPUT EXPOSURES'
    p=tmp_path/'drizzled.fits'
    fits.PrimaryHDU(np.ones((40,40),dtype='float32'),header=hdr).writeto(p)
    report=audit_epoch(p,tmp_path/'epoch.json')
    assert report['explicit_parsed_mjd_count']>=2
    assert report['min_header_mjd'] is not None
    assert report['max_header_mjd'] is not None
    assert report['header_mjd_span_days']>0
    assert report['effective_mosaic_epoch_verified'] is False
    assert report['gaia_epoch_can_be_inferred_automatically'] is False
    assert (tmp_path/'epoch.json').is_file()


def test_missing_epoch_stays_unknown(tmp_path):
    pytest.importorskip('astropy')
    from astropy.io import fits
    import numpy as np
    p=tmp_path/'nodate.fits'
    fits.PrimaryHDU(np.zeros((32,32),dtype='float32')).writeto(p)
    result=audit_epoch(p)
    assert result['min_header_mjd'] is None
    assert result['max_header_mjd'] is None
    assert result['effective_mosaic_epoch_verified'] is False
