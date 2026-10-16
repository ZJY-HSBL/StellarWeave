from pathlib import Path
import numpy as np
import pytest

from stellarweave.catalog import FIELDS, read_acsggct_catalog, quality_mask, project_catalog
from stellarweave.alignment import estimate_translation
from stellarweave.dataset import spatial_split, write_tiles

def make_catalog(path, rows=4, *, bad_ra=False):
    result = []
    for i in range(rows):
        a = np.zeros(23)
        a[0] = i+1
        a[1:3] = [100 + 10*i, 130+7*i]
        a[3] = 22. + i * .2
        a[4] = .03
        a[5] = .9
        a[6] = .03
        a[7] = 21. + i * .2
        a[8] = .04
        a[11:13] = [3,3]
        a[13:15] = [1,1]
        a[15:17] = [.1,.2]
        a[17:19] = [.1,.1]
        a[19:21] = [.1, .1]
        a[21:23] = [999. if bad_ra else 100. + i*.01, -30.]
        result.append(a)
    np.savetxt(path, result, fmt='%.10f')
    return np.array(result)

def test_read_schema_and_filter(tmp_path):
    sample = tmp_path/'catalog.txt'
    make_catalog(sample)
    cat = read_acsggct_catalog(sample)
    assert len(FIELDS) == 23
    assert cat['ra_deg'][0] == 100
    assert cat['mag_f606w'][0] == 22
    assert quality_mask(cat, 'f606w').sum() == 4
    cat['flag_f606w'][0] = 3
    assert quality_mask(cat, 'f606w').sum() == 3

def test_bad_catalog_rejected(tmp_path):
    path = tmp_path/'invalid.txt'
    make_catalog(path, bad_ra=True)
    with pytest.raises(ValueError, match='RA/Dec'):
        read_acsggct_catalog(path)

def test_missing_columns_rejected(tmp_path):
    path = tmp_path/'invalid.txt'
    path.write_text('1 2 3 4\n')
    with pytest.raises(ValueError, match='23'):
        read_acsggct_catalog(path)

def test_catalog_wcs_zero_based(tmp_path):
    path = tmp_path/'catalog.txt'
    make_catalog(path)
    cat = read_acsggct_catalog(path)
    class FakeWCS:
        def all_world2pix(self, ra, dec, origin):
            assert origin == 0
            return ra-100., dec+30.
    coords = project_catalog(cat, FakeWCS())
    assert abs(coords[0,0]) < 1e-12
    assert abs(coords[0,1]) < 1e-12

def test_translation_with_outliers():
    rng = np.random.default_rng(9)
    p = rng.uniform(10, 450, size=(180,2))
    true_shift = np.array([7.4, -6.2])
    detections = p + true_shift + rng.normal(0, .13, (180,2))
    extra = rng.uniform(0, 500, size=(80,2))
    outcome = estimate_translation(p, np.vstack([detections,extra]))
    np.testing.assert_allclose([outcome['dx_px'], outcome['dy_px']], true_shift, atol=.12)
    assert outcome['inliers'] >= 170
    assert outcome['rms_px'] < .3

def test_registration_fails_on_unrelated_fields():
    p = np.column_stack((np.arange(25)*10, np.zeros(25)))
    q = np.column_stack((np.arange(25)*10, np.ones(25)*500))
    with pytest.raises(ValueError):
        estimate_translation(p, q)

def test_spatial_splits_and_tiles(tmp_path):
    width = 1024
    for x in range(0,width-64+1,64):
        s = spatial_split(x, 64, width, guard=16)
        if s=='train': assert x+64 <= int(width*.6)-16
        if s=='val':
            assert x >= int(width*.6)+16
            assert x+64 <= int(width*.8)-16
        if s=='test': assert x >= int(width*.8)+16
    img=np.ones((256,1024),dtype=np.float32)
    valid=np.ones_like(img,dtype=bool)
    coords=np.array([[10., 10.],[680.,70.],[940.,10.]])
    fields={'source_id': np.array([1,2,3]),'mag':np.array([20.,21.,22.]),
            'mag_error':np.array([.1,.1,.1])}
    idx=write_tiles(img,valid,coords,fields,tmp_path,tile_size=64,guard=16)
    assert all(Path(tmp_path, rec['path']).is_file() for rec in idx)
    assert sum(rec['star_count'] for rec in idx)==3
    item=np.load(Path(tmp_path,idx[0]['path']))
    assert 'image' in item.files and 'labels' in item.files

def test_astropy_fits_pipeline(tmp_path):
    pytest.importorskip('astropy')
    from astropy.io import fits
    from astropy.wcs import WCS
    from stellarweave.simulation import simulated_field
    from stellarweave.dataset import prepare_dataset
    image, injected = simulated_field(size=512,n_stars=90,seed=99,min_flux=3000,max_flux=12000)
    w = WCS(naxis=2)
    w.wcs.crpix=[1,1]
    w.wcs.crval=[100.,-30.]
    w.wcs.cdelt=[-0.000014,0.000014]
    w.wcs.ctype=['RA---TAN','DEC--TAN']
    fits_path=tmp_path/'test.fits'
    fits.PrimaryHDU(image.astype(np.float32),header=w.to_header()).writeto(fits_path)
    ra,dec=w.all_pix2world(injected[:,0],injected[:,1],0)
    values=[]
    for i,(star,ar,decl) in enumerate(zip(injected,ra,dec)):
        row = np.zeros(23)
        row[0]=i+1
        row[1:3]=star[:2]
        row[3]=30-2.5*np.log10(star[2])
        row[4]=.03
        row[5]=.9
        row[6]=.03
        row[7]=row[3]-.9
        row[8]=.04
        row[11:13]=[2,2]
        row[13:15]=[1,1]
        row[15:17]=[.1,.1]
        row[17:19]=[.1,.1]
        row[19:21]=[.1,.1]
        row[21:23]=[ar,decl]
        values.append(row)
    cat=tmp_path/'test_catalog.txt'
    np.savetxt(cat,np.array(values),fmt='%.12f')
    dest=tmp_path/'tiles'
    result=prepare_dataset(fits_path,cat,dest,tile_size=64,guard=8)
    assert result['registration_passed_qc']
    assert result['tiles']['train']>0
    assert result['star_references_in_tiles']>0
    assert (dest/'manifest.json').is_file()
