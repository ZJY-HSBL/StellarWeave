"""Synthetic unit tests for PRF code; NOT telescope performance results."""
import json
import numpy as np
import pytest

from stellarweave.spatial_prf import (
    extract_stamps, template_from_training, split_training_test,
    compare_tiles, bootstrap_tile_median_ci,
)


def _stamp(sigma, rng, size=15):
    r=size//2
    y,x=np.mgrid[-r:r+1,-r:r+1]
    m=np.exp(-(x*x+y*y)/(2*sigma*sigma))
    m=m/m.sum()
    m=np.maximum(m+rng.normal(0,0.00007,size=m.shape),0)
    return m/m.sum()


def _nine_tiles(seed=3):
    rng=np.random.default_rng(seed)
    result=[]
    for idx in range(9):
        sigma=1.08 if idx==4 else 1.12+idx*.085
        stamps=np.stack([_stamp(sigma,rng) for _ in range(36)])
        xy=np.column_stack((np.arange(36)*2,np.arange(36)*3))
        result.append((stamps,xy,{'usable_stamps':36,'isolated':36,'candidates':36}))
    return result


def test_normalized_training_template():
    x=_nine_tiles()[0][0]
    model=template_from_training(x)
    assert model.shape==(15,15)
    assert np.isclose(model.sum(),1)
    assert (model>=0).all()


def test_local_prf_can_improve_holdout_on_spatially_varying_synthetic_psf():
    tiles,stars,models=compare_tiles(_nine_tiles(),seed=42)
    assert len(tiles)==9
    assert len(stars)==9*9
    assert set(models)==set(range(9))
    assert tiles[4]['median_paired_gain_pct']==pytest.approx(0,abs=1e-8)
    gains=[t['median_paired_gain_pct'] for t in tiles if t['index']!=4]
    assert sum(g>0 for g in gains)>=6
    ci=bootstrap_tile_median_ci(tiles)
    assert ci is not None and ci[0]>0


def test_train_test_split_has_no_leakage():
    stamps,xy,_=_nine_tiles()[0]
    a=split_training_test(stamps,xy,seed=123)
    b=split_training_test(stamps,xy,seed=123)
    assert np.array_equal(a['template'],b['template'])
    assert a['train_count']+a['test_count']==len(stamps)
    assert a['test_count']>=3
    assert a['holdout_stamps'].shape[0]==a['test_count']


def test_fail_closed_on_sparse_tile_and_missing_anchor():
    tiles=_nine_tiles()
    tiles[0]=(tiles[0][0][:4],tiles[0][1][:4],
              {'usable_stamps':4,'isolated':4,'candidates':4})
    rows,stars,models=compare_tiles(tiles)
    assert rows[0]['status']=='INSUFFICIENT_ISOLATED_STARS'
    assert 0 not in models
    assert all(s['tile_index']!=0 for s in stars)
    assert bootstrap_tile_median_ci([rows[0]]) is None
    tiles[4]=(tiles[4][0][:3],tiles[4][1][:3],
              {'usable_stamps':3,'isolated':3,'candidates':3})
    with pytest.raises(ValueError,match='center-tile'):
        compare_tiles(tiles)


def test_extract_stamps_from_synthetic_gaussian_stars():
    rng=np.random.default_rng(3)
    image=rng.normal(15,1,size=(260,260))
    yy,xx=np.mgrid[-7:8,-7:8]
    for y in range(30,231,35):
        for x in range(30,231,35):
            image[y-7:y+8,x-7:x+8]+=2000*np.exp(
                -(xx*xx+yy*yy)/(2*1.2**2))
    stamps,coords,meta=extract_stamps(image,np.ones_like(image,bool),cap=100)
    assert len(stamps)>=20
    assert len(stamps)==len(coords)
    assert np.allclose(stamps.sum(axis=(1,2)),1,atol=1e-8)
    assert meta['oversampling']==1
    assert meta['instrument_ePSf_calibrated'] is False


def test_mosaic_provenance_rejects_wrong_sha(tmp_path):
    pytest.importorskip('astropy')
    from astropy.io import fits
    from stellarweave.spatial_prf import audit_spatial_prf
    f=tmp_path/'wrong.fits'
    fits.PrimaryHDU(np.zeros((1536,1536),dtype='float32')).writeto(f)
    plate={'tile_size':512,'grid':3,'anchor_index':4,'tile_count':9,
           'heldout_tiles_physically_nonoverlapping':True,
           'anchor_status':'HELDOUT_QC_PASS',
           'fits_sha256':'wrong','science_hdu_index':0}
    p=tmp_path/'plate.json'
    p.write_text(json.dumps(plate))
    with pytest.raises(ValueError,match='checksum'):
        audit_spatial_prf(f,p,tmp_path/'out')
