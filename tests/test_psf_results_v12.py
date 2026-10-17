"""v0.12 figure pipeline: source provenance and incomplete-tile retention."""
import json
import pytest
from stellarweave.psf_results import summarize_spatial_psfs,render_spatial_psf_dashboard,CLUSTERS
from stellarweave.real_hst import sha256_file


def _reports(tmp_path,*,failed=False):
    psfs=[];plates=[]
    for c in CLUSTERS:
        folder=tmp_path/c
        folder.mkdir()
        pp=folder/'cross_region.json'
        plate={'cluster':c,'tile_size':512,'fits_sha256':'sha-validated-in-test',
               'tiles':[{'index':i,'tile_xy':[(i%3)*512,(i//3)*512]}
                        for i in range(9)]}
        pp.write_text(json.dumps(plate))
        rows=[]
        for i in range(9):
            bad=failed and c==CLUSTERS[1] and i==0
            rows.append({'tile_index':i,'tile_xy':[(i%3)*512,(i//3)*512],
                         'status':'INSUFFICIENT_ISOLATED_STARS' if bad
                                  else 'USABLE_NATIVE_PSF',
                         'used_stars':0 if bad else 20,
                         'metrics':None if bad else {'r50_px':1.25+i*.01,
                            'r80_px':2.4,'core_3x3_fraction':.6,
                            'moment_ellipticity_inner5px':.02},
                         'distance_to_anchor_total_variation':None if bad else i*.015})
        pr={'cluster':c,'kind':'REAL_HST_NATIVE_PIXEL_PSF_SPATIAL_DIAGNOSTIC',
            'native_pixel_psf_not_oversampled_ePSF':True,
            'not_calibrated_photometry':True,
            'not_independent_astrometry':True,
            'tile_count':9,'anchor_index':4,'tile_size':512,
            'input_fits_sha256':plate['fits_sha256'],
            'plate_report_sha256':sha256_file(pp),
            'tiles':rows,'usable_tiles':8 if failed and c==CLUSTERS[1] else 9,
            'failed_tiles':1 if failed and c==CLUSTERS[1] else 0,
            'anchor_usable':True,'median_r50_px':1.29,
            'r50_range_px':[1.25,1.33],
            'median_heldout_total_variation':.06}
        p=folder/'psf_spatial.json'
        p.write_text(json.dumps(pr))
        psfs.append(p);plates.append(pp)
    return psfs,plates


def test_psf_dashboard_retains_failed_tiles(tmp_path):
    psfs,plates=_reports(tmp_path,failed=True)
    out=tmp_path/'figs'
    d=render_spatial_psf_dashboard(psfs,plates,out)
    assert d['total_tiles']==27
    assert d['usable_tiles']==26 and d['failed_tiles']==1
    assert sum(r['status']=='INSUFFICIENT_ISOLATED_STARS'
               for r in d['per_tile'])==1
    for f in ('psf_spatial_r50.png','psf_spatial_mismatch.png'):
        assert (out/f).stat().st_size>8000
    assert (out/'per_tile.csv').exists()
    assert len(list((out/'source_psf_reports').glob('*.json')))==3


def test_psf_provenance_rejects_forged_fits(tmp_path):
    psfs,plates=_reports(tmp_path)
    d=json.loads(psfs[0].read_text())
    d['input_fits_sha256']='other-image'
    psfs[0].write_text(json.dumps(d))
    with pytest.raises(ValueError,match='provenance'):
        summarize_spatial_psfs(psfs,plates)


def test_missing_cluster_rejected(tmp_path):
    psfs,plates=_reports(tmp_path)
    with pytest.raises(ValueError,match='all three'):
        summarize_spatial_psfs(psfs[:2],plates)
