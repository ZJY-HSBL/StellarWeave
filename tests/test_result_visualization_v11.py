"""Reproducible scientific README plots; synthetic fixture is NOT HST data."""
import json
import pytest
from stellarweave.result_visualization import summarize,render_report,CLUSTERS


def _fixture(path,name,*,failed=False):
    tiles=[]
    for i in range(9):
        status=('ANCHOR_NOT_EVALUATED_AS_TRANSFER' if i==4
                else 'TRANSFER_QC_REJECTED' if failed and i==7
                else 'TRANSFER_QC_PASS')
        tile={'index':i,'tile_xy':[i%3*512,i//3*512],
              'fixed_anchor_transfer':{
                  'status':status,'refitted':False,'evaluated':100,
                  'matches':90 if status=='TRANSFER_QC_PASS' else 0,
                  'match_fraction':.9 if status=='TRANSFER_QC_PASS' else 0.,
                  'median_residual_px':.41 if status=='TRANSFER_QC_PASS' else None,
                  'p95_residual_px':.71 if status=='TRANSFER_QC_PASS' else None,
                  'match_over_null_ratio':12.}}
        tiles.append(tile)
    report={'cluster':name,'version':'0.9.0','grid':3,'tile_size':512,
        'tile_count':9,'anchor_index':4,'tiles':tiles,
        'anchor_status':'HELDOUT_QC_PASS',
        'catalog_rows':1000,'quality_selected_catalog':850,
        'transfer_tiles_pass':7 if failed else 8,
        'transfer_tiles_total':8,
        'heldout_tiles_physically_nonoverlapping':True,
        'independent_external_reference':False,
        'science_calibrated':False}
    path.write_text(json.dumps(report))


def test_result_figures_and_failure_retention(tmp_path):
    paths=[tmp_path/f'{name}.json' for name in CLUSTERS]
    for i,(name,path) in enumerate(zip(CLUSTERS,paths)):
        _fixture(path,name,failed=(i==2))
    out=tmp_path/'figures'
    result=render_report(paths,out,source_run_id='synthetic-test')
    assert result['heldout_regions_total']==24
    assert result['heldout_regions_passed']==23
    assert result['not_external_gaia_validation']
    assert len(result['per_tile'])==24
    assert sum(r['status']=='TRANSFER_QC_REJECTED' for r in result['per_tile'])==1
    assert len(list((out/'source_reports').glob('*.json')))==3
    for figure in ('spatial_residuals.png','spatial_match_fraction.png',
                   'catalog_source_counts.png'):
        assert (out/figure).stat().st_size>8000
    assert (out/'per_tile.csv').stat().st_size>100


def test_missing_cluster_fails_closed(tmp_path):
    p=tmp_path/'one.json'
    _fixture(p,CLUSTERS[0])
    with pytest.raises(ValueError,match='Need exactly'):
        summarize([p])


def test_duplicate_cluster_fails_closed(tmp_path):
    a,b=tmp_path/'a.json',tmp_path/'b.json'
    _fixture(a,CLUSTERS[0])
    _fixture(b,CLUSTERS[0])
    with pytest.raises(ValueError,match='Duplicate'):
        summarize([a,b])


def test_bad_science_claims_rejected(tmp_path):
    p=tmp_path/'bad.json'
    _fixture(p,CLUSTERS[0])
    d=json.loads(p.read_text())
    d['science_calibrated']=True
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='science-grade'):
        summarize([p])


def test_header_epoch_audit_retained_without_claiming_verified(tmp_path):
    paths=[]
    for i,name in enumerate(CLUSTERS):
        directory=tmp_path/name
        directory.mkdir()
        p=directory/'cross_region.json'
        _fixture(p,name)
        (directory/'epoch_audit.json').write_text(json.dumps({
            'effective_mosaic_epoch_verified':False,
            'min_header_mjd':53800.+i,
            'max_header_mjd':53800.+i,
            'gaia_epoch_can_be_inferred_automatically':False}))
        paths.append(p)
    out=tmp_path/'result'
    d=render_report(paths,out)
    assert len(d['header_epoch_hints'])==3
    assert all(not x['effective_mosaic_epoch_verified']
               for x in d['header_epoch_hints'].values())
    assert len(list((out/'epoch_audits').glob('*.json')))==3
    saved=json.loads((out/'summary.json').read_text())
    assert saved['header_epoch_hints']==d['header_epoch_hints']
