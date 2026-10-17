"""Reproducible real-HST PSF spatial evidence, never synthetic observations."""
from __future__ import annotations
import argparse
import csv
import json
import shutil
from pathlib import Path
import numpy as np

CLUSTERS=('ngc6397','ngc6752','ngc0104')
NAMES={'ngc6397':'NGC 6397','ngc6752':'NGC 6752','ngc0104':'NGC 104'}


def summarize_spatial_psfs(psf_paths,plate_paths):
    """Require a complete 3-cluster, 9-tile PSF analysis with matching FITS."""
    from .real_hst import sha256_file
    reports={};plates={};files={}
    for path in psf_paths:
        d=json.loads(Path(path).read_text(encoding='utf-8'))
        c=d.get('cluster')
        if c not in CLUSTERS or c in reports:
            raise ValueError('Unexpected/duplicate PSF cluster')
        reports[c]=d
        files[c]=Path(path)
    for path in plate_paths:
        d=json.loads(Path(path).read_text(encoding='utf-8'))
        c=d.get('cluster')
        if c not in CLUSTERS or c in plates:
            raise ValueError('Unexpected/duplicate plate cluster')
        plates[c]=d
    if set(reports)!=set(CLUSTERS) or set(plates)!=set(CLUSTERS):
        raise ValueError('Require all three PSF reports AND all three plate reports')
    clusters=[];tiles=[]
    for c in CLUSTERS:
        d=reports[c];p=plates[c]
        if (d.get('kind')!='REAL_HST_NATIVE_PIXEL_PSF_SPATIAL_DIAGNOSTIC'
            or d.get('native_pixel_psf_not_oversampled_ePSF') is not True
            or d.get('not_calibrated_photometry') is not True
            or d.get('not_independent_astrometry') is not True
            or d.get('tile_count')!=9 or d.get('anchor_index')!=4
            or d.get('tile_size')!=p.get('tile_size')
            or d.get('input_fits_sha256')!=p.get('fits_sha256')
            or d.get('plate_report_sha256')!=sha256_file(next(
                x for x in plate_paths if json.loads(Path(x).read_text())['cluster']==c))):
            raise ValueError(f'PSF provenance/geometry/QC mismatch: {c}')
        actual=sorted(d['tiles'],key=lambda x:x['tile_index'])
        ref=sorted(p['tiles'],key=lambda x:x['index'])
        if ([x['tile_index'] for x in actual]!=list(range(9))
            or [x['tile_xy'] for x in actual]!=[x['tile_xy'] for x in ref]):
            raise ValueError(f'PSF tile origins do not match frozen plate: {c}')
        usable=0
        for row in actual:
            status=row['status']
            if status not in ('USABLE_NATIVE_PSF','INSUFFICIENT_ISOLATED_STARS'):
                raise ValueError(f'Unexpected PSF tile status: {status}')
            good=status=='USABLE_NATIVE_PSF'
            if good!=(row['metrics'] is not None):
                raise ValueError('Missing or unearned PSF metrics')
            usable+=int(good)
            m=row['metrics'] or {}
            tiles.append({'cluster':c,'tile_index':row['tile_index'],
                          'status':status,'used_stars':int(row['used_stars']),
                          'r50_px':m.get('r50_px'),
                          'r80_px':m.get('r80_px'),
                          'core_3x3_fraction':m.get('core_3x3_fraction'),
                          'moment_ellipticity_inner5px':m.get('moment_ellipticity_inner5px'),
                          'total_variation_from_center':row['distance_to_anchor_total_variation']})
        if usable!=d['usable_tiles'] or 9-usable!=d['failed_tiles']:
            raise ValueError(f'Inconsistent PSF usable/failed counts: {c}')
        clusters.append({'cluster':c,'usable_tiles':usable,'failed_tiles':9-usable,
                         'anchor_usable':bool(d['anchor_usable']),
                         'median_r50_px':d['median_r50_px'],
                         'r50_range_px':d['r50_range_px'],
                         'median_heldout_total_variation':d['median_heldout_total_variation'],
                         'input_fits_sha256':d['input_fits_sha256'],
                         'psf_report_sha256':sha256_file(files[c])})
    return {'version':'0.12.0','kind':'REAL_HST_NATIVE_PIXEL_SPATIAL_PSF_SUMMARY',
            'clusters':clusters,'per_tile':tiles,'total_tiles':27,
            'usable_tiles':sum(c['usable_tiles'] for c in clusters),
            'failed_tiles':sum(c['failed_tiles'] for c in clusters),
            'not_oversampled_ePSF':True,
            'not_calibrated_photometry':True,
            'not_external_astrometry':True,
            'limitations':'Native-pixel stacks of peak-selected bright isolated sources. '
                          'No subpixel oversampling, CTE correction, field-dependent '
                          'photometric calibration, correlated-noise uncertainty or '
                          'independent astrometric validation.'}


def render_spatial_psf_dashboard(psf_paths,plate_paths,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .real_hst import sha256_file
    result=summarize_spatial_psfs(psf_paths,plate_paths)
    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    raw=out/'source_psf_reports'
    raw.mkdir(exist_ok=True)
    for path in psf_paths:
        c=json.loads(Path(path).read_text())['cluster']
        shutil.copyfile(path,raw/f'{c}.json')
        parent=Path(path).parent
        for f in parent.glob('tile_*_native_psf.npy'):
            shutil.copyfile(f,out/f'{c}_{f.name}')
        for name in ('psf_montage.png','psf_spatial.csv'):
            if (parent/name).exists():
                shutil.copyfile(parent/name,out/f'{c}_{name}')
    (out/'summary.json').write_text(
        json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    with (out/'per_tile.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(result['per_tile'][0]))
        writer.writeheader();writer.writerows(result['per_tile'])
    x=np.arange(9)
    fig,ax=plt.subplots(figsize=(10,4.8))
    for c in CLUSTERS:
        rows=[r for r in result['per_tile'] if r['cluster']==c]
        ax.plot(x,[r['r50_px'] if r['r50_px'] is not None else np.nan
                   for r in rows],marker='o',linewidth=1.8,label=NAMES[c])
    ax.set(xlabel='Predefined spatial tile index (center anchor = 4)',
           ylabel='Empirical native-pixel half-light radius r50 (px)',
           title='Observed HST PSF shape variation · native-pixel diagnostic')
    ax.set_xticks(x);ax.grid(alpha=.2);ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out/'psf_spatial_r50.png',dpi=165)
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(10,4.8))
    for c in CLUSTERS:
        rows=[r for r in result['per_tile'] if r['cluster']==c]
        ax.plot(x,[r['total_variation_from_center']
                   if r['total_variation_from_center'] is not None else np.nan
                   for r in rows],marker='o',linewidth=1.8,label=NAMES[c])
    ax.set(xlabel='Predefined spatial tile index (center anchor = 4)',
           ylabel='Total-variation distance to center-tile native PSF',
           title='Real HST empirical PSF mismatch (descriptive, uncalibrated)')
    ax.set_xticks(x);ax.set_ylim(bottom=0);ax.grid(alpha=.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out/'psf_spatial_mismatch.png',dpi=165)
    plt.close(fig)
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description='Render real HST PSF metrics/plots')
    p.add_argument('--psf-reports',nargs=3,required=True)
    p.add_argument('--plate-reports',nargs=3,required=True)
    p.add_argument('--out',required=True)
    a=p.parse_args(argv)
    result=render_spatial_psf_dashboard(a.psf_reports,a.plate_reports,a.out)
    print(json.dumps({'usable':result['usable_tiles'],'failed':result['failed_tiles'],
                      'out':str(Path(a.out).resolve())},indent=2))


if __name__=='__main__':
    main()
