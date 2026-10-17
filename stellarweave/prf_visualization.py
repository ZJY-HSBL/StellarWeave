"""Aggregate three genuine HST spatial-PRF holdout reports into README figures.

This code never fabricates measurements or drops failed tiles. All input
reports and per-star scores are retained. PRF shape L1 is NOT flux accuracy.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np

CLUSTERS=('ngc6397','ngc6752','ngc0104')
LABELS={'ngc6397':'NGC 6397','ngc6752':'NGC 6752','ngc0104':'NGC 104'}
TILE_INDICES=(0,1,2,3,5,6,7,8)


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for buf in iter(lambda:f.read(1024*1024),b''):
            h.update(buf)
    return h.hexdigest()


def validate(path):
    p=Path(path)
    d=json.loads(p.read_text(encoding='utf-8'))
    if d.get('cluster') not in CLUSTERS:
        raise ValueError(f'Unrecognized cluster: {p}')
    if (d.get('kind')!='REAL_HST_SPATIAL_NATIVE_PIXEL_PRF_HOLDOUT'
            or d.get('science_calibrated') is not False
            or d.get('oversampled_instrument_epsf') is not False
            or d.get('calibrated_photometric_accuracy') is not False
            or d.get('no_heldout_star_used_to_fit_a_scored_template') is not True):
        raise ValueError(f'Invalid PRF provenance or scientific claim: {p}')
    tiles=d.get('tiles',[])
    if [r.get('index') for r in tiles]!=list(range(9)):
        raise ValueError('All nine tiles, including failures, must be reported')
    n=sum(r.get('status')=='EVALUATED' for r in tiles if r['index']!=4)
    if n!=d.get('evaluated_heldout_tiles'):
        raise ValueError('Tile evaluation count disagrees with report')
    if len(d.get('heldout_star_scores',[]))!=d.get('total_heldout_stars'):
        raise ValueError('Per-star holdout evidence is incomplete')
    return d


def render_prf_summary(paths,out_dir,*,run_id=None,commit=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    reports={}
    paths_by_name={}
    for path in paths:
        d=validate(path)
        if d['cluster'] in reports:
            raise ValueError(f"Duplicate PRF cluster {d['cluster']}")
        reports[d['cluster']]=d
        paths_by_name[d['cluster']]=Path(path)
    if set(reports)!=set(CLUSTERS):
        raise ValueError(f'Expected three clusters, got {sorted(reports)}')
    out=Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    raw=out/'source_reports'
    raw.mkdir(exist_ok=True)
    rows=[]
    stars=[]
    groups=[]
    for name in CLUSTERS:
        d=reports[name]
        p=paths_by_name[name]
        shutil.copyfile(p,raw/f'{name}.json')
        for optional,suffix in [('heldout_stars.csv','heldout_stars.csv'),
                                ('spatial_prf_templates.png','prf_templates.png')]:
            src=p.with_name(optional)
            if not src.is_file() and p.parent.name=='source_reports':
                # Archived evidence keeps reports under source_reports/ and
                # matching CSV/PNG sidecars under the parent snapshot.
                src=p.parent.parent/f'{name}_{suffix}'
            if not src.is_file():
                raise FileNotFoundError(f'Missing required real evidence: {src}')
            shutil.copyfile(src,out/f'{name}_{suffix}')
        for tile in d['tiles']:
            rows.append({'cluster':name,'tile_index':tile['index'],
                         'status':tile['status'],
                         'train_stars':tile['train_stars'],
                         'heldout_stars':tile['heldout_stars'],
                         'local_l1':tile['median_local_l1'],
                         'anchor_l1':tile['median_anchor_l1'],
                         'gaussian_l1':tile['median_gaussian_l1'],
                         'paired_gain_pct':tile['median_paired_gain_pct']})
        for item in d['heldout_star_scores']:
            if item['tile_index']!=4:
                stars.append(dict(cluster=name,**item))
        subset=[s for s in stars if s['cluster']==name]
        groups.append({'cluster':name,
                       'evaluated_heldout_tiles':d['evaluated_heldout_tiles'],
                       'possible_heldout_tiles':8,
                       'local_prf_better_tiles':d['local_prf_better_tiles'],
                       'median_tile_paired_gain_pct':d['median_tile_paired_gain_pct'],
                       'bootstrap_95ci_median_tile_gain_pct':d['bootstrap_95ci_median_tile_gain_pct'],
                       'evaluated_heldout_stars':len(subset),
                       'median_star_local_l1':float(np.median([s['local_l1'] for s in subset])) if subset else None,
                       'median_star_anchor_l1':float(np.median([s['anchor_l1'] for s in subset])) if subset else None,
                       'median_star_gaussian_l1':float(np.median([s['gaussian_l1'] for s in subset])) if subset else None,
                       'report_sha256':sha(p),
                       'fits_sha256':d['fits_sha256'],
                       'plate_sha256':d['plate_sha256']})
    summary={'version':'0.12.0',
             'kind':'REAL_HST_NATIVE_PIXEL_PRF_CROSSVALIDATION',
             'source_github_actions_run_id':str(run_id or ''),
             'source_commit':str(commit or ''),
             'clusters':groups,
             'heldout_tile_slots':24,
             'evaluated_heldout_tiles':sum(g['evaluated_heldout_tiles'] for g in groups),
             'all_failed_tiles_retained':True,
             'external_gaia_validation':False,
             'calibrated_photometry':False,
             'oversampled_instrument_epsf':False,
             'metric':'L1 difference between normalized star stamp and train-only PRF; lower is better',
             'limitation':'Native-pixel empirical PRF and bright isolated candidate stamps. Not a measured photometric improvement or science-grade ePSF.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False),
                                    encoding='utf-8')
    with (out/'per_tile.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    if stars:
        with (out/'per_star.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(stars[0]))
            w.writeheader();w.writerows(stars)

    # Plot 1: tile-level paired median improvement, including negative results.
    fig,ax=plt.subplots(figsize=(9.6,4.8))
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
    for j,name in enumerate(CLUSTERS):
        rr=[r for r in rows if r['cluster']==name and r['tile_index']!=4]
        vals=[r['paired_gain_pct'] if r['status']=='EVALUATED' else np.nan for r in rr]
        ax.plot(np.arange(8),vals,'o-',label=LABELS[name],
                linewidth=1.8,markersize=5,color=colors[j])
    ax.axhline(0,color='black',linestyle='--',linewidth=1,alpha=.7)
    ax.set_xticks(range(8),[str(t) for t in TILE_INDICES])
    ax.set(xlabel='Held-out tile index (center 4 excluded)',
           ylabel='Median paired shape-error gain vs center PRF (%)',
           title='Real HST spatial PRF transfer: positive = local lower L1')
    ax.grid(alpha=.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out/'prf_tile_gains.png',dpi=165)
    plt.close(fig)

    # Plot 2: median of per-star heldout shape errors, same test stars.
    fig,ax=plt.subplots(figsize=(9.2,4.9))
    x=np.arange(3)
    for k,(field,label) in enumerate((('median_star_local_l1','Local PRF'),
                                     ('median_star_anchor_l1','Frozen center PRF'),
                                     ('median_star_gaussian_l1','Gaussian FWHM 2.6'))):
        vals=[g[field] if g[field] is not None else np.nan for g in groups]
        ax.bar(x+(k-1)*.24,vals,width=.23,label=label)
    ax.set_xticks(x,[LABELS[c] for c in CLUSTERS])
    ax.set(ylabel='Median held-out normalized-stamp L1 (lower is better)',
           title='Real HST PRF model comparison: bright isolated candidates')
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.grid(axis='y',alpha=.2)
    ax.legend(frameon=False,fontsize=9)
    fig.tight_layout()
    fig.savefig(out/'prf_heldout_comparison.png',dpi=165)
    plt.close(fig)
    return summary


def main(argv=None):
    p=argparse.ArgumentParser(description='Plot real HST PRF holdout experiments')
    p.add_argument('reports',nargs='+')
    p.add_argument('--out',required=True)
    p.add_argument('--run-id')
    p.add_argument('--source-commit')
    a=p.parse_args(argv)
    result=render_prf_summary(a.reports,a.out,run_id=a.run_id,
                              commit=a.source_commit)
    print(json.dumps({'evaluated_heldout_tiles':result['evaluated_heldout_tiles'],
         'heldout_tile_slots':24,
         'clusters':result['clusters'],
         'science_calibrated':False},indent=2))


if __name__=='__main__':
    main()
