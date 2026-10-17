"""Reproducible README plots from real HST cross-region GitHub artifacts.

This module never synthesizes observations. Inputs MUST be the three actual
cross_region.json reports produced by survey_cli; failed tiles are retained.
Plots summarize internal ACSGGCT image/catalog agreement, NOT Gaia validation.
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
NAMES={'ngc6397':'NGC 6397','ngc6752':'NGC 6752','ngc0104':'NGC 104'}
TILES=(0,1,2,3,5,6,7,8)


def _sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def _read_report(path):
    p=Path(path)
    d=json.loads(p.read_text(encoding='utf-8'))
    name=d.get('cluster')
    if name not in CLUSTERS:
        raise ValueError(f'Unknown cluster in {p}: {name!r}')
    if d.get('tile_count')!=9 or d.get('anchor_index')!=4:
        raise ValueError(f'Unexpected evaluation design in {p}')
    if d.get('tile_size')!=512 or d.get('grid')!=3:
        raise ValueError(f'Inconsistent 512-pixel 3x3 protocol in {p}')
    if d.get('independent_external_reference') is not False:
        raise ValueError(f'Unexpected external-reference provenance in {p}')
    if d.get('science_calibrated') is not False:
        raise ValueError('These plots cannot certify science-grade astrometry')
    tiles=d.get('tiles',[])
    if sorted(t['index'] for t in tiles)!=list(range(9)):
        raise ValueError(f'Incomplete or duplicate tiles in {p}')
    if d.get('transfer_tiles_total')!=8:
        raise ValueError(f'Expected 8 held-out regions in {p}')
    if d.get('anchor_status')!='HELDOUT_QC_PASS':
        raise ValueError(f'Anchor failed: {name}')
    if not d.get('heldout_tiles_physically_nonoverlapping',False):
        raise ValueError(f'Tiles not confirmed disjoint: {name}')
    return d


def summarize(paths,*,source_run_id=None,source_commit=None):
    """Read all regions, including failures; never silently drop bad tiles."""
    reports={}
    sources={}
    for path in paths:
        report=_read_report(path)
        name=report['cluster']
        if name in reports:
            raise ValueError(f'Duplicate report for {name}')
        reports[name]=report
        sources[name]={'sha256':_sha256(path)}
    if set(reports)!=set(CLUSTERS):
        raise ValueError(f'Need exactly {CLUSTERS}; got {sorted(reports)}')
    rows=[]
    groups=[]
    for name in CLUSTERS:
        d=reports[name]
        for tile in sorted(d['tiles'],key=lambda t:t['index']):
            if tile['index']==4:
                continue
            q=tile['fixed_anchor_transfer']
            if q.get('refitted') is not False:
                raise ValueError('Held-out transform was refitted')
            row={'cluster':name,'tile_index':tile['index'],
                 'status':q.get('status','UNKNOWN'),
                 'evaluated':int(q.get('evaluated',0)),
                 'matches':int(q.get('matches',0)),
                 'match_fraction':q.get('match_fraction'),
                 'median_residual_px':q.get('median_residual_px'),
                 'p95_residual_px':q.get('p95_residual_px'),
                 'match_over_null_ratio':q.get('match_over_null_ratio')}
            rows.append(row)
        subset=[r for r in rows if r['cluster']==name]
        if sum(r['status']=='TRANSFER_QC_PASS' for r in subset)!=d['transfer_tiles_pass']:
            raise ValueError(f'Pass-count disagreement for {name}')
        total=sum(r['evaluated'] for r in subset)
        matched=sum(r['matches'] for r in subset)
        med=[r['median_residual_px'] for r in subset
             if r['median_residual_px'] is not None]
        groups.append({'cluster':name,'catalog_rows':int(d['catalog_rows']),
                       'quality_selected_catalog':int(d['quality_selected_catalog']),
                       'anchor_status':d['anchor_status'],
                       'passed_tiles':d['transfer_tiles_pass'],
                       'evaluated_tiles':8,
                       'evaluated_reference_candidates':total,
                       'matched_reference_candidates':matched,
                       'pooled_match_fraction':matched/total if total else None,
                       'median_of_tile_median_residual_px':float(np.median(med)) if med else None,
                       'input_sha256':sources[name]['sha256']})
    return {'version':'0.11.0','dataset':'STScI ACSGGCT F606W mosaics',
            'kind':'INTERNAL_HST_CATALOG_SPATIAL_TRANSFER',
            'source_github_actions_run_id':str(source_run_id or ''),
            'source_commit':str(source_commit or ''),
            'clusters':groups,'per_tile':rows,
            'heldout_regions_total':24,
            'heldout_regions_passed':sum(r['status']=='TRANSFER_QC_PASS' for r in rows),
            'not_external_gaia_validation':True,
            'not_absolute_astrometry':True,
            'not_photometric_completeness':True,
            'caveat':'Same-survey reference catalog and HST image; bright-source subset, internal one-to-one matching. Independent Gaia and calibrated flux results NOT reported.'}


def render_report(paths,out,*,source_run_id=None,source_commit=None):
    """Create a fixed evidence snapshot plus three README-ready PNG figures."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    result=summarize(paths,source_run_id=source_run_id,
                     source_commit=source_commit)
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False),
                                    encoding='utf-8')
    fields=list(result['per_tile'][0])
    with (out/'per_tile.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader();w.writerows(result['per_tile'])
    raw=out/'source_reports'
    raw.mkdir(exist_ok=True)
    epochs=out/'epoch_audits'
    epochs.mkdir(exist_ok=True)
    header_epoch_hints={}
    for path in paths:
        d=_read_report(path)
        shutil.copyfile(path,raw/f"{d['cluster']}.json")
        # Epoch provenance is an independent FITS-header audit, NOT a
        # verified effective mosaic observation epoch for Gaia matching.
        candidate=Path(path).with_name('epoch_audit.json')
        if candidate.is_file():
            ep=json.loads(candidate.read_text(encoding='utf-8'))
            if ep.get('effective_mosaic_epoch_verified') is not False:
                raise ValueError('Unverified FITS date cannot be marked verified')
            shutil.copyfile(candidate,epochs/f"{d['cluster']}.json")
            header_epoch_hints[d['cluster']]={
                'min_header_mjd':ep.get('min_header_mjd'),
                'max_header_mjd':ep.get('max_header_mjd'),
                'effective_mosaic_epoch_verified':False}
    result['header_epoch_hints']=header_epoch_hints
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False),
                                    encoding='utf-8')
    labels=[NAMES[c] for c in CLUSTERS]
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
    x=np.arange(8)
    def common(ax,title,ylabel):
        ax.set(title=title,xlabel='Disjoint held-out tile index',ylabel=ylabel)
        ax.set_xticks(x,labels=[str(t) for t in TILES])
        ax.grid(alpha=.22)
        ax.legend(loc='best',frameon=False,fontsize=9)

    # Figure 1: actual per-tile median matched-source positional residual.
    fig,ax=plt.subplots(figsize=(9.5,4.6))
    for j,c in enumerate(CLUSTERS):
        r=[a for a in result['per_tile'] if a['cluster']==c]
        y=[a['median_residual_px'] for a in r]
        ax.plot(x,y,'o-',label=NAMES[c],linewidth=1.8,markersize=5,
                color=colors[j])
    common(ax,'HST/ACSGGCT frozen-transform transfer (internal)',
           'Median positional residual (pixels)')
    ax.set_ylim(.35,.45)  # Explicit zoom: emphasize real 0.38-0.43 px differences.
    ax.set_title('HST/ACSGGCT frozen-transform transfer (zoomed y-axis)')
    fig.tight_layout()
    fig.savefig(out/'spatial_residuals.png',dpi=165)
    plt.close(fig)

    # Figure 2: matched bright-reference fraction, including failed regions.
    fig,ax=plt.subplots(figsize=(9.5,4.6))
    for j,c in enumerate(CLUSTERS):
        r=[a for a in result['per_tile'] if a['cluster']==c]
        y=[a['match_fraction'] if a['match_fraction'] is not None else np.nan for a in r]
        ax.plot(x,y,'o-',label=NAMES[c],linewidth=1.8,markersize=5,
                color=colors[j])
    common(ax,'One-to-one matching in eight held-out HST regions',
           'Matched selected reference fraction')
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylim(.965,1.005)  # Explicit zoom: all observed rates are >=97%.
    ax.set_title('One-to-one HST reference matching (zoomed y-axis)')
    fig.tight_layout()
    fig.savefig(out/'spatial_match_fraction.png',dpi=165)
    plt.close(fig)

    # Figure 3: survey source counts and explicit quality-selection counts.
    fig,ax=plt.subplots(figsize=(8.7,4.5))
    xp=np.arange(3)
    catalog=np.array([g['catalog_rows'] for g in result['clusters']])
    quality=np.array([g['quality_selected_catalog'] for g in result['clusters']])
    ax.bar(xp-.18,catalog,width=.35,label='All survey catalog rows')
    ax.bar(xp+.18,quality,width=.35,label='Quality-selected sources')
    ax.set_xticks(xp,labels=labels)
    ax.set(ylabel='Number of sources',
           title='Real ACSGGCT survey catalogs (not detector predictions)')
    ax.legend(frameon=False,fontsize=9)
    ax.grid(axis='y',alpha=.2)
    for xx,yy in zip(xp-.18,catalog):
        ax.text(xx,yy,f'{yy:,}',ha='center',va='bottom',fontsize=8)
    for xx,yy in zip(xp+.18,quality):
        ax.text(xx,yy,f'{yy:,}',ha='center',va='bottom',fontsize=8)
    ax.set_ylim(0,max(catalog)*1.17)
    fig.tight_layout()
    fig.savefig(out/'catalog_source_counts.png',dpi=165)
    plt.close(fig)
    return result


def main(argv=None):
    p=argparse.ArgumentParser(
        description='Create verifiable HST experiment plots and JSON/CSV')
    p.add_argument('reports',nargs='+',help='Three cross_region.json reports')
    p.add_argument('--out',required=True)
    p.add_argument('--source-run-id')
    p.add_argument('--source-commit')
    args=p.parse_args(argv)
    result=render_report(args.reports,args.out,
                         source_run_id=args.source_run_id,
                         source_commit=args.source_commit)
    print(json.dumps({'clusters':[g['cluster'] for g in result['clusters']],
         'passed':result['heldout_regions_passed'],
         'total':result['heldout_regions_total'],
         'out':str(Path(args.out).resolve())},indent=2))


if __name__=='__main__':
    main()
