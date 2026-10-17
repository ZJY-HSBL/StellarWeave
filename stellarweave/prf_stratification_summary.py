"""Real HST PRF post-hoc stratification: full evidence and noncausal plots.

Predeclared features and two equal-frequency groups per cluster. Group
thresholds are medians of frozen TEST sources, computed after v0.12 model
scores were frozen. This is exploratory description, NOT a new test-set
optimization, causal explanation, or evidence of calibrated photometry.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import numpy as np

from .prf_stratification import FEATURES, CLUSTERS, OUTER
from .real_hst import sha256_file

NAMES = {'ngc6397':'NGC 6397','ngc6752':'NGC 6752','ngc0104':'NGC 104'}
MIN_PER_TILE_PER_GROUP = 2


def validate_frozen_reports(strat_paths, prf_paths):
    """Fail closed if any v0.12 source is missing, changed, or cherry-picked."""
    strat, prf, originals = {}, {}, {}
    for path in strat_paths:
        p = Path(path)
        d = json.loads(p.read_text(encoding='utf-8'))
        name = d.get('cluster')
        if name not in CLUSTERS or name in strat:
            raise ValueError('Missing/duplicate/unrecognized cluster in diagnostics')
        if (d.get('kind')!='FROZEN_HST_PRF_CROWDING_BRIGHTNESS_EXPLORATORY'
                or d.get('models_retrained') is not False
                or d.get('stars_removed_for_stratification') is not False
                or d.get('frozen_prf_scores_unchanged') is not True
                or d.get('analysis_is_posthoc_exploratory') is not True
                or d.get('science_calibrated') is not False
                or d.get('tile_count') != 8):
            raise ValueError('Invalid frozen exploratory diagnostic provenance')
        strat[name] = d
        originals[name] = p
    for path in prf_paths:
        p = Path(path)
        d = json.loads(p.read_text(encoding='utf-8'))
        name = d.get('cluster')
        if name not in CLUSTERS or name in prf:
            raise ValueError('Missing/duplicate/unrecognized original PRF cluster')
        prf[name] = (d,p)
    if set(strat) != set(CLUSTERS) or set(prf) != set(CLUSTERS):
        raise ValueError('Require all three HST clusters and original v0.12 reports')
    for name in CLUSTERS:
        d = strat[name]
        old, path = prf[name]
        if (d.get('original_prf_report_sha256') != sha256_file(path)
                or d.get('fits_sha256') != old.get('fits_sha256')
                or d.get('plate_sha256') != old.get('plate_sha256')
                or d.get('original_prf_kind') != old.get('kind')):
            raise ValueError(f'Frozen v0.12 input provenance mismatch: {name}')
        expected = {(int(s['tile_index']),int(s['holdout_index'])):s
                    for s in old['heldout_star_scores']
                    if int(s['tile_index']) in OUTER}
        rows = d['heldout_star_scores_with_proxies']
        if (len(expected) != len(rows)
                or len(rows) != d.get('all_original_outer_heldout_stars')):
            raise ValueError(f'Frozen heldout star count mismatch: {name}')
        seen = set()
        for s in rows:
            key = (int(s['tile_index']),int(s['holdout_index']))
            if key in seen or key not in expected:
                raise ValueError(f'Duplicate or unexpected test star: {name}')
            seen.add(key)
            old_star = expected[key]
            for field in ('tile_x','tile_y','local_l1','anchor_l1',
                          'gaussian_l1','paired_gain_pct'):
                if (not np.isfinite(float(s[field]))
                        or abs(float(s[field])-float(old_star[field]))>1e-10):
                    raise ValueError(f'Frozen PRF score modified: {name}/{key}')
            for feature in FEATURES:
                v = s[feature]
                if v is not None and not np.isfinite(float(v)):
                    raise ValueError('Non-finite diagnostic feature')
        if seen != set(expected):
            raise ValueError(f'Original frozen test set was cherry-picked: {name}')
        if {int(t['tile_index']) for t in d['tiles']} != set(OUTER):
            raise ValueError('All eight outer tiles must be present')
    return strat, originals


def _tile_block_ci(tile_values, *, seed=2027, samples=2000):
    vals = np.asarray(tile_values, dtype=float)
    if len(vals) < 4:
        return None
    rng = np.random.default_rng(seed)
    draw = rng.choice(vals, size=(samples,len(vals)),replace=True)
    return [float(x) for x in np.percentile(
        np.median(draw,axis=1),[2.5,97.5])]


def summarize_groups(strat):
    """Predeclared four image proxies, median split within each cluster.

    All thresholds are post-hoc descriptive and may vary by cluster. No
    model or feature is chosen for best performance on held-out sources.
    """
    groups, thresholds = [], []
    for cluster in CLUSTERS:
        stars = strat[cluster]['heldout_star_scores_with_proxies']
        for feature in FEATURES:
            valid = [s for s in stars if s[feature] is not None]
            threshold = float(np.median([float(s[feature]) for s in valid])) if valid else None
            thresholds.append({
                'cluster':cluster,'feature':feature,'median_split_threshold':threshold,
                'nonmissing':len(valid),'missing':len(stars)-len(valid),
                'threshold_is_posthoc_from_frozen_test_stars':True})
            for group in ('low','high'):
                subset = [s for s in valid if (
                    float(s[feature]) <= threshold if group=='low'
                    else float(s[feature]) > threshold)] if threshold is not None else []
                gains = [float(s['paired_gain_pct']) for s in subset]
                delta = [float(s['local_l1'])-float(s['anchor_l1']) for s in subset]
                tile_medians = []
                for tile in OUTER:
                    one = [float(s['paired_gain_pct']) for s in subset
                           if int(s['tile_index'])==tile]
                    if len(one) >= MIN_PER_TILE_PER_GROUP:
                        tile_medians.append(float(np.median(one)))
                groups.append({
                    'cluster':cluster,'feature':feature,'group':group,
                    'n':len(subset),
                    'tiles_with_at_least_2_stars':len(tile_medians),
                    'median_star_paired_gain_pct':float(np.median(gains)) if gains else None,
                    'median_star_local_minus_anchor_l1':float(np.median(delta)) if delta else None,
                    'fraction_stars_local_better':float(np.mean(np.asarray(delta)<0)) if delta else None,
                    'median_of_tile_median_gains_pct':float(np.median(tile_medians)) if tile_medians else None,
                    'tile_bootstrap_95ci_median_gain_pct':_tile_block_ci(
                        tile_medians,seed=2027+CLUSTERS.index(cluster)*31+
                        FEATURES.index(feature)*3+(group=='high')),
                })
    return groups, thresholds


def render_diagnostics(strat_paths, prf_paths, out_dir, *,
                       run_id=None, source_commit=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    reports, source_files = validate_frozen_reports(strat_paths,prf_paths)
    out = Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    source = out/'source_reports'
    source.mkdir(exist_ok=True)
    for name in CLUSTERS:
        shutil.copyfile(source_files[name],source/f'{name}.json')
    groups, thresholds = summarize_groups(reports)
    stars = [dict(cluster=name,**s)
             for name in CLUSTERS
             for s in reports[name]['heldout_star_scores_with_proxies']]
    summary = {
        'version':'0.13.0',
        'kind':'THREE_CLUSTER_FROZEN_HST_PRF_STRATIFIED_EXPLORATORY',
        'source_github_actions_run_id':str(run_id or ''),
        'source_commit':str(source_commit or ''),
        'clusters':list(CLUSTERS),
        'total_original_outer_heldout_stars':len(stars),
        'all_frozen_test_stars_retained':True,
        'all_24_outer_tiles_retained':True,
        'models_retrained':False,
        'science_calibrated':False,
        'causal_inference_supported':False,
        'stratification':'Post-hoc median splits of four predeclared pixel-derived proxies within each cluster, on unchanged v0.12 test stars',
        'primary_metric':'v0.12 held-out native-pixel PRF shape L1, not flux accuracy',
        'group_statistics':'star-pooled median and median of within-tile medians; exploratory tile-bootstrap CI',
        'group_definitions':thresholds,
        'groups':groups,
        'input_provenance':[
            {'cluster':name,'stratification_sha256':sha256_file(source_files[name]),
             'original_prf_sha256':reports[name]['original_prf_report_sha256'],
             'fits_sha256':reports[name]['fits_sha256']}
            for name in CLUSTERS],
        'limitations':[
            'Post-hoc subgroup associations cannot explain causality',
            'Six-sigma peak-neighbor distances miss unresolved and faint sources',
            'Candidates are deliberately bright/isolated, not a completeness sample',
            'Peak-over-tile-RMS is not calibrated SNR or flux',
            'Subgroups and image backgrounds differ across clusters',
            'No Gaia absolute astrometry or calibrated photometry',
        ],
    }
    (out/'summary.json').write_text(
        json.dumps(summary,indent=2,allow_nan=False),encoding='utf-8')
    for filename, data in [('per_star.csv',stars),('groups.csv',groups),
                           ('thresholds.csv',thresholds)]:
        with (out/filename).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)

    # Four feature panels, all clusters and both strata, including negative results.
    fig,axs=plt.subplots(2,2,figsize=(12,9))
    short={
        'peak_over_tile_rms':'Peak / robust tile RMS (not source SNR)',
        'nearest_detected_peak_px':'Distance to nearest 6-sigma detected peak (px)',
        'neighbors_within_25px':'Other 6-sigma peaks within 25 px',
        'annulus_mad_over_tile_rms':'5–10 px annulus texture / tile RMS'}
    x=np.arange(3)
    for ax,feature in zip(axs.flat,FEATURES):
        for i,group in enumerate(('low','high')):
            vals=[]
            for cluster in CLUSTERS:
                r=next(g for g in groups if g['cluster']==cluster
                       and g['feature']==feature and g['group']==group)
                vals.append(r['median_star_paired_gain_pct'] if
                            r['median_star_paired_gain_pct'] is not None else np.nan)
            ax.bar(x+(i-.5)*.34,vals,width=.33,label=group+' half')
        ax.axhline(0,color='black',linewidth=.9,linestyle='--')
        ax.set_title(short[feature],fontsize=10)
        ax.set_xticks(x,[NAMES[c] for c in CLUSTERS])
        ax.set_ylabel('Median held-out PRF paired shape gain (%)')
        ax.grid(axis='y',alpha=.2)
    axs[0,0].legend(frameon=False)
    fig.suptitle('Post-hoc frozen-test diagnostics — noncausal, not photometric accuracy')
    fig.tight_layout()
    fig.savefig(out/'stratified_prf_gains.png',dpi=165)
    plt.close(fig)

    # Spatial failure map for NGC 104: the same held-out sources, no refitting.
    ngc = reports['ngc0104']
    fig,axs=plt.subplots(3,3,figsize=(9.5,9))
    cmap='coolwarm_r'
    norm=TwoSlopeNorm(vmin=-50,vcenter=0,vmax=50)
    for idx,ax in enumerate(axs.flat):
        if idx==4:
            ax.text(.5,.5,'Frozen center PRF\n(training anchor)',
                    ha='center',va='center',transform=ax.transAxes,fontsize=10)
        else:
            subset=[s for s in ngc['heldout_star_scores_with_proxies']
                    if int(s['tile_index'])==idx]
            if subset:
                ax.scatter([s['tile_x'] for s in subset],
                           [s['tile_y'] for s in subset],
                           c=[s['paired_gain_pct'] for s in subset],
                           cmap=cmap,norm=norm,s=23,edgecolors='none')
        ax.set_xlim(0,512)
        ax.set_ylim(0,512)
        ax.set_aspect('equal')
        ax.set_title(f'Tile {idx}',fontsize=10)
        ax.tick_params(labelsize=7)
    sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap)
    fig.colorbar(sm,ax=axs.ravel().tolist(),shrink=.7,
                 label='Local PRF shape gain vs center (%) [clipped to ±50]')
    fig.suptitle('NGC 104: every frozen held-out candidate by tile (not a sky map)')
    fig.savefig(out/'ngc0104_failure_map.png',dpi=165,bbox_inches='tight')
    plt.close(fig)

    # Raw NGC 104 scatter, including sign reversals and possible confounding.
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.5))
    ns=ngc['heldout_star_scores_with_proxies']
    for ax,feature in zip(axs,('peak_over_tile_rms','nearest_detected_peak_px')):
        values=np.array([s[feature] for s in ns],float)
        gain=np.array([s['paired_gain_pct'] for s in ns],float)
        ax.scatter(values,gain,s=16,alpha=.65)
        ax.axhline(0,color='black',linestyle='--',linewidth=1)
        ax.set_xscale('log')
        ax.set_xlabel(short[feature])
        ax.set_ylabel('Paired shape gain vs center (%)')
        ax.grid(alpha=.2)
    fig.suptitle('NGC 104: frozen-test brightness and detected-neighbor proxies')
    fig.tight_layout()
    fig.savefig(out/'ngc0104_proxy_scatter.png',dpi=165)
    plt.close(fig)
    return summary


def main(argv=None):
    p=argparse.ArgumentParser(description='Plot full frozen HST PRF strata')
    p.add_argument('--strat-reports',nargs=3,required=True)
    p.add_argument('--prf-reports',nargs=3,required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--run-id')
    p.add_argument('--source-commit')
    a=p.parse_args(argv)
    result=render_diagnostics(
        a.strat_reports,a.prf_reports,a.out,
        run_id=a.run_id,source_commit=a.source_commit)
    print(json.dumps({
        'stars':result['total_original_outer_heldout_stars'],
        'strata':len(result['groups']),
        'frozen':True,'calibrated_photometry':False},indent=2))


if __name__=='__main__':
    main()
