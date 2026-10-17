"""Three-cluster, real HST expanded artificial-star results and README plots.

Every designed band and tile is present even when selected sites are sparse.
Counts are of planted trials, while three fluxes reuse the same pixel sites.
These are descriptive performance results, NOT true astrophysical completeness.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import numpy as np
from .crowding_injections import BANDS, CLUSTERS, TILES, FLUX_RMS_MULTIPLIERS
from .real_hst import sha256_file

LABELS={'ngc6397':'NGC 6397','ngc6752':'NGC 6752','ngc0104':'NGC 104'}
SHORT={'near_4_8':'Near 4–8 px','middle_8_12':'Middle 8–12 px',
       'outer_12_16':'Outer 12–16 px'}


def validate(path):
    path=Path(path)
    d=json.loads(path.read_text(encoding='utf-8'))
    if (d.get('cluster') not in CLUSTERS
        or d.get('kind')!='REAL_HST_EXPANDED_PREDECLARED_CROWDING_INJECTIONS'
        or d.get('all_injection_sites_chosen_before_outcomes') is not True
        or d.get('full_test_set_retained') is not True
        or d.get('all_models_fitted_on_noninjected_training_stars_only') is not True
        or d.get('stellar_population_completeness_estimated') is not False
        or d.get('injected_poisson_noise') is not False
        or d.get('mosaic_noise_physically_calibrated') is not False
        or d.get('forced_flux_uses_true_injection_position') is not True
        or d.get('bands')!={k:list(v) for k,v in BANDS.items()}):
        raise ValueError(f'Unverified/incompatible injection provenance: {path}')
    trials=d.get('trials',[])
    if len(trials)!=d['injection_attempts']:
        raise ValueError('Incomplete injected trial records')
    if sum(bool(r['detected']) for r in trials)!=d['recovered']:
        raise ValueError('Reported recovery count mismatch')
    seen=set();sites={}
    for r in trials:
        if (r['cluster']!=d['cluster'] or r['tile_index'] not in TILES
            or r['band'] not in BANDS
            or r['flux_over_tile_robust_pixel_rms'] not in FLUX_RMS_MULTIPLIERS):
            raise ValueError('Unexpected trial group')
        key=(r['tile_index'],r['band'],r['site_index'],
             r['flux_over_tile_robust_pixel_rms'])
        if key in seen:raise ValueError('Duplicate experiment trial')
        seen.add(key)
        site=key[:3]
        xy=(r['x'],r['y'])
        if site in sites and xy!=sites[site]:
            raise ValueError('Position moved between injected flux levels')
        sites[site]=xy
        lo,hi=BANDS[r['band']]
        v=float(r['nearest_original_detected_peak_px'])
        if not lo<=v<hi:
            raise ValueError('Trial not in declared detected-peak-distance band')
    if (len(sites)!=d['distinct_selected_injection_positions']
        or len(trials)!=len(sites)*len(FLUX_RMS_MULTIPLIERS)):
        raise ValueError('Some original sites or flux levels were cherry-picked')
    if {t['tile'] for t in d['tiles']}!=set(TILES):
        raise ValueError('Must retain all eight spatial tiles, including failures')
    return d


def aggregate(paths,out_dir,*,run_id=None,source_commit=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    reports={};input_files={}
    for path in paths:
        d=validate(path)
        c=d['cluster']
        if c in reports:raise ValueError('Duplicate HST cluster')
        reports[c]=d;input_files[c]=Path(path)
    if set(reports)!=set(CLUSTERS):
        raise ValueError('All three genuine telescope fields required')
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    raw=out/'source_reports';raw.mkdir(exist_ok=True)
    tile_rows=[];full=[];groups=[];clusters=[]
    for cluster in CLUSTERS:
        d=reports[cluster]
        shutil.copyfile(input_files[cluster],raw/f'{cluster}.json')
        full.extend(d['trials'])
        for tile in d['tiles']:
            tile_rows.append({'cluster':cluster,'tile_index':tile['tile'],
                              'status':tile['status'],
                              **{k+'_site_count':int(tile['bins'][k]) for k in BANDS}})
        for band in BANDS:
            for f in FLUX_RMS_MULTIPLIERS:
                trials=[r for r in d['trials'] if r['band']==band and
                        r['flux_over_tile_robust_pixel_rms']==f]
                n=len(trials)
                recovered=sum(r['detected'] for r in trials)
                groups.append({'cluster':cluster,'band':band,'flux_over_tile_rms':f,
                    'site_count':n,'attempted':n,'recovered':recovered,
                    'recovery_fraction':recovered/n if n else None,
                    'median_forced_local_abs_relative_error':float(np.median([
                        r['forced_local_abs_rel_error'] for r in trials])) if n else None,
                    'median_forced_center_abs_relative_error':float(np.median([
                        r['forced_center_abs_rel_error'] for r in trials])) if n else None,
                    'median_of_paired_local_minus_center_error':float(np.median([
                        r['forced_local_abs_rel_error']-r['forced_center_abs_rel_error']
                        for r in trials])) if n else None,
                    'median_detected_peak_position_error_px':float(np.median([
                        r['detector_peak_position_error_px'] for r in trials
                        if r['detected']])) if recovered else None})
        clusters.append({'cluster':cluster,'injection_trials':d['injection_attempts'],
                         'distinct_injection_sites':d['distinct_selected_injection_positions'],
                         'recovered':d['recovered'],
                         'recovery_fraction':d['recovered']/d['injection_attempts'] if d['injection_attempts'] else None,
                         'source_report_sha256':sha256_file(input_files[cluster]),
                         'original_fits_sha256':d['fits_sha256'],
                         'original_prf_report_sha256':d['prf_sha256'],
                         'lineage_audit':d['lineage_audit']})
    report={'version':'0.17.0','kind':'THREE_CLUSTER_EXPANDED_REAL_HST_PLANTED_STAR_BENCHMARK',
        'source_github_actions_run_id':str(run_id or ''),
        'source_commit':str(source_commit or ''),
        'clusters':clusters,'groups':groups,'original_tiles':tile_rows,
        'total_trials':len(full),
        'total_unique_pixel_sites':sum(x['distinct_injection_sites'] for x in clusters),
        'total_recovered':sum(x['recovered'] for x in clusters),
        'all_eight_tiles_per_cluster_retained':True,
        'site_geometry_chosen_before_injection_outcomes':True,
        'same_three_flux_levels_per_site':True,
        'instrument_noise_calibrated':False,
        'photometric_zeropoint_calibrated':False,
        'population_level_completeness':False,
        'mosaic_exposure_lineage_verified':False,
        'notes':'Crowding proxy = distance to detected 5-sigma peaks in real HST image, NOT unresolved-star density. Forced amplitude uses known planted site.'}
    (out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False),
                                    encoding='utf-8')
    for name,data in [('per_trial.csv',full),('groups.csv',groups),('tile_sites.csv',tile_rows)]:
        with (out/name).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]))
            w.writeheader();w.writerows(data)
    # Plot 1: full per-class recovery vs flux for all 3 HST fields.
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig,axs=plt.subplots(1,3,figsize=(14,4.5),sharey=True)
    for i,c in enumerate(CLUSTERS):
        ax=axs[i]
        for j,band in enumerate(BANDS):
            rr=[g for g in groups if g['cluster']==c and g['band']==band]
            ax.plot(FLUX_RMS_MULTIPLIERS,
                    [r['recovery_fraction'] if r['recovery_fraction'] is not None
                     else np.nan for r in rr],'o-',color=colors[j],
                    label=SHORT[band])
        ax.set(title=LABELS[c],xlabel='Injected flux / original tile robust RMS')
        ax.set_ylim(0,1.05)
        ax.grid(alpha=.22)
    axs[0].set_ylabel('Blind planted-peak recovery fraction')
    axs[0].legend(frameon=False,fontsize=8)
    fig.suptitle('Genuine HST F606W backgrounds · original detected-peak-distance bins')
    fig.tight_layout()
    fig.savefig(out/'recovery_crowding_distance.png',dpi=160)
    plt.close(fig)
    # Plot 2: error difference: positive means local is worse.
    fig,axs=plt.subplots(1,3,figsize=(14,4.5))
    for i,c in enumerate(CLUSTERS):
        ax=axs[i]
        for j,band in enumerate(BANDS):
            rr=[g for g in groups if g['cluster']==c and g['band']==band]
            ax.plot(FLUX_RMS_MULTIPLIERS,
                    [r['median_of_paired_local_minus_center_error']
                     if r['median_of_paired_local_minus_center_error'] is not None
                     else np.nan for r in rr],'o-',label=SHORT[band],color=colors[j])
        ax.axhline(0,color='black',linestyle='--',linewidth=1)
        ax.set(title=LABELS[c],xlabel='Injected flux / original tile robust RMS')
        ax.grid(alpha=.22)
    axs[0].set_ylabel('Paired |local error| − |center error| (negative favors local)')
    axs[0].legend(frameon=False,fontsize=8)
    fig.suptitle('Known-position forced native-pixel PRF errors, no flux calibration')
    fig.tight_layout()
    fig.savefig(out/'paired_flux_crowding.png',dpi=160)
    plt.close(fig)
    # Plot 3: unambiguous distinct site counts, not repeated flux experiments.
    fig,ax=plt.subplots(figsize=(9,4.8))
    x=np.arange(3)
    for j,b in enumerate(BANDS):
        vals=[sum(t[b+'_site_count'] for t in tile_rows if t['cluster']==c)
              for c in CLUSTERS]
        ax.bar(x+(j-1)*.25,vals,width=.24,label=SHORT[b],color=colors[j])
    ax.set_xticks(x,[LABELS[c] for c in CLUSTERS])
    ax.set(ylabel='Distinct predeclared physical injection sites',
           title='Actual valid positions per HST cluster (not threefold flux replicates)')
    ax.legend(frameon=False)
    ax.grid(axis='y',alpha=.2)
    fig.tight_layout()
    fig.savefig(out/'injection_site_coverage.png',dpi=160)
    plt.close(fig)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description='Plot expanded real-HST planted-star recoveries')
    p.add_argument('reports',nargs='+')
    p.add_argument('--out',required=True)
    p.add_argument('--run-id')
    p.add_argument('--source-commit')
    a=p.parse_args(argv)
    d=aggregate(a.reports,a.out,run_id=a.run_id,source_commit=a.source_commit)
    print(json.dumps({'trials':d['total_trials'],
        'distinct_sites':d['total_unique_pixel_sites'],
        'recovered':d['total_recovered'],
        'clusters':d['clusters'],
        'calibrated_photometry':False},indent=2))


if __name__=='__main__':
    main()
