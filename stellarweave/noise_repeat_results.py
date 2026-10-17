"""Aggregate true three-cluster HST background repeated-noise robustness.

Do not pool trial repeats as independent observations. Report between-seed
min/max and paired per-seed flux differences rather than spurious 2520-star
confidence intervals. All individual outcomes, including nondetections,
are archived for reproduction.
"""
from __future__ import annotations
import argparse
import csv
import json
import shutil
from pathlib import Path
import numpy as np

from .real_hst import sha256_file
from .noise_repeats import CLUSTERS, REPLICATES, ADDED_MODES, _trial_key

NAMES={'ngc6397':'NGC 6397','ngc6752':'NGC 6752','ngc0104':'NGC 104'}


def validate_report(path):
    path=Path(path)
    d=json.loads(path.read_text(encoding='utf-8'))
    if (d.get('kind')!='REAL_HST_FROZEN_PLANTED_STARS_REPEATED_NOISE_SENSITIVITY'
        or d.get('cluster') not in CLUSTERS
        or d.get('replicate_seeds')!=list(REPLICATES)
        or d.get('modes')!=list(ADDED_MODES)
        or d.get('calibrated_instrument_noise') is not False
        or d.get('physical_poisson_source_noise_added') is not False
        or d.get('star_positions_flux_shapes_frozen') is not True
        or d.get('no_source_model_retrained_on_injected_stars') is not True):
        raise ValueError('Unverified noise repetition science provenance')
    rows=d.get('trials',[])
    n=d.get('original_sites')
    if not isinstance(n,int) or n<=0 or len(rows)!=n*len(REPLICATES)*len(ADDED_MODES):
        raise ValueError('Missing original artificial-star noise-repeat outcomes')
    source_ids={};groups={}
    for row in rows:
        if row['cluster']!=d['cluster'] or row['mode'] not in ADDED_MODES or row['replicate_seed_id'] not in REPLICATES:
            raise ValueError('Unrecognized trial mode/seed/cluster')
        key=_trial_key(row)
        stable=tuple(row[k] for k in (
            'x','y','injected_flux_image_units','source_holdout_index'))
        if key in source_ids and source_ids[key]!=stable:
            raise ValueError('Frozen planted-star identity or flux changed')
        source_ids[key]=stable
        pair=(row['mode'],row['replicate_seed_id'])
        by=groups.setdefault(pair,set())
        if key in by:raise ValueError('Duplicated artificial star in repeat')
        by.add(key)
        for f in ('local_forced_abs_rel_error','center_forced_abs_rel_error'):
            if not np.isfinite(float(row[f])) or float(row[f])<0:
                raise ValueError('Invalid flux estimation diagnostic')
        if row['preexisting_noise_peak_at_truth'] and row['detected']:
            raise ValueError('A preexisting noisy peak falsely counted as recovered')
    if len(source_ids)!=n or set(groups)!={
        (mode,seed) for mode in ADDED_MODES for seed in REPLICATES}:
        raise ValueError('Some frozen sites or noise modes were discarded')
    if any(keys!=set(source_ids) for keys in groups.values()):
        raise ValueError('Noise mode/seed used different injection locations')
    audits=d.get('noise_audits',[])
    if len(audits)!=d.get('repeated_noise_audits'):
        raise ValueError('Noise audits are missing')
    tile_ids={row['tile_index'] for row in rows}
    target={(t,s,m) for t in tile_ids for s in REPLICATES for m in ADDED_MODES}
    found=[(x['tile_index'],x['replicate_seed_id'],x['mode']) for x in audits]
    if set(found)!=target or len(found)!=len(target):
        raise ValueError('Missing or duplicated noise-field provenance')
    for audit in audits:
        actual=float(audit['added_noise_pixel_rms'])
        expected=float(audit['expected_pixel_rms'])
        if not np.isfinite(actual) or abs(actual-expected)>1e-6*max(1.,expected):
            raise ValueError('Mismatch in equal-variance noise control')
        if audit['mode']=='correlated_gaussian' and (
                audit['lag1_added_field_correlation'] is None or
                audit['lag1_added_field_correlation']<.2):
            raise ValueError('Noise labeled correlated is not spatially correlated')
    return d


def aggregate(paths,out_dir,*,run_id=None,source_commit=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    reports={};sources={}
    for path in paths:
        r=validate_report(path)
        c=r['cluster']
        if c in reports:raise ValueError('Duplicate cluster')
        reports[c]=r;sources[c]=Path(path)
    if set(reports)!=set(CLUSTERS):
        raise ValueError('Require all three original HST clusters')
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    raw=out/'source_reports';raw.mkdir(exist_ok=True)
    rows=[];groups=[];summary_clusters=[]
    for cluster in CLUSTERS:
        r=reports[cluster]
        shutil.copyfile(sources[cluster],raw/f'{cluster}.json')
        rows.extend(r['trials'])
        details=[]
        for seed in REPLICATES:
            for mode in ADDED_MODES:
                one=[x for x in r['trials'] if x['mode']==mode
                     and x['replicate_seed_id']==seed]
                n=len(one);recovered=sum(bool(x['detected']) for x in one)
                group={
                    'cluster':cluster,'seed':seed,'mode':mode,
                    'attempted':n,'recovered':recovered,
                    'recovery_fraction':recovered/n,
                    'preexisting_noise_peaks':sum(x['preexisting_noise_peak_at_truth']
                                                  for x in one),
                    'median_local_forced_abs_relative_error':float(np.median(
                        [x['local_forced_abs_rel_error'] for x in one])),
                    'median_center_forced_abs_relative_error':float(np.median(
                        [x['center_forced_abs_rel_error'] for x in one])),
                    'median_paired_local_minus_center_abs_relative_error':float(np.median(
                        [x['local_forced_abs_rel_error']-x['center_forced_abs_rel_error']
                         for x in one]))}
                groups.append(group);details.append(group)
        aggregate_per_mode=[]
        for mode in ADDED_MODES:
            g=[x for x in details if x['mode']==mode]
            rates=[x['recovery_fraction'] for x in g]
            local=[x['median_local_forced_abs_relative_error'] for x in g]
            center=[x['median_center_forced_abs_relative_error'] for x in g]
            aggregate_per_mode.append({
                'mode':mode,'seed_repetitions':len(g),
                'mean_seed_recovery_fraction':float(np.mean(rates)),
                'min_seed_recovery_fraction':float(min(rates)),
                'max_seed_recovery_fraction':float(max(rates)),
                'median_of_seed_local_forced_errors':float(np.median(local)),
                'median_of_seed_center_forced_errors':float(np.median(center)),
                'seed_median_local_flux_error_range':[float(min(local)),float(max(local))],
                'not_an_independent_source_population_confidence_interval':True})
        summary_clusters.append({
            'cluster':cluster,'frozen_injection_sites':r['original_sites'],
            'repeated_noise_trials':r['repeated_noise_trials'],
            'source_report_sha256':sha256_file(sources[cluster]),
            'fits_sha256':r['fits_sha256'],
            'v015_baseline_sha256':r['v015_sha256'],
            'per_mode':aggregate_per_mode,
            'header_noise_audit':r['observation_metadata']})
    report={
        'version':'0.16.0',
        'kind':'THREE_CLUSTER_REPEATED_FROZEN_HST_NOISE_SENSITIVITY',
        'github_actions_run_id':str(run_id or ''),
        'source_commit':str(source_commit or ''),
        'replicate_seeds':list(REPLICATES),
        'modes':list(ADDED_MODES),
        'clusters':summary_clusters,
        'per_seed_groups':groups,
        'total_original_sites':sum(c['frozen_injection_sites'] for c in summary_clusters),
        'total_added_noise_trials':len(rows),
        'all_failed_injections_retained':True,
        'single_fixed_original_HST_background_per_cluster':True,
        'noise_strength_equal_between_modes':True,
        'HST_noise_physically_calibrated':False,
        'source_shot_noise_simulated':False,
        'instrumental_noise_covariance_measured':False,
        'flux_error_uses_known_truth_coordinate':True,
        'seed_ranges_are_descriptive_not_confidence_intervals':True,
        'warning':'Five random seeds cannot establish observational source-population or detector-noise confidence intervals'}
    (out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False),
                                    encoding='utf-8')
    audits=[dict(cluster=cluster,**a)
            for cluster in CLUSTERS for a in reports[cluster]['noise_audits']]
    for name,values in [('per_trial.csv',rows),('per_seed.csv',groups),('noise_audits.csv',audits)]:
        with (out/name).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(values[0]))
            w.writeheader();w.writerows(values)

    colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
    for filename,field,ylabel,top in [
        ('recovery_across_seeds.png','recovery_fraction',
         'Fraction of frozen planted stars detected',1.06),
        ('forced_flux_across_seeds.png','median_local_forced_abs_relative_error',
         'Per-seed median abs. known-coordinate local PRF error',None)]:
        fig,axs=plt.subplots(1,3,figsize=(13.8,4.5),sharey=top is not None)
        for ci,cluster in enumerate(CLUSTERS):
            ax=axs[ci]
            for mi,mode in enumerate(ADDED_MODES):
                vals=[next(g[field] for g in groups if
                           g['cluster']==cluster and g['seed']==seed and g['mode']==mode)
                      for seed in REPLICATES]
                ax.plot(range(1,len(REPLICATES)+1),vals,'o-',color=colors[mi],
                        label='Independent Gaussian' if mi==0 else 'Correlated Gaussian')
            ax.set(title=NAMES[cluster],xlabel='Predeclared repeat index')
            ax.set_xticks(range(1,len(REPLICATES)+1))
            if top is not None:ax.set_ylim(0,top)
            ax.grid(alpha=.2)
        axs[0].set_ylabel(ylabel)
        axs[0].legend(frameon=False,fontsize=8)
        fig.suptitle('True HST backgrounds, identical planting sites / varied added-noise seeds')
        fig.tight_layout()
        fig.savefig(out/filename,dpi=155)
        plt.close(fig)

    fig,ax=plt.subplots(figsize=(9.8,4.8))
    for ci,cluster in enumerate(CLUSTERS):
        for mi,mode in enumerate(ADDED_MODES):
            values=[g['median_paired_local_minus_center_abs_relative_error']
                    for g in groups if g['cluster']==cluster and g['mode']==mode]
            x=ci*2.6+mi*.8
            ax.scatter([x]*len(values),values,color=colors[mi],alpha=.7)
            ax.plot([x-.20,x+.20],[np.median(values)]*2,color=colors[mi],linewidth=3)
    ax.axhline(0,color='black',linestyle='--',linewidth=1)
    ax.set_xticks([ci*2.6+.4 for ci in range(3)],[NAMES[c] for c in CLUSTERS])
    ax.set_ylabel('Per-seed median paired |local| minus |center| amplitude error')
    ax.set_title('NGC 104 and other clusters: no post-hoc preferred seed')
    ax.grid(axis='y',alpha=.2)
    fig.tight_layout()
    fig.savefig(out/'paired_model_seed_spread.png',dpi=155)
    plt.close(fig)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description='Make real-HST repeated noise measurements and figures')
    p.add_argument('reports',nargs='+')
    p.add_argument('--out',required=True)
    p.add_argument('--run-id')
    p.add_argument('--source-commit')
    a=p.parse_args(argv)
    r=aggregate(a.reports,a.out,run_id=a.run_id,source_commit=a.source_commit)
    print(json.dumps({'original_sites':r['total_original_sites'],
                      'new_noise_trials':r['total_added_noise_trials'],
                      'groups':r['per_seed_groups'],
                      'calibrated_HST_detector_noise':False},indent=2))


if __name__=='__main__':
    main()
