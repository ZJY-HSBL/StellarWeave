"""Aggregate all real HST noise-sensitivity trials, failures, and plots.

Each noise level/mode has exactly the SAME frozen v0.14 injected sources.
Uncertainty intervals are Wilson over preselected sites and descriptive,
not a general astrophysical-population confidence interval.
"""
from __future__ import annotations
import argparse
import csv
import json
import shutil
from pathlib import Path
import numpy as np
from .noise_stress import MODES, CLUSTERS
from .real_hst import sha256_file
from .injection_results import wilson

NAMES={'ngc6397':'NGC 6397','ngc6752':'NGC 6752','ngc0104':'NGC 104'}
LABELS={'baseline':'Original HST',
        'iid_gaussian':'HST + independent Gaussian',
        'correlated_gaussian':'HST + correlated Gaussian'}


def _validate(path):
    p=Path(path)
    d=json.loads(p.read_text(encoding='utf-8'))
    if (d.get('cluster') not in CLUSTERS or
        d.get('kind')!='REAL_HST_FROZEN_INJECTION_CONTROLLED_NOISE_STRESS' or
        d.get('modes')!=list(MODES) or
        d.get('true_HST_detector_noise_model') is not False or
        d.get('no_heldout_source_in_training') is not True or
        d.get('photon_shot_noise_applied_to_real_HST') is not False or
        d.get('all_original_injection_positions_and_shapes_frozen') is not True):
        raise ValueError(f'Missing provenance or false science claim: {p}')
    trials=d.get('trials',[])
    if len(trials)!=d.get('expected_original_trials',-1)*len(MODES):
        raise ValueError('Incomplete frozen injections or noise controls')
    keys={mode:set() for mode in MODES}
    metadata={}
    for row in trials:
        mode=row['mode']
        if mode not in keys: raise ValueError('Unexpected noise mode')
        key=(row['tile_index'],row['class'],row['location_number'],
             row['flux_over_robust_tile_rms'])
        if key in keys[mode]: raise ValueError('Repeated trial')
        keys[mode].add(key)
        stable=(row['x'],row['y'],row['source_holdout_index'],
                row['injected_flux_image_units'])
        if key in metadata and metadata[key]!=stable:
            raise ValueError('Noise mode changed source location, shape, or injected flux')
        metadata[key]=stable
        for field in ('local_forced_abs_rel_error','center_forced_abs_rel_error'):
            if not np.isfinite(float(row[field])) or float(row[field])<0:
                raise ValueError('Nonfinite or negative forced flux error')
    if any(keys[mode]!=keys[MODES[0]] for mode in MODES):
        raise ValueError('A mode omitted or altered a trial')
    expected_audits={(key[0],mode) for key in keys[MODES[0]] for mode in MODES}
    seen_audits=[(int(a['tile_index']),a['mode']) for a in d.get('noise_audits',[])]
    if set(seen_audits)!=expected_audits or len(seen_audits)!=len(expected_audits):
        raise ValueError('Incomplete or duplicate noise auditing for original tiles')
    for row in d['noise_audits']:
        expected=float(row['target_added_pixel_rms'])
        measured=float(row['actual_added_pixel_rms'])
        if not np.isfinite(measured) or abs(measured-expected)>1e-6*max(1.,expected):
            raise ValueError('Noise pixel RMS does not match declared target')
        corr=row.get('added_field_lag1_correlation')
        if row['mode']=='correlated_gaussian' and (
                corr is None or float(corr)<.15):
            raise ValueError('Declared spatially correlated noise is uncorrelated')
    return d


def aggregate(paths,out_dir,*,run_id=None,source_commit=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    reports={};source={}
    for path in paths:
        d=_validate(path)
        name=d['cluster']
        if name in reports: raise ValueError('Duplicate HST field')
        reports[name]=d
        source[name]=Path(path)
    if set(reports)!=set(CLUSTERS):
        raise ValueError('Require all three observed HST fields')
    out=Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    raw=out/'source_reports';raw.mkdir(exist_ok=True)
    trials=[];groups=[];meta=[]
    for cluster in CLUSTERS:
        d=reports[cluster]
        shutil.copyfile(source[cluster],raw/f'{cluster}.json')
        for a in d['noise_audits']:
            meta.append(dict(cluster=cluster,**a))
        for r in d['trials']:
            trials.append(r)
        for mode in MODES:
            selection=[r for r in d['trials'] if r['mode']==mode]
            n=len(selection);recovered=sum(r['detected'] for r in selection)
            groups.append({
                'cluster':cluster,'mode':mode,
                'attempted':n,'recovered':recovered,
                'recovery_fraction':recovered/n if n else None,
                'wilson_95ci':wilson(recovered,n),
                'ambiguous_existing_peaks':sum(
                    r['preexisting_noise_peak_at_truth'] for r in selection),
                'median_abs_local_flux_error':float(np.median(
                    [r['local_forced_abs_rel_error'] for r in selection])) if n else None,
                'median_abs_center_flux_error':float(np.median(
                    [r['center_forced_abs_rel_error'] for r in selection])) if n else None,
                'median_paired_local_minus_center_flux_error':float(np.median(
                    [r['local_forced_abs_rel_error']-r['center_forced_abs_rel_error']
                     for r in selection])) if n else None})
    summary={'version':'0.15.0',
             'kind':'THREE_CLUSTER_REAL_HST_NOISE_SENSITIVITY_EXPERIMENT',
             'source_run_id':str(run_id or ''),
             'source_commit':str(source_commit or ''),
             'input_provenance':[
                 {'cluster':cluster,
                  'noise_report_sha256':sha256_file(source[cluster]),
                  'original_hst_fits_sha256':reports[cluster]['fits_sha256'],
                  'frozen_v014_injection_report_sha256':reports[cluster]['v014_baseline_sha256']}
                 for cluster in CLUSTERS],
             'clusters':list(CLUSTERS),'modes':list(MODES),
             'groups':groups,'noise_audits':meta,
             'total_mode_trials':len(trials),
             'total_original_sites':sum(d['expected_original_trials'] for d in reports.values()),
             'no_original_frozen_sites_removed':True,
             'added_noise_fraction_original_robust_tile_rms':.5,
             'correlation_sigma_px':1.3,
             'calibrated_HST_detector_noise':False,
             'physical_source_poisson_simulated':False,
             'source_flux_physical_units_calibrated':False,
             'full_catalog_completeness':False,
             'unweighted_known_position_forced_flux':True,
             'interpretation':'A controlled relative-noise stress test on genuine telescope backgrounds; relative performance only. Correlated Gaussian is NOT measured drizzle covariance.',
             'negative_results_included':True}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False),
                                    encoding='utf-8')
    for filename,values in [('per_trial.csv',trials),('groups.csv',groups),
                            ('noise_audits.csv',meta)]:
        with (out/filename).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(values[0]))
            w.writeheader();w.writerows(values)
    x=np.arange(len(CLUSTERS))
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig,ax=plt.subplots(figsize=(9.6,4.8))
    for k,mode in enumerate(MODES):
        values=[next(g['recovery_fraction'] for g in groups if
                     g['cluster']==c and g['mode']==mode) for c in CLUSTERS]
        ax.bar(x+(k-1)*.24,values,width=.23,color=colors[k],label=LABELS[mode])
    ax.set_xticks(x,[NAMES[c] for c in CLUSTERS])
    ax.set_ylim(0,1.06)
    ax.set_ylabel('Matched planted-source fraction')
    ax.set_title('Real HST backgrounds: unchanged injected sources, added noise controls')
    ax.legend(frameon=False,fontsize=9)
    ax.grid(axis='y',alpha=.2)
    fig.tight_layout()
    fig.savefig(out/'recovery_noise_sensitivity.png',dpi=160)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(9.6,4.8))
    for k,mode in enumerate(MODES):
        vals=[next(g['median_abs_local_flux_error'] for g in groups if
                   g['cluster']==c and g['mode']==mode) for c in CLUSTERS]
        ax.plot(x,vals,'o-',color=colors[k],label=LABELS[mode],markersize=6)
    ax.set_xticks(x,[NAMES[c] for c in CLUSTERS])
    ax.set_ylabel('Median |relative forced LOCAL amplitude error|')
    ax.set_title('Known-position PRF flux under controlled added noise (NOT calibrated)')
    ax.grid(alpha=.2)
    ax.legend(frameon=False,fontsize=9)
    fig.tight_layout()
    fig.savefig(out/'local_forced_flux_noise.png',dpi=160)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(9.6,4.8))
    for k,mode in enumerate(MODES):
        vals=[next(g['median_paired_local_minus_center_flux_error'] for g in groups if
                   g['cluster']==c and g['mode']==mode) for c in CLUSTERS]
        ax.plot(x,vals,'o-',color=colors[k],label=LABELS[mode],markersize=6)
    ax.axhline(0,color='black',linestyle='--',linewidth=1)
    ax.set_xticks(x,[NAMES[c] for c in CLUSTERS])
    ax.set_ylabel('Median paired |local error| minus |center error|')
    ax.set_title('Real HST forced PRF contrast: negative favors local model')
    ax.grid(alpha=.2)
    ax.legend(frameon=False,fontsize=9)
    fig.tight_layout()
    fig.savefig(out/'paired_prf_noise_robustness.png',dpi=160)
    plt.close(fig)
    return summary


def main(argv=None):
    p=argparse.ArgumentParser(description='Plot genuine three-cluster noise robustness')
    p.add_argument('reports',nargs='+')
    p.add_argument('--out',required=True)
    p.add_argument('--run-id')
    p.add_argument('--source-commit')
    a=p.parse_args(argv)
    r=aggregate(a.reports,a.out,run_id=a.run_id,source_commit=a.source_commit)
    print(json.dumps({'total_trials':r['total_mode_trials'],
                      'unique_frozen_sites':r['total_original_sites'],
                      'groups':r['groups'],
                      'HST_noise_calibrated':False},indent=2))


if __name__=='__main__':
    main()
