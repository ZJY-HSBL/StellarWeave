"""Reproducible 3-cluster real-HST injected-star benchmark figures.

Every trial, including non-detections, is counted. All training sources
and injection templates are frozen within each run; scientific flux units
and photon noise are not calibrated.
"""
from __future__ import annotations
import argparse
import csv
import json
import shutil
from pathlib import Path
import numpy as np
from .real_injection import CLUSTERS,NEIGHBOR_CLASSES,FLUX_RMS_MULTIPLIERS
from .real_hst import sha256_file

NAMES={'ngc6397':'NGC 6397','ngc6752':'NGC 6752','ngc0104':'NGC 104'}


def wilson(recovered,n):
    """95% Wilson binomial interval for injection recovery on fixed sites."""
    if n==0:
        return None
    z=1.959963984540054
    p=recovered/n
    denom=1+z*z/n
    ctr=(p+z*z/(2*n))/denom
    ext=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/denom
    return [float(max(0,ctr-ext)),float(min(1,ctr+ext))]


def validate(path):
    path=Path(path)
    r=json.loads(path.read_text(encoding='utf-8'))
    if (r.get('cluster') not in CLUSTERS
        or r.get('kind')!='REAL_HST_PLANTED_HELDOUT_STAR_RECOVERY'
        or r.get('science_calibrated') is not False
        or r.get('real_background_images') is not True
        or r.get('models_fitted_on_injection_sources') is not False
        or r.get('matched_source_position_known_for_forced_flux') is not True):
        raise ValueError(f'Invalid experimental provenance in {path}')
    trials=r.get('trials',[])
    if len(trials)!=r.get('total_attempted'):
        raise ValueError('Missing injected-star trials')
    if sum(s['detected'] for s in trials)!=r.get('total_recovered'):
        raise ValueError('Recovery count mismatch')
    for item in trials:
        if (item['cluster']!=r['cluster'] or
            item['class'] not in NEIGHBOR_CLASSES or
            item['flux_over_robust_tile_rms'] not in FLUX_RMS_MULTIPLIERS):
            raise ValueError('Unexpected star identity or experimental bin')
        if not np.isfinite([item['forced_local_abs_rel_error'],
                            item['forced_center_abs_rel_error']]).all():
            raise ValueError('Non-finite photometric diagnostic')
    if [t['tile'] for t in r['tiles']]!=[0,1,2,3,5,6,7,8]:
        raise ValueError('Missing or cherry-picked tiles')
    return r


def aggregate(paths,out_dir,*,run_id=None,source_commit=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    reports={}
    source={}
    for p in paths:
        r=validate(p)
        c=r['cluster']
        if c in reports:
            raise ValueError('Duplicate cluster')
        reports[c]=r;source[c]=Path(p)
    if set(reports)!=set(CLUSTERS):
        raise ValueError('Three original HST clusters required')
    out=Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    raw=out/'source_reports'
    raw.mkdir(exist_ok=True)
    for c in CLUSTERS: shutil.copyfile(source[c],raw/f'{c}.json')
    rows=[];groups=[];clusters=[]
    for c in CLUSTERS:
        r=reports[c]
        trials=r['trials']
        rows.extend(trials)
        clusters.append({
            'cluster':c,'attempted':len(trials),
            'recovered':r['total_recovered'],
            'recovery_fraction':r['total_recovered']/len(trials) if trials else None,
            'source_report_sha256':sha256_file(source[c]),
            'fits_sha256':r['fits_sha256'],
            'original_prf_sha256':r['prf_sha256']})
        for cls in NEIGHBOR_CLASSES:
            for flux in FLUX_RMS_MULTIPLIERS:
                subset=[s for s in trials if s['class']==cls
                        and s['flux_over_robust_tile_rms']==flux]
                n=len(subset)
                rec=sum(s['detected'] for s in subset)
                groups.append({
                    'cluster':c,'class':cls,'flux_over_tile_rms':flux,
                    'attempted':n,'recovered':rec,
                    'recovery_fraction':rec/n if n else None,
                    'recovery_wilson95':wilson(rec,n),
                    'median_detector_position_error_px':float(np.median([
                        x['detector_peak_position_error_px'] for x in subset if x['detected']
                    ])) if rec else None,
                    'median_forced_local_abs_rel_error':float(np.median([
                        x['forced_local_abs_rel_error'] for x in subset])) if n else None,
                    'median_forced_center_abs_rel_error':float(np.median([
                        x['forced_center_abs_rel_error'] for x in subset])) if n else None,
                    'median_local_minus_center_abs_rel_error':float(np.median([
                        x['forced_local_abs_rel_error']-x['forced_center_abs_rel_error']
                        for x in subset])) if n else None})
    summary={'version':'0.14.0',
             'kind':'REAL_HST_ARTIFICIAL_STAR_RECOVERY_THREE_CLUSTERS',
             'source_run_id':str(run_id or ''),
             'source_commit':str(source_commit or ''),
             'clusters':clusters,'groups':groups,
             'total_trials':len(rows),
             'total_recovered':sum(r['detected'] for r in rows),
             'source_noise_is_real_HST':True,
             'injection_photon_noise_simulated':False,
             'physical_flux_units_calibrated':False,
             'forced_flux_position_is_known':True,
             'new_astrometric_precision_proven':False,
             'stellar_catalog_completeness_proven':False,
             'instrument_epsf_calibrated':False,
             'all_nondetections_retained':True,
             'experiment':'Held-out HST star stamps injected one-at-a-time into original real HST image tiles, into positions selected by original peak separations',
             'limitations':['Known injection position supplied to forced model flux fits',
                            'Flux levels are relative to robust tile pixel RMS, not physical SNR',
                            'Relatively isolated positions may have unresolved neighbors',
                            'No injected photon shot noise or detector-noise re-simulation',
                            'No real-star population completeness or calibrated HST magnitude accuracy']}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False),encoding='utf-8')
    with (out/'per_trial.csv').open('w',newline='',encoding='utf-8') as f:
        if rows:
            w=csv.DictWriter(f,fieldnames=list(rows[0]))
            w.writeheader();w.writerows(rows)
    with (out/'groups.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(groups[0]))
        w.writeheader();w.writerows(groups)
    labels={'near_detected_peak':'Near original peaks (5–11 px)',
            'relatively_isolated':'Relatively isolated (≥16 px)'}
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig,axs=plt.subplots(1,3,figsize=(13,4.3),sharey=True)
    for i,c in enumerate(CLUSTERS):
        ax=axs[i]
        for j,cls in enumerate(NEIGHBOR_CLASSES):
            subset=[g for g in groups if g['cluster']==c and g['class']==cls]
            data=[g['recovery_fraction'] if g['recovery_fraction'] is not None
                  else np.nan for g in subset]
            ax.plot(FLUX_RMS_MULTIPLIERS,data,'o-',label=labels[cls],color=colors[j])
            for x,g in zip(FLUX_RMS_MULTIPLIERS,subset):
                if g['recovery_wilson95'] is None: continue
                lo,hi=g['recovery_wilson95']
                ax.plot([x,x],[lo,hi],color=colors[j],linewidth=2,alpha=.65)
        ax.set(title=NAMES[c],xlabel='Injected flux / robust tile RMS')
        ax.grid(alpha=.2)
        ax.set_ylim(-.035,1.06)
    axs[0].set_ylabel('Re-detected artificial-star fraction')
    axs[0].legend(frameon=False,fontsize=8)
    fig.suptitle('Real HST backgrounds · planted held-out PRF stamps (not catalog completeness)')
    fig.tight_layout()
    fig.savefig(out/'recovery_by_flux_crowding.png',dpi=165)
    plt.close(fig)

    fig,axs=plt.subplots(1,3,figsize=(13,4.3))
    for i,c in enumerate(CLUSTERS):
        ax=axs[i]
        for j,method in enumerate(('local','center')):
            vals=[]
            for f in FLUX_RMS_MULTIPLIERS:
                subset=[s for s in rows if s['cluster']==c and
                        s['flux_over_robust_tile_rms']==f]
                vals.append(float(np.median(
                    [s[f'forced_{method}_abs_rel_error'] for s in subset]))
                    if subset else np.nan)
            ax.plot(FLUX_RMS_MULTIPLIERS,vals,'o-',label=f'{method} PRF',
                    color=colors[j])
        ax.set(title=NAMES[c],xlabel='Injected flux / robust tile RMS',
               ylabel='Median absolute forced flux relative error')
        ax.grid(alpha=.2)
    axs[0].legend(frameon=False,fontsize=9)
    fig.suptitle('Known-coordinate forced native-pixel PRF fits · NOT blind calibrated photometry')
    fig.tight_layout()
    fig.savefig(out/'forced_flux_error.png',dpi=165)
    plt.close(fig)

    fig,axs=plt.subplots(1,3,figsize=(13,4.3),sharey=True)
    for i,c in enumerate(CLUSTERS):
        ax=axs[i]
        for j,cls in enumerate(NEIGHBOR_CLASSES):
            subset=[g for g in groups if g['cluster']==c and g['class']==cls]
            vals=[g['median_local_minus_center_abs_rel_error']
                  if g['median_local_minus_center_abs_rel_error'] is not None
                  else np.nan for g in subset]
            ax.plot(FLUX_RMS_MULTIPLIERS,vals,'o-',color=colors[j],
                    label=labels[cls])
        ax.axhline(0,color='black',linestyle='--',linewidth=1)
        ax.set(title=NAMES[c],xlabel='Injected flux / robust tile RMS')
        ax.grid(alpha=.2)
    axs[0].set_ylabel('Median |local flux error| − |center flux error|')
    axs[0].legend(frameon=False,fontsize=8)
    fig.suptitle('Paired forced-fit comparison · negative favors local · exploratory')
    fig.tight_layout()
    fig.savefig(out/'paired_forced_flux_difference.png',dpi=165)
    plt.close(fig)
    return summary


def main(argv=None):
    p=argparse.ArgumentParser(description='Generate real HST artificial star plots')
    p.add_argument('reports',nargs='+')
    p.add_argument('--out',required=True)
    p.add_argument('--run-id')
    p.add_argument('--source-commit')
    a=p.parse_args(argv)
    d=aggregate(a.reports,a.out,run_id=a.run_id,source_commit=a.source_commit)
    print(json.dumps({'trials':d['total_trials'],'recovered':d['total_recovered'],
                      'clusters':d['clusters'],'physical_flux_calibrated':False},indent=2))


if __name__=='__main__':
    main()
