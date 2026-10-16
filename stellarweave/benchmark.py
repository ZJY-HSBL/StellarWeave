"""Fair synthetic/guarded-FITS injection benchmark on identical trial images."""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .artificial_stars import (normalized_psf,inject_stars,sample_positions,
                                _greedy_unique_match,wilson_interval)
from .photometry import measure_sources
from .detection import scipy_peak_detector


def compare_photometry(image,psf,flux_levels,*,methods=('dao','psf','iterative'),
                       seed=42,trials=2,per_flux=4,valid_mask=None,
                       fwhm=2.6,threshold=5.,fit_shape=9,iterations=2,
                       match_radius=2.0,baseline_exclusion_radius=2.0):
    """Single injected test set per trial is reused across *all* methods.

    Recovery requires a new centroid not attributable to a baseline source.
    Flux error is only summarized for recovered injections isolated from
    baseline detections by >=3 FWHM; it remains an exploratory proxy.
    """
    base=np.asarray(image,float)
    if base.ndim!=2 or not np.isfinite(base).any():
        raise ValueError('Expected finite 2D image')
    valid=np.isfinite(base) if valid_mask is None else np.asarray(valid_mask,bool)&np.isfinite(base)
    if valid.shape!=base.shape:
        raise ValueError('Mask shape mismatch')
    if trials<1 or per_flux<1 or fwhm<=0 or match_radius<=0:
        raise ValueError('Invalid experiment parameters')
    flux=np.asarray(flux_levels,float)
    if flux.ndim!=1 or not len(flux) or not np.isfinite(flux).all() or (flux<=0).any():
        raise ValueError('Flux levels must be positive')
    if not methods or len(set(methods))!=len(methods):
        raise ValueError('Method names must be unique')
    for method in methods:
        if method not in ('scipy','dao','psf','iterative'):
            raise ValueError('Unknown method: '+method)
    kernel=normalized_psf(psf)
    radius=kernel.shape[0]//2
    rng=np.random.default_rng(seed)
    # Method-independent crowding reference.
    reference=scipy_peak_detector(np.where(valid,base,np.nan),
                                   fwhm=fwhm,nsigma=threshold)
    crowd_tree=cKDTree(reference) if len(reference) else None
    baseline={m:measure_sources(base,m,fwhm=fwhm,threshold=threshold,
                                fit_shape=fit_shape,psf_stamp=kernel if m in ('psf','iterative') else None,
                                iterations=iterations,valid_mask=valid)
              for m in methods}
    records=[]
    for trial in range(trials):
        n=int(per_flux*len(flux))
        positions=sample_positions(base.shape,valid,n,radius,rng=rng,
                                   min_separation=max(2*fwhm,kernel.shape[0]/2))
        assigned=np.repeat(flux,per_flux)[rng.permutation(n)]
        injected=inject_stars(base,np.column_stack((positions,assigned)),
                              kernel,valid_mask=valid)
        # Identical images, coordinates and assigned injected fluxes for all methods.
        distances=(crowd_tree.query(positions)[0] if crowd_tree is not None
                   else np.full(n,np.inf))
        for method in methods:
            after=measure_sources(injected,method,fwhm=fwhm,threshold=threshold,
                                  fit_shape=fit_shape,
                                  psf_stamp=kernel if method in ('psf','iterative') else None,
                                  iterations=iterations,valid_mask=valid)
            original=baseline[method]
            old=_greedy_unique_match(original[:,:2],after[:,:2],baseline_exclusion_radius)
            excluded={j for _,j,_ in old}
            keep=[i for i in range(len(after)) if i not in excluded]
            new=after[keep]
            links=_greedy_unique_match(positions,new[:,:2],match_radius)
            match={i:(j,d) for i,j,d in links}
            for i in range(n):
                distance=float(distances[i])
                group=('blend' if distance<1.5*fwhm else
                       'crowded' if distance<3*fwhm else 'isolated')
                hit=match.get(i)
                pos_error=float(hit[1]) if hit else None
                recovered_flux=float(new[hit[0],2]) if hit and np.isfinite(new[hit[0],2]) else None
                flux_error=(recovered_flux/float(assigned[i])-1
                            if recovered_flux is not None and group=='isolated' else None)
                records.append({
                    'method':method,'trial':int(trial),'injection_id':int(i),
                    'x_px':float(positions[i,0]),'y_px':float(positions[i,1]),
                    'injected_flux':float(assigned[i]),'crowding':group,
                    'nearest_reference_px':distance if np.isfinite(distance) else None,
                    'recovered':hit is not None,'position_error_px':pos_error,
                    'fitted_flux':recovered_flux,
                    'relative_flux_error_isolated':flux_error,
                    'formal_flux_error':(float(new[hit[0],3]) if hit and np.isfinite(new[hit[0],3]) else None)
                })
    groups=[]
    for method in methods:
        for level in flux:
            for crowd in ('all','isolated','crowded','blend'):
                rows=[r for r in records
                      if r['method']==method and r['injected_flux']==float(level)
                      and (crowd=='all' or r['crowding']==crowd)]
                total=len(rows)
                recovered=sum(r['recovered'] for r in rows)
                isolated_errors=[r['relative_flux_error_isolated'] for r in rows
                                  if r['relative_flux_error_isolated'] is not None]
                groups.append({'method':method,'flux':float(level),
                               'crowding':crowd,'injected':total,'recovered':recovered,
                               'completeness':recovered/total if total else None,
                               'wilson95':wilson_interval(recovered,total),
                               'median_relative_flux_error_isolated':
                               float(np.median(isolated_errors)) if isolated_errors else None})
    return {'metadata':{'version':'0.5.0','seed':int(seed),'methods':list(methods),
                        'trials':int(trials),'per_flux':int(per_flux),
                        'base_shape':list(base.shape),'fwhm':float(fwhm),
                        'threshold':float(threshold),'fit_shape':int(fit_shape),
                        'iterations':int(iterations),
                        'same_injected_images_all_methods':True,
                        'isolated_flux_only':True,'science_calibrated':False,
                        'injection_missing_extra_shot_noise':True},
            'records':records,'groups':groups}


def save_benchmark(result,out_dir):
    """Emit machine-readable data and visualization without masking missing bins."""
    out=Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    records=result['records']
    with (out/'source_results.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    (out/'benchmark.json').write_text(json.dumps(
        {k:v for k,v in result.items() if k!='records'},
        indent=2,allow_nan=False),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7.8,4.9))
    for method in result['metadata']['methods']:
        groups=[g for g in result['groups']
                if g['method']==method and g['crowding']=='all' and g['injected']>0]
        if not groups:
            continue
        x=np.array([g['flux'] for g in groups])
        y=np.array([g['completeness'] for g in groups])
        lows=np.array([g['wilson95'][0] for g in groups])
        highs=np.array([g['wilson95'][1] for g in groups])
        ax.errorbar(x,y,yerr=np.maximum(0,np.vstack((y-lows,highs-y))),
                    marker='o',capsize=2,label=method)
    ax.set_xscale('log')
    ax.set_ylim(-.04,1.06)
    ax.set(xlabel='Injected instrumental flux (same units as image)',
           ylabel='New-source recovery fraction',
           title='Exploratory artificial-star PSF baseline (not calibrated HST)')
    ax.legend()
    ax.grid(alpha=.15)
    fig.tight_layout()
    fig.savefig(out/'method_comparison.png',dpi=170)
    plt.close(fig)
    return str(out)
