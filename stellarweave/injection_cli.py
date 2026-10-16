"""Command-line artificial-star recovery, synthetic or guarded FITS."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .artificial_stars import gaussian_prf,normalized_psf,experiment
from .simulation import simulated_field
from .detection import scipy_peak_detector,dao_detector


def main(argv=None):
    p=argparse.ArgumentParser(description='StellarWeave v0.4.0 injection/recovery')
    p.add_argument('--fits',help='FITS science image; requires --psf-stamp and --image-unit')
    p.add_argument('--psf-stamp',help='Same-filter, pixel-integrated normalized PSF .npy')
    p.add_argument('--image-unit',help='Required for FITS; must agree with FITS BUNIT')
    p.add_argument('--size',type=int,default=320)
    p.add_argument('--crop-x',type=int)
    p.add_argument('--crop-y',type=int)
    p.add_argument('--detector',choices=['scipy','dao'],default='scipy')
    p.add_argument('--fwhm',type=float,default=2.6)
    p.add_argument('--threshold',type=float,default=4.)
    p.add_argument('--flux-levels',type=float,nargs='+',default=[350,800,1800,4000])
    p.add_argument('--trials',type=int,default=3)
    p.add_argument('--per-flux',type=int,default=8)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--match-radius',type=float,default=2.)
    p.add_argument('--baseline-exclusion-radius',type=float,default=2.)
    p.add_argument('--avoid-existing-radius',type=float,default=0.)
    p.add_argument('--out',default='outputs/injection_v04')
    a=p.parse_args(argv)
    if (a.size<48 or a.threshold<=0 or a.fwhm<=0 or a.match_radius<=0
            or a.baseline_exclusion_radius<=0):
        p.error('Size must be >=48; detector/matching parameters must be positive')
    if a.fits:
        if not a.psf_stamp or not a.image_unit:
            p.error('FITS mode REQUIRES both --psf-stamp and --image-unit')
        from .dataset import read_fits_science
        img,mask,_,meta=read_fits_science(a.fits)
        given=a.image_unit.strip().lower().replace(' ','')
        header=str(meta.get('bunit') or '').strip().lower().replace(' ','')
        if header and header!=given:
            p.error(f'BUNIT mismatch: FITS={meta["bunit"]!r}, argument={a.image_unit!r}')
        h,w=img.shape
        if min(h,w)<a.size:
            p.error('Image smaller than requested crop')
        x=max(0,min(w-a.size,(w-a.size)//2 if a.crop_x is None else a.crop_x))
        y=max(0,min(h-a.size,(h-a.size)//2 if a.crop_y is None else a.crop_y))
        image=img[y:y+a.size,x:x+a.size]
        valid=mask[y:y+a.size,x:x+a.size]
        kernel=normalized_psf(np.load(a.psf_stamp,allow_pickle=False))
        mode='FITS_EMPIRICAL_PSF_UNVERIFIED_SCIENCE'
        metadata={'source_fits':a.fits,'source_psf':a.psf_stamp,
                  'fits_header':meta,'crop_xy':[x,y],'image_unit':a.image_unit,
                  'psf_provenance_requires_manual_validation':True}
    else:
        if a.psf_stamp or a.image_unit:
            p.error('--psf-stamp/--image-unit require --fits')
        image,_=simulated_field(size=a.size,n_stars=max(10,a.size//4),
                                fwhm=a.fwhm,seed=a.seed)
        valid=np.isfinite(image)
        kernel=gaussian_prf(a.fwhm,size=max(9,2*int(np.ceil(3*a.fwhm))+1))
        mode='SYNTHETIC_GAUSSIAN_PRF'
        metadata={'synthetic_seed':a.seed,'unit':'synthetic electrons'}
    def detect(z):
        if a.detector=='dao':
            return dao_detector(z,fwhm=a.fwhm,nsigma=a.threshold)
        return scipy_peak_detector(z,fwhm=a.fwhm,nsigma=a.threshold)
    summary=experiment(image,kernel,a.flux_levels,detector=detect,valid_mask=valid,
                       trials=a.trials,stars_per_flux=a.per_flux,seed=a.seed,
                       match_radius=a.match_radius,
                       baseline_exclusion_radius=a.baseline_exclusion_radius,
                       fwhm_px=a.fwhm,avoid_existing_radius=a.avoid_existing_radius)
    rows=summary.pop('rows')
    out=Path(a.out)
    out.mkdir(parents=True,exist_ok=True)
    with (out/'injections.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    manifest={'version':'0.4.0','mode':mode,
              'parameters':{k:v for k,v in vars(a).items() if k!='out'},
              'source_metadata':metadata,'protocol':{
                 'injection':'additive pixel-integrated PRF / PSF, expectation-only',
                 'trials':'reinject into unmodified base image each time',
                 'recovery':'new centroids not matched to base detections',
                 'reference_catalog_is_not_exhaustive':True,
                 'unmatched_new_detections_are_not_certified_false_positives':True,
                 'real_hst_science_validated':False,
                 'instrumental_flux_not_calibrated_magnitude':True},
              **summary}
    (out/'manifest.json').write_text(
        json.dumps(manifest,indent=2,allow_nan=False),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7.5,4.8))
    for label in ('all','isolated','crowded','blend'):
        group=[r for r in summary['groups']
               if r['crowding']==label and r['injected']>0]
        if not group:
            continue
        x=np.asarray([r['flux'] for r in group])
        y=np.asarray([r['completeness'] for r in group])
        lo=np.asarray([r['wilson95'][0] for r in group])
        hi=np.asarray([r['wilson95'][1] for r in group])
        err=np.maximum(0,np.vstack((y-lo,hi-y)))
        ax.errorbar(x,y,yerr=err,capsize=3,marker='o',label=label)
    ax.set_xscale('log')
    ax.set_ylim(-.06,1.08)
    ax.set_xlabel('Injected instrumental flux (image units)')
    ax.set_ylabel('Recovery fraction (Wilson 95% CI)')
    ax.set_title(f'{mode}: {a.detector}, not science-calibrated')
    ax.grid(alpha=.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out/'completeness.png',dpi=175)
    plt.close(fig)
    print(json.dumps({'mode':mode,'baseline_detections':summary['baseline_detections'],
                      'injections':len(rows),'trials':len(summary['trials']),
                      'out':str(out.resolve())},indent=2))


if __name__=='__main__':
    main()
