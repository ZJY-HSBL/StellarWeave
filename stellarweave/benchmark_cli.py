"""Run fair PSFPhotometry vs iterative PSF recovery on shared injected images."""
import argparse
import json
from pathlib import Path
import numpy as np
from .simulation import simulated_field
from .artificial_stars import gaussian_prf,normalized_psf
from .benchmark import compare_photometry,save_benchmark


def main(argv=None):
    p=argparse.ArgumentParser(description='StellarWeave v0.5.0 PSF baseline comparison')
    p.add_argument('--fits')
    p.add_argument('--psf-stamp',help='Measured, pixel-integrated, same-band ePSF .npy required with --fits')
    p.add_argument('--image-unit',help='Instrumental FITS unit, compared to BUNIT')
    p.add_argument('--size',type=int,default=160)
    p.add_argument('--crop-x',type=int)
    p.add_argument('--crop-y',type=int)
    p.add_argument('--fwhm',type=float,default=2.6)
    p.add_argument('--threshold',type=float,default=5.)
    p.add_argument('--flux-levels',nargs='+',type=float,default=[700,2000,5000])
    p.add_argument('--methods',nargs='+',choices=['scipy','dao','psf','iterative'],
                   default=['dao','psf','iterative'])
    p.add_argument('--fit-shape',type=int,default=9)
    p.add_argument('--iterations',type=int,default=2)
    p.add_argument('--trials',type=int,default=2)
    p.add_argument('--per-flux',type=int,default=4)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--out',default='outputs/photometry_v05')
    a=p.parse_args(argv)
    if a.size<48 or a.fwhm<=0:
        p.error('size must be >=48 and FWHM > 0')
    if a.fits:
        if not a.psf_stamp or not a.image_unit:
            p.error('Real FITS requires a verified --psf-stamp and matching --image-unit')
        from .dataset import read_fits_science
        full,valid,_,info=read_fits_science(a.fits)
        requested=a.image_unit.strip().replace(' ','').lower()
        bunit=str(info.get('bunit') or '').strip().replace(' ','').lower()
        if not bunit:
            p.error('Missing FITS BUNIT; do not infer instrumental flux units')
        if bunit!=requested:
            p.error('Supplied unit does not match FITS BUNIT')
        h,w=full.shape
        if min(h,w)<a.size:
            p.error('FITS image smaller than specified crop')
        x=(w-a.size)//2 if a.crop_x is None else a.crop_x
        y=(h-a.size)//2 if a.crop_y is None else a.crop_y
        if not (0<=x<=w-a.size and 0<=y<=h-a.size):
            p.error('Crop offset outside image')
        image=full[y:y+a.size,x:x+a.size]
        valid=valid[y:y+a.size,x:x+a.size]
        kernel=normalized_psf(np.load(a.psf_stamp,allow_pickle=False))
        mode='REAL_FITS_UNVERIFIED_SCIENCE'
        notes={'source':a.fits,'psf_stamp':a.psf_stamp,
               'image_unit':info.get('bunit'),'crop':[x,y],
               'psf_position_and_normalization_must_be_validated':True}
    else:
        if a.psf_stamp or a.image_unit:
            p.error('--psf-stamp and --image-unit require --fits')
        image,_=simulated_field(size=a.size,n_stars=max(8,a.size//10),
                                fwhm=a.fwhm,seed=a.seed,min_flux=1500)
        valid=np.isfinite(image)
        kernel=gaussian_prf(fwhm=a.fwhm,size=max(9,2*int(np.ceil(3*a.fwhm))+1))
        mode='SYNTHETIC_GAUSSIAN_PRF'
        notes={'synthetic_seed':a.seed,'flux_unit':'toy electrons'}
    result=compare_photometry(
        image,kernel,a.flux_levels,methods=tuple(a.methods),seed=a.seed,
        trials=a.trials,per_flux=a.per_flux,valid_mask=valid,
        fwhm=a.fwhm,threshold=a.threshold,fit_shape=a.fit_shape,
        iterations=a.iterations)
    result['metadata']['mode']=mode
    result['metadata']['source']=notes
    result['metadata']['real_hst_scientific_validation']=False
    path=save_benchmark(result,a.out)
    print(json.dumps({'mode':mode,'methods':a.methods,'injections_per_method':
                      a.trials*a.per_flux*len(a.flux_levels),
                      'out':str(Path(path).resolve())},indent=2))


if __name__=='__main__':
    main()
