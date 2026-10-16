"""Photutils 3.x reference photometry; no science-calibrated uncertainties.

Use a pixel-integrated Gaussian PRF for toy images, or a verified ePSF
ImagePSF for real images. A background-subtracted image is fit in all cases.
"""
from __future__ import annotations
import numpy as np
from .detection import robust_background, scipy_peak_detector, dao_detector
from .artificial_stars import normalized_psf


def measure_sources(data, method, *, fwhm=2.6, threshold=5.,
                    fit_shape=9, psf_stamp=None, iterations=2, valid_mask=None):
    """Return an Nx4 array [x, y, fitted_flux, formal_flux_error].

    'scipy' and 'dao' are centroid-only controls, with NaN flux columns.
    'psf' and 'iterative' are native Photutils PSFPhotometry baselines.
    """
    image=np.asarray(data,float)
    if image.ndim!=2 or not np.isfinite(image).any():
        raise ValueError('Expected finite 2D image')
    if method not in ('scipy','dao','psf','iterative'):
        raise ValueError(f'Unsupported method: {method}')
    if fwhm<=0 or threshold<=0 or fit_shape<3 or fit_shape%2!=1 or iterations<1:
        raise ValueError('Invalid PSF fitting parameters')
    valid=np.isfinite(image)
    if valid_mask is not None:
        if np.shape(valid_mask)!=image.shape:
            raise ValueError('Mask must match image')
        valid &= np.asarray(valid_mask,bool)
    if not valid.any():
        raise ValueError('Image contains no valid pixels')

    background, rms=robust_background(np.where(valid,image,np.nan))
    data0=np.where(valid,image-background,0.)
    if method in ('scipy','dao'):
        # Both controls see the same sky-subtracted, invalid-pixel-masked image.
        # Detector background estimation is repeated for these quicklook controls.
        arr=np.where(valid,data0,np.nan)
        p=(scipy_peak_detector(arr,fwhm=fwhm,nsigma=threshold)
           if method=='scipy' else
           dao_detector(arr,fwhm=fwhm,nsigma=threshold))
        if len(p):
            good=valid[np.clip(np.rint(p[:,1]).astype(int),0,len(valid)-1),
                       np.clip(np.rint(p[:,0]).astype(int),0,valid.shape[1]-1)]
            p=p[good]
        return np.column_stack((p,np.full((len(p),2),np.nan)))

    from photutils.detection import DAOStarFinder
    from photutils.psf import (CircularGaussianPRF, ImagePSF,
                               SourceGrouper, PSFPhotometry,
                               IterativePSFPhotometry)
    if psf_stamp is None:
        model=CircularGaussianPRF(fwhm=fwhm)
        model.fwhm.fixed=True
    else:
        model=ImagePSF(normalized_psf(psf_stamp),oversampling=1)
    finder=DAOStarFinder(threshold=threshold*rms,fwhm=fwhm,
                         exclude_border=True)
    common=dict(psf_model=model,fit_shape=fit_shape,finder=finder,
                grouper=SourceGrouper(min_separation=2*fwhm),
                aperture_radius=max(1.5*fwhm,2.))
    if method=='psf':
        phot=PSFPhotometry(**common)
    else:
        phot=IterativePSFPhotometry(**common,maxiters=iterations,mode='new')
    # No physically justified variance map is available at this stage:
    # do not invent one or claim accurate reduced chi-square / confidence.
    result=phot(data0,mask=~valid)
    if result is None or not len(result):
        return np.empty((0,4))
    x=np.asarray(result['x_fit'],float)
    y=np.asarray(result['y_fit'],float)
    flux=np.asarray(result['flux_fit'],float)
    err=(np.asarray(result['flux_err'],float)
         if 'flux_err' in result.colnames else np.full(len(result),np.nan))
    flags=(np.asarray(result['flags'],int)
           if 'flags' in result.colnames else np.zeros(len(result),int))
    good=(np.isfinite(x)&np.isfinite(y)&np.isfinite(flux)&(flux>0)&(flags==0))
    return np.column_stack((x[good],y[good],flux[good],err[good]))
