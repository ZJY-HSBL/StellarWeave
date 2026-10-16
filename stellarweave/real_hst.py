"""FITS/ACS catalog preflight and local empirical PSF stamp diagnostics.

This is an image-quality and provenance audit, NOT an astrometric/photometric
science validation. The PSF is a native-pixel robust stack, not an oversampled,
spatially varying HST ePSF. Outputs explicitly carry this limitation.
"""
from __future__ import annotations
import csv
import json
import hashlib
from pathlib import Path
import numpy as np
from scipy.ndimage import shift as shift_image
from scipy.spatial import cKDTree
from .catalog import read_acsggct_catalog, quality_mask
from .detection import scipy_peak_detector, robust_background


def sha256_file(path):
    sha=hashlib.sha256()
    with open(path,'rb') as f:
        for buf in iter(lambda:f.read(1024*1024),b''):
            sha.update(buf)
    return sha.hexdigest()


def filter_tokens(headers):
    found=set()
    for hdr in headers:
        for key in ('FILTER','FILTER1','FILTER2'):
            value=hdr.get(key)
            if value is not None:
                found.add(str(value).strip().upper())
    return sorted(found)


def read_science_crop(path,size=512,x0=None,y0=None):
    """Inspect FITS SCI image with an optional matched ERR/DQ plane crop."""
    if not Path(path).is_file():
        raise FileNotFoundError(path)
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(path,memmap=True) as hdus:
        science=[(i,h) for i,h in enumerate(hdus)
                 if h.data is not None and getattr(h.data,'ndim',0)==2]
        if not science:
            raise ValueError('FITS contains no two-dimensional image')
        science.sort(key=lambda it:(it[1].header.get('EXTNAME','').upper()=='SCI',
                                   it[1].data.size),reverse=True)
        i,h=science[0]
        height,width=h.data.shape
        size=min(size,height,width)
        if size<33:
            raise ValueError('Science FITS image too small for local PSF estimation')
        if x0 is None: x0=(width-size)//2
        if y0 is None: y0=(height-size)//2
        if not (0<=x0<=width-size and 0<=y0<=height-size):
            raise ValueError('Crop origin outside FITS image')
        yslice=slice(y0,y0+size)
        xslice=slice(x0,x0+size)
        image=np.asarray(h.data[yslice,xslice],dtype=np.float64).copy()
        err=None; dq=None
        sci_ver=h.header.get('EXTVER',1)
        for extra in hdus:
            if extra is h or extra.data is None or extra.data.shape!=(height,width):
                continue
            name=extra.header.get('EXTNAME','').upper()
            ver=extra.header.get('EXTVER',1)
            if ver != sci_ver: continue
            if name=='DQ' and dq is None:
                dq=np.array(extra.data[yslice,xslice],copy=True)
            if name=='ERR' and err is None:
                err=np.array(extra.data[yslice,xslice],dtype=np.float64,copy=True)
        good=np.isfinite(image)
        if dq is not None: good &= (dq==0)
        if err is not None: good &= np.isfinite(err)&(err>0)
        headers=[hdus[0].header,h.header]
        try:
            wcs=WCS(h.header,hdus).celestial
            has_wcs=bool(wcs.has_celestial and wcs.pixel_n_dim==2)
        except Exception:
            has_wcs=False
        metadata={'file':str(path),'sha256':sha256_file(path),
                  'sci_hdu_index':i,'shape':[int(height),int(width)],
                  'crop_xy':[int(x0),int(y0)],'crop_shape':list(image.shape),
                  'valid_fraction':float(good.mean()),
                  'has_celestial_wcs':has_wcs,'filters':filter_tokens(headers),
                  'bunit':str(h.header.get('BUNIT',hdus[0].header.get('BUNIT',''))).strip(),
                  'exptime':h.header.get('EXPTIME',hdus[0].header.get('EXPTIME')),
                  'photflam':h.header.get('PHOTFLAM',hdus[0].header.get('PHOTFLAM')),
                  'has_dq_plane':dq is not None,'has_err_plane':err is not None,
                  'err_unit_not_assumed':True}
    return image,good,err,metadata


def _psf_candidates(image,valid,stamp_size=15,fwhm=2.6,
                    min_neighbors=4.0,cap=80):
    """Find isolated relative bright peaks on one tile only."""
    if stamp_size<7 or stamp_size%2!=1:
        raise ValueError('stamp_size must be odd >= 7')
    background,rms=robust_background(np.where(valid,image,np.nan))
    xy=scipy_peak_detector(np.where(valid,image,np.nan),fwhm=fwhm,nsigma=6)
    if not len(xy): return np.empty((0,2)),{'candidates':0,'isolated':0}
    tree=cKDTree(xy)
    dist=tree.query(xy,k=2)[0][:,1] if len(xy)>1 else np.array([np.inf])
    rad=stamp_size//2
    bright=[]
    for i,(x,y) in enumerate(xy):
        ix,iy=int(round(x)),int(round(y))
        if ix-rad<0 or iy-rad<0 or ix+rad>=image.shape[1] or iy+rad>=image.shape[0]:
            continue
        patch=image[iy-rad:iy+rad+1,ix-rad:ix+rad+1]
        patch_good=valid[iy-rad:iy+rad+1,ix-rad:ix+rad+1]
        if not patch_good.all() or dist[i]<min_neighbors*fwhm:
            continue
        local_peak=image[iy,ix]-background
        if local_peak<10*rms: continue
        bright.append((local_peak,float(x),float(y)))
    bright.sort(reverse=True)
    return np.asarray([[x,y] for _,x,y in bright[:cap]],float).reshape(-1,2),{
        'candidates':int(len(xy)),'isolated':len(bright),
        'image_background_median':float(background),
        'approx_background_rms':float(rms),
        'background_rms_not_instrument_noise':True}


def empirical_psf(image,valid,*,stamp_size=15,fwhm=2.6,min_stars=4,max_stars=80):
    """Build a centered, unit-integral, native-pixel PSF via robust star stacking.

    Source centroids, background and neighbor isolation are approximate. The
    returned stamp should not be used for scientific photometry without
    instrument-specific validation and radial/position-dependence checks.
    """
    image=np.asarray(image,float)
    valid=np.asarray(valid,bool)
    if valid.shape!=image.shape or image.ndim!=2:
        raise ValueError('Invalid image/mask shape')
    xy,meta=_psf_candidates(image,valid,stamp_size,fwhm,cap=max_stars)
    radius=stamp_size//2
    yy,xx=np.mgrid[-radius:radius+1,-radius:radius+1]
    ring=(np.hypot(xx,yy)>=radius*.75)
    stamps=[];accepted=[]
    for x,y in xy:
        ix,iy=int(round(x)),int(round(y))
        raw=image[iy-radius:iy+radius+1,ix-radius:ix+radius+1]
        local_sky=float(np.median(raw[ring]))
        patch=np.maximum(raw-local_sky,0)
        if patch.sum()<=0:continue
        cy,cx=np.indices(patch.shape)
        center_x=float((cx*patch).sum()/patch.sum())
        center_y=float((cy*patch).sum()/patch.sum())
        if max(abs(center_x-radius),abs(center_y-radius))>1.5:
            continue
        centered=shift_image(patch,(radius-center_y,radius-center_x),
                             order=3,mode='constant',cval=0,prefilter=True)
        centered=np.maximum(centered,0)
        if centered.sum()<=0: continue
        stamps.append(centered/centered.sum())
        accepted.append([float(x),float(y)])
    if len(stamps)<min_stars:
        raise ValueError(f'Only {len(stamps)} isolated stars usable; '
                         f'minimum={min_stars}. Select another crop.')
    psf=np.median(np.stack(stamps),axis=0)
    psf=np.maximum(psf,0)
    psf/=psf.sum()
    profile=np.sum(psf,axis=0)
    half=psf.max()/2
    inner=psf[psf>=half]
    meta.update({'used_stars':len(stamps),'stamp_size':int(stamp_size),
                 'normalization_sum':float(psf.sum()),
                 'oversampling':1,'pixel_integrated_assumption':'approximate',
                 'no_instrument_calibration':True})
    return psf,np.asarray(accepted),meta


def audit_product(fits_path,catalog_path,out_dir,*,size=512,x0=None,y0=None,
                  stamp_size=15,fwhm=2.6,min_stars=4):
    """Produce QA files and PSF estimate; not a validated HST result."""
    image,good,err,info=read_science_crop(fits_path,size=size,x0=x0,y0=y0)
    cat=read_acsggct_catalog(catalog_path)
    count=len(cat['source_id'])
    good606=int(quality_mask(cat,band='f606w').sum())
    good814=int(quality_mask(cat,band='f814w').sum())
    psf,stars,psfinfo=empirical_psf(image,good,stamp_size=stamp_size,
                                     fwhm=fwhm,min_stars=min_stars)
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    np.save(out/'empirical_prf_native.npy',psf)
    with (out/'selected_psf_stars.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f);writer.writerow(['crop_x_px','crop_y_px']);writer.writerows(stars)
    if err is not None:
        errvals=err[good]
        e={'finite_err_pixel_fraction':float(np.isfinite(err).mean()),
           'median_err_plane_value':float(np.median(errvals)) if len(errvals) else None,
           'err_units_unverified':True}
    else:
        e={'has_error':False,'warning':'No ERR plane; do not infer calibrated uncertainty'}
    report={'version':'0.6.0','kind':'PRELIMINARY_REAL_DATA_AUDIT',
            'input':info,'catalog':{
                'file':str(catalog_path),'sha256':sha256_file(catalog_path),
                'row_count':int(count),'quality_f606w':good606,'quality_f814w':good814,
                'coordinate_reference_not_verified':True},
            'crop_psf':psfinfo,'error_diagnostics':e,
            'warnings':['Image/catalog WCS alignment NOT validated by this audit',
                        'Real HST star positions and PSF model must be validated separately',
                        'PSF native stamp is not an oversampled instrument ePSF',
                        'Photometric zeropoints, noise correlations and PSF variations unverified']}
    (out/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False,default=str),
                                   encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(5,4))
    ax.imshow(psf,origin='lower',cmap='magma')
    ax.set(title=f'Native-pixel empirical PRF — {len(stars)} stars',
           xlabel='x px',ylabel='y px')
    fig.tight_layout()
    fig.savefig(out/'psf_preview.png',dpi=150)
    plt.close(fig)
    return report
