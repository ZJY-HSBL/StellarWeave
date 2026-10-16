"""Spatial PSF and HST catalog-image astrometry *diagnostics*.

Scientific guardrails:
  - A raw WCS nearest-neighbor residual is not a verified astrometric fit.
  - A low residual in a crowded stellar field can be a chance association.
  - A native-pixel median stellar stamp is NOT an oversampled instrumental ePSF.
  - Everything is diagnostic, with separate statuses and reproducible source hashes.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .catalog import read_acsggct_catalog, project_catalog, quality_mask
from .detection import scipy_peak_detector
from .real_hst import empirical_psf, read_science_crop
from .astrometry import fit_affine, apply_affine
from .robust_plate import plate_diagnostics


def match_unique(reference_xy, detected_xy, radius=2.5):
    """Greedy globally sorted one-to-one associations below a distance cut."""
    r=np.asarray(reference_xy,float).reshape(-1,2)
    d=np.asarray(detected_xy,float).reshape(-1,2)
    if radius<=0:
        raise ValueError('Match radius must be positive')
    if not len(r) or not len(d):
        return np.empty((0,3),float)
    nearby=cKDTree(d).query_ball_point(r,radius)
    edges=sorted((float(np.linalg.norm(r[i]-d[j])),i,j)
                 for i,indices in enumerate(nearby) for j in indices)
    taken_r,taken_d=set(),set()
    found=[]
    for dist,i,j in edges:
        if i not in taken_r and j not in taken_d:
            found.append([float(i),float(j),dist])
            taken_r.add(i)
            taken_d.add(j)
    return np.asarray(found,float).reshape(-1,3)


def correspondence_diagnostics(image,valid,projected_xy,mag,*,radius=3.0,
                               max_reference=350,min_inliers=12):
    """Report raw WCS matches and a conservatively QC-filtered affine proposal.

    Result remains only a candidate astrometric correspondence even when QC
    passes. This routine does NOT authorize export of science-quality labels.
    """
    image=np.asarray(image,float)
    valid=np.asarray(valid,bool)
    points=np.asarray(projected_xy,float).reshape(-1,2)
    mags=np.asarray(mag,float).reshape(-1)
    if image.ndim!=2 or valid.shape!=image.shape or len(points)!=len(mags):
        raise ValueError('Image, mask and reference positions mismatch')
    h,w=image.shape
    margin=max(radius+4,8)
    inside=(np.isfinite(points).all(axis=1)&np.isfinite(mags)
            &(points[:,0]>=margin)&(points[:,0]<w-margin)
            &(points[:,1]>=margin)&(points[:,1]<h-margin))
    ref=points[inside]
    mags=mags[inside]
    if len(ref):
        select=np.argsort(mags)[:max_reference]
        ref=ref[select]
    detected=scipy_peak_detector(np.where(valid,image,np.nan),nsigma=5)
    r={'status':'INSUFFICIENT_DATA','reference_in_crop':int(inside.sum()),
       'reference_used':int(len(ref)),'detected_peaks':int(len(detected)),
       'raw_wcs_match_radius_px':float(radius),
       'astrometric_science_validated':False}
    matched=match_unique(ref,detected,radius=radius)
    r['raw_wcs_match_count']=int(len(matched))
    r['raw_wcs_match_fraction']=float(len(matched)/len(ref)) if len(ref) else None
    r['raw_wcs_median_nearest_px']=(float(np.median(matched[:,2]))
                                     if len(matched) else None)
    r['raw_wcs_p95_px']=(float(np.percentile(matched[:,2],95))
                           if len(matched) else None)
    r['affine_proposal']=None
    if len(ref)<min_inliers or len(detected)<min_inliers:
        r['reason']='Too few candidate catalog stars or detected peaks'
        return r,ref,detected,np.array([[1,0,0],[0,1,0]],float)
    try:
        # WCS registration should ordinarily start close to identity. This is
        # not an unrestricted blind-astrometry plate solver.
        a=fit_affine(ref,detected,initial_shift=[0.,0.],
                     match_radius=radius,clip_radius=min(radius,1.6),
                     min_inliers=min_inliers,min_fraction=.20,
                     max_rms=.9,max_p95=1.5)
        r['affine_proposal']=a
        r['status']='CANDIDATE_AFFINE_QC_PASSED'
        return r,ref,detected,np.asarray(a['matrix'])
    except (ValueError,IndexError,np.linalg.LinAlgError) as exc:
        r['status']='AFFINE_QC_REJECTED'
        r['reason']=str(exc)
        return r,ref,detected,np.array([[1,0,0],[0,1,0]],float)


def spatial_psf_grid(image,valid,*,cells=2,stamp_size=15,fwhm=2.6,
                     min_stars=3):
    """Local 2x2 or 3x3 PSF *candidate* stamps; no interpolation claimed."""
    a=np.asarray(image,float); mask=np.asarray(valid,bool)
    if a.ndim!=2 or mask.shape!=a.shape or cells not in (2,3):
        raise ValueError('Expected 2D image/mask and cells=2 or 3')
    h,w=a.shape
    out=[]
    psfs=[]
    for iy in range(cells):
        for ix in range(cells):
            x0=int(np.floor(ix*w/cells));x1=int(np.floor((ix+1)*w/cells))
            y0=int(np.floor(iy*h/cells));y1=int(np.floor((iy+1)*h/cells))
            crop=a[y0:y1,x0:x1]
            good=mask[y0:y1,x0:x1]
            row={'x_cell':ix,'y_cell':iy,'pixel_bounds':[x0,y0,x1,y1],
                 'valid_fraction':float(good.mean()),'status':'REJECTED'}
            try:
                psf,stars,details=empirical_psf(crop,good,
                        stamp_size=stamp_size,fwhm=fwhm,
                        min_stars=min_stars,max_stars=70)
                xx,yy=np.meshgrid(np.arange(psf.shape[1]),np.arange(psf.shape[0]))
                cx=(psf*xx).sum();cy=(psf*yy).sum()
                var_x=float((psf*(xx-cx)**2).sum())
                var_y=float((psf*(yy-cy)**2).sum())
                row.update(status='CANDIDATE',used_stars=int(len(stars)),
                           psf_second_moment_x_px2=var_x,
                           psf_second_moment_y_px2=var_y)
                psfs.append(psf)
            except ValueError as exc:
                psfs.append(None)
                row.update(reason=str(exc),used_stars=0)
            out.append(row)
    passed=[p for p in psfs if p is not None]
    if passed:
        reference=np.median(np.stack(passed),axis=0)
        reference=np.maximum(reference,0)
        reference/=reference.sum()
        for entry,psf in zip(out,psfs):
            if psf is not None:
                entry['l1_distance_to_median_psf']=float(np.abs(psf-reference).sum())
    return out,psfs


def projected_crop_positions(fits_path,table,meta,selected):
    """WCS sky coordinates to crop-local 0-based pixels, without assumed shifts."""
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(fits_path,memmap=True) as h:
        i=int(meta['sci_hdu_index'])
        wcs=WCS(h[i].header,h).celestial
        if not wcs.has_celestial or wcs.pixel_n_dim!=2:
            raise ValueError('Selected FITS extension has no 2-D celestial WCS')
        positions=project_catalog(table,wcs,selected)
    return positions-np.asarray(meta['crop_xy'],float)


def write_spatial_report(fits_path,catalog_path,out_dir,*,size=512,x0=None,
                         y0=None,cells=2,stamp_size=15,fwhm=2.6,
                         min_stars=3,matching_radius=3.0):
    """Generate versioned exploratory registration and local-PSF QC artifacts."""
    image,valid,err,meta=read_science_crop(fits_path,size,x0,y0)
    catalog=read_acsggct_catalog(catalog_path)
    quality=quality_mask(catalog,'f606w')
    positions=projected_crop_positions(fits_path,catalog,meta,quality)
    mag=np.asarray(catalog['mag_f606w'][quality])
    registration,refs,detected,matrix=correspondence_diagnostics(
        image,valid,positions,mag,radius=matching_radius)
    # Separate wide-shift search and held-out source validation, independent
    # from the legacy narrow-WCS nearest-neighbor diagnostic.
    h,w=image.shape
    pos=np.asarray(positions,float).reshape(-1,2)
    onfield=(np.isfinite(pos).all(axis=1)&np.isfinite(mag)
             &(pos[:,0]>8)&(pos[:,0]<w-8)
             &(pos[:,1]>8)&(pos[:,1]<h-8))
    candidates=pos[onfield]
    candidatemag=mag[onfield]
    if len(candidates)>500:
        candidates=candidates[np.argsort(candidatemag)[:500]]
    wide_plate=plate_diagnostics(
        candidates,detected,shape=image.shape,
        max_shift_px=192.,train_limit=260,holdout_limit=150,seed=42,
        initial_radius=3.5,match_radius=1.4)
    grid,stamps=spatial_psf_grid(image,valid,cells=cells,
                                stamp_size=stamp_size,fwhm=fwhm,
                                min_stars=min_stars)
    out=Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    for r,p in zip(grid,stamps):
        if p is not None:
            np.save(out/f"psf_y{r['y_cell']}_x{r['x_cell']}.npy",p)
    report={'version':'0.7.0','kind':'EXPLORATORY_SPATIAL_HST_DIAGNOSTIC',
            'image_source':meta['file'],'image_sha256':meta['sha256'],
            'catalog_source':str(catalog_path),
            'catalog_rows':len(catalog['source_id']),
            'quality_f606w':int(quality.sum()),
            'crop_origin_xy':meta['crop_xy'],
            'catalog_wcs_registration':registration,
            'wide_shift_heldout_plate_search':wide_plate,
            'psf_grid':grid,'psf_candidate_cells':sum(x['status']=='CANDIDATE' for x in grid),
            'scientific_astrometry_validated':False,
            'instrument_grade_epsf':False,
            'science_photometry_claim':False,
            'uncertainty_model_calibrated':False,
            'notes':['High stellar density allows chance nearest-neighbor matches',
                     'WCS residual matching is not independent external astrometric calibration',
                     'PSF variations are crude native-pixel candidate differences',
                     'ERR extension availability is not sufficient to model drizzle covariance']}
    (out/'spatial_diagnostics.json').write_text(
        json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(cells,cells,figsize=(6.4,6.4),squeeze=False)
    for ax,cell,stamp in zip(axs.ravel(),grid,stamps):
        if stamp is not None:
            ax.imshow(stamp,origin='lower',cmap='magma')
            ax.set_title(f"({cell['x_cell']},{cell['y_cell']}) n={cell['used_stars']}")
        else:
            ax.text(.5,.5,"REJECTED",ha='center',va='center',transform=ax.transAxes)
            ax.set_title(f"({cell['x_cell']},{cell['y_cell']})")
        ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle('Spatial native-pixel PSF candidates (NOT calibrated ePSF)')
    fig.tight_layout()
    fig.savefig(out/'spatial_psf_grid.png',dpi=150)
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(7,7))
    finite=image[valid]
    low,high=np.percentile(finite,[5,99.8]) if len(finite) else (0,1)
    ax.imshow(image,cmap='gray',origin='lower',vmin=low,
              vmax=max(high,low+1e-8))
    if len(refs):
        # Always show raw catalog WCS independently of any proposed fit.
        ax.scatter(refs[:,0],refs[:,1],s=18,facecolors='none',
                   edgecolors='cyan',linewidths=.45,label='Raw WCS reference')
        corrected=wide_plate.get('training_affine')
        if wide_plate['status']=='HELDOUT_QC_PASS' and corrected is not None:
            proposed=apply_affine(refs,corrected['matrix'])
            ax.scatter(proposed[:,0],proposed[:,1],s=27,facecolors='none',
                       edgecolors='lime',linewidths=.65,
                       label='Wide-shift affine (held-out QC passed)')
        elif registration['status']=='CANDIDATE_AFFINE_QC_PASSED':
            transformed=apply_affine(refs,matrix)
            ax.scatter(transformed[:,0],transformed[:,1],s=27,facecolors='none',
                       edgecolors='yellow',linewidths=.65,
                       label='Narrow-match affine candidate')
    if len(detected):
        sample=detected[:400]
        ax.scatter(sample[:,0],sample[:,1],s=5,color='orange',
                   alpha=.35,label='Image peaks (first 400)')
    ax.legend(loc='upper right',fontsize=7)
    ax.set_title('Raw WCS / candidate affine / image peaks — not absolute calibration')
    fig.tight_layout()
    fig.savefig(out/'astrometry_overlay.png',dpi=155)
    plt.close(fig)
    return report
