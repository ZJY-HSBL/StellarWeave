"""Reproducible artificial-star injection/recovery (instrumental flux only).

This is a methodological baseline, NOT calibrated HST completeness. Real
images require a compatible pixel-integrated PSF stamp. Additive injection
does not reproduce extra Poisson shot noise or drizzle noise correlations.
"""
from __future__ import annotations
import numpy as np
from scipy.ndimage import shift as image_shift
from scipy.spatial import cKDTree
from scipy.special import ndtr


def gaussian_prf(fwhm=2.6, size=17):
    """Pixel-integrated, unit-normalized toy Gaussian response."""
    if fwhm <= 0 or size < 5 or size % 2 != 1:
        raise ValueError('Positive FWHM and odd size >= 5 required')
    sig=fwhm/2.354820045
    x=np.arange(size)-size//2
    part=ndtr((x+.5)/sig)-ndtr((x-.5)/sig)
    psf=np.outer(part,part)
    return psf/psf.sum()


def normalized_psf(stamp):
    k=np.asarray(stamp,dtype=float)
    if (k.ndim!=2 or k.shape[0]!=k.shape[1] or k.shape[0]%2!=1
            or k.shape[0]<5 or not np.isfinite(k).all()
            or np.any(k<0) or k.sum()<=0):
        raise ValueError('PSF must be finite, nonnegative, odd square >=5 and nonzero')
    k=k/k.sum()
    y,x=np.indices(k.shape)
    midpoint=(k.shape[0]-1)/2
    center=np.array([(x*k).sum(),(y*k).sum()])
    if np.max(np.abs(center-midpoint))>1:
        raise ValueError('PSF centroid differs from stamp center by >1 pixel')
    return k


def inject_stars(image,stars,psf,valid_mask=None):
    """Inject rows (x,y,total_flux) without clipping, conserving added flux."""
    original=np.asarray(image,dtype=float)
    if original.ndim!=2:
        raise ValueError('Expected two-dimensional image')
    valid=np.isfinite(original) if valid_mask is None else np.asarray(valid_mask,bool)&np.isfinite(original)
    if valid.shape!=original.shape:
        raise ValueError('Invalid mask shape')
    kernel=normalized_psf(psf)
    rows=np.asarray(stars,float).reshape(-1,3)
    if not np.isfinite(rows).all() or np.any(rows[:,2]<=0):
        raise ValueError('Finite coordinates and positive fluxes required')
    out=original.copy()
    radius=kernel.shape[0]//2
    for x,y,flux in rows:
        ix,iy=int(np.floor(x+.5)),int(np.floor(y+.5))
        x0,x1=ix-radius,ix+radius+1
        y0,y1=iy-radius,iy+radius+1
        if x0<0 or y0<0 or x1>out.shape[1] or y1>out.shape[0]:
            raise ValueError('Injection footprint extends beyond image border')
        if not valid[y0:y1,x0:x1].all():
            raise ValueError('Injection overlaps masked/nonfinite pixels')
        shifted=image_shift(kernel,(y-iy,x-ix),order=3,
                            mode='constant',cval=0,prefilter=True)
        shifted=np.maximum(shifted,0)
        if shifted.sum()<=0:
            raise ValueError('PSF interpolation failed')
        out[y0:y1,x0:x1]+=flux*shifted/shifted.sum()
    return out


def sample_positions(shape,valid_mask,count,radius,*,rng,
                     min_separation=12.,avoid_positions=None,
                     avoid_radius=0.,max_attempts=250000):
    """Reproducibly sample positions with complete valid injection footprints."""
    h,w=map(int,shape)
    valid=np.asarray(valid_mask,bool)
    if valid.shape!=(h,w) or count<0 or radius<1:
        raise ValueError('Invalid sample geometry')
    if count==0:
        return np.empty((0,2))
    if h<=2*radius+2 or w<=2*radius+2:
        raise ValueError('Image too small for PSF')
    avoid=np.asarray(avoid_positions if avoid_positions is not None else [],float).reshape(-1,2)
    tree=cKDTree(avoid) if len(avoid) else None
    positions=[]
    for _ in range(max_attempts):
        x=rng.uniform(radius+1,w-radius-2)
        y=rng.uniform(radius+1,h-radius-2)
        ix,iy=int(np.floor(x+.5)),int(np.floor(y+.5))
        if not valid[iy-radius:iy+radius+1,ix-radius:ix+radius+1].all():
            continue
        if tree is not None and tree.query([x,y])[0]<avoid_radius:
            continue
        if positions and np.min(np.sum((np.asarray(positions)-[x,y])**2,axis=1))<min_separation**2:
            continue
        positions.append([float(x),float(y)])
        if len(positions)==count:
            return np.asarray(positions)
    raise ValueError('Cannot place stars with requested mask/separation')


def _greedy_unique_match(a,b,radius):
    a=np.asarray(a,float).reshape(-1,2)
    b=np.asarray(b,float).reshape(-1,2)
    if not len(a) or not len(b):
        return []
    lists=cKDTree(b).query_ball_point(a,radius)
    edges=sorted((float(np.linalg.norm(a[i]-b[j])),i,j)
                 for i,js in enumerate(lists) for j in js)
    used_a,used_b,found=set(),set(),[]
    for dist,i,j in edges:
        if i in used_a or j in used_b:
            continue
        used_a.add(i)
        used_b.add(j)
        found.append((i,j,dist))
    return found


def recover_injected(injected_xy,before_xy,after_xy,radius=2.,
                     baseline_exclusion_radius=2.):
    """Credit only NEW detections, excluding candidates matching baseline stars.

    This peak-based protocol is intentionally conservative for blends and
    does not evaluate PSF-model deblending performance.
    """
    injection=np.asarray(injected_xy,float).reshape(-1,2)
    before=np.asarray(before_xy,float).reshape(-1,2)
    after=np.asarray(after_xy,float).reshape(-1,2)
    matched_old=_greedy_unique_match(before,after,baseline_exclusion_radius)
    old_ids={j for _,j,_ in matched_old}
    new=after[[i for i in range(len(after)) if i not in old_ids]]
    found=_greedy_unique_match(injection,new,radius)
    success=np.zeros(len(injection),bool)
    offsets=np.full(len(injection),np.nan)
    for i,_,d in found:
        success[i]=True
        offsets[i]=d
    return success,offsets,len(new)-len(found)


def wilson_interval(success,total,z=1.95996398454):
    """Approximate 95% Wilson interval for independently sampled injection sites."""
    if total<=0:
        return [None,None]
    p=success/total
    denom=1+z*z/total
    center=(p+z*z/(2*total))/denom
    half=z*np.sqrt(p*(1-p)/total+z*z/(4*total*total))/denom
    return [float(max(0,center-half)),float(min(1,center+half))]


def experiment(base_image,psf,flux_levels,*,detector,valid_mask=None,
               trials=3,stars_per_flux=12,seed=42,match_radius=2.,
               baseline_exclusion_radius=2.,fwhm_px=2.6,
               injection_min_separation=None,avoid_existing_radius=0.):
    """Fresh independent injections into fixed base image; summaries by flux/density."""
    image=np.asarray(base_image,float)
    valid=np.isfinite(image) if valid_mask is None else np.asarray(valid_mask,bool)&np.isfinite(image)
    kernel=normalized_psf(psf)
    if image.ndim!=2 or valid.shape!=image.shape or not np.isfinite(image[valid]).all():
        raise ValueError('Invalid image/mask')
    levels=np.asarray(flux_levels,float)
    if levels.ndim!=1 or not len(levels) or not np.isfinite(levels).all() or np.any(levels<=0):
        raise ValueError('Flux levels must be positive')
    if trials<1 or stars_per_flux<1 or fwhm_px<=0:
        raise ValueError('Trials, stars per flux and FWHM must be positive')
    baseline=np.asarray(detector(image),float).reshape(-1,2)
    base_tree=cKDTree(baseline) if len(baseline) else None
    rng=np.random.default_rng(seed)
    records,runs=[],[]
    separation=injection_min_separation or max(2*fwhm_px,kernel.shape[0]/2)
    for trial in range(trials):
        n=stars_per_flux*len(levels)
        pos=sample_positions(image.shape,valid,n,kernel.shape[0]//2,
                             rng=rng,min_separation=separation,
                             avoid_positions=baseline,avoid_radius=avoid_existing_radius)
        flux=np.repeat(levels,stars_per_flux)[rng.permutation(n)]
        modified=inject_stars(image,np.column_stack((pos,flux)),kernel,valid_mask=valid)
        after=np.asarray(detector(modified),float).reshape(-1,2)
        success,offsets,unmatched=recover_injected(
            pos,baseline,after,match_radius,baseline_exclusion_radius)
        nearest=base_tree.query(pos)[0] if base_tree is not None else np.full(n,np.inf)
        for i,(coord,f) in enumerate(zip(pos,flux)):
            d=float(nearest[i])
            crowd='blend' if d<1.5*fwhm_px else 'crowded' if d<3*fwhm_px else 'isolated'
            records.append({'trial':trial,'injection_id':int(i),'x_px':float(coord[0]),
                            'y_px':float(coord[1]),'flux':float(f),
                            'nearest_baseline_px':d if np.isfinite(d) else None,
                            'crowding':crowd,'recovered':bool(success[i]),
                            'position_error_px':float(offsets[i]) if success[i] else None})
        runs.append({'trial':trial,'baseline_detections':int(len(baseline)),
                     'after_detections':int(len(after)),
                     'unmatched_new_candidates_not_false_positives':int(unmatched),
                     'injected':int(n),'recovered':int(success.sum()),
                     'injected_flux_total':float(flux.sum()),
                     'image_flux_difference':float((modified-image).sum())})
    groups=[]
    for level in levels:
        for crowd in ('all','isolated','crowded','blend'):
            selected=[r for r in records if r['flux']==float(level)
                      and (crowd=='all' or r['crowding']==crowd)]
            n=len(selected)
            recovered=sum(r['recovered'] for r in selected)
            groups.append({'flux':float(level),'crowding':crowd,
                           'injected':n,'recovered':recovered,
                           'completeness':recovered/n if n else None,
                           'wilson95':wilson_interval(recovered,n)})
    return {'baseline_detections':int(len(baseline)),'rows':records,
            'groups':groups,'trials':runs,'seed':int(seed)}
