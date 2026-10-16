"""Conservative WCS-guided plate registration with held-out verification.

This is not blind astrometry: the initial WCS should have approximately
correct orientation/scale and shift by at most max_shift_px. Search for
translation hypotheses via sparse displacement voting, fit an affine model
on TRAIN sources only, then evaluate distinct withheld reference sources.

A dense field can create spurious associations: low residual alone is not
proof of correspondence. The null-shift controls are empirical diagnostics,
not statistically calibrated false-match probabilities.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree
from .astrometry import apply_affine


def unique_links(reference_xy, detected_xy, radius):
    """Global nearest-edge greedy one-to-one pairs; indices stay tied to input."""
    a=np.asarray(reference_xy,dtype=float).reshape(-1,2)
    b=np.asarray(detected_xy,dtype=float).reshape(-1,2)
    if radius<=0:
        raise ValueError('radius must be positive')
    if not len(a) or not len(b):
        return np.empty((0,3),dtype=float)
    pairs=cKDTree(b).query_ball_point(a,radius)
    edges=sorted((float(np.linalg.norm(a[i]-b[j])),i,j)
                 for i,js in enumerate(pairs) for j in js)
    seen_a=set()
    seen_b=set()
    matched=[]
    for d,i,j in edges:
        if i in seen_a or j in seen_b:
            continue
        seen_a.add(i)
        seen_b.add(j)
        matched.append([float(i),float(j),d])
    return np.asarray(matched,dtype=float).reshape(-1,3)


def translation_hypotheses(reference_xy,detected_xy,*,max_shift_px=192.,
                           bin_width=2.,max_hypotheses=12):
    """Sparse 2-D displacement voting: several coarse shifts, no label fitting."""
    a=np.asarray(reference_xy,float).reshape(-1,2)
    b=np.asarray(detected_xy,float).reshape(-1,2)
    if not len(a) or not len(b):
        return []
    if max_shift_px<=0 or bin_width<=0 or max_hypotheses<1:
        raise ValueError('Bad shift voting parameters')
    # Bound search complexity and retain input order for reproducibility.
    if len(a)*len(b)>2_000_000:
        raise ValueError('Too many displacement pairs; reduce the source subsets')
    deltas=(b[None,:,:]-a[:,None,:]).reshape(-1,2)
    allowed=np.all(np.abs(deltas)<=max_shift_px,axis=1)
    votes=deltas[allowed]
    if not len(votes):
        return []
    from scipy.ndimage import maximum_filter
    extent=max_shift_px
    bins=int(np.ceil((2*extent)/bin_width))
    hist,edges_x,edges_y=np.histogram2d(
        votes[:,0],votes[:,1],bins=bins,
        range=[[-extent,extent],[-extent,extent]])
    localmax=maximum_filter(hist,size=5,mode='constant')
    peaks=np.argwhere((hist==localmax)&(hist>0))
    peaks=sorted(peaks,key=lambda ij:(-hist[tuple(ij)],int(ij[0]),int(ij[1])))
    hypotheses=[]
    for ix,iy in peaks:
        center=np.array([.5*(edges_x[ix]+edges_x[ix+1]),
                         .5*(edges_y[iy]+edges_y[iy+1])])
        # Pool neighboring bins to recover sub-bin translation.
        nearby=np.linalg.norm(votes-center,axis=1)<=max(2*bin_width,3.0)
        if not nearby.any():
            continue
        center=np.median(votes[nearby],axis=0)
        if any(np.linalg.norm(center-p['shift'])<2.5*bin_width
               for p in hypotheses):
            continue
        hypotheses.append({'shift':center.tolist(),
                           'bin_votes':int(hist[ix,iy]),
                           'support_votes':int(nearby.sum())})
        if len(hypotheses)>=max_hypotheses:
            break
    return hypotheses


def _affine_from_pairs(a,b):
    if len(a)<3:
        raise ValueError('Insufficient affine pairs')
    design=np.column_stack((a,np.ones(len(a))))
    if np.linalg.matrix_rank(design)<3:
        raise ValueError('Degenerate affine geometry')
    coefficients,*_=np.linalg.lstsq(design,b,rcond=None)
    return coefficients.T


def _fit_from_seed(reference,detected,shift,*,initial_radius=3.5,
                   residual_radius=1.4,max_scale_error=.07):
    """Iterative mutual-uniqueness constrained fit; no withheld sources used."""
    transform=np.array([[1.,0.,shift[0]],[0.,1.,shift[1]]])
    best=None
    for t in range(8):
        r=initial_radius if t<2 else max(2.0,residual_radius*1.6)
        pairs=unique_links(apply_affine(reference,transform),detected,r)
        if len(pairs)<8:
            break
        i=pairs[:,0].astype(int)
        j=pairs[:,1].astype(int)
        try:
            candidate=_affine_from_pairs(reference[i],detected[j])
        except ValueError:
            break
        # Reject geometrically impossible fits before re-associating.
        matrix=candidate[:,:2]
        singular=np.linalg.svd(matrix,compute_uv=False)
        if np.linalg.det(matrix)<=0 or np.max(np.abs(singular-1.))>max_scale_error:
            break
        residual=np.linalg.norm(apply_affine(reference[i],candidate)-detected[j],axis=1)
        clipped=residual<max(residual_radius*1.5,np.median(residual)*2.)
        if clipped.sum()<8:
            break
        try:
            new_transform=_affine_from_pairs(reference[i[clipped]],detected[j[clipped]])
        except ValueError:
            break
        transform=new_transform
        best=transform.copy()
    if best is None:
        return None
    pairs=unique_links(apply_affine(reference,best),detected,residual_radius)
    if len(pairs)<8:
        return None
    i=pairs[:,0].astype(int)
    j=pairs[:,1].astype(int)
    best=_affine_from_pairs(reference[i],detected[j])
    singular=np.linalg.svd(best[:,:2],compute_uv=False)
    if np.linalg.det(best[:,:2])<=0 or np.max(np.abs(singular-1.))>max_scale_error:
        return None
    matched=unique_links(apply_affine(reference,best),detected,residual_radius)
    if len(matched)<8:
        return None
    xy=reference[matched[:,0].astype(int)]
    span=np.ptp(xy,axis=0)
    allspan=np.maximum(np.ptp(reference,axis=0),1.)
    coverage=span/allspan
    error=matched[:,2]
    return {'matrix':best.tolist(),'matched':int(len(matched)),
            'match_fraction':float(len(matched)/len(reference)),
            'rms_px':float(np.sqrt(np.mean(error**2))),
            'p95_px':float(np.percentile(error,95)),
            'median_px':float(np.median(error)),
            'coverage_fraction_xy':coverage.tolist(),
            'scale_singular_values':singular.tolist()}


def plate_diagnostics(reference_xy,detected_xy,*,shape=None,
                      max_shift_px=192.,train_limit=320,
                      holdout_limit=180,seed=42,
                      initial_radius=3.5,match_radius=1.4,
                      min_train_matches=20,min_holdout_matches=12,
                      min_holdout_fraction=.12,min_null_advantage=3.):
    """Find a proposed affine using train only, independently verify on holdout.

    Inputs are WCS-projected references and detected image peaks. Both may
    have unmatched objects. The output status is conservative and a pass is
    not equivalent to absolute astrometric calibration.
    """
    a=np.asarray(reference_xy,float).reshape(-1,2)
    b=np.asarray(detected_xy,float).reshape(-1,2)
    if train_limit<min_train_matches or holdout_limit<min_holdout_matches:
        raise ValueError('Reference limits smaller than minimum match requirements')
    if len(a)<train_limit//4 or len(b)<min_train_matches:
        return {'status':'INSUFFICIENT_SOURCES','reference_count':len(a),
                'detected_count':len(b),'science_calibrated':False}
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Sources must have finite coordinates')
    rng=np.random.default_rng(seed)
    # Caller preselects brighter sources. Randomize *within that selection*
    # for disjoint fit/validation sets.
    n=min(len(a),train_limit+holdout_limit)
    sample=rng.permutation(len(a))[:n]
    # For sparse tiles, an unconditional 150-source holdout previously left
    # only one training reference. Preserve a meaningful train/holdout split.
    hold_target=min(holdout_limit,max(min_holdout_matches,n//3))
    n_train=min(train_limit,n-hold_target)
    if n_train<min_train_matches:
        return {'status':'INSUFFICIENT_TRAIN','reference_count':len(a),
                'training_reference_count':n_train,'holdout_count':n-n_train,
                'science_calibrated':False}
    train_indices=sample[:n_train]
    hold_indices=sample[n_train:]
    if len(hold_indices)<min_holdout_matches:
        return {'status':'INSUFFICIENT_HOLDOUT','reference_count':len(a),
                'holdout_count':len(hold_indices),'science_calibrated':False}
    train=a[train_indices]
    hold=a[hold_indices]
    hypotheses=translation_hypotheses(train,b,max_shift_px=max_shift_px)
    models=[]
    # Include WCS identity in comparisons, but do not use heldout results
    # to select the fitted matrix.
    seeds=[{'shift':[0.,0.],'bin_votes':0,'support_votes':0}]+hypotheses
    for hypothesis in seeds:
        candidate=_fit_from_seed(train,b,hypothesis['shift'],
                  initial_radius=initial_radius,residual_radius=match_radius)
        if candidate:
            models.append((candidate,hypothesis))
    if not models:
        return {'status':'NO_RELIABLE_TRAIN_FIT',
                'train_count':len(train),'holdout_count':len(hold),
                'candidate_shifts':hypotheses[:7],
                'science_calibrated':False}
    # All model selection uses TRAIN statistics only.
    selected,hypothesis=max(models,key=lambda x:(x[0]['matched'],-x[0]['rms_px']))
    matrix=np.asarray(selected['matrix'])
    held=unique_links(apply_affine(hold,matrix),b,match_radius)
    # Disjoint reference IDs, but detected peaks might be shared across the
    # fit and holdout sets in a crowd. Record this limitation.
    null_counts=[]
    for dx,dy in ((39,0),(-39,0),(0,39),(0,-39),(27,31),(-27,-31)):
        shifted=apply_affine(hold,matrix)+np.array([dx,dy])
        if shape is not None:
            height,width=shape
            inside=((shifted[:,0]>=0)&(shifted[:,0]<width)&
                    (shifted[:,1]>=0)&(shifted[:,1]<height))
            shifted=shifted[inside]
        count=len(unique_links(shifted,b,match_radius))
        null_counts.append(int(count))
    null_median=float(np.median(null_counts))
    held_fraction=len(held)/len(hold)
    held_median=float(np.median(held[:,2])) if len(held) else None
    passed=(
        selected['matched']>=min_train_matches
        and min(selected['coverage_fraction_xy'])>=.40
        and selected['rms_px']<=.85 and selected['p95_px']<=1.35
        and len(held)>=min_holdout_matches
        and held_fraction>=min_holdout_fraction
        and len(held)>=min_null_advantage*max(1.,null_median)
        and held_median is not None and held_median<=1.1
    )
    return {'status':'HELDOUT_QC_PASS' if passed else 'HELDOUT_QC_REJECTED',
            'science_calibrated':False,'training_reference_count':int(len(train)),
            'holdout_reference_count':int(len(hold)),
            'detected_count':int(len(b)),
            'selected_seed':hypothesis,'training_affine':selected,
            'holdout_matches':int(len(held)),
            'holdout_match_fraction':float(held_fraction),
            'holdout_median_residual_px':held_median,
            'holdout_p95_residual_px':float(np.percentile(held[:,2],95)) if len(held) else None,
            'null_shift_counts':null_counts,'null_shift_median':null_median,
            'holdout_vs_null_ratio':float(len(held)/max(1.,null_median)),
            'candidate_seed_count':len(seeds),
            'note':'No independent external astrometric calibration; holdout matches are internal corroboration, not ground truth.'}
