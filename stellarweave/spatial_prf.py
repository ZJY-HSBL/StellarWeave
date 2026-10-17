"""Cross-validated spatial native-pixel PRF diagnostics on real HST mosaics.

An empirical PRF (pixel response function) is NOT an oversampled instrument
ePSF. This experiment measures normalized held-out star-stamp shape error,
NOT calibrated flux accuracy, astrometry, completeness, or detection recall.

The nine disjoint 512px tiles are taken from an independently generated,
predefined cross_region.json geometry. Each tile splits its candidate stars
into train/test ONCE with a fixed seed. The anchor PRF is fitted only to
center-tile training stars; a local PRF uses only that tile's training stars.
Every comparison is scored on held-out stars, never fitting test stamps.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import shift as shift_image

from .real_hst import _psf_candidates, sha256_file

STAMP_SIZE = 15
ANCHOR_INDEX = 4
TILE_SIZE = 512


def extract_stamps(image, valid, *, stamp_size=STAMP_SIZE,
                   fwhm=2.6, min_neighbors=3.0, cap=160):
    """Find and centroid-normalize isolated native-pixel star cutouts.

    Stars may still be blended or saturated: this is a candidate-based
    diagnostic, not a validated astronomical ePSF calibration.
    """
    image=np.asarray(image,float)
    valid=np.asarray(valid,bool)
    if image.ndim!=2 or valid.shape!=image.shape:
        raise ValueError('Image and valid mask must be matching 2D arrays')
    xy,meta=_psf_candidates(image,valid,stamp_size=stamp_size,
                             fwhm=fwhm,min_neighbors=min_neighbors,cap=cap)
    r=stamp_size//2
    yy,xx=np.mgrid[-r:r+1,-r:r+1]
    ring=(np.hypot(xx,yy)>=.75*r)
    gx,gy=np.meshgrid(np.arange(stamp_size),np.arange(stamp_size))
    stamps=[]
    kept=[]
    for x,y in xy:
        ix,iy=int(round(x)),int(round(y))
        raw=image[iy-r:iy+r+1,ix-r:ix+r+1]
        if raw.shape!=(stamp_size,stamp_size):
            continue
        sky=float(np.median(raw[ring]))
        patch=np.maximum(raw-sky,0.)
        total=float(patch.sum())
        if total<=0 or not np.isfinite(total):
            continue
        cx=float(np.sum(gx*patch)/total)
        cy=float(np.sum(gy*patch)/total)
        if max(abs(cx-r),abs(cy-r))>1.5:
            continue
        centered=shift_image(patch,(r-cy,r-cx),order=3,
                             mode='constant',cval=0,prefilter=True)
        centered=np.maximum(centered,0)
        norm=float(centered.sum())
        if norm<=0 or not np.isfinite(norm):
            continue
        stamps.append(centered/norm)
        kept.append((float(x),float(y)))
    meta=dict(meta,usable_stamps=len(stamps),
              stamp_size=stamp_size,
              min_neighbors_fwhm=min_neighbors,
              star_selection='image-detected isolated bright candidates',
              oversampling=1,
              instrument_ePSf_calibrated=False)
    return np.asarray(stamps,dtype=float).reshape(-1,stamp_size,stamp_size), \
           np.asarray(kept,dtype=float).reshape(-1,2),meta


def template_from_training(stamps):
    """Robust native-pixel stack; never accepts held-out stamps."""
    x=np.asarray(stamps,float)
    if x.ndim!=3 or len(x)<2 or x.shape[1]!=x.shape[2]:
        raise ValueError('Need >=2 square training stamps')
    if not np.isfinite(x).all() or np.any(x<0):
        raise ValueError('Training stamps must be finite and nonnegative')
    template=np.median(x,axis=0)
    template=np.maximum(template,0)
    if template.sum()<=0:
        raise ValueError('Empty PRF template')
    return template/template.sum()


def gaussian_reference(stamp_size=STAMP_SIZE,fwhm=2.6):
    """Fixed isotropic Gaussian baseline (not tuned on test stamps)."""
    r=stamp_size//2
    y,x=np.mgrid[-r:r+1,-r:r+1]
    sigma=fwhm/(2*np.sqrt(2*np.log(2)))
    model=np.exp(-.5*(x*x+y*y)/(sigma*sigma))
    return model/model.sum()


def split_training_test(stamps,xy,*,seed,min_train=8,min_test=3):
    """Predeclared 75/25 random split, no test-set model selection."""
    n=len(stamps)
    if n<min_train+min_test:
        return None
    rng=np.random.default_rng(seed)
    order=rng.permutation(n)
    n_test=min(n-min_train,max(min_test,int(round(n*.25))))
    train=order[n_test:]
    test=order[:n_test]
    return {'template':template_from_training(stamps[train]),
            'holdout_stamps':stamps[test],
            'holdout_xy':xy[test],
            'train_count':int(len(train)),
            'test_count':int(len(test))}


def _l1(stamps,template):
    """Per-star total variation-like L1 of unit-sum normalized pixel arrays."""
    return np.sum(np.abs(stamps-template),axis=(1,2))


def compare_tiles(tile_inputs,*,seed=42,stamp_size=STAMP_SIZE):
    """Pure, testable comparison of pre-extracted disjoint tile stamp arrays.

    tile_inputs: nine (stamps, coordinates, metadata) tuples in tile order.
    Returns report rows, held-out star scores and per-tile training PRFs.
    """
    if len(tile_inputs)!=9:
        raise ValueError('Exactly nine predefined spatial tiles are required')
    prepared=[]
    for idx,(stamps,xy,meta) in enumerate(tile_inputs):
        split=split_training_test(stamps,xy,seed=seed+idx)
        prepared.append((split,meta))
    anchor=prepared[ANCHOR_INDEX][0]
    if anchor is None:
        raise ValueError('Insufficient isolated center-tile stars for frozen anchor PRF')
    anchor_model=anchor['template']
    gaussian=gaussian_reference(stamp_size)
    rows=[]
    per_star=[]
    templates={}
    for idx,(split,meta) in enumerate(prepared):
        row={'index':idx,'role':'ANCHOR' if idx==ANCHOR_INDEX else 'HELDOUT_REGION',
             'detected_candidates':int(meta.get('candidates',0)),
             'isolated_candidates':int(meta.get('isolated',0)),
             'usable_stamps':int(meta.get('usable_stamps',0)),
             'train_stars':0,'heldout_stars':0,
             'status':'INSUFFICIENT_ISOLATED_STARS',
             'median_local_l1':None,'median_anchor_l1':None,
             'median_gaussian_l1':None,
             'median_paired_gain_pct':None}
        if split is None:
            rows.append(row)
            continue
        template=split['template']
        templates[idx]=template
        hold=split['holdout_stamps']
        loc=_l1(hold,template)
        anc=_l1(hold,anchor_model)
        gauss=_l1(hold,gaussian)
        gain=100*(anc-loc)/np.maximum(anc,1e-12)
        row.update({'status':'EVALUATED','train_stars':split['train_count'],
                    'heldout_stars':split['test_count'],
                    'median_local_l1':float(np.median(loc)),
                    'median_anchor_l1':float(np.median(anc)),
                    'median_gaussian_l1':float(np.median(gauss)),
                    'median_paired_gain_pct':float(np.median(gain))})
        for j,(xy,l,a,g,z) in enumerate(zip(split['holdout_xy'],loc,anc,gauss,gain)):
            per_star.append({'tile_index':idx,'holdout_index':j,
                             'tile_x':float(xy[0]),'tile_y':float(xy[1]),
                             'local_l1':float(l),'anchor_l1':float(a),
                             'gaussian_l1':float(g),
                             'paired_gain_pct':float(z)})
        rows.append(row)
    return rows,per_star,templates


def bootstrap_tile_median_ci(tile_rows,*,seed=2026,n_resamples=2000):
    """Bootstrap over held-out TILES, not correlated stars within tiles."""
    values=np.array([r['median_paired_gain_pct'] for r in tile_rows
                     if r['index']!=ANCHOR_INDEX and r['status']=='EVALUATED'],
                    dtype=float)
    if len(values)<4:
        return None
    rng=np.random.default_rng(seed)
    samples=rng.choice(values,size=(n_resamples,len(values)),replace=True)
    bounds=np.percentile(np.median(samples,axis=1),[2.5,97.5])
    return [float(bounds[0]),float(bounds[1])]


def audit_spatial_prf(fits_path,plate_path,out_dir,*,cluster=None,seed=42):
    """Run one real FITS mosaic and write full, non-cherry-picked evidence."""
    from astropy.io import fits
    fits_path=Path(fits_path)
    plate_path=Path(plate_path)
    with plate_path.open(encoding='utf-8') as f:
        plate=json.load(f)
    if (plate.get('tile_size')!=TILE_SIZE or plate.get('grid')!=3
            or plate.get('anchor_index')!=ANCHOR_INDEX
            or plate.get('tile_count')!=9
            or plate.get('heldout_tiles_physically_nonoverlapping') is not True):
        raise ValueError('Requires audited, predefined, disjoint 3x3 512px tiles')
    if plate.get('anchor_status')!='HELDOUT_QC_PASS':
        raise ValueError('Spatial geometry anchor has not passed QC')
    sha=sha256_file(fits_path)
    if sha!=plate.get('fits_sha256'):
        raise ValueError('FITS checksum differs from geometry provenance')
    inputs=[]
    valid_fractions=[]
    err_present=False
    dq_present=False
    with fits.open(fits_path,memmap=True) as hdus:
        sci=hdus[int(plate['science_hdu_index'])]
        if sci.data is None or getattr(sci.data,'ndim',0)!=2:
            raise ValueError('Selected science plane missing')
        ver=sci.header.get('EXTVER',1)
        err=dq=None
        for h in hdus:
            if h is sci or h.data is None or h.data.shape!=sci.data.shape:
                continue
            if h.header.get('EXTVER',1)!=ver:
                continue
            ext=h.header.get('EXTNAME','').upper()
            if ext=='ERR' and err is None:
                err=h.data
            if ext=='DQ' and dq is None:
                dq=h.data
        err_present=err is not None
        dq_present=dq is not None
        tiles=sorted(plate['tiles'],key=lambda t:t['index'])
        if [t['index'] for t in tiles]!=list(range(9)):
            raise ValueError('Missing or duplicate tile geometry')
        origins=[tuple(map(int,t['tile_xy'])) for t in tiles]
        for i,(x,y) in enumerate(origins):
            if x<0 or y<0 or x+TILE_SIZE>sci.data.shape[1] or y+TILE_SIZE>sci.data.shape[0]:
                raise ValueError('Tile outside science image')
            for x2,y2 in origins[:i]:
                if abs(x-x2)<TILE_SIZE and abs(y-y2)<TILE_SIZE:
                    raise ValueError('Overlapping tile geometry')
        for x,y in origins:
            ys=slice(y,y+TILE_SIZE)
            xs=slice(x,x+TILE_SIZE)
            patch=np.asarray(sci.data[ys,xs],dtype=float)
            valid=np.isfinite(patch)
            if dq is not None:
                valid&=(np.asarray(dq[ys,xs])==0)
            if err is not None:
                local_err=np.asarray(err[ys,xs],float)
                valid&=np.isfinite(local_err)&(local_err>0)
            valid_fractions.append(float(valid.mean()))
            inputs.append(extract_stamps(patch,valid))
    rows,per_star,templates=compare_tiles(inputs,seed=seed)
    for i,row in enumerate(rows):
        row['tile_xy']=list(origins[i])
        row['valid_pixel_fraction']=valid_fractions[i]
    evals=[r for r in rows if r['index']!=ANCHOR_INDEX and r['status']=='EVALUATED']
    gains=[r['median_paired_gain_pct'] for r in evals]
    result={'version':'0.12.0','kind':'REAL_HST_SPATIAL_NATIVE_PIXEL_PRF_HOLDOUT',
            'cluster':cluster or plate.get('cluster',fits_path.stem),
            'fits_sha256':sha,'plate_sha256':sha256_file(plate_path),
            'input_geometry':'same predeclared disjoint 3x3 tiles as internal astrometry',
            'tile_size':TILE_SIZE,'grid':3,'anchor_index':ANCHOR_INDEX,
            'seed':seed,'stamp_size':STAMP_SIZE,'candidate_fwhm_px':2.6,
            'isolation_min_neighbors_fwhm':3.,
            'model_training':'median stack of 75% local bright isolated candidates',
            'holdout':'25% disjoint candidate stars per tile; deterministic split',
            'metric':'sum(abs(unit-normalized heldout stamp - PRF template))',
            'anchor_model_fitted_on_center_train_only':True,
            'no_heldout_star_used_to_fit_a_scored_template':True,
            'anchor_and_local_comparison_same_test_stars':True,
            'gaussian_baseline_fwhm_px':2.6,
            'science_calibrated':False,
            'oversampled_instrument_epsf':False,
            'independent_gaia_validation':False,
            'calibrated_photometric_accuracy':False,
            'error_plane_available':err_present,'dq_plane_available':dq_present,
            'error_plane_units_verified':False,
            'evaluated_heldout_tiles':len(evals),
            'possible_heldout_tiles':8,
            'local_prf_better_tiles':sum(g>0 for g in gains),
            'median_tile_paired_gain_pct':float(np.median(gains)) if gains else None,
            'bootstrap_95ci_median_tile_gain_pct':bootstrap_tile_median_ci(rows,seed=seed+777),
            'total_heldout_stars':sum(r['heldout_stars'] for r in rows),
            'tiles':rows,'heldout_star_scores':per_star,
            'limitations':['PRF is native-pixel, not oversampled ePSF',
                           'Isolation/centroiding are approximate; saturated and blended sources possible',
                           'Held-out stamp shape error is not flux or positional accuracy',
                           'Each star is scored against local and anchor PRF without fitting its flux',
                           'Bootstrap resamples tiles, not independent stars',
                           'No Gaia absolute astrometry or completeness claim']}
    out=Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    (out/'spatial_prf.json').write_text(json.dumps(result,indent=2,allow_nan=False),
                                        encoding='utf-8')
    with (out/'heldout_stars.csv').open('w',newline='',encoding='utf-8') as f:
        fields=['tile_index','holdout_index','tile_x','tile_y',
                'local_l1','anchor_l1','gaussian_l1','paired_gain_pct']
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        w.writerows(per_star)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(3,3,figsize=(8.5,8.5))
    scale=max((float(m.max()) for m in templates.values()),default=1.)
    for idx,ax in enumerate(axs.ravel()):
        row=rows[idx]
        if idx in templates:
            ax.imshow(templates[idx],origin='lower',cmap='magma',
                      vmin=0,vmax=scale,interpolation='nearest')
        else:
            ax.text(.5,.5,'Insufficient stars',transform=ax.transAxes,
                    ha='center',va='center',fontsize=9)
        ax.set(title=f"Tile {idx}: {row['train_stars']} train / {row['heldout_stars']} test")
        ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle('Real HST native-pixel PRFs (train-only; not ePSF)',fontsize=12)
    fig.tight_layout()
    fig.savefig(out/'spatial_prf_templates.png',dpi=160)
    plt.close(fig)
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description='Cross-validate spatial native-pixel PRF on real HST')
    p.add_argument('--fits',required=True)
    p.add_argument('--plate-report',required=True)
    p.add_argument('--cluster')
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--out',required=True)
    a=p.parse_args(argv)
    d=audit_spatial_prf(a.fits,a.plate_report,a.out,
                        cluster=a.cluster,seed=a.seed)
    print(json.dumps({'cluster':d['cluster'],
        'evaluated_heldout_tiles':d['evaluated_heldout_tiles'],
        'local_prf_better_tiles':d['local_prf_better_tiles'],
        'median_tile_paired_gain_pct':d['median_tile_paired_gain_pct'],
        'bootstrap_95ci':d['bootstrap_95ci_median_tile_gain_pct'],
        'total_heldout_stars':d['total_heldout_stars'],
        'calibrated_photometry':False},indent=2))


if __name__=='__main__':
    main()
