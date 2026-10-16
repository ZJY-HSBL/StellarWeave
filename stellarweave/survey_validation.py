"""Cross-region and cross-cluster ACSGGCT plate-solution diagnostics (v0.9).

Unlike the v0.8 center-only test, this benchmark reserves physically disjoint
FITS tiles. An affine transform is chosen using *only the predefined central
tile*. Its fixed parameters are then evaluated on the other eight tiles,
WITHOUT re-fitting. Independent per-tile fits are separately labeled
exploratory and must never be confused with zero-shot transfer.

This is internal FITS/catalog agreement, not independent absolute astrometry.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import numpy as np

from .catalog import read_acsggct_catalog, quality_mask, project_catalog
from .detection import scipy_peak_detector
from .robust_plate import plate_diagnostics, unique_links
from .astrometry import apply_affine
from .real_hst import sha256_file


NULL_OFFSETS=((39.,0.),(-39.,0.),(0.,39.),(0.,-39.),(27.,31.),(-27.,-31.))


def tile_origins(shape,size,grid=3):
    """Generate a fixed 3x3 disjoint tiling, anchored to outer image margins."""
    h,w=map(int,shape)
    if size<64 or grid not in (2,3) or h<grid*size or w<grid*size:
        raise ValueError('Need image >= grid*size in both dimensions and size >=64')
    xs=[int(round(i*(w-size)/(grid-1))) for i in range(grid)]
    ys=[int(round(i*(h-size)/(grid-1))) for i in range(grid)]
    points=[(x,y) for y in ys for x in xs]
    for i,(x1,y1) in enumerate(points):
        for x2,y2 in points[i+1:]:
            if abs(x1-x2)<size and abs(y1-y2)<size:
                raise AssertionError('A proposed evaluation tiling overlaps')
    return points


def global_affine(crop_affine,origin):
    """Lift pixel crop affine to the full mosaic pixel coordinate system."""
    matrix=np.asarray(crop_affine,float)
    if matrix.shape!=(2,3):
        raise ValueError('Expected a 2x3 affine matrix')
    o=np.asarray(origin,float)
    a=matrix[:,:2]
    lifted=matrix.copy()
    lifted[:,2]=matrix[:,2]+o-a@o
    return lifted


def _clean_references(positions,magnitudes):
    p=np.asarray(positions,float).reshape(-1,2)
    m=np.asarray(magnitudes,float).reshape(-1)
    if len(p)!=len(m):
        raise ValueError('Mismatched WCS positions and magnitudes')
    valid=np.isfinite(p).all(axis=1)&np.isfinite(m)
    return p[valid],m[valid]


def _eligible(positions,magnitudes,origin,size,margin=8,max_stars=480):
    local=positions-np.asarray(origin,float)
    keep=((local[:,0]>=margin)&(local[:,0]<size-margin)
          &(local[:,1]>=margin)&(local[:,1]<size-margin))
    loc=local[keep]
    mags=magnitudes[keep]
    if len(loc)>max_stars:
        loc=loc[np.argsort(mags,kind='stable')[:max_stars]]
    return loc


def _offset_null_counts(points,peaks,size,radius):
    counts=[]
    for d in NULL_OFFSETS:
        p=points+np.asarray(d,float)
        inside=((p[:,0]>=0)&(p[:,0]<size)&(p[:,1]>=0)&(p[:,1]<size))
        counts.append(len(unique_links(p[inside],peaks,radius)))
    return counts


def transfer_score(global_reference, magnitudes, full_affine, origin,
                   detected, size, *, match_radius=1.4,max_stars=300):
    """Assess frozen center-tile transform on a non-overlapping spatial tile.

    Reference coordinates and source-selection use transformed target positions,
    never test-tile fit optimization. The stellar list is from the same
    pipeline catalog; not an independent Gaia/HST external catalog.
    """
    if full_affine is None:
        return {'status':'NO_ANCHOR_FIT','evaluated':0,'matches':0}
    transformed=apply_affine(global_reference,full_affine)
    local=_eligible(transformed,magnitudes,origin,size,margin=8,
                    max_stars=max_stars)
    if len(local)<12:
        return {'status':'INSUFFICIENT_TEST_REFERENCES',
                'evaluated':len(local),'matches':0}
    matched=unique_links(local,detected,match_radius)
    null=_offset_null_counts(local,detected,size,match_radius)
    null_mid=float(np.median(null))
    ratio=len(matched)/max(1.,null_mid)
    median=float(np.median(matched[:,2])) if len(matched) else None
    p95=float(np.percentile(matched[:,2],95)) if len(matched) else None
    fraction=len(matched)/len(local)
    # This is a strictly non-adaptive acceptance rule; any failure is retained.
    status=('TRANSFER_QC_PASS'
            if len(matched)>=12 and fraction>=.25
               and ratio>=3. and median is not None and median<=1.1
               and p95 is not None and p95<=1.4
            else 'TRANSFER_QC_REJECTED')
    return {'status':status,'evaluated':int(len(local)),
            'matches':int(len(matched)),'match_fraction':float(fraction),
            'median_residual_px':median,'p95_residual_px':p95,
            'null_counts':null,'null_median':null_mid,
            'match_over_null_ratio':float(ratio),'refitted':False}


def assess_mosaic(image,valid,reference_xy,magnitudes,*,tile_size=512,grid=3,
                  seed=42,match_radius=1.4,max_shift=192.):
    """Spatially disjoint center-anchor plus frozen-transfer tests."""
    im=np.asarray(image,float)
    good=np.asarray(valid,bool)
    if im.ndim!=2 or good.shape!=im.shape:
        raise ValueError('Invalid image/valid mask geometry')
    positions,mags=_clean_references(reference_xy,magnitudes)
    origins=tile_origins(im.shape,tile_size,grid=grid)
    center=(grid*grid)//2
    tiles=[]
    # Per-tile local fits are diagnostics only; *never* use them to
    # re-parameterize the frozen center anchor tested on other tiles.
    for index,(x,y) in enumerate(origins):
        patch=im[y:y+tile_size,x:x+tile_size]
        mask=good[y:y+tile_size,x:x+tile_size]
        peaks=scipy_peak_detector(np.where(mask,patch,np.nan),
                                  fwhm=2.6,nsigma=5.)
        local=_eligible(positions,mags,(x,y),tile_size,max_stars=480)
        proposal=plate_diagnostics(local,peaks,shape=patch.shape,
                     seed=seed+index,max_shift_px=max_shift,
                     train_limit=260,holdout_limit=150,
                     initial_radius=3.5,match_radius=match_radius,
                     min_train_matches=20,min_holdout_matches=12)
        tiles.append({'index':index,'role':'ANCHOR_TRAIN' if index==center else 'HELDOUT_REGION',
                      'tile_xy':[x,y],'detected_peaks':len(peaks),
                      'reference_candidates':len(local),
                      'valid_fraction':float(mask.mean()),
                      'local_candidate_fit':proposal,
                      '_peaks':peaks})
    anchor=tiles[center]['local_candidate_fit']
    model=(global_affine(anchor['training_affine']['matrix'],
                         origins[center])
           if anchor.get('status')=='HELDOUT_QC_PASS' else None)
    for tile in tiles:
        idx=tile['index']
        if idx==center:
            tile['fixed_anchor_transfer']={'status':'ANCHOR_NOT_EVALUATED_AS_TRANSFER',
                                           'evaluated':0,'matches':0}
        else:
            tile['fixed_anchor_transfer']=transfer_score(
                positions,mags,model,origins[idx],tile['_peaks'],tile_size,
                match_radius=match_radius)
        del tile['_peaks']
    held=[t['fixed_anchor_transfer'] for t in tiles if t['role']=='HELDOUT_REGION']
    passes=sum(v['status']=='TRANSFER_QC_PASS' for v in held)
    failures=[{'tile_xy':t['tile_xy'],'reason':t['fixed_anchor_transfer']['status']}
              for t in tiles if t['role']=='HELDOUT_REGION'
              and t['fixed_anchor_transfer']['status']!='TRANSFER_QC_PASS']
    return {'version':'0.9.0','kind':'INTERNAL_ASTROMETRY_SPATIAL_GENERALIZATION',
            'image_shape':list(im.shape),'tile_size':int(tile_size),
            'grid':grid,'tile_count':len(tiles),'anchor_index':center,
            'anchor_status':anchor.get('status'),
            'anchor_affine_global':model.tolist() if model is not None else None,
            'transfer_tiles_pass':passes,'transfer_tiles_total':len(held),
            'failed_transfer_tiles':failures,'tiles':tiles,
            'science_calibrated':False,
            'independent_external_reference':False,
            'heldout_tiles_physically_nonoverlapping':True,
            'transform_frozen_on_heldout_regions':True,
            'warning':'Candidate FITS-to-catalog consistency, NOT absolute astrometry or independently calibrated photometry'}


def validate_cluster(fits_path,catalog_path,out_dir,*,cluster=None,
                     tile_size=512,grid=3,seed=42):
    """Read an official public FITS/catalog once, evaluate disjoint mosaic tiles."""
    from astropy.io import fits
    from astropy.wcs import WCS
    fits_path=Path(fits_path)
    catalog_path=Path(catalog_path)
    if not fits_path.is_file() or not catalog_path.is_file():
        raise FileNotFoundError('The FITS and corresponding catalog must exist')
    stars=read_acsggct_catalog(catalog_path)
    mask=quality_mask(stars,band='f606w')
    if not mask.any():
        raise ValueError('No quality-selected stars; verify target catalog schema')
    with fits.open(fits_path,memmap=True) as hdus:
        options=[(i,h) for i,h in enumerate(hdus)
                 if h.data is not None and getattr(h.data,'ndim',0)==2]
        options.sort(key=lambda item:(
            item[1].header.get('EXTNAME','').upper()=='SCI',
            item[1].data.size),reverse=True)
        if not options:
            raise ValueError('No science image in FITS')
        i,science=options[0]
        if science.data.shape[0]<grid*tile_size or science.data.shape[1]<grid*tile_size:
            raise ValueError('FITS too small for non-overlapping requested grid')
        world=WCS(science.header,hdus).celestial
        if not world.has_celestial or world.pixel_n_dim!=2:
            raise ValueError('No celestial WCS in selected science image')
        projected=project_catalog(stars,world,mask)
        img=np.asarray(science.data,dtype=float)  # one copy for safety
        valid=np.isfinite(img)
        for h in hdus:
            if h is not science and h.data is not None and h.data.shape==img.shape:
                if h.header.get('EXTNAME','').upper()=='DQ':
                    valid &= (np.asarray(h.data)==0)
        report=assess_mosaic(img,valid,projected,stars['mag_f606w'][mask],
                             tile_size=tile_size,grid=grid,seed=seed)
        report['science_hdu_index']=i
        report['fits_bunit']=science.header.get('BUNIT')
    report['cluster']=cluster or fits_path.stem
    report['fits_sha256']=sha256_file(fits_path)
    report['catalog_sha256']=sha256_file(catalog_path)
    report['catalog_rows']=len(stars['source_id'])
    report['quality_selected_catalog']=int(mask.sum())
    out=Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    (out/'cross_region.json').write_text(
        json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    columns=('tile_id','role','x','y','n_refs','n_peaks',
             'local_status','local_matched','local_rms',
             'transfer_status','transfer_n','transfer_fraction',
             'transfer_median','transfer_over_null')
    with (out/'cross_region.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f)
        writer.writerow(columns)
        for row in report['tiles']:
            fit=row['local_candidate_fit']
            tf=row['fixed_anchor_transfer']
            lf=fit.get('training_affine') or {}
            writer.writerow([row['index'],row['role'],*row['tile_xy'],
                   row['reference_candidates'],row['detected_peaks'],
                   fit.get('status'),lf.get('matched'),lf.get('rms_px'),
                   tf.get('status'),tf.get('matches'),tf.get('match_fraction'),
                   tf.get('median_residual_px'),tf.get('match_over_null_ratio')])
    return report
