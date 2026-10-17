"""Spatial native-pixel HST PSF diagnostics from real ACSGGCT FITS tiles.

Exploratory instrumental shape diagnostics only. Native-pixel stacked PRFs
are NOT supersampled HST ePSFs, a calibrated photometric PSF, or proof of
faint-star detection completeness. No tile is omitted on PSF failure.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from .real_hst import empirical_psf, sha256_file


def psf_shape_metrics(psf):
    """Nonparametric encircled-energy radii and moment ellipticity in pixels."""
    a=np.asarray(psf,dtype=float)
    if a.ndim!=2 or a.shape[0]!=a.shape[1] or a.shape[0]%2!=1:
        raise ValueError('PSF must be a square odd-sized native-pixel stamp')
    if not np.isfinite(a).all() or np.min(a)<0 or a.sum()<=0:
        raise ValueError('PSF must be finite, nonnegative and nonempty')
    a=a/a.sum()
    n=a.shape[0]
    yy,xx=np.mgrid[:n,:n].astype(float)
    c=(n-1)/2.
    rr=np.hypot(xx-c,yy-c)
    order=np.argsort(rr.ravel(),kind='stable')
    cum=np.cumsum(a.ravel()[order])
    radii=rr.ravel()[order]
    # Piecewise interpolation of the sampled radial encircled-energy curve.
    r50=float(np.interp(.5,cum,radii,left=0.,right=radii[-1]))
    r80=float(np.interp(.8,cum,radii,left=0.,right=radii[-1]))
    core=float(a[int(c)-1:int(c)+2,int(c)-1:int(c)+2].sum())
    inner=rr<=min(5.,c)
    w=np.where(inner,a,0.)
    w=w/w.sum()
    mx=float((xx*w).sum())
    my=float((yy*w).sum())
    dx=xx-mx;dy=yy-my
    cov=np.array([[(w*dx*dx).sum(),(w*dx*dy).sum()],
                  [(w*dx*dy).sum(),(w*dy*dy).sum()]])
    eigen=np.linalg.eigvalsh(cov)
    ellipticity=float((eigen[-1]-eigen[0])/eigen.sum()) if eigen.sum()>0 else 0.
    return {'r50_px':r50,'r80_px':r80,
            'half_light_diameter_px':2*r50,
            'core_3x3_fraction':core,
            'moment_ellipticity_inner5px':ellipticity}


def assess_psf_grid(image,valid,origins,*,tile_size,cluster,
                    stamp_size=15,min_stars=8,max_stars=80):
    """Nine disjoint predeclared tiles; no adaptive search for good PSF regions."""
    a=np.asarray(image)
    mask=np.asarray(valid,bool)
    if a.ndim!=2 or mask.shape!=a.shape or len(origins)!=9:
        raise ValueError('Need 2D image/mask and exactly nine predefined origins')
    if stamp_size<7 or stamp_size%2!=1 or min_stars<2 or max_stars<min_stars:
        raise ValueError('Invalid PSF stamp/star selection')
    if tile_size<64:
        raise ValueError('Tile size too small')
    origins=[tuple(map(int,p)) for p in origins]
    if len(set(origins))!=9:
        raise ValueError('Repeated PSF tile origin')
    for i,(x,y) in enumerate(origins):
        if x<0 or y<0 or x+tile_size>a.shape[1] or y+tile_size>a.shape[0]:
            raise ValueError('PSF tile exceeds science image')
        for u,v in origins[i+1:]:
            if abs(x-u)<tile_size and abs(y-v)<tile_size:
                raise ValueError('PSF tiles overlap')
    rows=[];psfs={}
    for index,(x,y) in enumerate(origins):
        patch=np.asarray(a[y:y+tile_size,x:x+tile_size],float)
        good=mask[y:y+tile_size,x:x+tile_size]
        row={'tile_index':index,'tile_xy':[x,y],
             'valid_fraction':float(good.mean()),
             'status':'NOT_EVALUATED','used_stars':0,
             'metrics':None,'distance_to_anchor_total_variation':None}
        try:
            psf,stars,meta=empirical_psf(
                patch,good,stamp_size=stamp_size,min_stars=min_stars,
                max_stars=max_stars)
            row.update({'status':'USABLE_NATIVE_PSF',
                        'used_stars':int(meta['used_stars']),
                        'isolated_candidates':int(meta['isolated']),
                        'metrics':psf_shape_metrics(psf)})
            psfs[index]=psf
        except ValueError as exc:
            if ('isolated stars usable' not in str(exc)
                    and 'No finite image pixels' not in str(exc)):
                raise
            row.update({'status':'INSUFFICIENT_ISOLATED_STARS',
                        'failure_reason':str(exc)})
        rows.append(row)
    if 4 in psfs:
        for row in rows:
            i=row['tile_index']
            if i in psfs:
                row['distance_to_anchor_total_variation']=float(
                    .5*np.abs(psfs[i]-psfs[4]).sum())
    usable=[r for r in rows if r['status']=='USABLE_NATIVE_PSF']
    r50=[r['metrics']['r50_px'] for r in usable]
    tv=[r['distance_to_anchor_total_variation'] for r in usable
        if r['tile_index']!=4 and r['distance_to_anchor_total_variation'] is not None]
    report={'version':'0.12.0','kind':'REAL_HST_NATIVE_PIXEL_PSF_SPATIAL_DIAGNOSTIC',
            'cluster':cluster,'tile_size':tile_size,'stamp_size':stamp_size,
            'min_usable_stars':min_stars,'max_usable_stars':max_stars,
            'anchor_index':4,'tile_count':9,'tiles':rows,
            'usable_tiles':len(usable),'failed_tiles':9-len(usable),
            'anchor_usable':4 in psfs,
            'r50_range_px':[float(min(r50)),float(max(r50))] if r50 else None,
            'median_r50_px':float(np.median(r50)) if r50 else None,
            'median_heldout_total_variation':float(np.median(tv)) if tv else None,
            'same_tile_origins_as_astrometric_benchmark':True,
            'native_pixel_psf_not_oversampled_ePSF':True,
            'not_calibrated_photometry':True,
            'not_independent_astrometry':True,
            'limitations':'Peak-selected isolated bright sources; background and centroid '
                          'approximations, unknown sampling/CTE/drizzle correlation. '
                          'Spatial differences are descriptive, not calibrated PSF error.'}
    return report,psfs


def analyze_fits_psf(fits_path,plate_report,out_dir,*,min_stars=8,
                     stamp_size=15):
    """Read FITS tiles using EXACT existing cross-region benchmark origins."""
    from astropy.io import fits
    fits_path=Path(fits_path)
    plate=json.loads(Path(plate_report).read_text(encoding='utf-8'))
    if plate.get('tile_count')!=9 or plate.get('anchor_index')!=4:
        raise ValueError('Expected nine-tile frozen-plate manifest')
    if plate.get('heldout_tiles_physically_nonoverlapping') is not True:
        raise ValueError('Plate manifest lacks nonoverlap QC')
    sha=sha256_file(fits_path)
    if sha!=plate.get('fits_sha256'):
        raise ValueError('PSF FITS SHA-256 differs from frozen-plate FITS')
    tile_size=int(plate['tile_size'])
    origins=[None]*9
    for t in plate['tiles']:
        idx=int(t['index'])
        if idx not in range(9) or origins[idx] is not None:
            raise ValueError('Bad tile indices in plate manifest')
        origins[idx]=t['tile_xy']
    if any(x is None for x in origins):
        raise ValueError('Missing plate tile')
    with fits.open(fits_path,memmap=True) as hdus:
        i=int(plate['science_hdu_index'])
        sci=hdus[i]
        shape=sci.data.shape
        # Copy only nine 512x512 regions; no full 100MB image allocation.
        # Build a minimal bounding region that preserves original coordinates.
        # The tile origins are 3x3 contiguous but not assumed to be at 0,0.
        xmin=min(x for x,y in origins)
        ymin=min(y for x,y in origins)
        xmax=max(x for x,y in origins)+tile_size
        ymax=max(y for x,y in origins)+tile_size
        if not (0<=xmin<xmax<=shape[1] and 0<=ymin<ymax<=shape[0]):
            raise ValueError('Plate tile outside original FITS image')
        crop=np.asarray(sci.data[ymin:ymax,xmin:xmax],float)
        good=np.isfinite(crop)
        for h in hdus:
            if h is sci or h.data is None or h.data.shape!=shape:
                continue
            if h.header.get('EXTVER',1)!=sci.header.get('EXTVER',1):
                continue
            name=h.header.get('EXTNAME','').upper()
            if name=='DQ':
                good &= (h.data[ymin:ymax,xmin:xmax]==0)
            elif name=='ERR':
                err=np.asarray(h.data[ymin:ymax,xmin:xmax])
                good &= np.isfinite(err)&(err>0)
    local_origins=[(x-xmin,y-ymin) for x,y in origins]
    report,psfs=assess_psf_grid(crop,good,local_origins,
        tile_size=tile_size,cluster=plate.get('cluster'),
        min_stars=min_stars,stamp_size=stamp_size)
    for row,(x,y) in zip(report['tiles'],origins):
        row['tile_xy']=[int(x),int(y)]
    report['input_fits_sha256']=sha
    report['plate_report_sha256']=sha256_file(plate_report)
    out=Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    (out/'psf_spatial.json').write_text(
        json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    with (out/'psf_spatial.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f)
        w.writerow(['tile_index','status','used_stars','r50_px','r80_px',
                    'ellipticity','distance_to_anchor_total_variation'])
        for row in report['tiles']:
            m=row['metrics'] or {}
            w.writerow([row['tile_index'],row['status'],row['used_stars'],
                        m.get('r50_px'),m.get('r80_px'),
                        m.get('moment_ellipticity_inner5px'),
                        row['distance_to_anchor_total_variation']])
    for index,psf in psfs.items():
        np.save(out/f'tile_{index}_native_psf.npy',psf)
    # Genuine image-derived 3x3 PSF mosaic; failed tiles remain visible.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,3,figsize=(8.6,8.4))
    vmax=max((p.max() for p in psfs.values()),default=1.)
    for row,ax in zip(report['tiles'],axes.flat):
        i=row['tile_index']
        if i in psfs:
            ax.imshow(psfs[i],origin='lower',cmap='magma',vmin=0,vmax=vmax)
            ax.set_title(f"Tile {i} · n={row['used_stars']} · "
                         f"r50={row['metrics']['r50_px']:.2f}px",fontsize=9)
        else:
            ax.set_facecolor('#e8e8e8')
            ax.text(.5,.5,'INSUFFICIENT\nISOLATED STARS',
                    transform=ax.transAxes,ha='center',va='center',fontsize=9)
            ax.set_title(f'Tile {i} · rejected',fontsize=9)
        ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle(f"{plate.get('cluster','HST')} · native-pixel empirical PSFs "
                 "(NOT oversampled ePSF)",fontsize=11)
    fig.tight_layout()
    fig.savefig(out/'psf_montage.png',dpi=160)
    plt.close(fig)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description='Real HST spatial native-pixel PSF audit')
    p.add_argument('--fits',required=True)
    p.add_argument('--plate-report',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--min-stars',type=int,default=8)
    p.add_argument('--stamp-size',type=int,default=15)
    args=p.parse_args(argv)
    d=analyze_fits_psf(args.fits,args.plate_report,args.out,
        min_stars=args.min_stars,stamp_size=args.stamp_size)
    print(json.dumps({'cluster':d['cluster'],'usable_tiles':d['usable_tiles'],
                      'failed_tiles':d['failed_tiles'],
                      'median_r50_px':d['median_r50_px'],
                      'median_heldout_total_variation':d['median_heldout_total_variation'],
                      'calibrated':False},indent=2))


if __name__=='__main__':
    main()
