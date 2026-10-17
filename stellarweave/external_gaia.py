"""Epoch-aware external Gaia DR3 validation of FROZEN HST plate solutions.

Gaia sources are never used to fit the affine or select the HST anchor. The
validation is an independent catalog check ONLY if the Gaia file really
comes from the ESA Gaia Archive and the observation epoch is defensible.
An unverified Gaia-shaped CSV does not establish observational accuracy.

Gaia proper motions are pmra = mu_alpha*cos(delta), pmdec (mas/yr).
Astropy SkyCoord handles propagation; no default HST epoch is invented.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .astrometry import apply_affine
from .robust_plate import unique_links

# Explicitly define null offsets; do not rely on private unrelated constants.
NULL_OFFSETS=((39.,0.),(-39.,0.),(0.,39.),(0.,-39.),(27.,31.),(-27.,-31.))
FIELDS=('source_id','ra','dec','pmra','pmdec','ref_epoch',
        'phot_g_mean_mag','ruwe','visibility_periods_used')


def read_gaia_dr3_csv(path,*,max_ruwe=1.4,mag_max=19.5,
                      min_visibility=8):
    """Strict Gaia DR3 csv reader with documented astrometric-quality cuts."""
    f=Path(path)
    if not f.is_file():
        raise FileNotFoundError(str(f))
    rows=[]
    with f.open('r',encoding='utf-8-sig',newline='') as stream:
        reader=csv.DictReader(stream)
        if reader.fieldnames is None or not set(FIELDS).issubset(reader.fieldnames):
            raise ValueError('Gaia CSV needs source_id, ra, dec, pmra, pmdec, '
                             'ref_epoch, G magnitude, ruwe, visibility_periods_used')
        for row in reader:
            try:
                source_id=int(row['source_id'])
                vals=[float(row[name]) for name in FIELDS[1:]]
            except (ValueError,TypeError,OverflowError) as exc:
                # Unmeasured Gaia DR3 fields can be empty and are excluded.
                continue
            ra,dec,pmra,pmdec,epoch,g,ruwe,visible=vals
            if (source_id<=0 or not all(np.isfinite(vals))
                    or not (0<=ra<=360 and -90<=dec<=90)
                    or not (2000<=epoch<=2035)
                    or ruwe<0 or visible<0):
                continue
            if ruwe>max_ruwe or g>mag_max or visible<min_visibility:
                continue
            rows.append((source_id,*vals))
    if not rows:
        raise ValueError('No Gaia DR3 sources survive the required quality cuts')
    if len(set(r[0] for r in rows))!=len(rows):
        raise ValueError('Duplicate Gaia source_id detected')
    return {key:np.asarray([r[i] for r in rows],
            dtype=np.int64 if i==0 else float) for i,key in enumerate(FIELDS)}


def propagate_to_hst_epoch(table,hst_mjd):
    """Propagate all Gaia positions to the explicitly declared HST MJD."""
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    import astropy.units as u
    if not np.isfinite(hst_mjd) or hst_mjd<40000 or hst_mjd>90000:
        raise ValueError('Provide real HST observation MJD; no guessed epoch')
    sc=SkyCoord(ra=table['ra']*u.deg,dec=table['dec']*u.deg,
                pm_ra_cosdec=table['pmra']*u.mas/u.yr,
                pm_dec=table['pmdec']*u.mas/u.yr,
                obstime=Time(table['ref_epoch'],format='jyear'))
    projected=sc.apply_space_motion(new_obstime=Time(hst_mjd,format='mjd'))
    return np.column_stack((projected.ra.deg,projected.dec.deg))


def _checksum(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def evaluate_frozen_plate(fits_path,gaia_csv,plate_report,out_path, *,
                          hst_mjd,eval_tile_index,matching_px=2.,
                          min_sources=5,provenance='GAIA_DR3_USER_SUPPLIED'):
    """Use existing v0.9 frozen plate to match Gaia stars on a held-out tile.

    Does not optimize the transform on Gaia data. A pass is evidence of
    cross-survey agreement, not proof of full science astrometry.
    """
    from astropy.io import fits
    from astropy.wcs import WCS
    from .detection import scipy_peak_detector
    with Path(plate_report).open() as f:
        plate=json.load(f)
    if plate.get('anchor_status')!='HELDOUT_QC_PASS':
        raise ValueError('Frozen plate anchor has not passed its own QC')
    if plate.get('anchor_affine_global') is None:
        raise ValueError('No fixed anchor affine in the plate manifest')
    anchor=int(plate['anchor_index'])
    tiles=plate['tiles']
    selected=[t for t in tiles if t['index']==eval_tile_index]
    if len(selected)!=1 or eval_tile_index==anchor:
        raise ValueError('External Gaia evaluation tile MUST be withheld from anchor fit')
    tile=selected[0]
    if tile.get('fixed_anchor_transfer',{}).get('status')!='TRANSFER_QC_PASS':
        raise ValueError('Gaia tile has no validated internal frozen-transfer frame')
    if matching_px<=0 or min_sources<1:
        raise ValueError('Invalid match thresholds')
    gaia=read_gaia_dr3_csv(gaia_csv)
    sky=propagate_to_hst_epoch(gaia,float(hst_mjd))
    with fits.open(fits_path,memmap=True) as hdus:
        i=int(plate['science_hdu_index'])
        image_hdu=hdus[i]
        wcs=WCS(image_hdu.header,hdus).celestial
        if not wcs.has_celestial or wcs.pixel_n_dim!=2:
            raise ValueError('Missing celestial WCS')
        x,y=wcs.all_world2pix(sky[:,0],sky[:,1],0)
        points=apply_affine(np.column_stack((x,y)),
                            np.asarray(plate['anchor_affine_global'],float))
        origin=np.asarray(tile['tile_xy'],float)
        size=int(plate['tile_size'])
        local=points-origin
        inside=(np.isfinite(local).all(axis=1)
                &(local[:,0]>=5)&(local[:,0]<size-5)
                &(local[:,1]>=5)&(local[:,1]<size-5))
        local=local[inside]
        if len(local)<min_sources:
            return _write({'status':'INSUFFICIENT_GAIA_STARS',
                   'gaia_eligible':int(len(local)),
                   'observation_mjd':float(hst_mjd),
                   'absolute_astrometry_certified':False},out_path)
        x0,y0=map(int,origin)
        patch=np.asarray(image_hdu.data[y0:y0+size,x0:x0+size],
                         dtype=float)
        if patch.shape!=(size,size):
            raise ValueError('Expected exact image crop')
        peaks=scipy_peak_detector(patch,fwhm=2.6,nsigma=5.)
    matched=unique_links(local,peaks,matching_px)
    offsets=[]
    for dx,dy in NULL_OFFSETS:
        shifted=local+np.array([dx,dy])
        ok=((shifted[:,0]>=0)&(shifted[:,0]<size)
             &(shifted[:,1]>=0)&(shifted[:,1]<size))
        offsets.append(int(len(unique_links(shifted[ok],peaks,matching_px))))
    ratio=len(matched)/max(1.,float(np.median(offsets)))
    report={'version':'0.10.0','kind':'FROZEN_PLATE_GAIA_DR3_EXTERNAL_AUDIT',
            'status':'MATCHED' if len(matched)>=min_sources else 'INCONCLUSIVE',
            'gaia_source':'ESA Gaia DR3 (file provenance must be independently verified)',
            'gaia_provenance_label':provenance,
            'gaia_file_sha256':_checksum(gaia_csv),
            'plate_file_sha256':_checksum(plate_report),
            'fits_file_sha256':_checksum(fits_path),
            'observation_mjd':float(hst_mjd),
            'gaia_epoch_propagation':'Astropy SkyCoord.apply_space_motion',
            'fitting_uses_gaia':False,
            'evaluated_tile_index':int(eval_tile_index),
            'anchor_tile_index':int(anchor),'gaia_eligible':int(len(local)),
            'gaia_matched':int(len(matched)),
            'gaia_matched_fraction':float(len(matched)/len(local)),
            'median_position_residual_px':float(np.median(matched[:,2])) if len(matched) else None,
            'p95_position_residual_px':float(np.percentile(matched[:,2],95)) if len(matched) else None,
            'null_shift_counts':offsets,
            'match_over_null_ratio':float(ratio),
            'external_observations_verified_by_runtime':False,
            'absolute_astrometry_certified':False,
            'warning':'Valid Gaia source provenance, mosaic observation epoch, '
                      'proper motion and correlated systematics must be checked before any science claim.'}
    return _write(report,out_path)


def _write(result,path):
    out=Path(path)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    return result
