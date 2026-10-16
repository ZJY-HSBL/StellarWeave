"""Strict reader for the archived ACSGGCT 23-column RDVIQ catalog.

Verified against the *actual* STScI NGC6397 public header on 2026-10-17:
id x y Vvega err VIvega err Ivega err Vground Iground Nv Ni
wV wI xsig ysig othv othi qfitV qfitI RA Dec

Older synthetic tests erroneously assumed 21 columns and mixed up the
final two sky-coordinate columns. Reject that format rather than guessing.
The quality cuts are conservative reference-selection heuristics and are
not validated astronomical completeness criteria.
"""
from pathlib import Path
import numpy as np

HEADER = ('id','x','y','Vvega','err','VIvega','err','Ivega','err',
          'Vground','Iground','Nv','Ni','wV','wI','xsig','ysig',
          'othv','othi','qfitV','qfitI','RA','Dec')
FIELDS = (
    'source_id','x_ref','y_ref','mag_f606w','magerr_f606w',
    'color_f606w_f814w','color_error','mag_f814w','magerr_f814w',
    'ground_v','ground_i','n_f606w','n_f814w',
    'flag_f606w','flag_f814w','x_sigma','y_sigma',
    'crowding_f606w','crowding_f814w','psfq_f606w','psfq_f814w',
    'ra_deg','dec_deg',
)


def read_acsggct_catalog(path):
    """Read only known 23-column rows and verify celestial coordinate fields."""
    rows=[]
    seen_header=False
    with Path(path).open('r',encoding='utf-8',errors='replace') as f:
        for n,line in enumerate(f,1):
            s=line.strip()
            if not s or s.startswith(('#','!','%')):
                continue
            values=s.split()
            if values[0].lower()=='id':
                if seen_header or rows or tuple(z.lower() for z in values)!=tuple(z.lower() for z in HEADER):
                    raise ValueError(f'Catalog line {n}: unexpected 23-column header')
                seen_header=True
                continue
            if len(values)!=len(FIELDS):
                raise ValueError(f'Catalog line {n}: expected exactly 23 columns; got {len(values)}')
            try:
                row=[float(token.replace('D','E').replace('d','e')) for token in values]
            except ValueError as exc:
                raise ValueError(f'Catalog line {n}: invalid numeric column') from exc
            rows.append(row)
    if not rows:
        raise ValueError('Catalog contains no numeric stars')
    arr=np.asarray(rows,dtype=float)
    result={name:arr[:,i].copy() for i,name in enumerate(FIELDS)}
    good=np.isfinite(result['ra_deg'])&np.isfinite(result['dec_deg'])
    if not good.any():
        raise ValueError('No finite RA/Dec values found in catalog')
    if not np.all((result['ra_deg'][good]>=0)&(result['ra_deg'][good]<=360)
                  &(result['dec_deg'][good]>=-90)&(result['dec_deg'][good]<=90)):
        raise ValueError('Catalog RA/Dec columns fail sky-coordinate range validation')
    ids=result['source_id']
    if not np.isfinite(ids).all() or np.any(ids<=0) or np.any(ids!=np.floor(ids)):
        raise ValueError('Catalog source IDs are not positive integers')
    return result


def quality_mask(catalog,band='f606w',max_mag_error=0.2,
                 max_crowding=0.5,max_psf_quality=0.5,
                 min_exposures=1,require_unsaturated=True):
    """Heuristic quality selector using wV/wI, othV/othI and qfitV/qfitI.

    The wV/wI code values 1 and 2 are retained for continuity with the
    published reference selection. Do not interpret these codes as verified
    saturation bits without consulting the survey documentation.
    """
    if band not in ('f606w','f814w'):
        raise ValueError('band must be f606w or f814w')
    mag=catalog['mag_'+band]
    error=catalog['magerr_'+band]
    crowd=catalog['crowding_'+band]
    fit=catalog['psfq_'+band]
    exposures=catalog['n_'+band]
    flag=catalog['flag_'+band]
    m=(np.isfinite(mag)&(mag>-5)&(mag<40)
       &np.isfinite(error)&(error>=0)&(error<=max_mag_error)
       &np.isfinite(crowd)&(crowd>=0)&(crowd<=max_crowding)
       &np.isfinite(fit)&(fit>=0)&(fit<=max_psf_quality)
       &(exposures>=min_exposures)
       &np.isfinite(catalog['ra_deg'])&np.isfinite(catalog['dec_deg']))
    if require_unsaturated:
        m &= np.isin(flag,(1,2))
    return m


def project_catalog(catalog,wcs,mask=None):
    """WCS world coordinates to 0-based pixel positions (no assumed x/y frame)."""
    if mask is None:
        mask=np.ones(len(catalog['source_id']),bool)
    ra=catalog['ra_deg'][mask]
    dec=catalog['dec_deg'][mask]
    if not len(ra):
        return np.empty((0,2),dtype=float)
    x,y=wcs.all_world2pix(ra,dec,0)
    return np.column_stack((np.asarray(x,float),np.asarray(y,float)))
