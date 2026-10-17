"""Fetch a bounded Gaia DR3 source table via the official ESA Gaia TAP archive.

Optional network dependency: astroquery. Queries and metadata are saved to
establish provenance. This is public-source retrieval, not science validation.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np


COLUMNS=('source_id','ra','dec','pmra','pmdec','ref_epoch',
         'phot_g_mean_mag','ruwe','visibility_periods_used')


def query_from_fits(fits_path,*,max_rows=3000,mag_max=19.5):
    """Derive a sky cone from WCS image corners without assuming cluster coords."""
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    if not isinstance(max_rows,int) or not 100<=max_rows<=10000:
        raise ValueError('max_rows must be 100..10000')
    if not 8<=mag_max<=20.5:
        raise ValueError('mag_max outside supported selection range')
    with fits.open(fits_path,memmap=True) as h:
        science=[z for z in h if z.data is not None and getattr(z.data,'ndim',0)==2]
        science.sort(key=lambda z: (z.header.get('EXTNAME','').upper()=='SCI',
                                   z.data.size),reverse=True)
        if not science:
            raise ValueError('FITS has no 2-D image')
        im=science[0]
        hgt,wid=im.data.shape
        w=WCS(im.header,h).celestial
        if not w.has_celestial or w.pixel_n_dim!=2:
            raise ValueError('Selected FITS lacks 2D sky WCS')
        cx,cy=w.all_pix2world((wid-1)/2.,(hgt-1)/2.,0)
        locs=np.array([[0,0],[wid-1,0],[0,hgt-1],[wid-1,hgt-1],
                       [(wid-1)/2,0],[(wid-1)/2,hgt-1],
                       [0,(hgt-1)/2],[wid-1,(hgt-1)/2]],float)
        ra,dec=w.all_pix2world(locs[:,0],locs[:,1],0)
        center=SkyCoord(ra=cx*u.deg,dec=cy*u.deg)
        edge=SkyCoord(ra=ra*u.deg,dec=dec*u.deg)
        radius=float(center.separation(edge).deg.max()+.003)
    if not np.isfinite([cx,cy,radius]).all() or radius<=0 or radius>.20:
        raise ValueError('Sky-cone radius outside safe limits: investigate FITS WCS')
    adql=(f"SELECT TOP {max_rows} "+
          ", ".join(COLUMNS)+
          " FROM gaiadr3.gaia_source"+
          f" WHERE CONTAINS(POINT('ICRS',ra,dec),"+
          f" CIRCLE('ICRS',{cx:.9f},{cy:.9f},{radius:.9f}))=1"+
          f" AND phot_g_mean_mag < {mag_max:.2f}"+
          " AND ruwe < 1.4 AND visibility_periods_used >= 8"+
          " AND pmra IS NOT NULL AND pmdec IS NOT NULL"+
          " ORDER BY phot_g_mean_mag ASC")
    return adql,{'center_ra':float(cx),'center_dec':float(cy),
                 'radius_deg':radius,'max_rows':max_rows,
                 'gmag_max':mag_max,'source':'ESA gaiadr3.gaia_source'}


def fetch_gaia_csv(fits_path,out_csv,*,max_rows=3000,mag_max=19.5):
    from astroquery.gaia import Gaia
    from .real_hst import sha256_file
    query,region=query_from_fits(fits_path,max_rows=max_rows,mag_max=mag_max)
    result=Gaia.launch_job_async(query,output_format='csv').get_results()
    if len(result)>=max_rows:
        raise RuntimeError('Gaia TOP limit reached: results incomplete. '
                           'Use a smaller region or brighter threshold.')
    if not len(result):
        raise ValueError('No Gaia quality-selected DR3 sources returned')
    out=Path(out_csv)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',encoding='utf-8',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(COLUMNS)
        for row in result:
            writer.writerow([str(row[k]) for k in COLUMNS])
    provenance={'gaia_release':'DR3','query':query,'region':region,
                'rows':len(result),'input_fits_sha256':sha256_file(fits_path),
                'gaia_csv_sha256':sha256_file(out),
                'user_query_must_be_rechecked_for_science':True}
    (out.parent/(out.stem+'_provenance.json')).write_text(
        json.dumps(provenance,indent=2),encoding='utf-8')
    return provenance


def main(argv=None):
    p=argparse.ArgumentParser(description='Query ESA Gaia DR3 source subset')
    p.add_argument('--fits',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--max-rows',type=int,default=3000)
    p.add_argument('--gmag-max',type=float,default=19.5)
    p.add_argument('--dry-run',action='store_true')
    a=p.parse_args(argv)
    if a.dry_run:
        q,info=query_from_fits(a.fits,max_rows=a.max_rows,mag_max=a.gmag_max)
        print(json.dumps({'query':q,'region':info},indent=2))
        return
    print(json.dumps(fetch_gaia_csv(
        a.fits,a.out,max_rows=a.max_rows,mag_max=a.gmag_max),indent=2))


if __name__=='__main__':
    main()
