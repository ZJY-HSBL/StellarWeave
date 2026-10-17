"""Conservative HST FITS observation-epoch provenance audit.

This tool inventories explicit FITS timing keywords. It does NOT assert a
unique science-valid mosaic epoch, because a drizzled stack can combine
exposures across dates and instruments. No Gaia validation may use an
unverified single timestamp silently.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np

TIME_KEYS=('EXPSTART','EXPEND','MJD-OBS','MJDOBS','DATE-OBS',
           'DATE-BEG','DATE-END','DATE-AVG','DATEOBS','TSTART','TSTOP')
MJD_KEYS={'EXPSTART','EXPEND','MJD-OBS','MJDOBS'}


def audit_epoch(fits_path,out_path=None):
    from astropy.io import fits
    from astropy.time import Time
    path=Path(fits_path)
    with fits.open(path,memmap=True) as hdus:
        headers=[]
        for i,h in enumerate(hdus):
            raw={}
            converted={}
            for key in TIME_KEYS:
                if key not in h.header:
                    continue
                val=h.header[key]
                raw[key]=str(val)
                try:
                    if key in MJD_KEYS:
                        mjd=float(val)
                        if 40000<=mjd<=90000 and np.isfinite(mjd):
                            converted[key]=mjd
                    elif key.startswith('DATE'):
                        mjd=float(Time(str(val),format='isot').mjd)
                        if 40000<=mjd<=90000 and np.isfinite(mjd):
                            converted[key]=mjd
                except (ValueError,TypeError):
                    try:
                        if key.startswith('DATE'):
                            mjd=float(Time(str(val),format='fits').mjd)
                            if 40000<=mjd<=90000:
                                converted[key]=mjd
                    except (ValueError,TypeError):
                        pass
            headers.append({'hdu':i,'extname':str(h.header.get('EXTNAME','PRIMARY')),
                            'naxis':int(h.header.get('NAXIS',0)),
                            'explicit_timing_keywords':raw,
                            'parsed_mjd':converted,
                            'history_lines':len(h.header.get('HISTORY',[]))
                            if 'HISTORY' in h.header else 0})
    mjds=[v for h in headers for v in h['parsed_mjd'].values()]
    result={'kind':'HST_MOSAIC_EPOCH_PROVENANCE_AUDIT',
            'file':path.name,
            'headers':headers,
            'explicit_parsed_mjd_count':len(mjds),
            'min_header_mjd':min(mjds) if mjds else None,
            'max_header_mjd':max(mjds) if mjds else None,
            'header_mjd_span_days':max(mjds)-min(mjds) if mjds else None,
            'effective_mosaic_epoch_verified':False,
            'gaia_epoch_can_be_inferred_automatically':False,
            'warning':'FITS mosaic header dates are not necessarily a weighted observation epoch. '
                      'Inspect exposure provenance before external Gaia astrometry.'}
    if out_path is not None:
        dest=Path(out_path)
        dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(json.dumps(result,indent=2,allow_nan=False),
                        encoding='utf-8')
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description='Audit FITS timing keywords, never guess HST epoch')
    p.add_argument('--fits',required=True)
    p.add_argument('--out')
    a=p.parse_args(argv)
    result=audit_epoch(a.fits,a.out)
    print(json.dumps({'file':result['file'],
         'explicit_parsed_mjd_count':result['explicit_parsed_mjd_count'],
         'min_header_mjd':result['min_header_mjd'],
         'max_header_mjd':result['max_header_mjd'],
         'header_mjd_span_days':result['header_mjd_span_days'],
         'effective_mosaic_epoch_verified':False},indent=2))


if __name__=='__main__':
    main()
