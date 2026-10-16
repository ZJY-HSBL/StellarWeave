"""Inspect a real ACSGGCT FITS + catalog pair and produce PSF QA artifacts."""
import argparse
import json
from pathlib import Path
from .real_hst import audit_product


def main(argv=None):
    p=argparse.ArgumentParser(description='StellarWeave v0.6 FITS+catalog preflight')
    p.add_argument('--fits',required=True)
    p.add_argument('--catalog',required=True)
    p.add_argument('--out',default='outputs/real_hst_audit')
    p.add_argument('--size',type=int,default=512)
    p.add_argument('--crop-x',type=int)
    p.add_argument('--crop-y',type=int)
    p.add_argument('--stamp-size',type=int,default=15)
    p.add_argument('--fwhm',type=float,default=2.6)
    p.add_argument('--min-stars',type=int,default=4)
    a=p.parse_args(argv)
    result=audit_product(a.fits,a.catalog,a.out,size=a.size,
                         x0=a.crop_x,y0=a.crop_y,stamp_size=a.stamp_size,
                         fwhm=a.fwhm,min_stars=a.min_stars)
    print(json.dumps({'kind':result['kind'],
                      'image':result['input']['file'],
                      'valid_pixel_fraction':result['input']['valid_fraction'],
                      'catalog_rows':result['catalog']['row_count'],
                      'psf_stars':result['crop_psf']['used_stars'],
                      'output':str(Path(a.out).resolve())},indent=2))


if __name__=='__main__':
    main()
