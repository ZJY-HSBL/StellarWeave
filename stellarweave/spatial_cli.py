"""Run exploratory spatial HST native-PSF and astrometry diagnostics."""
import argparse
import json
from pathlib import Path
from .spatial_validation import write_spatial_report


def main(argv=None):
    parser=argparse.ArgumentParser(description='StellarWeave v0.7.0 spatial HST QC')
    parser.add_argument('--fits',required=True)
    parser.add_argument('--catalog',required=True)
    parser.add_argument('--out',default='outputs/spatial_hst_v07')
    parser.add_argument('--size',type=int,default=512)
    parser.add_argument('--crop-x',type=int)
    parser.add_argument('--crop-y',type=int)
    parser.add_argument('--grid-cells',type=int,choices=(2,3),default=2)
    parser.add_argument('--stamp-size',type=int,default=15)
    parser.add_argument('--fwhm',type=float,default=2.6)
    parser.add_argument('--min-stars',type=int,default=3)
    parser.add_argument('--match-radius',type=float,default=3.)
    args=parser.parse_args(argv)
    result=write_spatial_report(args.fits,args.catalog,args.out,
        size=args.size,x0=args.crop_x,y0=args.crop_y,cells=args.grid_cells,
        stamp_size=args.stamp_size,fwhm=args.fwhm,min_stars=args.min_stars,
        matching_radius=args.match_radius)
    print(json.dumps({'status':result['kind'],
        'catalog_rows':result['catalog_rows'],
        'registration':result['catalog_wcs_registration']['status'],
        'raw_wcs_matches':result['catalog_wcs_registration']['raw_wcs_match_count'],
        'psf_candidate_cells':result['psf_candidate_cells'],
        'output':str(Path(args.out).resolve()),
        'science_grade':False},indent=2))


if __name__=='__main__':
    main()
