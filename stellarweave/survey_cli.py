"""Disjoint spatial validation of a registered ACSGGCT FITS/catalog pair."""
import argparse
import json
from pathlib import Path
from .survey_validation import validate_cluster


def main(argv=None):
    p=argparse.ArgumentParser(
        description='StellarWeave v0.9 non-overlapping region transfer (NOT calibrated)')
    p.add_argument('--fits',required=True)
    p.add_argument('--catalog',required=True)
    p.add_argument('--cluster',required=True)
    p.add_argument('--out',default='outputs/cross_region_v09')
    p.add_argument('--tile-size',type=int,default=512)
    p.add_argument('--grid',choices=[2,3],type=int,default=3)
    p.add_argument('--seed',type=int,default=42)
    args=p.parse_args(argv)
    report=validate_cluster(args.fits,args.catalog,args.out,
        cluster=args.cluster,tile_size=args.tile_size,grid=args.grid,seed=args.seed)
    print(json.dumps({'version':report['version'],'cluster':report['cluster'],
         'catalog_rows':report['catalog_rows'],
         'anchor_status':report['anchor_status'],
         'pass':report['transfer_tiles_pass'],
         'test_tiles':report['transfer_tiles_total'],
         'failures':report['failed_transfer_tiles'],
         'science_calibrated':report['science_calibrated'],
         'output':str(Path(args.out).resolve())},indent=2))


if __name__=='__main__':
    main()
