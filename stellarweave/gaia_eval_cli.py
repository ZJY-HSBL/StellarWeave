"""Evaluate a Gaia DR3 CSV on a HST tile NOT used for plate fitting."""
import argparse
import json
from pathlib import Path
from .external_gaia import evaluate_frozen_plate


def main(argv=None):
    p=argparse.ArgumentParser(description='Epoch-aware frozen-plate Gaia DR3 validation')
    p.add_argument('--fits',required=True)
    p.add_argument('--gaia-csv',required=True)
    p.add_argument('--plate-report',required=True,
                   help='v0.9 cross_region.json (frozen center-fit)')
    p.add_argument('--hst-mjd',required=True,type=float,
                   help='Verified HST observation MJD; NEVER assume Gaia epoch')
    p.add_argument('--eval-tile',required=True,type=int,
                   help='Must differ from center anchor index')
    p.add_argument('--out',default='outputs/external_gaia_v010.json')
    p.add_argument('--radius-px',default=2.,type=float)
    p.add_argument('--min-gaia-matches',default=5,type=int)
    a=p.parse_args(argv)
    result=evaluate_frozen_plate(
        a.fits,a.gaia_csv,a.plate_report,a.out,hst_mjd=a.hst_mjd,
        eval_tile_index=a.eval_tile,matching_px=a.radius_px,
        min_sources=a.min_gaia_matches)
    print(json.dumps({'status':result['status'],
        'eligible':result.get('gaia_eligible'),
        'matched':result.get('gaia_matched'),
        'median_px':result.get('median_position_residual_px'),
        'calibrated':False,'output':str(Path(a.out).resolve())},indent=2))


if __name__=='__main__':
    main()
