"""Convert a matching HST ACSGGCT FITS + ASCII catalog into training tiles."""
import argparse
import json
from .dataset import prepare_dataset

def main(argv=None):
    p = argparse.ArgumentParser(description='StellarWeave v0.3.0 HST dataset builder')
    p.add_argument('--fits', required=True, help='ACSGGCT reference image (.fits)')
    p.add_argument('--catalog', required=True, help='Corresponding RDVIQ text catalog')
    p.add_argument('--out', default='data/processed')
    p.add_argument('--band', choices=['f606w','f814w'], default='f606w')
    p.add_argument('--tile-size', type=int, default=256)
    p.add_argument('--guard', type=int, default=64)
    p.add_argument('--alignment', choices=['affine','wcs','manual'], default='affine')
    p.add_argument('--dx', type=float)
    p.add_argument('--dy', type=float)
    p.add_argument('--split', choices=['spatial','test-only'], default='spatial')
    p.add_argument('--max-mag-error', type=float, default=.2)
    p.add_argument('--max-psf-quality', type=float, default=.5)
    p.add_argument('--max-crowding', type=float, default=.5)
    args = p.parse_args(argv)
    kwargs = vars(args)
    kwargs['fits_path'] = kwargs.pop('fits')
    kwargs['catalog_path'] = kwargs.pop('catalog')
    kwargs['out_dir'] = kwargs.pop('out')
    result = prepare_dataset(**kwargs)
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
