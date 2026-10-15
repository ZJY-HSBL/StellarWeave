"""Associate F606W and F814W *through documented shared catalog IDs*.

The ACS survey catalog already contains paired magnitudes. We map that same
reference list independently onto each science image, rather than matching
two independently detected lists and calling the result ground truth.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .catalog import read_acsggct_catalog, quality_mask, project_catalog
from .dataset import read_fits_science, plot_alignment_qc
from .registration_pipeline import register_catalog
from .download_acsggct import checksum


def paired_reference(catalog, *, max_mag_error=.2, max_crowding=.5,
                     max_psf_quality=.5):
    """Require high-quality *both-band* photometry; derive paired magnitudes."""
    q606 = quality_mask(catalog, band='f606w', max_mag_error=max_mag_error,
                        max_crowding=max_crowding, max_psf_quality=max_psf_quality)
    q814 = quality_mask(catalog, band='f814w', max_mag_error=max_mag_error,
                        max_crowding=max_crowding, max_psf_quality=max_psf_quality)
    shared = q606 & q814
    # Avoid misinterpreting extreme values as physical colors.
    color = catalog['mag_f606w'] - catalog['mag_f814w']
    shared &= np.isfinite(color) & (color > -5) & (color < 10)
    return shared


def build_dualband(f606w_path, f814w_path, catalog_path, out_dir,
                   *, alignment='affine', max_mag_error=.2,
                   max_crowding=.5, max_psf_quality=.5):
    if Path(f606w_path).resolve() == Path(f814w_path).resolve():
        raise ValueError('F606W and F814W must use distinct image files')
    table = read_acsggct_catalog(catalog_path)
    mask = paired_reference(table, max_mag_error=max_mag_error,
                            max_crowding=max_crowding,
                            max_psf_quality=max_psf_quality)
    if int(mask.sum()) < 20:
        raise ValueError('Insufficient two-band catalog reference stars')
    catalog = {name: arr[mask] for name, arr in table.items()}
    images = {}
    positions = {}
    quality = {}
    info = {}
    for band, source in (('f606w',f606w_path), ('f814w',f814w_path)):
        img, valid, wcs, metadata = read_fits_science(source)
        projected = project_catalog(catalog, wcs)
        positions[band], quality[band], passed = register_catalog(
            img, projected, catalog['mag_'+band], mode=alignment)
        images[band] = {'shape': list(img.shape), 'valid_fraction': float(valid.mean()),
                        'source': str(source), 'sha256': checksum(source)}
        info[band] = passed
    if alignment == 'affine' and not all(info.values()):
        raise RuntimeError('One or more bands failed registration validation')
    coords606, coords814 = positions['f606w'], positions['f814w']
    def inside(pts, shape):
        h,w = shape
        return ((pts[:,0] >= 0) & (pts[:,0] < w)
                & (pts[:,1] >= 0) & (pts[:,1] < h)
                & np.isfinite(pts).all(axis=1))
    keep = inside(coords606, images['f606w']['shape']) & inside(
        coords814, images['f814w']['shape'])
    if int(keep.sum()) < 20:
        raise ValueError('Too few paired reference stars lie inside both images')
    out = Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    fields = ['catalog_id','ra_deg','dec_deg','f606w_x','f606w_y',
              'f814w_x','f814w_y','mag_f606w','mag_f814w',
              'color_f606w_minus_f814w','error_f606w','error_f814w']
    with (out/'paired_reference.csv').open('w',newline='',encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        for i in np.flatnonzero(keep):
            writer.writerow([int(catalog['source_id'][i]),
                             float(catalog['ra_deg'][i]),float(catalog['dec_deg'][i]),
                             *coords606[i].tolist(),*coords814[i].tolist(),
                             float(catalog['mag_f606w'][i]),
                             float(catalog['mag_f814w'][i]),
                             float(catalog['mag_f606w'][i]-catalog['mag_f814w'][i]),
                             float(catalog['magerr_f606w'][i]),
                             float(catalog['magerr_f814w'][i])])
    report = {'version':'0.3.0','catalog_source':str(catalog_path),
              'catalog_sha256':checksum(catalog_path),
              'registered': bool(all(info.values())),
              'requires_visual_inspection': True,
              'stellar_catalog_not_complete_ground_truth': True,
              'science_calibrated_photometry': False,
              'paired_records': int(keep.sum()), 'images': images,
              'astrometric_quality': quality,
              'quality_selection': {'max_mag_error':max_mag_error,
                                    'max_crowding':max_crowding,
                                    'max_psf_quality':max_psf_quality}}
    (out/'dualband_manifest.json').write_text(
        json.dumps(report,indent=2),encoding='utf-8')
    # No image resampling or flux transformations are performed.
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description='ACSGGCT F606W/F814W catalog association')
    p.add_argument('--f606w', required=True, help='F606W HST FITS file')
    p.add_argument('--f814w', required=True, help='F814W HST FITS file')
    p.add_argument('--catalog',required=True)
    p.add_argument('--out',default='data/processed/paired')
    p.add_argument('--alignment',choices=['affine','wcs'],default='affine')
    a=p.parse_args(argv)
    summary=build_dualband(a.f606w,a.f814w,a.catalog,a.out,alignment=a.alignment)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
