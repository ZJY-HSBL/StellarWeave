"""Convert a documented ACSGC FITS mosaic and star catalog into safe tiles.

Only uses external FITS WCS and a catalog-sky WCS plus affine validation; never assumes
RDVIQ xref/yref equal image pixels. No trained model or PSF-fitting is claimed.
"""
import csv
import json
from pathlib import Path
import numpy as np
from .catalog import read_acsggct_catalog, quality_mask, project_catalog
from .alignment import estimate_translation
from .download_acsggct import checksum

LABEL_COLUMNS = ('x_pixel', 'y_pixel', 'mag_vega', 'mag_error', 'catalog_id')

def read_fits_science(path):
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(path, memmap=True) as hdul:
        options = [h for h in hdul if h.data is not None and getattr(h.data, 'ndim', 0) == 2]
        options.sort(key=lambda h: (h.header.get('EXTNAME', '').upper() == 'SCI', h.data.size), reverse=True)
        if not options:
            raise ValueError('No two-dimensional FITS science image')
        hdu = options[0]
        wcs = WCS(hdu.header, hdul).celestial
        if not wcs.has_celestial or wcs.pixel_n_dim != 2:
            raise ValueError('FITS science HDU has no valid two-dimensional celestial WCS')
        image = np.asarray(hdu.data, dtype=np.float32).copy()
        dq = None
        for h in hdul:
            if h.header.get('EXTNAME', '').upper() == 'DQ' and h.data is not None and h.data.shape == image.shape:
                dq = np.asarray(h.data).copy()
                break
        valid = np.isfinite(image)
        if dq is not None:
            valid &= (dq == 0)
        info = {'fits_hdu': hdu.name, 'bunit': hdu.header.get('BUNIT', None),
                'image_shape': list(image.shape), 'dq_mask_applied': dq is not None}
        return image, valid, wcs, info

def spatial_split(x0, tile_size, width, guard=64):
    """Disjoint x-stripes with gap to reduce train/val/test spatial leakage."""
    if tile_size < 1 or guard < 0:
        raise ValueError('Invalid tile size or guard')
    b1 = int(width * .6)
    b2 = int(width * .8)
    x1 = x0 + tile_size
    if x1 <= b1 - guard:
        return 'train'
    if x0 >= b1 + guard and x1 <= b2 - guard:
        return 'val'
    if x0 >= b2 + guard and x1 <= width:
        return 'test'
    return None

def write_tiles(image, valid, catalog_positions, reference_rows, out_dir,
                tile_size=256, guard=64, mode='spatial'):
    """Export raw image tiles, validity masks, magnitude labels, and provenance.

    Split is geographic in pixel coordinates, never random overlapping tiles.
    All exported labels are reference detections, NOT complete ground truth.
    """
    if mode not in ('spatial', 'test-only'):
        raise ValueError('mode must be spatial or test-only')
    if tile_size <= 0 or image.ndim != 2 or valid.shape != image.shape:
        raise ValueError('Bad image or tile size')
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    h, w = image.shape
    pos = np.asarray(catalog_positions, dtype=float).reshape(-1, 2)
    if len(pos) != len(reference_rows['source_id']):
        raise ValueError('Catalog/position mismatch')
    files = []
    for y0 in range(0, h - tile_size + 1, tile_size):
        for x0 in range(0, w - tile_size + 1, tile_size):
            split = 'test' if mode == 'test-only' else spatial_split(x0, tile_size, w, guard=guard)
            if split is None:
                continue
            tile = image[y0:y0+tile_size, x0:x0+tile_size]
            good = valid[y0:y0+tile_size, x0:x0+tile_size]
            if np.mean(good) < .95:
                continue
            inside = ((pos[:, 0] >= x0) & (pos[:, 0] < x0+tile_size)
                      & (pos[:, 1] >= y0) & (pos[:, 1] < y0+tile_size)
                      & np.isfinite(pos).all(axis=1))
            star = np.column_stack((pos[inside, 0]-x0, pos[inside, 1]-y0,
                                    reference_rows['mag'][inside],
                                    reference_rows['mag_error'][inside],
                                    reference_rows['source_id'][inside]))
            name = f'{split}_x{x0:05d}_y{y0:05d}'
            dest = out / split
            dest.mkdir(exist_ok=True)
            np.savez_compressed(dest / f'{name}.npz', image=tile.astype(np.float32),
                                valid=good, labels=star.astype(np.float64))
            files.append({'path': str((Path(split) / (name+'.npz')).as_posix()),
                          'split': split, 'x0': x0, 'y0': y0,
                          'star_count': len(star), 'valid_fraction': float(good.mean())})
    if not files:
        raise ValueError('No eligible tiles were produced. Check image dimensions, mask, guard and tile size.')
    with (out / 'index.csv').open('w', newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, fieldnames=files[0].keys())
        writer.writeheader()
        writer.writerows(files)
    return files

def prepare_dataset(fits_path, catalog_path, out_dir, band='f606w',
                    tile_size=256, guard=64, alignment='affine',
                    dx=None, dy=None, split='spatial',
                    max_mag_error=.2, max_psf_quality=.5, max_crowding=.5):
    image, valid, wcs, header = read_fits_science(fits_path)
    table = read_acsggct_catalog(catalog_path)
    mask = quality_mask(table, band=band, max_mag_error=max_mag_error,
                        max_psf_quality=max_psf_quality, max_crowding=max_crowding)
    if mask.sum() == 0:
        raise ValueError('No stars passed conservative reference-quality selection')
    pos = project_catalog(table, wcs, mask)
    q = {k: v[mask] for k,v in table.items()}
    from .registration_pipeline import register_catalog
    pos, report, verified = register_catalog(
        image, pos, q['mag_'+band], mode=alignment, dx=dx, dy=dy)
    # WCS-only/manual registration is diagnostic; forbid normal train/val
    # exports so unverified labels cannot silently enter a training set.
    if not verified and split != 'test-only':
        raise ValueError('Unverified registration requires --split test-only')
    rows = {'source_id': q['source_id'], 'mag': q['mag_'+band],
            'mag_error': q['magerr_'+band]}
    outputs = write_tiles(image, valid, pos, rows, out_dir, tile_size,
                          guard=guard, mode=split)
    summary = {'version': '0.3.0', 'reference_is_incomplete': True,
               'science_grade_photometry': False,
               'catalog_schema': 'ACSGGCT public ASCII:23-column schema validated on NGC6397',
               'fits_path': str(fits_path), 'catalog_path': str(catalog_path),
               'band': band, 'quality_defaults': {'max_mag_error': max_mag_error,
                'max_psf_quality': max_psf_quality, 'max_crowding': max_crowding},
               'fits_metadata': header, 'alignment': report,
               'registration_passed_qc': verified,
               'registration_requires_visual_inspection': True,
               'input_sha256': {'image': checksum(fits_path), 'catalog': checksum(catalog_path)},
               'split_policy': split, 'guard_px': guard, 'tile_size': tile_size,
               'eligible_catalog_stars': int(mask.sum()),
               'tiles': {k: sum(f['split']==k for f in outputs)
                         for k in ('train','val','test')},
               'star_references_in_tiles': sum(f['star_count'] for f in outputs),
               'label_columns': LABEL_COLUMNS}
    plot_alignment_qc(image, pos, out_dir)
    Path(out_dir, 'manifest.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary

def plot_alignment_qc(image, positions, output_directory):
    """Quicklook for mandatory human inspection of candidate-label alignment."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    h, w = image.shape
    x0, y0 = max(0,(w-512)//2), max(0,(h-512)//2)
    crop = image[y0:min(y0+512,h), x0:min(x0+512,w)]
    finite = crop[np.isfinite(crop)]
    lo, hi = np.percentile(finite,[5,99.8]) if len(finite) else (0,1)
    pts = positions[(positions[:,0]>=x0) & (positions[:,0]<x0+crop.shape[1])
                    & (positions[:,1]>=y0) & (positions[:,1]<y0+crop.shape[0])]
    fig,ax=plt.subplots(figsize=(7,7))
    ax.imshow(crop,cmap='gray',origin='lower',vmin=lo,vmax=max(hi,lo+1e-8))
    if len(pts):
        ax.scatter(pts[:,0]-x0,pts[:,1]-y0,s=22,facecolors='none',
                   edgecolors='cyan',linewidths=.65)
    ax.set_title('Reference star WCS/shift overlay — inspect manually')
    ax.set_xlabel('Crop x pixel')
    ax.set_ylabel('Crop y pixel')
    fig.tight_layout()
    fig.savefig(Path(output_directory)/'alignment_qc.png',dpi=150)
    plt.close(fig)
