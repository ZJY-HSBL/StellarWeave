"""Anderson et al. (2008) 21-column ACSGC photometric catalog reader.

The source catalog is a reference measurement, NOT exhaustive detection truth.
Column order: Anderson et al. 2008 AJ 135 2055, Table 4.
"""
from pathlib import Path
import numpy as np

FIELDS = (
    'source_id', 'x_ref', 'y_ref', 'mag_f606w', 'magerr_f606w',
    'color_f606w_f814w', 'color_error', 'mag_f814w', 'magerr_f814w',
    'ground_v', 'ground_i', 'n_f606w', 'n_f814w', 'flag_f606w',
    'flag_f814w', 'crowding_f606w', 'crowding_f814w', 'psfq_f606w',
    'psfq_f814w', 'ra_deg', 'dec_deg',
)

def read_acsggct_catalog(path):
    """Load a documented 21-column RDVIQ catalog; reject unrecognized layouts.

    Additional columns, if any, are ignored, as *only* the core 21-column
    schema from the 2008 description is understood. Live .adj.zpt variant
    must be checked against FITS/WCS diagnostics before science use.
    """
    rows = []
    with Path(path).open('r', encoding='utf-8', errors='replace') as handle:
        for lineno, line in enumerate(handle, 1):
            content = line.strip()
            if not content or content.startswith(('#', '!', '%')):
                continue
            pieces = content.replace('D', 'E').replace('d', 'e').split()
            if len(pieces) < len(FIELDS):
                raise ValueError(f'Catalog line {lineno}: expected >=21 columns, got {len(pieces)}')
            try:
                values = [float(p) for p in pieces[:len(FIELDS)]]
            except ValueError as exc:
                raise ValueError(f'Catalog line {lineno}: invalid numeric column') from exc
            rows.append(values)
    if not rows:
        raise ValueError('Catalog is empty')
    arr = np.asarray(rows, dtype=np.float64)
    result = {key: arr[:, col].copy() for col, key in enumerate(FIELDS)}
    valid = np.isfinite(result['ra_deg']) & np.isfinite(result['dec_deg'])
    if np.count_nonzero(valid) == 0 or not np.all(
        (result['ra_deg'][valid] >= 0) & (result['ra_deg'][valid] <= 360)
        & (result['dec_deg'][valid] >= -90) & (result['dec_deg'][valid] <= 90)
    ):
        raise ValueError('Catalog RA/Dec columns fail validation; wrong schema?')
    return result

def quality_mask(catalog, band='f606w', max_mag_error=0.2,
                 max_crowding=0.5, max_psf_quality=0.5,
                 min_exposures=1, require_unsaturated=True):
    """Conservative training-label filter, not an astronomical completeness definition."""
    if band not in ('f606w', 'f814w'):
        raise ValueError('band must be f606w or f814w')
    mag = catalog['mag_' + band]
    error = catalog['magerr_' + band]
    crowd = catalog['crowding_' + band]
    q = catalog['psfq_' + band]
    n = catalog['n_' + band]
    flag = catalog['flag_' + band]
    m = (np.isfinite(mag) & (mag > -5) & (mag < 40)
         & np.isfinite(error) & (error >= 0) & (error <= max_mag_error)
         & np.isfinite(crowd) & (crowd >= 0) & (crowd <= max_crowding)
         & np.isfinite(q) & (q >= 0) & (q <= max_psf_quality)
         & (n >= min_exposures)
         & np.isfinite(catalog['ra_deg']) & np.isfinite(catalog['dec_deg']))
    if require_unsaturated:
        m &= np.isin(flag, (1, 2))
    return m

def project_catalog(catalog, wcs, mask=None):
    """World coordinates -> 0-based pixel centers using FITS WCS."""
    if mask is None:
        mask = np.ones(len(catalog['source_id']), dtype=bool)
    ra = catalog['ra_deg'][mask]
    dec = catalog['dec_deg'][mask]
    if not len(ra):
        return np.empty((0, 2), dtype=float)
    x, y = wcs.all_world2pix(ra, dec, 0)
    return np.column_stack((np.asarray(x, float), np.asarray(y, float)))
