"""Coordinate-frame alignment for the reference catalog and one FITS image."""
import numpy as np
from .astrometry import fit_affine, apply_affine
from .detection import scipy_peak_detector


def register_catalog(image, projected_xy, magnitudes, *, mode='affine',
                     dx=None, dy=None, max_bright=500):
    """Return corrected pixel coordinates, report, and verified-QC boolean.

    Conservative only: a pass proves internal centroid correspondence within
    the set thresholds, NOT absolute scientific astrometry.
    """
    pos = np.asarray(projected_xy, dtype=float).reshape(-1, 2)
    mag = np.asarray(magnitudes, dtype=float).reshape(-1)
    if len(mag) != len(pos):
        raise ValueError('Magnitude and position arrays have different lengths')
    height, width = image.shape
    onchip = (np.isfinite(pos).all(axis=1) & np.isfinite(mag)
              & (pos[:, 0] >= -60) & (pos[:, 0] < width + 60)
              & (pos[:, 1] >= -60) & (pos[:, 1] < height + 60))
    if mode in ('affine', 'auto'):
        if dx is not None or dy is not None:
            raise ValueError('Cannot supply manual offsets for automatic registration')
        idx = np.flatnonzero(onchip)
        if len(idx) < 20:
            raise ValueError('Fewer than 20 projected candidate stars in FITS field')
        bright = idx[np.argsort(mag[idx])[:max_bright]]
        detections = scipy_peak_detector(image, nsigma=6)
        report = fit_affine(pos[bright], detections)
        corrected = apply_affine(pos, report['matrix'])
        return corrected, report, True
    if mode == 'wcs':
        if dx is not None or dy is not None:
            raise ValueError('Manual offsets require manual registration mode')
        report = {'status': 'UNVERIFIED_WCS_ONLY', 'matrix': [[1,0,0],[0,1,0]],
                  'warning': 'WCS alone was not validated with detected stars'}
        return pos.copy(), report, False
    if mode == 'manual':
        if dx is None or dy is None:
            raise ValueError('Both dx and dy are required')
        report = {'status': 'UNVERIFIED_USER_SHIFT',
                  'matrix': [[1,0,float(dx)],[0,1,float(dy)]],
                  'warning': 'Operator-specified shift has no residual-based verification'}
        return apply_affine(pos, report['matrix']), report, False
    raise ValueError('mode must be affine, wcs, or manual')
