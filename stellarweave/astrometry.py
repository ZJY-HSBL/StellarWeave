"""Robust 2-D affine registration from WCS catalog pixels to star centroids.

Six parameters correct first-order scale, rotation, and shear. This is not
a complete geometric-distortion or absolute astrometric calibration.
"""
import numpy as np
from scipy.spatial import cKDTree
from .alignment import estimate_translation


def apply_affine(points, matrix):
    """Apply a 2x3 matrix to an N-by-2 array of 0-based pixel coordinates."""
    p = np.asarray(points, dtype=float).reshape(-1, 2)
    a = np.asarray(matrix, dtype=float)
    if a.shape != (2, 3):
        raise ValueError('Affine matrix must have shape (2,3)')
    return p @ a[:, :2].T + a[:, 2]


def _unique_nearest(source, detected, predicted, radius):
    """Compute greedy nearest one-to-one matches within a maximum distance."""
    distances, indices = cKDTree(detected).query(predicted, distance_upper_bound=radius)
    edges = sorted((float(distances[i]), int(i), int(indices[i]))
                   for i in np.where(np.isfinite(distances))[0])
    used_sources, used_detections = set(), set()
    pairs = []
    for distance, source_id, detection_id in edges:
        if source_id in used_sources or detection_id in used_detections:
            continue
        used_sources.add(source_id)
        used_detections.add(detection_id)
        pairs.append((source_id, detection_id, distance))
    return pairs


def _fit_least_squares(source, target):
    if len(source) < 3:
        raise ValueError('An affine fit requires at least 3 matched stars')
    x = np.column_stack((source, np.ones(len(source))))
    if np.linalg.matrix_rank(x) < 3:
        raise ValueError('Catalog stars are nearly collinear')
    coefficients, _, _, _ = np.linalg.lstsq(x, target, rcond=None)
    return coefficients.T


def _local_stats(source, residual, cells=3, minimum=4):
    """Residual diagnostics across a regular source-coordinate grid."""
    if len(source) == 0:
        return []
    lo = source.min(axis=0)
    width = np.maximum(np.ptp(source, axis=0), 1e-9)
    tile = np.minimum((cells * (source - lo) / width).astype(int), cells-1)
    output = []
    for y in range(cells):
        for x in range(cells):
            mask = (tile[:, 0] == x) & (tile[:, 1] == y)
            if mask.sum() < minimum:
                continue
            distances = np.linalg.norm(residual[mask], axis=1)
            output.append({'x_cell': x, 'y_cell': y, 'n': int(mask.sum()),
                           'rms_px': float(np.sqrt(np.mean(distances**2))),
                           'median_dx_px': float(np.median(residual[mask, 0])),
                           'median_dy_px': float(np.median(residual[mask, 1]))})
    return output


def fit_affine(source_xy, detected_xy, *, initial_shift=None, search_radius=35.0,
               match_radius=3.0, clip_radius=1.5, min_inliers=20,
               min_fraction=.30, max_rms=.65, max_p95=1.3, max_scale_error=.04):
    """Fit WCS-to-image affine transform; reject underconstrained solutions.

    Begin with a global-shift initialization, use iterated one-to-one matching
    and least squares, and robustly exclude large residuals. A visual
    inspection remains necessary for actual astronomical data.
    """
    source = np.asarray(source_xy, float).reshape(-1, 2)
    detected = np.asarray(detected_xy, float).reshape(-1, 2)
    source = source[np.isfinite(source).all(axis=1)]
    detected = detected[np.isfinite(detected).all(axis=1)]
    if len(source) < min_inliers or len(detected) < min_inliers:
        raise ValueError('Insufficient stars for affine registration')
    if min_inliers < 3 or match_radius <= 0 or clip_radius <= 0:
        raise ValueError('Invalid fitting thresholds')
    if initial_shift is None:
        start = estimate_translation(source, detected, search_radius=search_radius,
                                     tolerance=match_radius, min_inliers=min_inliers,
                                     min_fraction=min_fraction, max_reference=len(source))
        shift = [start['dx_px'], start['dy_px']]
    else:
        shift = np.asarray(initial_shift, dtype=float).reshape(2)
    matrix = np.array([[1., 0., shift[0]], [0., 1., shift[1]]], float)
    previous = None
    for iteration in range(8):
        pairs = _unique_nearest(source, detected, apply_affine(source, matrix),
                                match_radius if iteration == 0 else max(clip_radius * 2, 2.))
        if len(pairs) < min_inliers:
            raise ValueError('Insufficient unique matches for affine fitting')
        idx = np.array([p[0] for p in pairs], int)
        obs = np.array([p[1] for p in pairs], int)
        proposal = _fit_least_squares(source[idx], detected[obs])
        residual = detected[obs] - apply_affine(source[idx], proposal)
        norms = np.linalg.norm(residual, axis=1)
        keep = norms < max(clip_radius, float(np.median(norms)) * 2.5)
        if int(keep.sum()) < min_inliers:
            raise ValueError('Excessive astrometric outliers')
        matrix = _fit_least_squares(source[idx[keep]], detected[obs[keep]])
        if previous is not None and np.max(np.abs(previous - matrix)) < 1e-7:
            break
        previous = matrix.copy()

    pairs = _unique_nearest(source, detected, apply_affine(source, matrix),
                            clip_radius)
    if len(pairs) < min_inliers:
        raise ValueError('Final residual threshold rejected too many matches')
    idx = np.array([p[0] for p in pairs], int)
    obs = np.array([p[1] for p in pairs], int)
    matrix = _fit_least_squares(source[idx], detected[obs])
    residual = detected[obs] - apply_affine(source[idx], matrix)
    distance = np.linalg.norm(residual, axis=1)
    rms = float(np.sqrt(np.mean(distance**2)))
    p95 = float(np.percentile(distance, 95))
    frac = len(pairs) / len(source)
    a = matrix[:, :2]
    singular = np.linalg.svd(a, compute_uv=False)
    determinant = np.linalg.det(a)
    if (frac < min_fraction or rms > max_rms or p95 > max_p95 or determinant <= 0
            or np.max(np.abs(singular-1.)) > max_scale_error):
        raise ValueError('Affine fit failed QC: '
                         f'n={len(pairs)}, fraction={frac:.3f}, '
                         f'RMS={rms:.3f}px, p95={p95:.3f}px, scale={singular.tolist()}')
    bx = np.ptp(source[idx, 0]) / max(float(np.ptp(source[:, 0])), 1e-9)
    by = np.ptp(source[idx, 1]) / max(float(np.ptp(source[:, 1])), 1e-9)
    if bx < .45 or by < .45:
        raise ValueError('Affine fit covers too little of the reference field')
    return {'status': 'AFFINE_QC_PASSED', 'matrix': matrix.tolist(),
            'inliers': int(len(pairs)), 'reference_count': int(len(source)),
            'matched_fraction': float(frac), 'rms_px': rms, 'p95_px': p95,
            'coverage_x': float(bx), 'coverage_y': float(by),
            'singular_values': [float(s) for s in singular],
            'local_residual_grid': _local_stats(source[idx], residual),
            'mean_dx_px': float(residual[:, 0].mean()),
            'mean_dy_px': float(residual[:, 1].mean())}
