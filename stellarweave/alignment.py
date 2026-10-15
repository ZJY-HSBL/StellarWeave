"""Conservative translation fit between a WCS-projected star list and detections.

A fit accepted here is only a preliminary global registration, never evidence
that astrometric distortions or local residuals are calibrated.
"""
import numpy as np
from scipy.spatial import cKDTree

def estimate_translation(catalog_xy, detected_xy, search_radius=35.,
                         tolerance=1.25, min_inliers=12, min_fraction=.55,
                         bin_width=1.0, max_reference=500):
    refs = np.asarray(catalog_xy, dtype=float).reshape(-1, 2)
    found = np.asarray(detected_xy, dtype=float).reshape(-1, 2)
    refs = refs[np.isfinite(refs).all(axis=1)][:max_reference]
    found = found[np.isfinite(found).all(axis=1)]
    if len(refs) < min_inliers or len(found) < min_inliers:
        raise ValueError('Too few catalog positions or detections to register')
    tree = cKDTree(found)
    pairs = tree.query_ball_point(refs, search_radius)
    deltas = [found[np.asarray(js, int)] - ref for ref, js in zip(refs, pairs) if js]
    if not deltas:
        raise ValueError('No registration candidates near projected positions')
    votes = np.concatenate(deltas)
    grid = np.floor(votes / bin_width).astype(int)
    keys, counts = np.unique(grid, axis=0, return_counts=True)
    best_key = keys[np.argmax(counts)]
    approx = (best_key + .5) * bin_width
    neighbors = np.linalg.norm(votes - approx, axis=1) <= 2.0 * bin_width
    shift = np.median(votes[neighbors], axis=0)
    for _ in range(3):
        distances, indices = tree.query(refs + shift, distance_upper_bound=tolerance)
        keep = np.isfinite(distances)
        if not keep.any():
            raise ValueError('Translation did not produce any inlier')
        candidates = sorted(((distances[i], i, indices[i]) for i in np.where(keep)[0]))
        used = set()
        ref_ids = []
        det_ids = []
        for _, i, j in candidates:
            if int(j) not in used:
                used.add(int(j))
                ref_ids.append(i)
                det_ids.append(j)
        offset = found[det_ids] - refs[ref_ids]
        shift = np.median(offset, axis=0)
    residual = np.linalg.norm(found[det_ids] - (refs[ref_ids] + shift), axis=1)
    n = len(ref_ids)
    frac = n / len(refs)
    rms = float(np.sqrt(np.mean(residual ** 2)))
    if n < min_inliers or frac < min_fraction or rms > tolerance:
        raise ValueError(f'Unreliable astrometric registration: matches={n}, '
                         f'fraction={frac:.3f}, RMS={rms:.3f} px. '
                         'Inspect FITS WCS/catalog frame, do not export labels.')
    return {'dx_px': float(shift[0]), 'dy_px': float(shift[1]),
            'inliers': n, 'reference_count': len(refs),
            'matched_fraction': float(frac), 'rms_px': rms,
            'search_radius_px': float(search_radius)}
