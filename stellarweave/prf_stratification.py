"""Post-hoc, frozen-holdout HST PRF diagnostics for crowding and brightness.

Uses EXACT v0.12 held-out stars and their unchanged PRF model scores. No
model is refitted or selected using these features. These are exploratory
image-derived proxies, not physical stellar flux/SNR or proof of causation.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .detection import robust_background, scipy_peak_detector
from .real_hst import sha256_file

FEATURES = (
    'peak_over_tile_rms',
    'nearest_detected_peak_px',
    'neighbors_within_25px',
    'annulus_mad_over_tile_rms',
)
CLUSTERS = ('ngc6397', 'ngc6752', 'ngc0104')
OUTER = (0, 1, 2, 3, 5, 6, 7, 8)
ANNULUS_INNER_PX = 5
ANNULUS_OUTER_PX = 10
NEIGHBOR_RADIUS_PX = 25


def measure_tile_proxies(image, valid, holdout):
    """Attach image-derived proxies to already-scored held-out sources.

    A nearest neighbor is another 6-sigma *detected peak*, not necessarily
    another resolved star. The annulus MAD measures local image texture;
    both proxies are contaminated by blends and source wings.
    """
    image = np.asarray(image, dtype=float)
    valid = np.asarray(valid, dtype=bool)
    if image.ndim != 2 or image.shape != valid.shape:
        raise ValueError('Expected same-shaped 2D image and validity mask')
    sky, noise = robust_background(np.where(valid, image, np.nan))
    detections = scipy_peak_detector(
        np.where(valid, image, np.nan), fwhm=2.6, nsigma=6)
    if len(detections) < 2:
        raise ValueError('Insufficient independent detected peaks for neighbors')
    tree = cKDTree(detections)
    yy, xx = np.mgrid[-ANNULUS_OUTER_PX:ANNULUS_OUTER_PX+1,
                      -ANNULUS_OUTER_PX:ANNULUS_OUTER_PX+1]
    annulus = (xx*xx + yy*yy >= ANNULUS_INNER_PX**2) & (
        xx*xx + yy*yy <= ANNULUS_OUTER_PX**2)
    rows = []
    for star in holdout:
        x, y = float(star['tile_x']), float(star['tile_y'])
        if not np.isfinite([x, y]).all():
            raise ValueError('Non-finite frozen holdout coordinate')
        dist, nearest = tree.query([x, y], k=1)
        if dist > 0.25:
            raise ValueError('Frozen holdout source not present in image peaks')
        all_dists = tree.query(detections[int(nearest)], k=2)[0]
        neighbor_distance = float(all_dists[1])
        count = int(len(tree.query_ball_point([x, y], NEIGHBOR_RADIUS_PX)) - 1)
        ix, iy = int(round(x)), int(round(y))
        if not valid[iy, ix]:
            raise ValueError('Frozen holdout peak lies on invalid image pixel')
        peak_contrast = max(0.0, (float(image[iy, ix]) - sky) / noise)
        r = ANNULUS_OUTER_PX
        texture = None
        if (ix-r >= 0 and iy-r >= 0 and
                ix+r < image.shape[1] and iy+r < image.shape[0]):
            patch = image[iy-r:iy+r+1, ix-r:ix+r+1]
            good = valid[iy-r:iy+r+1, ix-r:ix+r+1]
            samples = patch[annulus & good]
            if len(samples) >= 80:
                med = float(np.median(samples))
                texture = float(1.4826*np.median(np.abs(samples-med))/noise)
        row = dict(star)
        row.update({
            'peak_over_tile_rms': float(peak_contrast),
            'nearest_detected_peak_px': neighbor_distance,
            'neighbors_within_25px': count,
            'annulus_mad_over_tile_rms': texture,
            'tile_robust_background': float(sky),
            'tile_robust_rms': float(noise),
            'peak_detector_threshold_sigma': 6.0,
        })
        rows.append(row)
    return rows, {'detected_peaks': int(len(detections)),
                  'tile_robust_background': float(sky),
                  'tile_robust_rms': float(noise)}


def run_frozen_stratification(fits_path, plate_path, prf_path, out_dir, *,
                              cluster=None):
    """Read exact v0.12 frozen test set and independently inspect HST pixels."""
    from astropy.io import fits

    fits_path, plate_path, prf_path = map(
        Path, (fits_path, plate_path, prf_path))
    plate = json.loads(plate_path.read_text(encoding='utf-8'))
    prf = json.loads(prf_path.read_text(encoding='utf-8'))
    cluster = cluster or prf.get('cluster')
    if cluster not in CLUSTERS or plate.get('cluster') != cluster or (
            prf.get('cluster') != cluster):
        raise ValueError('Cluster identity mismatch')
    if (prf.get('kind') != 'REAL_HST_SPATIAL_NATIVE_PIXEL_PRF_HOLDOUT'
            or prf.get('science_calibrated') is not False
            or prf.get('no_heldout_star_used_to_fit_a_scored_template') is not True
            or plate.get('tile_size') != 512 or plate.get('grid') != 3
            or plate.get('tile_count') != 9
            or plate.get('heldout_tiles_physically_nonoverlapping') is not True
            or plate.get('anchor_status') != 'HELDOUT_QC_PASS'
            or prf.get('seed') != 42):
        raise ValueError('Unverified frozen v0.12 experiment provenance')
    image_sha = sha256_file(fits_path)
    plate_sha = sha256_file(plate_path)
    if image_sha != plate.get('fits_sha256') or (
            image_sha != prf.get('fits_sha256') or
            plate_sha != prf.get('plate_sha256')):
        raise ValueError('Input FITS or frozen plate checksum mismatch')
    tiles = sorted(plate['tiles'], key=lambda x: x['index'])
    if [x['index'] for x in tiles] != list(range(9)):
        raise ValueError('Incomplete 3x3 tile grid')
    stars = prf['heldout_star_scores']
    if len(stars) != prf.get('total_heldout_stars'):
        raise ValueError('Incomplete original held-out star scores')
    seen = set()
    for star in stars:
        key = (int(star['tile_index']), int(star['holdout_index']))
        if key in seen:
            raise ValueError('Duplicate held-out star')
        seen.add(key)
    selected = {i: [x for x in stars if x['tile_index'] == i] for i in OUTER}
    if sum(map(len, selected.values())) != sum(
            t['heldout_stars'] for t in prf['tiles'] if t['index'] in OUTER):
        raise ValueError('Frozen held-out star counts disagree')
    rows = []
    tile_diagnostics = []
    with fits.open(fits_path, memmap=True) as hdus:
        sci = hdus[int(plate['science_hdu_index'])]
        if sci.data is None or sci.data.ndim != 2:
            raise ValueError('Science FITS plane unavailable')
        ver = sci.header.get('EXTVER', 1)
        err = dq = None
        for h in hdus:
            if h is sci or h.data is None or h.data.shape != sci.data.shape:
                continue
            if h.header.get('EXTVER', 1) != ver:
                continue
            if h.header.get('EXTNAME', '').upper() == 'ERR' and err is None:
                err = h.data
            if h.header.get('EXTNAME', '').upper() == 'DQ' and dq is None:
                dq = h.data
        for t in tiles:
            idx = t['index']
            if idx not in OUTER:
                continue
            x, y = map(int, t['tile_xy'])
            if (x < 0 or y < 0 or x+512 > sci.data.shape[1] or
                    y+512 > sci.data.shape[0]):
                raise ValueError('Tile outside FITS image')
            ys, xs = slice(y, y+512), slice(x, x+512)
            image = np.asarray(sci.data[ys, xs], dtype=float)
            valid = np.isfinite(image)
            if dq is not None:
                valid &= (np.asarray(dq[ys, xs]) == 0)
            if err is not None:
                e = np.asarray(err[ys, xs], dtype=float)
                valid &= np.isfinite(e) & (e > 0)
            subset = selected[idx]
            if not subset:
                tile_diagnostics.append({
                    'tile_index': idx, 'status': 'NO_HELDOUT_STARS',
                    'heldout_stars': 0})
                continue
            measured, info = measure_tile_proxies(image, valid, subset)
            for star in measured:
                star['tile_index'] = idx
                star['global_x'] = float(x+star['tile_x'])
                star['global_y'] = float(y+star['tile_y'])
            rows.extend(measured)
            tile_diagnostics.append({
                'tile_index': idx, 'status': 'EVALUATED',
                'heldout_stars': len(measured), 'tile_xy': [x, y], **info})
    if len(rows) != sum(map(len, selected.values())):
        raise ValueError('Not all original held-out stars were evaluated')
    report = {
        'version': '0.13.0',
        'kind': 'FROZEN_HST_PRF_CROWDING_BRIGHTNESS_EXPLORATORY',
        'cluster': cluster,
        'fits_sha256': image_sha,
        'plate_sha256': plate_sha,
        'original_prf_report_sha256': sha256_file(prf_path),
        'original_prf_kind': prf['kind'],
        'frozen_prf_scores_unchanged': True,
        'models_retrained': False,
        'stars_removed_for_stratification': False,
        'training_data_used_for_threshold_selection': False,
        'analysis_is_posthoc_exploratory': True,
        'science_calibrated': False,
        'physical_crowding_measured': False,
        'peak_snr_calibrated': False,
        'all_original_outer_heldout_stars': len(rows),
        'tile_count': len(OUTER),
        'tiles': tile_diagnostics,
        'features': list(FEATURES),
        'feature_definitions': {
            'peak_over_tile_rms': 'Peak minus robust tile sky, divided by sigma-clipped tile pixel RMS; not physical source SNR',
            'nearest_detected_peak_px': 'Distance to next detected 6-sigma peak; unresolved neighbors invisible',
            'neighbors_within_25px': 'Count of detected 6-sigma peaks within 25 native pixels, excluding self',
            'annulus_mad_over_tile_rms': 'MAD-based annulus (5 to 10 px) texture / robust tile pixel RMS; may include source wings',
        },
        'heldout_star_scores_with_proxies': rows,
        'limitations': [
            'Post-hoc descriptive stratification cannot identify causes',
            '6-sigma peak catalog misses unresolved or faint neighbors',
            'Bright isolated candidates only; not representative of faint crowded stars',
            'All source L1 values are from unchanged native-pixel v0.12 templates',
            'Image peak-over-tile-RMS is not a calibrated flux SNR',
            'No independent Gaia absolute astrometry or calibrated photometry',
        ],
    }
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out/'prf_stratification.json').write_text(
        json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    if rows:
        with (out/'prf_stratification.csv').open(
                'w', newline='', encoding='utf-8') as f:
            fields = list(rows[0])
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
    return report


def main(argv=None):
    p = argparse.ArgumentParser(
        description='Frozen v0.12 HST PRF brightness/neighbor diagnostics')
    p.add_argument('--fits', required=True)
    p.add_argument('--plate-report', required=True)
    p.add_argument('--prf-report', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--cluster')
    a = p.parse_args(argv)
    d = run_frozen_stratification(
        a.fits, a.plate_report, a.prf_report, a.out, cluster=a.cluster)
    print(json.dumps({
        'cluster': d['cluster'],
        'heldout_stars_unchanged': d['all_original_outer_heldout_stars'],
        'tiles': d['tile_count'],
        'models_retrained': False,
        'posthoc_exploratory': True,
    }, indent=2))


if __name__ == '__main__':
    main()
