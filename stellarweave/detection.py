"""Quicklook source detection, not PSF photometry."""
import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter

def robust_background(data, iterations=4):
    values = np.asarray(data, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        raise ValueError("No finite image pixels")
    for _ in range(iterations):
        med = np.median(values)
        sig = 1.4826 * np.median(np.abs(values - med))
        if sig <= 0:
            sig = max(float(np.std(values)), 1e-8)
        clipped = values[np.abs(values - med) < 3 * sig]
        if not len(clipped):
            break
        values = clipped
    return float(np.median(values)), max(float(np.std(values)), 1e-8)

def scipy_peak_detector(data, fwhm=2.6, nsigma=4.0, min_separation=3):
    arr = np.asarray(data, dtype=float)
    background, noise = robust_background(arr)
    signal = np.nan_to_num(arr - background, nan=0.0, posinf=0.0, neginf=0.0)
    sigma = fwhm / 2.354820045
    if sigma <= 0:
        raise ValueError("fwhm must be positive")
    filtered = gaussian_filter(signal, sigma=sigma)
    rad = int(np.ceil(4 * sigma))
    yy, xx = np.mgrid[-rad:rad+1, -rad:rad+1]
    kernel = np.exp(-(xx * xx + yy * yy) / (2 * sigma * sigma))
    kernel /= kernel.sum()
    sigma_filtered = noise * np.sqrt(np.sum(kernel ** 2))
    neighborhood = max(3, 2 * (min_separation // 2) + 1)
    peaks = filtered == maximum_filter(filtered, size=neighborhood)
    mask = peaks & (filtered > nsigma * sigma_filtered)
    margin = max(4, int(np.ceil(3 * sigma)))
    mask[:margin, :] = False
    mask[-margin:, :] = False
    mask[:, :margin] = False
    mask[:, -margin:] = False
    y, x = np.nonzero(mask)
    return np.column_stack((x.astype(float), y.astype(float)))

def dao_detector(data, fwhm=2.6, nsigma=5.0):
    from astropy.stats import sigma_clipped_stats
    from photutils.detection import DAOStarFinder
    arr = np.asarray(data, dtype=float)
    invalid = ~np.isfinite(arr)
    _, median, std = sigma_clipped_stats(arr, sigma=3, mask=invalid)
    finder = DAOStarFinder(threshold=nsigma * float(std), fwhm=fwhm,
                           sharpness_range=(0.1, 1.5), exclude_border=True)
    found = finder(np.where(invalid, median, arr) - float(median))
    if found is None:
        return np.empty((0, 2))
    return np.column_stack((np.asarray(found["xcentroid"], dtype=float),
                            np.asarray(found["ycentroid"], dtype=float)))
