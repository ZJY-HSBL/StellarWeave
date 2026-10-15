"""Synthetic Gaussian-PSF stars; not a calibrated HST instrument simulation."""
import numpy as np

def simulated_field(size=512, n_stars=300, fwhm=2.6, seed=42, background=75.0,
                    read_noise=3.0, min_flux=120.0, max_flux=15000.0):
    rng = np.random.default_rng(seed)
    sigma = fwhm / 2.354820045
    model = np.full((size, size), background, dtype=float)
    margin = max(10, int(np.ceil(5 * sigma)))
    if size <= 2 * margin:
        raise ValueError("Image too small for selected PSF")
    x = rng.uniform(margin, size - margin, n_stars)
    y = rng.uniform(margin, size - margin, n_stars)
    flux = np.exp(rng.uniform(np.log(min_flux), np.log(max_flux), n_stars))
    radius = int(np.ceil(5 * sigma))
    for cx, cy, f in zip(x, y, flux):
        x0 = max(0, int(np.floor(cx)) - radius)
        x1 = min(size, int(np.floor(cx)) + radius + 2)
        y0 = max(0, int(np.floor(cy)) - radius)
        y1 = min(size, int(np.floor(cy)) + radius + 2)
        yy, xx = np.mgrid[y0:y1, x0:x1]
        kernel = np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sigma ** 2))
        kernel /= kernel.sum()
        model[y0:y1, x0:x1] += f * kernel
    measured = rng.poisson(np.maximum(model, 0)).astype(float)
    measured += rng.normal(0, read_noise, model.shape)
    return measured, np.column_stack((x, y, flux))
