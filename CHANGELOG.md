# Changelog

## v0.2.0 — 2026-10-15

**Real-data pipeline foundation and conservative registration.**

- Added a documented 21-column ACSGGCT catalog parser, strict RA/Dec schema validation and per-band quality filtering.
- Added explicit NGC 6397, NGC 6752 and NGC 104 public product-pair downloader, size guard and SHA-256.
- Added FITS science-HDU loader with celestial WCS and basic DQ mask support.
- Added global astrometric translation fitting, crossmatch count/fraction and RMS fail-safe checks.
- Added nonoverlapping spatial tiles, quality labels, provenance manifests, and registration QC overlay.
- Added synthetic FITS/WCS integration tests (Astropy-enabled CI) and explicit caveats.

**Not demonstrated:** live archived file/schema compatibility, complete science astrometry, calibrated flux/PSF measurement, trustworthy real-HST detection completeness, or deep learning improvements.

## v0.1.0 — 2026-10-15

Synthetic Gaussian-PSF research prototype, basic detection baselines, metrics, FITS quicklook and automated prereleases.
