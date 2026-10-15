# StellarWeave

**StellarWeave** is an open-source research prototype for crowded-field stellar-source detection and future PSF-aware astronomical photometry.

**星织**：面向密集星场的恒星检测、混叠分离与可靠性约束测光研究项目。

> **Research status — v0.1.0 (prototype).** This release contains reproducible synthetic experiments, a SciPy peak baseline, an optional Photutils DAOStarFinder baseline, a FITS quicklook loader, MAST data-discovery/download helpers, and automated tests. **It does not include a validated PSF-fitting photometry pipeline, a trained deep-learning model, evaluated real-HST catalog matches, or a claim of scientific discovery.**

## Quick start / 快速开始

Requires Python 3.11 or 3.12.

```bash
python -m pip install -r requirements.txt
python -m pytest -q
python -m stellarweave.cli --detector scipy --size 512 --stars 300 --seed 42 --out outputs/demo
python -m stellarweave.cli --detector dao --size 512 --stars 300 --seed 42 --out outputs/dao_demo
```

Outputs: `detections.csv`, `synthetic_truth.csv`, `summary.json`, `detections.png`. Results are toy-Gaussian-PSF smoke tests **not real HST results**.

## Real data / 真实天文数据

The public **HST ACS Survey of Galactic Globular Clusters** offers FITS images and photometric star catalogs. Start from [NGC 6397](https://archive.stsci.edu/pub/hlsp/acsggct/ngc6397/), and consider NGC 6752 and NGC 104 for independent-cluster tests. Public catalogs are reference measurements, not perfectly complete ground truth.

Optional MAST query:

```bash
python -m stellarweave.download_hst --object "NGC 6397" --max-mb 250 --out data/hst
```

For a locally downloaded FITS mosaic:

```bash
python -m stellarweave.cli --input data/hst/your_image.fits --detector dao --size 512 --out outputs/hst_quicklook
```

This command produces **candidate positions only**; without matched truth it does **not** calculate accuracy or claim calibrated photometry. FITS science units (often electrons/s), ERR/DQ extensions, PSF variation, WCS, exposure time, and pixel-area corrections must be checked before scientific use.

## Research roadmap / 科研路线

1. Reproduce DAOStarFinder and iterative PSF photometry; inspect ACS image headers and quality masks.
2. Join matched HST catalog references using validated WCS and photometric calibration.
3. Perform artificial-star injection/recovery, reporting completeness versus magnitude, crowding, and separation in PSF FWHM units.
4. Test neural candidate detection coupled to local PSF joint fitting and false-positive / uncertainty calibration.
5. Validate across held-out sky regions and independent clusters. Compare against established pipelines; **do not use a random overlapping-patch split**.

Research hypotheses include improvements under spatially varying/mismatched PSF and reliability constraints. They are hypotheses, not demonstrated novel results.

## Releases / 版本管理

Stable file layout, automated smoke-test CI, and GitHub Releases are used to track versions. A successful push of a new package version on `main` creates a prerelease with a source ZIP and SHA-256 digest; tagged `vX.Y.Z` pushes also work. Only bump `stellarweave.__version__` when you intend a new release.

## Research integrity

Synthetic metrics refer to injected toy stars. Do not treat them as astronomical completeness measurements. The purpose of v0.1.0 is a reproducible, extensible baseline that will become a physically grounded scientific pipeline.

## Sources

- [HST MAST archive](https://archive.stsci.edu/)
- [ACS globular cluster survey](https://archive.stsci.edu/prepds/acsggct/)
- [Photutils documentation](https://photutils.readthedocs.io/)
- [Astroquery MAST documentation](https://astroquery.readthedocs.io/en/stable/mast/mast_obsquery.html)

License: MIT.