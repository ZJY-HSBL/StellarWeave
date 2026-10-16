# StellarWeave — 星织

**Astronomical Image Processing · Crowded-field Stellar Detection · Reliable Photometry**

[![CI](https://github.com/ZJY-HSBL/StellarWeave/actions/workflows/ci.yml/badge.svg)](https://github.com/ZJY-HSBL/StellarWeave/actions/workflows/ci.yml) · [GitHub Releases](https://github.com/ZJY-HSBL/StellarWeave/releases)

**v0.7.0 research prototype / 科研原型。** This repository does NOT yet contain a trained neural network, science-grade PSF photometry, or validated real-HST quantitative results. The catalog/image pipeline is tested against synthetic FITS/WCS fixtures, and any real ACS image/catalog pairing requires manual QC.

## Purpose / 研究目标

StellarWeave is designed to investigate faint stellar source detection, deblending, and measurement reliability in crowded telescope imagery. We are developing physically informed baselines, not merely applying YOLO to astronomy.

## v0.7.0 — Spatial PSF and Catalog-WCS Diagnostics / 空间 PSF 与星表配准诊断

**This release performs exploratory spatial diagnostics on genuine HST data. It does NOT deliver a validated astrometric solution or science-grade ePSF.**

First download the original STScI ACSGGCT NGC 6397 F606W image and joint two-filter catalog using `stellarweave.download_acsggct`, as described below. Then run:

```bash
python -m stellarweave.spatial_cli \
  --fits data/raw/ngc6397/hlsp_acsggct_hst_acs-wfc_ngc6397_f606w_v2_img.fits \
  --catalog data/raw/ngc6397/hlsp_acsggct_hst_acs-wfc_ngc6397_r.rdviq.cal.adj.zpt \
  --size 512 --grid-cells 2 --stamp-size 15 --min-stars 3 \
  --match-radius 3 --out outputs/ngc6397_spatial_v07
```

Outputs: `spatial_diagnostics.json` records reference counts, one-to-one raw WCS matches and quality-checked *candidate* affine fit; `astrometry_overlay.png` overlays catalog positions and detected peaks; `spatial_psf_grid.png` shows the local native-pixel PSF estimates, and `psf_yN_xM.npy` provides the individual diagnostic stamps.

### First real NGC 6397 findings

From the successful public-archive GitHub Actions run on October 17, 2026:

| Quantity | Diagnostic result |
|---|---:|
| Public photometric catalog | 16,100 sources |
| Reference stars projected into selected 512×512 image crop | 533 |
| Brighter reference stars selected | 350 |
| Pixel-space brightness peaks detected | 679 |
| Raw WCS-position pairs at <=3 pixels | 36 |
| Fraction of selected references matched | 10.3% |
| Candidate PSF cells with local stars | 4 / 4 |
| Affine registration proposal | **REJECTED BY QC** |

The raw matched positions have a **median separation of about 2.47 px**; with a crowded stellar field, many pairs may be accidental. The affine estimate failed the final inlier/residual criterion. This is a scientifically relevant **negative diagnostic result** and must not be reported as successful stellar-catalog registration.

**Consequences:** no trusted real-data training labels or completeness measurements are emitted by this diagnostic. The next step is to check the particular FITS WCS/reference-frame mapping, broaden and independently validate plate-solving initialization, and test matched stars against independent withheld measurements. The local native-pixel PSFs must also be checked for undersampling, crowding biases, spatial dependence and instrument calibration; they are not a full spatially adaptive oversampled HST ePSF.

Read the [real-HST workflow record](https://github.com/ZJY-HSBL/StellarWeave/actions/workflows/real-hst.yml) and download the `ngc6397-real-hst-audit` diagnostic artifact from the successful run for plots and machine-readable data.
## v0.6.0 — Genuine HST FITS & Catalog Audit / 真实天文观测数据审计

**Real official NGC 6397 was successfully downloaded and inspected in GitHub Actions.** In the tested release candidate, the official ASCII catalog contained **16,100 records**, and the selected 512×512 FITS crop had **100% valid pixels**. A robust native-pixel candidate PSF was estimated from **79 isolated bright source candidates**. These are raw data-quality figures, **not detection completeness or measurement accuracy**.

Official directory: https://archive.stsci.edu/pub/hlsp/acsggct/ngc6397/

Download the actual F606W image and paired catalog:

```bash
python -m stellarweave.download_acsggct --cluster ngc6397 --band f606w --out data/raw
```

Audit the files and create a native-pixel candidate PRF:

```bash
python -m stellarweave.real_hst_cli \
  --fits data/raw/ngc6397/hlsp_acsggct_hst_acs-wfc_ngc6397_f606w_v2_img.fits \
  --catalog data/raw/ngc6397/hlsp_acsggct_hst_acs-wfc_ngc6397_r.rdviq.cal.adj.zpt \
  --size 512 --stamp-size 15 --fwhm 2.6 \
  --out outputs/ngc6397_real_audit
```

Outputs: `audit.json` for FITS unit/exposure/filter/ERR/DQ/WCS presence, catalog row counts, and SHA-256; `selected_psf_stars.csv` for isolated candidate positions; `empirical_prf_native.npy` (unit-normalized local PSF stamp), and `psf_preview.png`.

The audit records FITS scientific pixels, not scientifically verified WCS star matches. The PSF is only a single-tile, native-pixel empirical estimate: **not an oversampled HST instrumental ePSF**, not photometrically calibrated, and not guaranteed transferable to another location. Real HST detection completeness, per-star flux accuracy, noise covariance and physical uncertainty propagation are **not measured**.

[View the Real HST Archive Audit workflow](https://github.com/ZJY-HSBL/StellarWeave/actions/workflows/real-hst.yml). GitHub artifacts contain only small QA products, not redistributed original FITS/catalog files.
## v0.5.0 — PSF photometry benchmark / 精密测光基线

This research baseline compares **DAOStarFinder centroids**, **PSFPhotometry**, and **IterativePSFPhotometry**, as implemented in Photutils 3.x, using the **same injected images, source positions, and input flux levels** in every trial. An optional SciPy centroid control is also available.

Install dependencies, then run a small reproducible synthetic benchmark:

```bash
python -m pip install -r requirements.txt
python -m stellarweave.benchmark_cli \
  --size 160 --trials 2 --per-flux 4 \
  --flux-levels 700 2000 5000 \
  --methods dao psf iterative --seed 42 \
  --out outputs/photometry_v05

```

Outputs: `source_results.csv` contains per-star detection/recovery, localization error, and fitted flux (where available); `benchmark.json` stores matched-sample recovery and Wilson intervals for all methods; `method_comparison.png` plots completeness-vs-instrumental-flux on identical injections.

**Fair comparison protocol:** we create one set of non-overlapping artificial stars per trial, inject each set into the same untouched image, and run each method on precisely those pixels. Baseline source detections are removed before scoring newly detected sources. For recovered isolated injected stars, fitted flux is compared with known instrumental injected flux. Near an original star, flux attribution is ambiguous and is **not used** in the flux-error summary.

**Science limitations:** all measurements here are *toy* or *exploratory*. A single analytic or image PRF can suffer PSF mismatch; reduced chi-square and reported formal flux errors cannot be interpreted as calibrated uncertainty without a physically justified per-pixel error array. Matched-source recovery does not by itself establish deblending completeness. Real-FITS flux and PSF are instrument/filter-specific, and this release does not model injected photon noise, drizzle correlations, detector gain/read noise, or spatial PSF variability. No quantitative HST results are claimed.

Example *only* for a real image (do not run until PSF and FITS units are verified):

```bash
python -m stellarweave.benchmark_cli \
  --fits data/raw/ngc6397/YOUR_VALID_F606W_IMAGE.fits \
  --psf-stamp data/psf/YOUR_VERIFIED_SAME_BAND_PRF.npy \
  --image-unit ELECTRONS/S --size 160 \
  --flux-levels 0.05 0.2 0.8 \
  --methods dao psf iterative --out outputs/hst_psf_benchmark

```

See [Photutils PSF photometry docs](https://photutils.readthedocs.io/en/3.1.0/user_guide/psf.html), [PSFPhotometry](https://photutils.readthedocs.io/en/3.1.0/api/photutils.psf.PSFPhotometry.html), and [IterativePSFPhotometry](https://photutils.readthedocs.io/en/3.1.0/api/photutils.psf.IterativePSFPhotometry.html).

## v0.4.0 — Artificial-star injection and recovery / 人工星注入恢复

**No manual star annotation is needed for artificial-star truth**: injection position and **instrumental** flux are saved exactly. This first implementation is a *methodological* completeness experiment, **not calibrated HST faint-star completeness**.

Quick synthetic experiment (no external telescope files required):

```bash
python -m stellarweave.injection_cli --size 320 --detector scipy \
  --flux-levels 350 800 1800 4000 --trials 3 --per-flux 8 \
  --seed 42 --out outputs/injection_v04
```

Run on a real FITS science image only after obtaining a **pixel-integrated, properly centered and matched-filter/instrument PSF** as a two-dimensional odd-sized NumPy \`.npy\` file. Its sum is normalized to one, and injected total flux must be in the original FITS image's \`BUNIT\`; the CLI rejects a mismatched supplied unit:

```bash
python -m stellarweave.injection_cli \
  --fits data/raw/ngc6397/example_science.fits \
  --psf-stamp data/psf/ngc6397_f606w_prf.npy \
  --image-unit ELECTRONS/S --size 320 --crop-x 512 --crop-y 512 \
  --fwhm 2.6 --detector scipy --trials 5 --per-flux 10 \
  --flux-levels 0.01 0.03 0.10 0.30 \
  --out outputs/ngc6397_injections
```

The above FITS filename, PSF filename and flux levels are **illustrative**; select real files and a sensible flux range based on image calibration. A PSF extracted from stars, without accurate background subtraction and verification that it represents a pixel-integrated PRF, is **not automatically valid**.

Outputs: \`injections.csv\` (per injected source; original x/y, flux, distance to original detections, crowding class, recovery flag and position error), \`manifest.json\` (settings, baseline counts, trial totals and Wilson intervals), and \`completeness.png\` (recovery curves). Each trial injects into the **same untouched base image** to avoid cumulative contamination. Injected-source spacing is controlled to avoid artificial crowding from other injected sources.

**Interpretation limitations**:

- The baseline must be run on the **unmodified** image; recovery is credited only for *new* detections, excluding previously existing centroids. This is intentionally conservative close to existing sources and is **not a PSF-deblending benchmark**.
- "Unmatched new candidates" are **not equivalent to false positives**, because the original sky star catalog is not exhaustive. The distance to the nearest original detection is only a proxy for actual crowding.
- The supplied PSF is treated as fixed across the selected image tile. The expectation-only addition doesn't propagate new photon shot noise or account for drizzled noise correlations, pixel-area correction, exposure weighting, or PSF variation.
- Recovery versus *instrumental flux* must not be reported as a calibrated completeness-versus-magnitude curve until image units and photometric zero points are validated.
- Wilson 95% intervals assume independent Bernoulli observations; repeated injections into the same background may be correlated. More robust later work will use image/tile and injection-batch bootstrapping.
- Baseline \`scipy\` and \`dao\` source detectors are available. Strong science claims require validation against DOLPHOT/Photutils iterative PSF photometry and real ACS image products.

Reference: [Photutils PSF photometry](https://photutils.readthedocs.io/en/stable/user_guide/psf.html) and [pixel-integrated PSF simulation](https://photutils.readthedocs.io/en/stable/reference/psf_api.html).


## Source data / 公开数据（不需人工逐星标注）

Use the [ACS Survey of Galactic Globular Clusters](https://archive.stsci.edu/prepds/acsggct/) from NASA/STScI's MAST archive. Its stellar photometric catalogs are **reference measurements**, not exhaustive ground truth. [Anderson et al. (2008)](https://arxiv.org/abs/0804.2025) Table 4 documents the actual 23-column RDVIQ catalog: ID, reference-frame x/y, F606W and F814W magnitudes/errors, number of exposures, unsaturated/saturated flags, neighbor contamination, PSF quality, and RA/Dec.

Default NGC 6397, then independent clusters NGC 6752 and NGC 104. Note that the actual `.rdviq.cal.adj.zpt` file header has been checked and verified as 23-column and alignment against the specific image.

## Install / 安装

Python 3.11–3.12:

```bash
python -m pip install -r requirements.txt
python -m pytest -q
```

Minimal synthetic test without astronomy dependencies:

```bash
python -m pip install -r requirements-smoke.txt
python -m stellarweave.cli --detector scipy --size 512 --stars 300 --seed 42 --out outputs/demo
```

Synthetic precision, recall and localization error must **never be reported as HST scientific performance**.

## Generate HST training tiles / 准备真实 FITS 与星表

**1. Download the public image/catalog pair (internet required):**

```bash
python -m stellarweave.download_acsggct --cluster ngc6397 --band f606w --out data/raw
```

Specific archive URLs are used instead of searching arbitrary observations. If a link changes, the downloader reports the URL and stops. Files total hundreds of MB. The MAST directory is [here](https://archive.stsci.edu/pub/hlsp/acsggct/ngc6397/).

**2. Join the documented catalog to the matching FITS via sky WCS:**

```bash
python -m stellarweave.prepare_hst \
  --fits data/raw/ngc6397/hlsp_acsggct_hst_acs-wfc_ngc6397_f606w_v2_img.fits \
  --catalog data/raw/ngc6397/hlsp_acsggct_hst_acs-wfc_ngc6397_r.rdviq.cal.adj.zpt \
  --band f606w --alignment affine --tile-size 256 --guard 64 \
  --out data/processed/ngc6397
```

**3. Inspect the output before training.** `alignment_qc.png` plots candidate label positions on the actual image, and `manifest.json` records reference-quality cuts, SHA-256 input hashes, estimated translation, crossmatch inliers/RMS, and train/validation/test tile counts.

The default alignment initializes from a global translation and fits a **six-parameter affine transformation** with one-to-one matching, residual clipping, RMS/p95/coverage QC, and a 3×3 local-residual grid. It is not a complete optical-distortion or absolute astrometric solution. If auto alignment fails, investigate the frame; do **not** silently use unverified labels.

For diagnostic-only projection:

```bash
python -m stellarweave.prepare_hst --fits PATH --catalog PATH --alignment wcs --out data/diagnostic
```

This explicitly records `UNVERIFIED_WCS_ONLY`. An operator-specified `--alignment manual --dx N --dy N` is also unverified. Neither output is benchmark ground truth.

**Export structure:**

```text
data/processed/ngc6397/
  manifest.json
  alignment_qc.png
  index.csv
  train/*.npz
  val/*.npz
  test/*.npz
```

Each `.npz` contains raw float32 `image`, a `valid` mask, and `labels` with tile-local, 0-based (x,y), calibrated Vega magnitude, its error, and source ID. These are **point labels, not bounding boxes**; no arbitrary magnitude-to-flux conversion is performed.

The default data split uses spatially separated vertical stripes and a guard region to limit data leakage. Boundary tiles are intentionally discarded, **not counted as missed detections**. To reserve an entire different cluster: add `--split test-only`, and never tune thresholds on it.

## Two-filter association / F606W + F814W 关联（v0.3）

The public ACSGGCT NGC 6397 / NGC 104 / NGC 6752 directories contain a single **F606W reference FITS image** and a two-filter photometry catalog; **no corresponding F814W FITS reference image** is present in those directories. Earlier README instructions claiming that an F814W image can be downloaded from the same directory were incorrect. The downloader now rejects such a request.

The joint catalog already has F606W and F814W magnitudes for color–magnitude analyses. The optional `stellarweave.dualband` program requires a **separately obtained, matched-field F814W FITS file with independently checked WCS**; it does not fetch or validate such a file automatically.

```bash
python -m stellarweave.download_acsggct --cluster ngc6397 --band f606w
```

Products: `paired_reference.csv` with each common source ID, RA/Dec, independently mapped F606W and F814W pixel coordinates, calibrated Vega magnitudes, color, and errors; `dualband_manifest.json` includes per-filter affine residual QC and input SHA-256 hashes. The two bands are **not resampled onto one pixel grid**, and those magnitudes originate from the published catalog; no new photometric inference is performed.

The affine fitting uses a conservative initialization based on a near-identity WCS projection. Large initial rotations, scales, and incorrect WCS may fail closed rather than producing misleading labels. For field calibration one must examine the accepted matches spatially and adjust image/catalog selection or use professional astrometric tools.

Unverified `--alignment wcs` or manual-shift positions are **diagnostic only**, and v0.3 prevents their export into regular train/validation sets unless the output split is explicitly `test-only`. This is still not scientific benchmark truth.

## Scientific limitations / 科学限制

1. A catalog built from multiple exposures is not a perfect label set for one mosaic; faint or blended stars can be missing in either.
2. The actual `.adj.zpt` archive contains 23 columns. v0.6 correctly maps wV/wI, xsig/ysig, othV/othI and qfitV/qfitI before RA/Dec. These are survey diagnostics; their interpretation needs domain validation.
3. Verified affine correspondence is not complete astronomical astrometry; a QC plot and the 3×3 residual grid still require human inspection.
4. ACS `BUNIT`, exposure, zeropoints, DQ masks, correlated noise and PSF vary. This release does not do calibrated PSF flux fitting.
5. The separate GitHub workflow has downloaded and audited the real NGC 6397 FITS and catalog. This does not yet establish their star-by-star astrometric registration or scientifically calibrated PSF photometry.

## Roadmap / 开发路线

- **v0.3**: affine/local astrometric residual maps and reference-catalog two-band association implemented on constructed data; **real archive verification remains outstanding**.
- **v0.4**: synthetic and guarded FITS injection/recovery protocol implemented; real scientific calibration and noise modeling remain outstanding.
- **v0.5**: common-injection PSF and iterative-PSF baselines implemented; science-grade uncertainty validation remains future work.
- **v0.6**: real archive ingest, fixed 23-column catalog, DQ/ERR/BUNIT/WCS metadata audit, native-pixel PSF diagnostic; instrument-calibrated physical models remain future work.
- **v0.7**: actual NGC 6397 raw WCS/catalog matching and spatial native-pixel PSF diagnostics, with rejected affine QC explicitly recorded.
- **v0.8 (planned)**: independently verified plate alignment, cross-matched reference stars, per-location instrumental ePSF, physically consistent uncertainty model and genuine recovery/completeness testing.

Cite the [survey archive](https://archive.stsci.edu/prepds/acsggct/), [Sarajedini et al. (2007)](https://arxiv.org/abs/astro-ph/0612598), and [Anderson et al. (2008)](https://arxiv.org/abs/0804.2025) when using the data.

MIT License.
