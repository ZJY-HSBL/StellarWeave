# StellarWeave — 星织

**Astronomical Image Processing · Crowded-field Stellar Detection · Reliable Photometry**

[![CI](https://github.com/ZJY-HSBL/StellarWeave/actions/workflows/ci.yml/badge.svg)](https://github.com/ZJY-HSBL/StellarWeave/actions/workflows/ci.yml) · [GitHub Releases](https://github.com/ZJY-HSBL/StellarWeave/releases)

**v0.4.0 research prototype / 科研原型。** This repository does NOT yet contain a trained neural network, science-grade PSF photometry, or validated real-HST quantitative results. The catalog/image pipeline is tested against synthetic FITS/WCS fixtures, and any real ACS image/catalog pairing requires manual QC.

## Purpose / 研究目标

StellarWeave is designed to investigate faint stellar source detection, deblending, and measurement reliability in crowded telescope imagery. We are developing physically informed baselines, not merely applying YOLO to astronomy.

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

Use the [ACS Survey of Galactic Globular Clusters](https://archive.stsci.edu/prepds/acsggct/) from NASA/STScI's MAST archive. Its stellar photometric catalogs are **reference measurements**, not exhaustive ground truth. [Anderson et al. (2008)](https://arxiv.org/abs/0804.2025) Table 4 documents the original 21-column RDVIQ catalog: ID, reference-frame x/y, F606W and F814W magnitudes/errors, number of exposures, unsaturated/saturated flags, neighbor contamination, PSF quality, and RA/Dec.

Default NGC 6397, then independent clusters NGC 6752 and NGC 104. Note that the archived `.rdviq.cal.adj.zpt` file variant still requires a live check of column compatibility and alignment against the specific image.

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

First download the F606W/F814W science mosaics and the accompanying **shared** ACS photometric catalog from the [same cluster directory](https://archive.stsci.edu/pub/hlsp/acsggct/ngc6397/). The official catalog contains the identity and photometry in both filters; **do not** infer identity from a positional nearest-neighbor search across images.

```bash
python -m stellarweave.download_acsggct --cluster ngc6397 --band f606w
python -m stellarweave.download_acsggct --cluster ngc6397 --band f814w
python -m stellarweave.dualband \
  --f606w data/raw/ngc6397/hlsp_acsggct_hst_acs-wfc_ngc6397_f606w_v2_img.fits \
  --f814w data/raw/ngc6397/hlsp_acsggct_hst_acs-wfc_ngc6397_f814w_v2_img.fits \
  --catalog data/raw/ngc6397/hlsp_acsggct_hst_acs-wfc_ngc6397_r.rdviq.cal.adj.zpt \
  --out data/processed/ngc6397/paired
```

Products: `paired_reference.csv` with each common source ID, RA/Dec, independently mapped F606W and F814W pixel coordinates, calibrated Vega magnitudes, color, and errors; `dualband_manifest.json` includes per-filter affine residual QC and input SHA-256 hashes. The two bands are **not resampled onto one pixel grid**, and those magnitudes originate from the published catalog; no new photometric inference is performed.

The affine fitting uses a conservative initialization based on a near-identity WCS projection. Large initial rotations, scales, and incorrect WCS may fail closed rather than producing misleading labels. For field calibration one must examine the accepted matches spatially and adjust image/catalog selection or use professional astrometric tools.

Unverified `--alignment wcs` or manual-shift positions are **diagnostic only**, and v0.3 prevents their export into regular train/validation sets unless the output split is explicitly `test-only`. This is still not scientific benchmark truth.

## Scientific limitations / 科学限制

1. A catalog built from multiple exposures is not a perfect label set for one mosaic; faint or blended stars can be missing in either.
2. The original published RDVIQ column definitions may need checking against the archived `.adj.zpt` variant.
3. Verified affine correspondence is not complete astronomical astrometry; a QC plot and the 3×3 residual grid still require human inspection.
4. ACS `BUNIT`, exposure, zeropoints, DQ masks, correlated noise and PSF vary. This release does not do calibrated PSF flux fitting.
5. Real HST download + catalog correspondence **could not be run in this execution environment**; end-to-end tests use constructed FITS fixtures.

## Roadmap / 开发路线

- **v0.3**: affine/local astrometric residual maps and reference-catalog two-band association implemented on constructed data; **real archive verification remains outstanding**.
- **v0.4**: synthetic and guarded FITS injection/recovery protocol implemented; real scientific calibration and noise modeling remain outstanding.
- **v0.5**: PSF photometry baselines and uncertainty validation; then evaluate physically informed neural improvements.

Cite the [survey archive](https://archive.stsci.edu/prepds/acsggct/), [Sarajedini et al. (2007)](https://arxiv.org/abs/astro-ph/0612598), and [Anderson et al. (2008)](https://arxiv.org/abs/0804.2025) when using the data.

MIT License.
