# StellarWeave — 星织

**Astronomical Image Processing · Crowded-field Stellar Detection · Reliable Photometry**

[![CI](https://github.com/ZJY-HSBL/StellarWeave/actions/workflows/ci.yml/badge.svg)](https://github.com/ZJY-HSBL/StellarWeave/actions/workflows/ci.yml) · [GitHub Releases](https://github.com/ZJY-HSBL/StellarWeave/releases)

**v0.2.0 research prototype / 科研原型。** This repository does NOT yet contain a trained neural network, science-grade PSF photometry, or validated real-HST quantitative results. The catalog/image pipeline is tested against synthetic FITS/WCS fixtures, and any real ACS image/catalog pairing requires manual QC.

## Purpose / 研究目标

StellarWeave is designed to investigate faint stellar source detection, deblending, and measurement reliability in crowded telescope imagery. We are developing physically informed baselines, not merely applying YOLO to astronomy.

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
  --band f606w --alignment auto --tile-size 256 --guard 64 \
  --out data/processed/ngc6397
```

**3. Inspect the output before training.** `alignment_qc.png` plots candidate label positions on the actual image, and `manifest.json` records reference-quality cuts, SHA-256 input hashes, estimated translation, crossmatch inliers/RMS, and train/validation/test tile counts.

The auto-registration attempts only a **global shift** and rejects matches with too few one-to-one inliers or excessive residual. It cannot guarantee that spatially varying distortions and instrument calibration are correct. If auto alignment fails, investigate the frame; do **not** silently use unverified labels.

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

## Scientific limitations / 科学限制

1. A catalog built from multiple exposures is not a perfect label set for one mosaic; faint or blended stars can be missing in either.
2. The original published RDVIQ column definitions may need checking against the archived `.adj.zpt` variant.
3. Verified *translation-only* fit is not a complete astrometric solution; visual QC and local residual checks are required.
4. ACS `BUNIT`, exposure, zeropoints, DQ masks, correlated noise and PSF vary. This release does not do calibrated PSF flux fitting.
5. Real HST download + catalog correspondence **could not be run in this execution environment**; end-to-end tests use constructed FITS fixtures.

## Roadmap / 开发路线

- **v0.3**: real-field catalog verification, affine/local astrometric residual map, two-band associations.
- **v0.4**: artificial-star injection and recovery, completeness vs. brightness and stellar density.
- **v0.5**: PSF photometry baselines and uncertainty validation; then evaluate physically informed neural improvements.

Cite the [survey archive](https://archive.stsci.edu/prepds/acsggct/), [Sarajedini et al. (2007)](https://arxiv.org/abs/astro-ph/0612598), and [Anderson et al. (2008)](https://arxiv.org/abs/0804.2025) when using the data.

MIT License.
