# Changelog

## v0.14.0 — 2026-10-17

**Real HST ACSGGCT F606W artificial-star injection-and-recovery controls.**

- Added `stellarweave/real_injection.py`: inject normalized **held-out observed star stamps** one at a time into original non-synthetic HST mosaic images, with predeclared source flux relative to measured robust tile pixel RMS (80, 200, 500), and predetermined detection-neighbor distance classes. Rerun source detection blindly, count all successes/failures and calculate forced, known-coordinate local-vs-center native-pixel PRF flux amplitude errors.
- Added `stellarweave/injection_results.py`: real three-cluster trial aggregation, Wilson 95% binomial intervals, comparison plots, all unsuccessful injections, original source JSON with SHA-256 provenance and CSV output. No star is injected into a model used to fit that star.
- Actual experiments: NGC6397 84/96; NGC6752 91/96; NGC0104 54/60; **229/252 total**. NGC104 only has four sufficiently isolated positions per flux level under the fixed class definition; incomplete cells are recorded rather than silently filled.
- Example NGC104 near-existing-peak, low-injection-flux known-position absolute relative amplitude errors: local 35.56%, frozen center 30.74% — a negative result for local PRF. Local PRF is **not uniformly superior**. These are not calibrated stellar photometric biases.
- Published three genuine scientific figures, every attempted injection, per-group Wilson confidence intervals, complete raw reports and hash-linked summaries under `docs/results/real_injections_v14/` and embedded plots in README.
- Added regression tests verifying source placement, no pre-existing peak accepted as recovery, known-flux forced fit, held-out source identities, deterministic selection and chart/source integrity. Real observational workflow plus Python 3.11/3.12 CI.

**Scientific restrictions:** real HST background and real measured star-shape templates, but injected stars are artificial and have no additional simulated source photon shot noise, subpixel centroiding or instrument calibration. Flux fits use truth positions; results are *not* real catalog completeness, calibrated photometric error, Gaia astrometric precision or validated oversampled HST ePSF.

## v0.13.0 — 2026-10-17

**Real HST PSF/PRF failure analysis without test-set refitting or cherry-picking.**

- Added `stellarweave/prf_stratification.py` to reopen original HST FITS mosaics and compute peak-to-tile-RMS contrast, nearest 6-sigma detected-peak distance, detected neighbor count within 25 px and 5–10 px annular image texture for the exact v0.12 frozen test stars.
- Added `stellarweave/prf_stratification_summary.py`: validates every original held-out star and its unchanged L1 scores, verifies FITS/plate/PRF provenance hashes, and generates three real-data diagnostic PNGs plus complete per-star, group and threshold CSVs.
- Results: NGC 104 has median 16 detected neighbors within 25 px (NGC 6752: 7; NGC 6397: 3) and peak/tile-RMS contrast 218 (vs 778 and 879). Within NGC 104, low vs high peak-contrast subgroups show median star gains −1.02% vs +2.99%; low vs high neighbor-count subgroups show +1.30% vs −0.93%. **These are exploratory associations, not causal evidence or photometric improvements.**
- Retained all 838 frozen outer-tile test stars; annulus-texture proxy is missing for 20 sources, explicitly accounted for without dropping the baseline test records.
- Archived all three exact original PRF baseline reports used in this run alongside three stratification source reports to permit hash-checked regeneration even when independent runs differ in last-digit floating-point output.
- Extended three-cluster GitHub Actions to rerun the real HST pixel-proxy analysis, archive raw data and plot genuine results directly into README; added tests for source integrity, score tampering, full test-set retention and snapshot reproduction.

**Not achieved:** a causal account of NGC 104, calibrated photometric accuracy, physical source SNR, independent Gaia astrometry or an instrument-calibrated oversampled ePSF.

## v0.12.0 — 2026-10-17

**Cross-validated native-pixel PRF shape diagnostics on real HST/ACS mosaics.**

- Added `stellarweave/psf_spatial.py` and `stellarweave/psf_results.py` for nine-tile real PSF shape diagnostics (r50, r80, inner moments, core fraction, and total-variation distance to center). 27/27 candidate PSF tiles usable across NGC 6397, NGC 6752 and NGC 104.
- Added `stellarweave/spatial_prf.py`: strict disjoint train/test star split in every tile, frozen center-PRF vs local PRF vs fixed Gaussian baseline, per-star heldout L1 shape error, tile-bootstrap exploratory intervals, FITS/plate SHA-256 cross-checks and fail-closed geometry validation.
- Added `stellarweave/prf_visualization.py`: cross-cluster plots, all 838 outer-tile held-out star scores, 27 tile summaries, original reports, per-cluster PRF template montage PNGs, source SHA-256 and GitHub Actions provenance.
- Actual NGC 6397 median tile paired PRF gain +18.59% (tile-bootstrap 95% interval +5.66 to +26.52); NGC 6752 +16.83% (+4.59 to +25.66); NGC 104 +0.35% (−0.51 to +2.25) and pooled per-star median local L1 worse than frozen center (0.4000 vs 0.3861). **No universal improvement claim.**
- Published source data and genuine measured PNGs under `docs/results/real_hst_v12/` and `docs/results/real_prf_v12/`, with documented reproduction from archived evidence, including negative results.
- Added synthetic fixture tests for no leakage, sparse tiles, provenance rejection, source integrity and full dashboard regeneration; all science claims explicitly limited to native-pixel stamp shape rather than photometric precision.

**Not achieved:** independent Gaia-vs-HST astrometry, oversampled instrument ePSF, calibrated flux error/completeness, physical noise covariance or confirmed weighted HST mosaic observation epoch.

## v0.11.0 — 2026-10-17

**Real HST evidence published as data, charts and provenance rather than claims.**

- Added `stellarweave/result_visualization.py`, rendering three publication-style README figures directly from the three real ACSGGCT cross-region JSON reports. Includes per-tile median residuals, matched bright-reference fractions and catalog source-count comparisons.
- Added `docs/results/real_hst_v09/` with all three unfiltered raw reports, 24-row per-tile CSV, derived summary JSON, three PNG charts, SHA-256 provenance, and source GitHub Actions run ID.
- Actual data: NGC6397 2239/2266 matches, 0.399 px median-of-tile-medians; NGC6752 2374/2400, 0.407 px; NGC0104 2366/2400, 0.405 px. All 24 held-out tiles pass predeclared internal QC. Same-survey references: NOT independent absolute astrometry or detection completeness.
- Extended GitHub Actions to reproduce all figures and commit the results snapshot to main automatically, with rejection/failure retention and schema checks. Added tests that fail on missing/duplicated cluster reports, science-accuracy mislabeling and missing plot output.
- Added `stellarweave/epoch_audit.py` to inventory FITS header timing metadata without guessing an effective stacked-image epoch. Actual HST archive mosaics contain day-level header hints: NGC6397 MJD 53884, NGC6752 MJD 53879, NGC0104 MJD 53807. These are NOT independently validated weighted mosaic epochs.
- Included all epoch-audit JSON records alongside the figures and tests that enforce `effective_mosaic_epoch_verified=false`.

**Outstanding:** independent Gaia-vs-HST external astrometric accuracy, trained stellar source detection, local ePSF calibration, correlated-noise model, calibrated faint-star completeness, and measured science-grade photometric uncertainty.
## v0.10.0 — 2026-10-17

**Opt-in Gaia DR3 cross-survey validation, epoch propagation and fail-closed scientific QC.**

- Added an optional Gaia DR3 ESA TAP query client that derives a bounded cone from the target FITS celestial WCS, selects bright/astrometrically quality-filtered sources and logs exact query/provenance hashes.
- Added strict DR3 CSV parsing, source-ID checks, RUWE/magnitude/visibility-period filters, and Astropy proper-motion propagation from Gaia reference epoch to explicitly supplied HST observation MJD.
- Added per-withheld-tile comparison to a previously fitted frozen HST plate solution; Gaia never participates in model training or candidate selection.
- Added one-to-one candidate matches, position residual percentiles, displaced-catalog null controls and conservative non-calibration flags.
- Added a manually dispatched GitHub Gaia workflow requiring cluster name, independently verified HST observation MJD and held-out tile index.
- Added regression tests for Gaia filtering, proper-motion sign and scale, FITS WCS projection, excluded anchor tile and bounded ADQL query generation.

**Scope restriction:** live Gaia-to-real-HST measurements have NOT been run or validated as part of this release; CI results are from constructed Gaia-shaped data and existing real HST ACSGGCT-only tests. Without a justified mosaic reference epoch, external astrometric accuracy must not be claimed. Two-dimensional proper-motion propagation omits parallax and radial velocity contributions.
## v0.9.0 — 2026-10-17

**Real HST cross-region and independent three-cluster validation.**

- Added `stellarweave/survey_validation.py` and `survey_cli.py` for predeclared non-overlapping mosaic regions, central anchor affine registration and frozen-transform evaluation in eight disjoint tiles.
- Added per-region independent local fits as diagnostic-only comparisons, never substituted into the frozen-anchor transfer score.
- Added transfer-match fractions, median/p95 residuals, six displaced-catalog null controls, rejection reasons and source counts.
- Corrected the first full-extent mosaic tiling: many corner tiles have no reference stars. The final footprint-centered tiling is determined by median WCS-projected catalog positions before peak detection, not by outcome selection.
- Corrected reference splitting when sparse tiles have fewer reference stars than the v0.8 fixed holdout target.
- Added generated-FITS/WCS regression tests and a three-cluster GitHub Actions matrix with machine-readable per-cluster and aggregate summaries; acceptance thresholds prevent silent performance regressions.
- Genuine HST ACSGGCT catalog rows: NGC6397=16,100, NGC6752=52,817, NGC0104=158,154. Center-fit internal QC passes 3/3; fixed-transfer tile QC passes 8/8 for each cluster (24/24 total).

**Not proven:** calibrated absolute sky coordinates, independent external validation, per-star completeness or photometric accuracy, an adaptive HST ePSF model, or novel neural-model superiority.
## v0.8.0 — 2026-10-17

**Wide-shift astrometric search recovers real NGC 6397 catalog-to-image association.**

- Added bounded sparse displacement voting with multiple large-translation initialization proposals and iterative, geometry-constrained affine fits.
- Split reference stars into disjoint training and withheld validation groups; withheld stars never participate in selecting the transform.
- Added six displaced-catalog null-control checks, affine sanity gates, footprint coverage and held-out match/residual diagnostics.
- Added a real HST workflow diagnostic report and an overlay that distinguishes raw WCS positions from the accepted candidate affine solution.
- Real NGC 6397 central 512x512 crop: train matched **255/260**, held-out matched **147/150**; affine translation **(+10.9713, -0.9434) px**, train RMS **0.4058 px**, held-out median residual **0.3759 px**.
- The previous 3-pixel raw-WCS detector only matched 36/350 reference stars; the correct displacement explains the failure of the narrow initial search.
- Automated tests include wide-shift recovery with outliers, independent held-out matching, incorrect-field rejection and generated FITS/WCS fixtures.

**Caution:** HELDOUT_QC_PASS is internally consistent astrometry, not an externally calibrated plate solution. No measured true HST completeness, instrument ePSF or calibrated photometric uncertainty is claimed.
## v0.7.0 — 2026-10-17

**Real NGC 6397 spatial PSF diagnostics and falsifiable astrometric QC.**

- Added per-crop WCS projection and one-to-one source matching, with conservative affine proposal validation and explicit rejection status.
- Added grid-based native-pixel candidate PSF extraction and position-dependent PSF-shape residual diagnostics; produces JSON, PNG overlays and NPY PSF stamps.
- Added CLI and synthetic/WCS end-to-end regression tests.
- Extended the real HST GitHub workflow to produce spatial PSF grids, astrometry overlay and provenance-rich diagnostics.
- Real NGC 6397 (512px crop): 16,100 catalog records, 533 projected into crop, 350 candidate reference stars, 679 image peaks, and just 36 raw WCS matches at <=3px.
- Affine proposal **failed quality control**. The program correctly rejected the source mapping rather than exporting unverified training labels.
- All four grid cells produced preliminary native-pixel PSF candidate stacks. These are not spatially varying calibrated ePSFs.

**Outstanding:** validated plate solution, astrometric distortion model, instrument-calibrated ePSF, noise correlations, true real-HST photometric accuracy and artificial-star completeness.
## v0.6.0 — 2026-10-17

**Real ACSGGCT FITS/catalog audit and a critical column-schema correction.**

- Official NGC 6397 directory verified: one F606W reference FITS and one joint F606W/F814W star catalog; no corresponding F814W image in this HLSP folder.
- Corrected severe 21-column catalog assumption. The actual archive header contains 23 columns; RA and Dec are columns 22/23, not 20/21.
- Added strict header-aware catalog reader; updated source quality diagnostics and all synthetic regression fixtures.
- Added FITS SCI/ERR/DQ crop auditing, BUNIT/exposure/WCS metadata and SHA-256 provenance.
- Added isolated-star native-pixel empirical PRF stamp estimator, summary JSON and QA preview.
- Added a separate real-archive GitHub Actions workflow that downloads genuine NASA/STScI NGC 6397 products and uploads small QA artifacts.
- The initial real-run yielded 16,100 catalog entries, a 512×512 fully valid crop and 79 candidate PSF stars. No science-grade photometric uncertainty or completeness result is claimed.

**Unfinished:** quantitative real-image/cross-catalog registration, calibrated HST ePSF or noise model, detector gain and drizzle covariance, real-star flux-bias and completeness benchmarking.
## v0.5.0 — 2026-10-17

**Common-sample PSF fitting baseline for artificial-star recovery.**

- Added native Photutils 3.x PSFPhotometry and IterativePSFPhotometry methods with a pixel-integrated Gaussian PRF or externally supplied, unit-normalized ImagePSF.
- Added DAOStarFinder and SciPy centroid baselines with a shared detection/measurement output schema.
- Added a fair benchmark that uses identical injected artificial-star positions and instrumental fluxes for each method, each trial.
- Added matched-new-source recovery, crowding strata, localization and *isolated-only* fitted-flux error records, Wilson completeness intervals, and plots.
- Added FITS-mode hard checks for an empirical PRF stamp and BUNIT-consistent units; synthetic and science results are explicitly distinguished.
- Added Photutils science integration tests and CI benchmark smoke execution.

**Still not delivered:** real-HST completeness/flux-bias validation, spatially variable instrumental ePSF, physical variance maps or trustworthy uncertainty coverage, end-to-end DOLPHOT benchmark, trained neural model.

## v0.4.0 — 2026-10-17

**Artificial-star injection/recovery baseline with reproducible experimental protocol.**

- Implemented a unit-normalized, pixel-integrated Gaussian PRF generator and external empirical \`.npy\` PSF validation.
- Implemented flux-conserving subpixel stamp placement that fails closed near invalid pixels/edges.
- Implemented deterministic independent trial placement, source separation, and nearest-original-detection crowding bins.
- Reported recovery of *new* centroids only (existing baseline stars cannot be counted as injected-source hits), using one-to-one matches.
- Added per-star CSV logs, trial-wise provenance/flux accounting, flux/crowding completeness, approximate Wilson 95% intervals, and PNG curves.
- Guarded FITS mode with mandatory PSF stamp and explicit BUNIT-compatible instrumental flux units; included CLI and regression tests.
- Results are exploratory, not validated HST photometric completeness. PSF spatial variation and extra shot noise are not yet modeled.

## v0.3.0 — 2026-10-15

**Six-parameter astrometric registration and shared-ID two-filter catalog pairing.**

- Added translation-initialized iterative affine least-squares mapping with one-to-one matches and outlier clipping.
- Added fail-closed inlier fraction, RMS, 95th percentile residual, source field coverage and scale/shear sanity gates; record a 3-by-3 local residual diagnostic.
- Updated FITS-to-catalog dataset preparation to store full affine quality results rather than assuming a shift-only fit is sufficient.
- Added high-quality F606W/F814W common-ID catalog export with independent per-band WCS and affine mapping, calibrated catalog magnitudes, colors and measurement errors.
- Added generated-data regression and Astropy FITS/WCS integration tests.
- Disabled ordinary train/validation exports for manual or unverified WCS alignment.
- Documented that archived ACSGGCT product compatibility and real-field residuals remain unverified.

**Unfinished:** live MAST paired-file verification, high-precision crowded-star PSF photometry, artificial-star recovery, calibrated measurement uncertainty, any trained deep model.

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
