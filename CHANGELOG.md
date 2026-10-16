# Changelog

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
