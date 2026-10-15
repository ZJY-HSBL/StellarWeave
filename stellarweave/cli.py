"""Synthetic smoke experiment or position-only quicklook of a FITS image."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .simulation import simulated_field
from .detection import dao_detector, scipy_peak_detector
from .evaluation import detection_metrics, flux_binned_recall
from .fits_io import load_fits_patch

def main(argv=None):
    p = argparse.ArgumentParser(description="StellarWeave v0.1.0")
    p.add_argument("--input", help="Path to real FITS (omit for synthetic)")
    p.add_argument("--detector", choices=["scipy", "dao"], default="scipy")
    p.add_argument("--size", type=int, default=512)
    p.add_argument("--stars", type=int, default=300)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--fwhm", type=float, default=2.6)
    p.add_argument("--sigma", type=float, default=4)
    p.add_argument("--x", type=int)
    p.add_argument("--y", type=int)
    p.add_argument("--out", default="outputs/first_run")
    args = p.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.input:
        img, metadata = load_fits_patch(args.input, size=args.size, x=args.x, y=args.y)
        truth, mode = None, "REAL_FITS_without_ground_truth"
    else:
        img, truth = simulated_field(size=args.size, n_stars=args.stars,
                                     seed=args.seed, fwhm=args.fwhm)
        metadata, mode = {"seed": args.seed, "toy_gaussian_psf": True}, "SYNTHETIC_GAUSSIAN_PSF"
    detector = dao_detector if args.detector == "dao" else scipy_peak_detector
    pred = detector(img, fwhm=args.fwhm, nsigma=args.sigma)
    with (out / "detections.csv").open("w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["x_pixel", "y_pixel"])
        wr.writerows(pred)
    summary = {"mode": mode, "detector": args.detector, "sigma_threshold": args.sigma,
               "fwhm_assumed_px": args.fwhm, "image_shape": list(img.shape),
               "detections": len(pred), "metadata": metadata}
    if truth is not None:
        summary["metrics"] = detection_metrics(truth[:, :2], pred)
        summary["flux_binned_recall"] = flux_binned_recall(truth, pred)
        with (out / "synthetic_truth.csv").open("w", newline="", encoding="utf-8") as f:
            wr = csv.writer(f)
            wr.writerow(["x_pixel", "y_pixel", "injected_flux_electrons"])
            wr.writerows(truth)
    else:
        summary["warning"] = "No reference catalog. No precision/recall or photometry is reported."
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    fig, ax = plt.subplots(figsize=(8, 8))
    finite = img[np.isfinite(img)]
    low, high = np.percentile(finite, [5, 99.75])
    ax.imshow(img, cmap="gray", origin="lower", vmin=low, vmax=high)
    if len(pred):
        ax.scatter(pred[:, 0], pred[:, 1], facecolors="none",
                   edgecolors="tomato", linewidths=.65, s=25)
    ax.set_title(f"{mode} | {args.detector} | count={len(pred)}")
    ax.set_xlabel("x (pixels)")
    ax.set_ylabel("y (pixels)")
    fig.tight_layout()
    fig.savefig(out / "detections.png", dpi=175)
    plt.close(fig)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("Artifacts:", out.absolute())

if __name__ == "__main__":
    main()
