"""Find and download one public HST calibrated mosaic from MAST."""
import argparse
from pathlib import Path

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--object", default="NGC 6397")
    p.add_argument("--radius", default="0.04 deg")
    p.add_argument("--max-mb", type=float, default=250)
    p.add_argument("--out", default="data/hst")
    args = p.parse_args(argv)
    from astroquery.mast import Observations
    print("Querying MAST; internet access required")
    obs = Observations.query_object(args.object, radius=args.radius)
    if not len(obs):
        raise RuntimeError("No HST observations found for this target")
    eligible = [i for i, r in enumerate(obs)
                if str(r["obs_collection"]).upper() == "HST"
                and str(r["dataproduct_type"]).lower() == "image"
                and str(r["dataRights"]).upper() == "PUBLIC"]
    if not eligible:
        raise RuntimeError("No public HST image observations")
    preferred = [i for i in eligible
                 if "ACS/WFC" in str(obs[i]["instrument_name"]).upper()
                 and any(f in str(obs[i]["filters"]).upper() for f in ("F606W", "F814W"))]
    selected = obs[(preferred or eligible)[:60]]
    products = Observations.get_product_list(selected)
    candidates = []
    for r in products:
        filename = str(r["productFilename"]).lower()
        if not filename.endswith(("_drc.fits", "_drz.fits")):
            continue
        size_mb = float(r["size"]) / (1024 ** 2)
        if 0 < size_mb <= args.max_mb:
            candidates.append((size_mb, r))
    if not candidates:
        raise RuntimeError("No suitable size-capped DRC/DRZ found; select manually at https://mast.stsci.edu")
    size_mb, item = min(candidates, key=lambda x: x[0])
    destination = Path(args.out)
    destination.mkdir(parents=True, exist_ok=True)
    path = destination / str(item["productFilename"])
    print("Selected", path.name, round(size_mb, 1), "MB")
    status = Observations.download_file(item["dataURI"], local_path=str(path))
    print("MAST status:", status)
    print("Inspect FITS units, ERR/DQ extensions and associated catalog before any science claims.")

if __name__ == "__main__":
    main()
