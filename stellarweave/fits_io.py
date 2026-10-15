"""HST FITS quicklook: preserve provenance; NOT calibrated photometry."""
import numpy as np

def load_fits_patch(path, size=512, x=None, y=None):
    from astropy.io import fits
    with fits.open(path, memmap=True) as hdul:
        choices = [h for h in hdul if h.data is not None
                   and getattr(h.data, "ndim", 0) == 2
                   and h.header.get("EXTNAME", "").upper() == "SCI"]
        if not choices:
            choices = [h for h in hdul if h.data is not None
                       and getattr(h.data, "ndim", 0) == 2]
        if not choices:
            raise ValueError("No 2D image HDU found")
        hdu = max(choices, key=lambda h: h.data.size)
        height, width = hdu.data.shape
        x0 = max(0, min(width - size, (width - size)//2 if x is None else x))
        y0 = max(0, min(height - size, (height - size)//2 if y is None else y))
        patch = np.array(hdu.data[y0:min(y0+size, height),
                                  x0:min(x0+size, width)], dtype=float)
        metadata = {"source_file": str(path), "hdu_name": hdu.name,
                    "offset_x": int(x0), "offset_y": int(y0),
                    "bunit": hdu.header.get("BUNIT", ""),
                    "exptime": hdu.header.get("EXPTIME", None)}
        return patch, metadata
