"""Download *specific* HST ACS Globular Cluster Survey reference products.

Source directory: https://archive.stsci.edu/pub/hlsp/acsggct/
URLs are explicit, and fail closed if the archive layout has changed.
"""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

TARGETS = ('ngc6397', 'ngc6752', 'ngc0104')
BASE = 'https://archive.stsci.edu/pub/hlsp/acsggct'

def filenames(target, band='f606w'):
    if target not in TARGETS:
        raise ValueError(f'Unsupported target: {target}. Allowed: {TARGETS}')
    if band != 'f606w':
        raise ValueError(
            'ACSGGCT HLSP reference directories archive only F606W FITS mosaics; '
            'F814W catalog magnitudes exist, but an F814W image must be obtained '
            'and validated separately from another HST observation.')
    image = f'hlsp_acsggct_hst_acs-wfc_{target}_{band}_v2_img.fits'
    catalog = f'hlsp_acsggct_hst_acs-wfc_{target}_r.rdviq.cal.adj.zpt'
    return image, catalog

def download_one(url, path, max_mb=500):
    """Stream data into .part, enforce a cap and reject HTML/short responses."""
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 512:
        if dest.suffix.lower() == '.fits':
            with dest.open('rb') as check:
                if check.read(8) != b'SIMPLE  ':
                    raise ValueError('Cached FITS file fails SIMPLE-header validation')
        return {'file': str(dest), 'bytes': dest.stat().st_size,
                'sha256': checksum(dest), 'cached': True}
    partial = dest.with_name(dest.name + '.part')
    maximum = int(max_mb * 1024 * 1024)
    sha = hashlib.sha256()
    count = 0
    try:
        request = Request(url, headers={'User-Agent': 'StellarWeave/0.6 (academic data client)'})
        with urlopen(request, timeout=60) as response, partial.open('wb') as out:
            length = response.headers.get('Content-Length')
            if length and int(length) > maximum:
                raise ValueError('Response exceeds size cap')
            first = response.read(8192)
            if first.lstrip().lower().startswith((b'<!doctype', b'<html')):
                raise ValueError('Archive returned HTML, not science data')
            if dest.suffix.lower() == '.fits' and not first.startswith(b'SIMPLE'):
                raise ValueError('Response is not a standard FITS primary header')
            chunks = [first]
            while True:
                for buf in chunks:
                    count += len(buf)
                    if count > maximum:
                        raise ValueError('Download exceeded size cap')
                    out.write(buf)
                    sha.update(buf)
                block = response.read(1024 * 1024)
                if not block:
                    break
                chunks = [block]
        if count < 512:
            raise ValueError('Unexpectedly small archive file')
        partial.replace(dest)
    except Exception:
        partial.unlink(missing_ok=True)
        raise
    return {'file': str(dest), 'bytes': count, 'sha256': sha.hexdigest(),
            'cached': False}

def checksum(path):
    h = hashlib.sha256()
    with open(path, 'rb') as inp:
        for chunk in iter(lambda: inp.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def main(argv=None):
    parser = argparse.ArgumentParser(description='Fetch public ACSGGCT image/catalog pair')
    parser.add_argument('--cluster', choices=TARGETS, default='ngc6397')
    parser.add_argument('--band', choices=['f606w','f814w'], default='f606w')
    parser.add_argument('--out', default='data/raw')
    parser.add_argument('--max-mb', type=float, default=500.)
    args = parser.parse_args(argv)
    image, catalog = filenames(args.cluster, args.band)
    directory = Path(args.out)/args.cluster
    result = {'cluster': args.cluster, 'band': args.band, 'files': {}}
    for name, fname in [('image',image), ('catalog',catalog)]:
        url = f'{BASE}/{args.cluster}/{fname}'
        try:
            info = download_one(url, directory/fname, max_mb=args.max_mb)
        except Exception as exc:
            raise RuntimeError(f'Failed to download {name} from {url}: {exc}. '
                               'Check the official data directory and product names.') from exc
        result['files'][name] = {'url': url, **info}
    (directory/'downloads.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
