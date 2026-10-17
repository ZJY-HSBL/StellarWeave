"""v0.15: controlled additive noise stress tests on real HST mosaics.

Original v0.14 stars, flux levels, coordinates and held-out PRFs are FROZEN.
The background is a real ACSGGCT HST science tile; two added perturbations
have equal pixel RMS: iid Gaussian and correlated Gaussian (sigma=1.3 px).

IMPORTANT: The added fields are not physically calibrated detector/Poisson
noise or an inferred Drizzle covariance. The original real HST noise remains.
No photometric or astronomical completeness certification is claimed.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from .detection import robust_background, scipy_peak_detector
from .real_hst import sha256_file
from .real_injection import STAMP_SIZE, MATCH_RADIUS_PX, CLUSTERS, TILES, \
    forced_amplitude
from .spatial_prf import extract_stamps, split_training_test

MODES = ('baseline', 'iid_gaussian', 'correlated_gaussian')
ADDED_RMS_FRACTION = 0.5
CORRELATION_SIGMA_PX = 1.3


def added_noise(shape,valid,*,mode,tile_rms,seed,
                rms_fraction=ADDED_RMS_FRACTION):
    """Equal per-pixel RMS modes, preserving predeclared seed and geometry."""
    if mode not in MODES:
        raise ValueError(f'Unknown noise stress mode {mode}')
    mask=np.asarray(valid,bool)
    if shape!=mask.shape or not np.any(mask):
        raise ValueError('Invalid noise mask')
    if not np.isfinite(tile_rms) or tile_rms<=0:
        raise ValueError('Invalid input image RMS')
    if not np.isfinite(rms_fraction) or not 0<rms_fraction<=2:
        raise ValueError('Invalid controlled noise strength')
    if mode=='baseline':
        return np.zeros(shape,dtype=float)
    rng=np.random.default_rng(int(seed))
    draw=rng.standard_normal(shape)
    if mode=='correlated_gaussian':
        draw=gaussian_filter(draw,CORRELATION_SIGMA_PX,mode='reflect')
    draw=draw-np.mean(draw[mask])
    std=float(np.std(draw[mask]))
    if std<=0 or not np.isfinite(std):
        raise ValueError('Degenerate controlled Gaussian realization')
    draw *= rms_fraction*tile_rms/std
    draw[~mask]=0.
    return draw


def noise_horizontal_correlation(noise,mask):
    """Lag-one Pearson product for input Gaussian field, not real HST noise."""
    a=np.asarray(noise,float)
    valid=np.asarray(mask,bool)
    good=valid[:,:-1]&valid[:,1:]
    if good.sum()<30:
        return None
    x=a[:,:-1][good]
    y=a[:,1:][good]
    sx=np.std(x);sy=np.std(y)
    return float(np.corrcoef(x,y)[0,1]) if sx>0 and sy>0 else None


def poisson_electron_fluctuation(expected_electrons,rng):
    """Unit-tested shot-noise primitive for VERIFIED electron-count inputs only.

    It is NOT automatically applied to drizzled HST science mosaics. A
    single output pixel in such data cannot in general be treated as an
    independent detector electron count.
    """
    value=np.asarray(expected_electrons,float)
    if np.any(~np.isfinite(value)) or np.any(value<0):
        raise ValueError('Expected nonnegative finite electrons')
    return rng.poisson(value).astype(float)-value


def audit_science_noise_metadata(header,*,err_present,dq_present):
    """Preserve raw metadata without declaring detector/electron calibration."""
    keys=('BUNIT','EXPTIME','TEXPTIME','CCDGAIN','GAIN','PHOTFLAM',
          'NCOMBINE','DRIZSCAL','D001SCAL','PIXFRAC','INSTRUME')
    return {'header_raw':{k:str(header[k]) for k in keys if k in header},
            'err_plane_available':bool(err_present),
            'dq_plane_available':bool(dq_present),
            'effective_exposure_and_gain_calibrated':False,
            'electron_counts_verified':False,
            'drizzle_covariance_calibrated':False,
            'poisson_shot_noise_applied_to_real_HST':False,
            'reason':'No validated electron-counts/exposure-weight/mosaic covariance transform; only controlled additive Gaussian noise permitted.'}


def one_augmented_trial(image,valid,trial,shape,local,anchor,baseline_peaks):
    """Same trial as v0.14; consider a preexisting noisy peak inconclusive."""
    x,y=int(trial['x']),int(trial['y'])
    flux=float(trial['injected_flux_image_units'])
    if shape.shape!=(STAMP_SIZE,STAMP_SIZE) or not np.isclose(shape.sum(),1,atol=1e-6):
        raise ValueError('Invalid holdout injection template')
    old=np.asarray(baseline_peaks,float).reshape(-1,2)
    preexisting=bool(len(old) and cKDTree(old).query([x,y],k=1)[0]<=MATCH_RADIUS_PX)
    r=STAMP_SIZE//2
    planted=np.asarray(image,float).copy()
    planted[y-r:y+r+1,x-r:x+r+1]+=flux*shape
    peaks=scipy_peak_detector(np.where(valid,planted,np.nan),fwhm=2.6,nsigma=5.)
    distance=float(cKDTree(peaks).query([x,y],k=1)[0]) if len(peaks) else None
    recovered=bool(not preexisting and distance is not None and distance<=MATCH_RADIUS_PX)
    amp_local=forced_amplitude(planted,(x,y),local)
    amp_anchor=forced_amplitude(planted,(x,y),anchor)
    return {
        'preexisting_noise_peak_at_truth':preexisting,
        'detected':recovered,
        'raw_nearest_postinject_peak_px':distance,
        'position_error_px':distance if recovered else None,
        'local_forced_abs_rel_error':abs(amp_local-flux)/flux,
        'center_forced_abs_rel_error':abs(amp_anchor-flux)/flux,
        'local_forced_signed_rel_error':(amp_local-flux)/flux,
        'center_forced_signed_rel_error':(amp_anchor-flux)/flux,
        'forced_amplitude_uses_true_coordinate':True}


def run_real_noise_stress(fits_path,plate_path,prf_path,base_path,out_dir,
                          *,cluster=None,seed=42):
    """Replay unchanged v0.14 trials in three paired noise modes."""
    from astropy.io import fits
    paths=list(map(Path,(fits_path,plate_path,prf_path,base_path)))
    fits_path,plate_path,prf_path,base_path=paths
    plate=json.loads(plate_path.read_text())
    prf=json.loads(prf_path.read_text())
    base=json.loads(base_path.read_text())
    cluster=cluster or base.get('cluster')
    if (cluster not in CLUSTERS or any(
            x.get('cluster')!=cluster for x in (base,plate,prf))
            or base.get('kind')!='REAL_HST_PLANTED_HELDOUT_STAR_RECOVERY'
            or prf.get('kind')!='REAL_HST_SPATIAL_NATIVE_PIXEL_PRF_HOLDOUT'
            or prf.get('seed')!=seed or base.get('seed')!=seed
            or plate.get('anchor_index')!=4 or plate.get('grid')!=3
            or plate.get('tile_size')!=512
            or prf.get('no_heldout_star_used_to_fit_a_scored_template') is not True
            or base.get('models_fitted_on_injection_sources') is not False
            or plate.get('heldout_tiles_physically_nonoverlapping') is not True):
        raise ValueError('Unverified frozen v0.14 source provenance')
    fits_sha=sha256_file(fits_path)
    if (fits_sha!=plate.get('fits_sha256')
            or fits_sha!=prf.get('fits_sha256')
            or fits_sha!=base.get('fits_sha256')
            or sha256_file(plate_path)!=prf.get('plate_sha256')
            or sha256_file(prf_path)!=base.get('prf_sha256')):
        raise ValueError('Frozen FITS/plate/PRF source checksum mismatch')
    originals=base.get('trials',[])
    if len(originals)!=base.get('total_attempted') or sum(
            r['detected'] for r in originals)!=base.get('total_recovered'):
        raise ValueError('Baseline injection report incomplete')
    inputs={}
    with fits.open(fits_path,memmap=True) as h:
        science=h[int(plate['science_hdu_index'])]
        v=science.header.get('EXTVER',1)
        err=dq=None
        for item in h:
            if item is science or item.data is None or item.data.shape!=science.data.shape:
                continue
            if item.header.get('EXTVER',1)!=v: continue
            ext=item.header.get('EXTNAME','').upper()
            if ext=='ERR': err=item.data
            if ext=='DQ': dq=item.data
        audit=audit_science_noise_metadata(science.header,err_present=err is not None,
                                           dq_present=dq is not None)
        for t in sorted(plate['tiles'],key=lambda t:t['index']):
            i=int(t['index'])
            x,y=map(int,t['tile_xy'])
            sl=(slice(y,y+512),slice(x,x+512))
            im=np.asarray(science.data[sl],float)
            if im.shape!=(512,512):
                raise ValueError('Tile crop shape mismatch')
            mask=np.isfinite(im)
            if dq is not None: mask&=np.asarray(dq[sl])==0
            if err is not None:
                er=np.asarray(err[sl],float)
                mask&=np.isfinite(er)&(er>0)
            stamps,xy,_=extract_stamps(im,mask)
            split=split_training_test(stamps,xy,seed=seed+i)
            if split is None:
                if i==4 or any(t['tile_index']==i for t in originals):
                    raise ValueError('Baseline included stars from an unfit PRF tile')
                continue
            inputs[i]=(im,mask,split)
    anchor=inputs[4][2]['template']
    trials=[];noise_audits=[]
    keys=set()
    for trial in originals:
        key=(int(trial['tile_index']),trial['class'],
             int(trial['location_number']),float(trial['flux_over_robust_tile_rms']))
        if key in keys: raise ValueError('Duplicate base trial key')
        keys.add(key)
        if trial['tile_index'] not in inputs:
            raise ValueError('Missing original PRF tile for injection')
    for i in TILES:
        group=[t for t in originals if int(t['tile_index'])==i]
        if not group: continue
        im,mask,split=inputs[i]
        orig_rms=robust_background(np.where(mask,im,np.nan))[1]
        for j,mode in enumerate(MODES):
            perturb=added_noise(im.shape,mask,mode=mode,tile_rms=orig_rms,
                                seed=seed+173*i+100003*j)
            noisy=im+perturb
            peaks=scipy_peak_detector(np.where(mask,noisy,np.nan),fwhm=2.6,nsigma=5.)
            noise_audits.append({
                'tile_index':i,'mode':mode,
                'added_noise_fraction_of_tile_rms':0 if j==0 else ADDED_RMS_FRACTION,
                'actual_added_pixel_rms':float(np.std(perturb[mask])),
                'target_added_pixel_rms':float(orig_rms)*ADDED_RMS_FRACTION if j else 0.,
                'added_field_lag1_correlation':noise_horizontal_correlation(perturb,mask),
                'baseline_peak_count_after_augmentation':len(peaks)})
            for trial in group:
                source_idx=int(trial['source_holdout_index'])
                if not 0<=source_idx<len(split['holdout_stamps']):
                    raise ValueError('Frozen heldout source index out of bounds')
                source=split['holdout_stamps'][source_idx]
                if mode=='baseline':
                    result={
                        'preexisting_noise_peak_at_truth':False,
                        'detected':bool(trial['detected']),
                        'position_error_px':trial['detector_peak_position_error_px'],
                        'local_forced_abs_rel_error':trial['forced_local_abs_rel_error'],
                        'center_forced_abs_rel_error':trial['forced_center_abs_rel_error'],
                        'local_forced_signed_rel_error':trial['forced_local_signed_rel_error'],
                        'center_forced_signed_rel_error':trial['forced_center_signed_rel_error'],
                        'raw_nearest_postinject_peak_px':trial['detector_peak_position_error_px'],
                        'forced_amplitude_uses_true_coordinate':True}
                else:
                    result=one_augmented_trial(noisy,mask,trial,source,
                                 split['template'],anchor,peaks)
                trials.append({
                    'cluster':cluster,'tile_index':i,
                    'class':trial['class'],'location_number':trial['location_number'],
                    'flux_over_robust_tile_rms':trial['flux_over_robust_tile_rms'],
                    'source_holdout_index':source_idx,
                    'x':trial['x'],'y':trial['y'],
                    'injected_flux_image_units':trial['injected_flux_image_units'],
                    'mode':mode,**result})
    # Every original trial must appear in every stress condition (including failures).
    if len(trials)!=len(originals)*len(MODES):
        raise ValueError('Some frozen v0.14 trials disappeared')
    result={
        'version':'0.15.0',
        'kind':'REAL_HST_FROZEN_INJECTION_CONTROLLED_NOISE_STRESS',
        'cluster':cluster,'fits_sha256':fits_sha,
        'plate_sha256':sha256_file(plate_path),
        'prf_sha256':sha256_file(prf_path),
        'v014_baseline_sha256':sha256_file(base_path),
        'seed':seed,
        'modes':list(MODES),
        'additional_gaussian_std_fraction_tile_rms':ADDED_RMS_FRACTION,
        'correlated_gaussian_filter_sigma_px':CORRELATION_SIGMA_PX,
        'noise_units':'relative to robust observed pixel RMS in image units',
        'true_HST_detector_noise_model':False,
        'real_HST_background':True,
        'no_heldout_source_in_training':True,
        'all_original_injection_positions_and_shapes_frozen':True,
        'added_noise_has_same_pixel_variance_across_modes':True,
        'photon_shot_noise_applied_to_real_HST':False,
        'instrument_psf_calibrated':False,
        'absolute_photometry_calibrated':False,
        'expected_original_trials':len(originals),
        'attempts':len(trials),
        'noise_audits':noise_audits,
        'header_noise_audit':audit,
        'trials':trials,
        'limitations':[
            'Correlated mode is Gaussian-filtered IID noise with matched pixel RMS, NOT a calibrated Drizzle noise covariance',
            'Original background contains the original telescope noise in every mode',
            'No source Poisson shot noise because mosaic electron/exposure mapping is unverified',
            'No completeness claims for unobserved underlying star population',
            'Known-source coordinates used for forced PRF amplitudes',
            'Extra-noise baseline peaks are checked to avoid false recovery at existing sources',
            'Exactly the same original v0.14 planting positions, shape IDs and flux values are used for every mode']
    }
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    (out/'noise_stress.json').write_text(
        json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    with (out/'noise_stress.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(trials[0]))
        w.writeheader();w.writerows(trials)
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description='Real HST frozen-artificial-star noise controls')
    p.add_argument('--fits',required=True)
    p.add_argument('--plate-report',required=True)
    p.add_argument('--prf-report',required=True)
    p.add_argument('--base-report',required=True)
    p.add_argument('--cluster',required=True)
    p.add_argument('--out',required=True)
    a=p.parse_args(argv)
    d=run_real_noise_stress(a.fits,a.plate_report,a.prf_report,a.base_report,
                            a.out,cluster=a.cluster)
    stats={}
    for mode in MODES:
        rows=[r for r in d['trials'] if r['mode']==mode]
        stats[mode]={'attempted':len(rows),'recovered':sum(r['detected'] for r in rows),
                     'already_present_peaks':sum(r['preexisting_noise_peak_at_truth']
                                                 for r in rows)}
    print(json.dumps({'cluster':d['cluster'],'modes':stats,
                      'HST_instrument_noise_calibrated':False},indent=2))


if __name__=='__main__':
    main()
