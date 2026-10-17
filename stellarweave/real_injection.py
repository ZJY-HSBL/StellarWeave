"""Artificial-star injections into real HST F606W mosaics (v0.14).

Real observed backgrounds + held-out observed star shapes, not a synthetic
noise-free astronomy benchmark. Inputs are frozen 3x3 disjoint HST tiles and
the same deterministic 75/25 train/test PRF split as v0.12.

Each injection is added ONE AT A TIME. Both comparison PRFs are trained on
non-injected sources, never on the injected held-out shape. Detection is
re-run blind on the augmented image. Forced flux fits deliberately use the
KNOWN injection coordinate, and are reported separately from detection.
Physical flux calibration, stellar completeness and noise covariance are
not established by this experiment.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

from .detection import robust_background, scipy_peak_detector
from .real_hst import sha256_file
from .spatial_prf import extract_stamps, split_training_test

CLUSTERS=('ngc6397','ngc6752','ngc0104')
TILES=(0,1,2,3,5,6,7,8)
FLUX_RMS_MULTIPLIERS=(80.,200.,500.)
NEIGHBOR_CLASSES=('near_detected_peak','relatively_isolated')
INJECTIONS_PER_CELL=2
STAMP_SIZE=15
MATCH_RADIUS_PX=2.5


def predeclare_positions(image, valid, detections, *, seed, per_class=2):
    """Choose positions before examining any injected-image outcomes.

    Classes use separation from OBSERVED 5-sigma peaks; faint unresolved
    stars can still contaminate the "isolated" class. If a class cannot
    supply enough locations, preserve an explicit absent stratum.
    """
    good=np.asarray(valid,bool)
    if good.ndim!=2 or np.asarray(image).shape!=good.shape:
        raise ValueError('Mask/image geometry differs')
    if per_class<1:
        raise ValueError('Need >=1 location per crowding class')
    peaks=np.asarray(detections,float).reshape(-1,2)
    if not len(peaks):
        raise ValueError('Cannot define detection-based classes without peaks')
    # Fixed grid, deterministically shuffled. Tested sources and injected
    # detections never influence this geometry.
    yy,xx=np.mgrid[20:good.shape[0]-20:5,20:good.shape[1]-20:5]
    locations=np.column_stack([xx.ravel(),yy.ravel()]).astype(int)
    rng=np.random.default_rng(seed)
    locations=locations[rng.permutation(len(locations))]
    distances=cKDTree(peaks).query(locations,k=1)[0]
    positions={}
    for cls in NEIGHBOR_CLASSES:
        if cls=='near_detected_peak':
            selection=(distances>=5.)&(distances<11.)
        else:
            selection=distances>=16.
        coords=[]
        for (x,y),eligible in zip(locations,selection):
            if not eligible:
                continue
            if not good[y-9:y+10,x-9:x+10].all():
                continue
            coords.append([int(x),int(y)])
            if len(coords)==per_class:
                break
        positions[cls]=coords
    return positions


def forced_amplitude(image,position,psf):
    """One-source forced native-pixel PRF flux (unweighted, annular sky).

    This is NOT blinded source discovery; no physical ERR-plane weights or
    aperture/CTE/calibration corrections have been applied.
    """
    x,y=map(int,position)
    r=psf.shape[0]//2
    cut=np.asarray(image[y-r:y+r+1,x-r:x+r+1],float)
    if cut.shape!=psf.shape:
        raise ValueError('Forced fit outside image')
    yy,xx=np.mgrid[-r:r+1,-r:r+1]
    sky=np.median(cut[(xx*xx+yy*yy)>=(r*.88)**2])
    residual=cut-sky
    return float(np.sum(residual*psf)/np.sum(psf*psf))


def inject_and_measure(image,valid,position,source_shape,flux,baseline_peaks,
                       local_model,anchor_model):
    """Add one known-flux star; detect blind and fit two FROZEN PRFs."""
    if not np.isfinite(flux) or flux<=0:
        raise ValueError('Injected flux must be finite and positive')
    x,y=map(int,position)
    shape=np.asarray(source_shape,float)
    if shape.shape!=(STAMP_SIZE,STAMP_SIZE) or not np.isclose(shape.sum(),1,atol=1e-6):
        raise ValueError('Source must be a normalized independent 15x15 stamp')
    r=STAMP_SIZE//2
    new=np.asarray(image,float).copy()
    new[y-r:y+r+1,x-r:x+r+1]+=flux*shape
    before=np.asarray(baseline_peaks,float).reshape(-1,2)
    was_present=bool(len(before) and np.min(np.hypot(before[:,0]-x,before[:,1]-y))<=MATCH_RADIUS_PX)
    if was_present:
        raise ValueError('Injection location already has a baseline detection')
    found=scipy_peak_detector(np.where(valid,new,np.nan),fwhm=2.6,nsigma=5.)
    distances=np.linalg.norm(found-np.array([x,y]),axis=1) if len(found) else np.array([])
    matched=bool(len(distances) and float(distances.min())<=MATCH_RADIUS_PX)
    poserr=float(distances.min()) if matched else None
    fit_local=forced_amplitude(new,(x,y),local_model)
    fit_center=forced_amplitude(new,(x,y),anchor_model)
    return {'detected':matched,
            'detector_peak_position_error_px':poserr,
            'forced_local_flux':fit_local,
            'forced_center_flux':fit_center,
            'forced_local_signed_rel_error':(fit_local-flux)/flux,
            'forced_center_signed_rel_error':(fit_center-flux)/flux,
            'forced_local_abs_rel_error':abs(fit_local-flux)/flux,
            'forced_center_abs_rel_error':abs(fit_center-flux)/flux,
            'forcing_uses_known_position':True,
            'blind_detector':'scipy_peak_detector_fwhm2.6_nsigma5'}


def run_injection_tile(image,valid,stamps,xy,anchor_model,*,
                       tile_index,seed=42,per_class=INJECTIONS_PER_CELL):
    """Use EXACT deterministic training+held-out split as v0.12."""
    sp=split_training_test(stamps,xy,seed=seed+tile_index)
    if sp is None:
        return [],{'tile':tile_index,'status':'INSUFFICIENT_PRF_STARS',
                  'positions':{},'injections':0}
    local_model=sp['template']
    source_shapes=sp['holdout_stamps']
    noise=robust_background(np.where(valid,image,np.nan))[1]
    baseline=scipy_peak_detector(np.where(valid,image,np.nan),fwhm=2.6,nsigma=5.)
    positions=predeclare_positions(image,valid,baseline,
                                   seed=seed+500+tile_index,per_class=per_class)
    rows=[]
    for cls in NEIGHBOR_CLASSES:
        coords=positions[cls]
        for f_idx,scale in enumerate(FLUX_RMS_MULTIPLIERS):
            for k,position in enumerate(coords):
                # The injected PRF is a HELDOUT observed source stamp and
                # cannot appear in either training template.
                shape=source_shapes[(f_idx*per_class+k)%len(source_shapes)]
                flux=float(scale*noise)
                d=inject_and_measure(image,valid,position,shape,flux,baseline,
                                     local_model,anchor_model)
                rows.append({'tile_index':tile_index,'class':cls,
                    'location_number':k,'x':position[0],'y':position[1],
                    'flux_over_robust_tile_rms':float(scale),
                    'injected_flux_image_units':flux,
                    'tile_robust_pixel_rms':float(noise),
                    'source_holdout_index':(f_idx*per_class+k)%len(source_shapes),
                    **d})
    return rows,{'tile':tile_index,'status':'EVALUATED',
                 'positions':positions,'injections':len(rows),
                 'local_prf_training_stars':sp['train_count'],
                 'independent_injection_source_stamps':sp['test_count'],
                 'baseline_detected_peaks':len(baseline),
                 'tile_robust_pixel_rms':float(noise)}


def run_real_hst_injections(fits_path,plate_path,prf_path,out_dir, *,
                            cluster=None,seed=42):
    """Real HST pixels, known planted stars; save ALL negative results."""
    from astropy.io import fits
    fits_path,plate_path,prf_path=map(Path,[fits_path,plate_path,prf_path])
    plate=json.loads(plate_path.read_text(encoding='utf-8'))
    prf=json.loads(prf_path.read_text(encoding='utf-8'))
    cluster=cluster or plate.get('cluster')
    if (cluster not in CLUSTERS or plate.get('cluster')!=cluster
        or prf.get('cluster')!=cluster
        or plate.get('grid')!=3 or plate.get('tile_size')!=512
        or plate.get('anchor_index')!=4
        or plate.get('heldout_tiles_physically_nonoverlapping') is not True
        or plate.get('anchor_status')!='HELDOUT_QC_PASS'
        or prf.get('kind')!='REAL_HST_SPATIAL_NATIVE_PIXEL_PRF_HOLDOUT'
        or prf.get('seed')!=seed
        or prf.get('no_heldout_star_used_to_fit_a_scored_template') is not True):
        raise ValueError('Unverified HST frozen-PRF benchmark context')
    image_sha=sha256_file(fits_path)
    if (image_sha!=plate.get('fits_sha256')
        or image_sha!=prf.get('fits_sha256')
        or sha256_file(plate_path)!=prf.get('plate_sha256')):
        raise ValueError('FITS/plate/PRF SHA256 mismatch')
    data={}
    with fits.open(fits_path,memmap=True) as h:
        sci=h[int(plate['science_hdu_index'])]
        if sci.data is None or sci.data.ndim!=2:
            raise ValueError('No science plane')
        ver=sci.header.get('EXTVER',1)
        dq=err=None
        for x in h:
            if x is sci or x.data is None or x.data.shape!=sci.data.shape:
                continue
            if x.header.get('EXTVER',1)!=ver:
                continue
            if x.header.get('EXTNAME','').upper()=='DQ':
                dq=x.data
            if x.header.get('EXTNAME','').upper()=='ERR':
                err=x.data
        for t in plate['tiles']:
            ix=int(t['index'])
            x,y=map(int,t['tile_xy'])
            patch=np.asarray(sci.data[y:y+512,x:x+512],float)
            mask=np.isfinite(patch)
            if dq is not None: mask&=(np.asarray(dq[y:y+512,x:x+512])==0)
            if err is not None:
                e=np.asarray(err[y:y+512,x:x+512],float)
                mask&=np.isfinite(e)&(e>0)
            if patch.shape!=(512,512):
                raise ValueError('Incomplete tile image')
            stamps,xy,meta=extract_stamps(patch,mask)
            data[ix]=(patch,mask,stamps,xy)
    if set(data)!=set(range(9)):
        raise ValueError('Missing tile')
    anchor=split_training_test(data[4][2],data[4][3],seed=seed+4)
    if anchor is None:
        raise ValueError('No independent anchor PRF')
    all_rows=[];tiles=[]
    for idx in TILES:
        patch,mask,stamps,xy=data[idx]
        group,status=run_injection_tile(patch,mask,stamps,xy,anchor['template'],
                   tile_index=idx,seed=seed)
        tiles.append(status)
        for row in group:
            all_rows.append({'cluster':cluster,**row})
    groups=[]
    for cls in NEIGHBOR_CLASSES:
        for scale in FLUX_RMS_MULTIPLIERS:
            subset=[s for s in all_rows if s['class']==cls and
                    s['flux_over_robust_tile_rms']==scale]
            det=sum(s['detected'] for s in subset)
            groups.append({'class':cls,'flux_over_tile_rms':scale,
                           'attempted':len(subset),'recovered':det,
                           'recovery_fraction':det/len(subset) if subset else None,
                           'median_detection_error_px':float(np.median(
                               [s['detector_peak_position_error_px'] for s in subset
                                if s['detected']])) if det else None,
                           'median_abs_forced_local_relative_flux_error':float(
                               np.median([s['forced_local_abs_rel_error'] for s in subset]))
                               if subset else None,
                           'median_abs_forced_center_relative_flux_error':float(
                               np.median([s['forced_center_abs_rel_error'] for s in subset]))
                               if subset else None})
    result={'version':'0.14.0','kind':'REAL_HST_PLANTED_HELDOUT_STAR_RECOVERY',
            'cluster':cluster,'fits_sha256':image_sha,
            'plate_sha256':sha256_file(plate_path),
            'prf_sha256':sha256_file(prf_path),
            'seed':seed,'tile_ids':list(TILES),'tile_size':512,
            'crowding_bins':{
                'near_detected_peak':'[5,11) px from nearest original 5-sigma peak',
                'relatively_isolated':'>=16 px from nearest original 5-sigma peak'},
            'flux_levels':list(FLUX_RMS_MULTIPLIERS),
            'per_class_location_target':INJECTIONS_PER_CELL,
            'source_model':'frozen nontraining held-out observed HST star stamp',
            'injection':'one at a time, exact same real HST pixels for every trial',
            'detector':'SciPy smoothed-peak, 5sigma, distance <=2.5px',
            'forced_flux':'known-position unweighted one-source PRF amplitude; local annulus median sky',
            'metric_recovery':'incremental planted-source detection, NOT catalog completeness',
            'metric_flux':'relative injected flux-unit amplitude error, NOT physical/calibrated photometric bias',
            'matched_source_position_known_for_forced_flux':True,
            'science_calibrated':False,
            'instrument_epsf_calibrated':False,
            'flux_unit_calibrated':False,
            'real_background_images':True,
            'injected_sources_are_natural_or_observed_stars':False,
            'models_fitted_on_injection_sources':False,
            'total_attempted':len(all_rows),'total_recovered':sum(x['detected'] for x in all_rows),
            'tiles':tiles,'groups':groups,'trials':all_rows,
            'limitations':[
                'Injection flux in raw FITS mosaic units (not instrumental magnitude)',
                'Injected shapes cloned from held-out bright observed native-pixel stamps',
                'Position selection is a predeclared detected-peak-distance proxy, not true crowding',
                'Real backgrounds can contain unresolved or saturated stars',
                'No Poisson shot noise or detector covariance re-simulation added to injected stars',
                'Forced flux uses the true injected source location, unlike blind photometry',
                'No completeness in unknown real underlying star population is claimed']}
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    (out/'artificial_stars.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    with (out/'artificial_stars.csv').open('w',newline='',encoding='utf-8') as f:
        if all_rows:
            writer=csv.DictWriter(f,fieldnames=list(all_rows[0]))
            writer.writeheader();writer.writerows(all_rows)
        else:
            f.write('cluster,tile_index,class,detected\n')
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description='Real HST artificial star injection with frozen HST PRF')
    p.add_argument('--fits',required=True)
    p.add_argument('--plate-report',required=True)
    p.add_argument('--prf-report',required=True)
    p.add_argument('--cluster',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--seed',type=int,default=42)
    a=p.parse_args(argv)
    r=run_real_hst_injections(a.fits,a.plate_report,a.prf_report,a.out,
                              cluster=a.cluster,seed=a.seed)
    print(json.dumps({'cluster':r['cluster'],'attempted':r['total_attempted'],
                      'recovered':r['total_recovered'],
                      'groups':r['groups'],'calibrated':False},indent=2))


if __name__=='__main__':
    main()
