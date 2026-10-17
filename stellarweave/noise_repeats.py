"""Repeat *unchanged* v0.14 artificial stars in five independently seeded noise fields.

The images/backgrounds and non-injected PRFs are real archived ACSGGCT HST
products; the added fields are artificial equal-pixel-RMS Gaussian controls.
Header BUNIT=ELECTRONS alone cannot prove that resampled/drizzled image
pixels are independent Poisson electron draws. No detector physics is claimed.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .noise_stress import (
    added_noise, one_augmented_trial, noise_horizontal_correlation,
    audit_science_noise_metadata, MODES, ADDED_RMS_FRACTION,
)
from .real_injection import TILES, CLUSTERS
from .spatial_prf import extract_stamps, split_training_test
from .detection import robust_background, scipy_peak_detector
from .real_hst import sha256_file

REPLICATES=(11,23,37,53,71)
ADDED_MODES=('iid_gaussian','correlated_gaussian')
VERSION='0.16.0'


def _trial_key(row):
    return (int(row['tile_index']),str(row['class']),
            int(row['location_number']),float(row['flux_over_robust_tile_rms']))


def validate_inputs(plate,prf,original,v15,*,cluster,seed,fits_sha,plate_sha,prf_sha,original_sha):
    if cluster not in CLUSTERS or any(r.get('cluster')!=cluster for r in
                                      (plate,prf,original,v15)):
        raise ValueError('Cluster mismatch in HST source reports')
    if (plate.get('grid')!=3 or plate.get('tile_size')!=512
        or plate.get('anchor_index')!=4
        or plate.get('heldout_tiles_physically_nonoverlapping') is not True
        or plate.get('anchor_status')!='HELDOUT_QC_PASS'
        or prf.get('kind')!='REAL_HST_SPATIAL_NATIVE_PIXEL_PRF_HOLDOUT'
        or prf.get('seed')!=seed
        or prf.get('no_heldout_star_used_to_fit_a_scored_template') is not True
        or original.get('kind')!='REAL_HST_PLANTED_HELDOUT_STAR_RECOVERY'
        or original.get('seed')!=seed
        or original.get('models_fitted_on_injection_sources') is not False
        or v15.get('kind')!='REAL_HST_FROZEN_INJECTION_CONTROLLED_NOISE_STRESS'
        or v15.get('true_HST_detector_noise_model') is not False):
        raise ValueError('Unverified original experiment protocol')
    if (plate.get('fits_sha256')!=fits_sha or prf.get('fits_sha256')!=fits_sha
        or original.get('fits_sha256')!=fits_sha or v15.get('fits_sha256')!=fits_sha
        or prf.get('plate_sha256')!=plate_sha
        or original.get('prf_sha256')!=prf_sha
        or v15.get('v014_baseline_sha256')!=original_sha):
        raise ValueError('Source checksum drift in original v0.14/v0.15 experiments')
    rows=original.get('trials',[])
    if len(rows)!=original.get('total_attempted') or sum(bool(x['detected']) for x in rows)!=original.get('total_recovered'):
        raise ValueError('Original detection trials missing')
    original_map={_trial_key(x):x for x in rows}
    if len(original_map)!=len(rows):
        raise ValueError('Duplicate original star identity')
    v15_baseline={_trial_key(x):x for x in v15.get('trials',[]) if x['mode']=='baseline'}
    if original_map.keys()!=v15_baseline.keys():
        raise ValueError('v0.15 baseline changed original star selection')
    for key,row in original_map.items():
        if bool(row['detected'])!=bool(v15_baseline[key]['detected']):
            raise ValueError('v0.15 baseline changed detection success')
        if int(row['source_holdout_index'])!=int(v15_baseline[key]['source_holdout_index']):
            raise ValueError('Original heldout shape changed between snapshots')
    return rows


def run_replicates(fits_path,plate_path,prf_path,original_path,v15_path,
                   out_dir,*,cluster,seed=42):
    from astropy.io import fits
    fits_path,plate_path,prf_path,original_path,v15_path=map(
        Path,(fits_path,plate_path,prf_path,original_path,v15_path))
    plate=json.loads(plate_path.read_text(encoding='utf-8'))
    prf=json.loads(prf_path.read_text(encoding='utf-8'))
    original=json.loads(original_path.read_text(encoding='utf-8'))
    v15=json.loads(v15_path.read_text(encoding='utf-8'))
    fits_sha=sha256_file(fits_path)
    old=validate_inputs(
        plate,prf,original,v15,cluster=cluster,seed=seed,
        fits_sha=fits_sha,plate_sha=sha256_file(plate_path),
        prf_sha=sha256_file(prf_path),original_sha=sha256_file(original_path))
    prepared={}
    with fits.open(fits_path,memmap=True) as hdus:
        sci=hdus[int(plate['science_hdu_index'])]
        if sci.data is None or sci.data.ndim!=2:
            raise ValueError('Missing two-dimensional science image')
        extver=sci.header.get('EXTVER',1)
        err=dq=None
        for h in hdus:
            if h is sci or h.data is None or h.data.shape!=sci.data.shape:
                continue
            if h.header.get('EXTVER',1)!=extver:
                continue
            if h.header.get('EXTNAME','').upper()=='ERR':err=h.data
            if h.header.get('EXTNAME','').upper()=='DQ':dq=h.data
        audit=audit_science_noise_metadata(sci.header,err_present=err is not None,
                                           dq_present=dq is not None)
        for item in sorted(plate['tiles'],key=lambda x:int(x['index'])):
            i=int(item['index'])
            x,y=map(int,item['tile_xy'])
            section=(slice(y,y+512),slice(x,x+512))
            im=np.asarray(sci.data[section],float)
            if im.shape!=(512,512):raise ValueError('HST mosaic tile outside science image')
            valid=np.isfinite(im)
            if dq is not None:valid &= np.asarray(dq[section])==0
            if err is not None:
                e=np.asarray(err[section],float)
                valid &= np.isfinite(e)&(e>0)
            stamps,xy,_=extract_stamps(im,valid)
            split=split_training_test(stamps,xy,seed=seed+i)
            if split is None and (i==4 or any(int(t['tile_index'])==i for t in old)):
                raise ValueError('Source model no longer fits the original selected region')
            prepared[i]=(im,valid,split)
    anchor=prepared[4][2]['template']
    details=[];noise_audits=[]
    for i in TILES:
        originals=[x for x in old if int(x['tile_index'])==i]
        if not originals:
            continue
        im,valid,split=prepared[i]
        if split is None:raise ValueError('Missing original local PRF')
        rms=float(robust_background(np.where(valid,im,np.nan))[1])
        for rep in REPLICATES:
            for mode_idx,mode in enumerate(ADDED_MODES):
                # The seed is stable, nonoverlapping across cluster, tile,
                # replicate and noise mode, and is explicitly saved.
                noise_seed=(seed+101*CLUSTERS.index(cluster)+
                            1009*i+104729*rep+1000003*mode_idx)
                noise=added_noise(im.shape,valid,mode=mode,
                                  tile_rms=rms,seed=noise_seed)
                augmented=im+noise
                before=scipy_peak_detector(np.where(valid,augmented,np.nan),
                                           fwhm=2.6,nsigma=5.)
                noise_audits.append({
                    'tile_index':i,'replicate_seed_id':rep,'mode':mode,
                    'generation_seed':noise_seed,
                    'added_noise_pixel_rms':float(np.std(noise[valid])),
                    'expected_pixel_rms':rms*ADDED_RMS_FRACTION,
                    'lag1_added_field_correlation':noise_horizontal_correlation(noise,valid)})
                for row in originals:
                    shape_idx=int(row['source_holdout_index'])
                    if not 0<=shape_idx<len(split['holdout_stamps']):
                        raise ValueError('Frozen injected shape ID no longer valid')
                    measured=one_augmented_trial(
                        augmented,valid,row,split['holdout_stamps'][shape_idx],
                        split['template'],anchor,before)
                    details.append({
                        'cluster':cluster,'tile_index':i,
                        'class':row['class'],
                        'location_number':row['location_number'],
                        'flux_over_robust_tile_rms':row['flux_over_robust_tile_rms'],
                        'x':row['x'],'y':row['y'],
                        'source_holdout_index':shape_idx,
                        'injected_flux_image_units':row['injected_flux_image_units'],
                        'replicate_seed_id':rep,
                        'mode':mode,
                        'generation_seed':noise_seed,**measured})
    expected=len(old)*len(REPLICATES)*len(ADDED_MODES)
    if len(details)!=expected:
        raise ValueError('Some original planted sources were omitted from repeats')
    report={
        'version':VERSION,
        'kind':'REAL_HST_FROZEN_PLANTED_STARS_REPEATED_NOISE_SENSITIVITY',
        'cluster':cluster,'fits_sha256':fits_sha,
        'plate_sha256':sha256_file(plate_path),
        'prf_sha256':sha256_file(prf_path),
        'v014_sha256':sha256_file(original_path),
        'v015_sha256':sha256_file(v15_path),
        'replicate_seeds':list(REPLICATES),
        'modes':list(ADDED_MODES),
        'original_sites':len(old),
        'repeated_noise_trials':len(details),
        'repeated_noise_audits':len(noise_audits),
        'added_rms_fraction_observed_tile_rms':ADDED_RMS_FRACTION,
        'gaussian_smoothing_sigma_px':1.3,
        'no_original_sites_removed':True,
        'no_source_model_retrained_on_injected_stars':True,
        'star_positions_flux_shapes_frozen':True,
        'real_hst_background':True,
        'calibrated_instrument_noise':False,
        'mosaic_electron_counts_verified':False,
        'physical_poisson_source_noise_added':False,
        'drizzle_covariance_fitted':False,
        'known_coordinate_forced_amplitude_not_blind_photometry':True,
        'observation_metadata':audit,
        'trials':details,
        'noise_audits':noise_audits,
        'caveats':[
            'Repeated noise realizations are paired on the same non-randomly selected injection sites; the between-seed range is not a population-level confidence interval',
            'Correlated noise is constructed by Gaussian smoothing, not derived from HST Drizzle weights',
            'Header electron BUNIT alone cannot validate a physical Poisson mosaic-noise model',
            'Only controlled extra-noise fluctuations are repeated; actual telescope image stays fixed',
            'Source shape templates are bright isolated observed star cutouts, not calibrated HST ePSFs',
            'Detection and known-truth-location forced photometry remain distinct metrics',
        ]}
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    (out/'noise_repeats.json').write_text(
        json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    with (out/'noise_repeats.csv').open('w',newline='',encoding='utf-8') as f:
        if details:
            w=csv.DictWriter(f,fieldnames=list(details[0]))
            w.writeheader();w.writerows(details)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description='Repeat frozen real-HST artificial stars under added noise')
    for x in ('fits','plate-report','prf-report','original-report','v15-report','cluster','out'):
        p.add_argument('--'+x,required=True)
    a=p.parse_args(argv)
    result=run_replicates(a.fits,a.plate_report,a.prf_report,
                          a.original_report,a.v15_report,a.out,cluster=a.cluster)
    print(json.dumps({'cluster':result['cluster'],'original_sites':result['original_sites'],
                      'noise_trials':result['repeated_noise_trials'],
                      'noise_audits':result['repeated_noise_audits'],
                      'physical_noise_calibrated':False},indent=2))


if __name__=='__main__':
    main()
