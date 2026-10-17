"""Expanded crowding-stratified planted stars on real HST/ACS mosaics.

This is a *new, predeclared* injection-site protocol; it does not change
v0.14-v0.16 results. Sites are selected using only original 5-sigma peak
geometry, before observing injection recovery or fitted flux. Three finite
nearest-detected-peak bands replace the sparsely populated >=16px class in
NGC 104. All outcomes and missing bin capacities are retained.

Planted source shapes are actual held-out HST bright-star cutouts; local and
frozen-center native-pixel PRFs are derived only from non-injected train stars.
No physical Poisson noise, ePSF or HST zeropoint calibration is implied.
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
from .real_injection import (
    CLUSTERS, TILES, FLUX_RMS_MULTIPLIERS, inject_and_measure)
from .spatial_prf import extract_stamps, split_training_test
from .noise_stress import audit_science_noise_metadata

BANDS={'near_4_8':(4.,8.),'middle_8_12':(8.,12.),
       'outer_12_16':(12.,16.)}
PER_TILE_BAND=4
GRID_SPACING=4
SITES_MIN_SEPARATION=12.
SEED=42


def select_distance_bands(shape,valid,original_peaks,*,seed,
                          per_band=PER_TILE_BAND,spacing=GRID_SPACING,
                          minimum_inter_site_px=SITES_MIN_SEPARATION):
    """Fixed seeded site selection on ORIGINAL image, before injection.

    Not true astrophysical crowding. The 'outer' band is 12-16 pixels from
    another *detected* peak, and is emphatically NOT an isolated-star sample.
    """
    valid=np.asarray(valid,bool)
    if len(shape)!=2 or tuple(shape)!=valid.shape or spacing<1 or per_band<1:
        raise ValueError('Invalid geometry or selection protocol')
    original_peaks=np.asarray(original_peaks,float).reshape(-1,2)
    if not len(original_peaks):
        raise ValueError('Cannot define peak-distance bands without original peaks')
    yy,xx=np.mgrid[20:shape[0]-20:spacing,20:shape[1]-20:spacing]
    coordinates=np.column_stack((xx.ravel(),yy.ravel())).astype(int)
    rng=np.random.default_rng(seed)
    coordinates=coordinates[rng.permutation(len(coordinates))]
    distances=cKDTree(original_peaks).query(coordinates,k=1)[0]
    groups={}
    selected_overall=[]
    for name,(lo,hi) in BANDS.items():
        sites=[]
        for (x,y),distance in zip(coordinates,distances):
            if distance<lo or distance>=hi:continue
            if not valid[y-9:y+10,x-9:x+10].all():continue
            if selected_overall and min(
                    np.hypot(np.asarray(selected_overall)[:,0]-x,
                             np.asarray(selected_overall)[:,1]-y))<minimum_inter_site_px:
                continue
            selected_overall.append([int(x),int(y)])
            sites.append({'x':int(x),'y':int(y),
                          'nearest_original_detected_peak_px':float(distance)})
            if len(sites)>=per_band:break
        groups[name]=sites
    return groups


def tile_experiment(image,valid,stamps,xy,anchor_model,*,
                    cluster,tile_index,seed=SEED):
    split=split_training_test(stamps,xy,seed=seed+tile_index)
    if split is None:
        return [],{'tile':tile_index,'status':'INSUFFICIENT_PRF_TRAIN',
                   'bins':{k:0 for k in BANDS}}
    image=np.asarray(image,float)
    noise=float(robust_background(np.where(valid,image,np.nan))[1])
    peaks=scipy_peak_detector(np.where(valid,image,np.nan),fwhm=2.6,nsigma=5.)
    site_groups=select_distance_bands(image.shape,valid,peaks,
                                     seed=seed+1700+tile_index)
    rows=[]
    for band,sites in site_groups.items():
        for k,s in enumerate(sites):
            for fi,scale in enumerate(FLUX_RMS_MULTIPLIERS):
                # All three fluxes share one physical pixel site; they are
                # NOT independent injection positions.
                # Source-shape choice is fixed before inspecting outcomes.
                shape_index=(k+fi*PER_TILE_BAND+
                             list(BANDS).index(band)*PER_TILE_BAND*len(FLUX_RMS_MULTIPLIERS)) % len(split['holdout_stamps'])
                planted_flux=float(noise*scale)
                fit=inject_and_measure(image,valid,(s['x'],s['y']),
                          split['holdout_stamps'][shape_index],planted_flux,
                          peaks,split['template'],anchor_model)
                rows.append({
                    'cluster':cluster,'tile_index':tile_index,'band':band,
                    'site_index':k,'x':s['x'],'y':s['y'],
                    'nearest_original_detected_peak_px':s['nearest_original_detected_peak_px'],
                    'flux_over_tile_robust_pixel_rms':float(scale),
                    'injected_flux_image_units':planted_flux,
                    'heldout_stamp_index':shape_index,
                    'tile_robust_pixel_rms':noise,
                    **fit})
    meta={'tile':tile_index,'status':'EVALUATED','bins':{
          band:len(sites) for band,sites in site_groups.items()},
          'base_detected_peaks':len(peaks),'train_stars':split['train_count'],
          'heldout_stars':split['test_count'],'attempts':len(rows)}
    return rows,meta


def audit_fits_lineage(hdus,science_index):
    """Inspect exposure ancestry without claiming recovered input frames."""
    header_keys=('INSTRUME','DETECTOR','FILTER1','FILTER2','FILTER','BUNIT',
                 'EXPTIME','TEXPTIME','NCOMBINE','CCDGAIN','GAIN','DRIZSCAL',
                 'D001SCAL','PIXFRAC','EXPSTART','EXPEND','DATE-OBS','DATE-AVG')
    ext=[]
    for i,h in enumerate(hdus):
        hdr=h.header
        ext.append({'index':i,'extname':str(hdr.get('EXTNAME','PRIMARY')),
                    'pixel_dimensions':list(h.data.shape) if h.data is not None else None,
                    'header_values':{k:str(hdr[k]) for k in header_keys if k in hdr},
                    'history_record_count':len(hdr.get('HISTORY',[])) if 'HISTORY' in hdr else 0})
    selected=ext[science_index]
    return {
        'source_hdu_index':science_index,
        'hdu_inventory':ext,
        'selected_science_header':selected['header_values'],
        'underlying_exposure_ids_verified':False,
        'effective_exposure_map_verified':False,
        'pixelwise_noise_covariance_measured':False,
        'electrons_per_output_pixel_verified':False,
        'physical_poisson_noise_applied':False,
        'warning':'A stitched ACS science HDU and BUNIT=ELECTRONS do not establish independent detector-electron observations per mosaic pixel. Exposure lineage/weight maps must be independently verified.'}


def run_extended_real_injections(fits_path,plate_path,prf_path,out_dir,*,
                                 cluster,seed=SEED):
    from astropy.io import fits
    fits_path,plate_path,prf_path=map(Path,(fits_path,plate_path,prf_path))
    plate=json.loads(plate_path.read_text())
    prf=json.loads(prf_path.read_text())
    if (cluster not in CLUSTERS or plate.get('cluster')!=cluster
        or prf.get('cluster')!=cluster or plate.get('tile_size')!=512
        or plate.get('grid')!=3 or plate.get('anchor_index')!=4
        or plate.get('heldout_tiles_physically_nonoverlapping') is not True
        or plate.get('anchor_status')!='HELDOUT_QC_PASS'
        or prf.get('kind')!='REAL_HST_SPATIAL_NATIVE_PIXEL_PRF_HOLDOUT'
        or prf.get('seed')!=seed
        or prf.get('no_heldout_star_used_to_fit_a_scored_template') is not True):
        raise ValueError('Unverified HST source models or geometry')
    fits_sha=sha256_file(fits_path)
    if (fits_sha!=plate.get('fits_sha256')
        or fits_sha!=prf.get('fits_sha256')
        or sha256_file(plate_path)!=prf.get('plate_sha256')):
        raise ValueError('Original FITS/plate/PRF SHA-256 changed')
    patches={}
    with fits.open(fits_path,memmap=True) as hdus:
        sci=hdus[int(plate['science_hdu_index'])]
        if sci.data is None or sci.data.ndim!=2:
            raise ValueError('Missing HST science plane')
        lineage=audit_fits_lineage(hdus,int(plate['science_hdu_index']))
        extver=sci.header.get('EXTVER',1)
        dq=err=None
        for h in hdus:
            if h is sci or h.data is None or h.data.shape!=sci.data.shape:
                continue
            if h.header.get('EXTVER',1)!=extver:continue
            if h.header.get('EXTNAME','').upper()=='DQ':dq=h.data
            if h.header.get('EXTNAME','').upper()=='ERR':err=h.data
        noise_meta=audit_science_noise_metadata(
            sci.header,err_present=err is not None,dq_present=dq is not None)
        for entry in plate['tiles']:
            tile=int(entry['index'])
            x,y=map(int,entry['tile_xy'])
            bounds=(slice(y,y+512),slice(x,x+512))
            im=np.asarray(sci.data[bounds],float)
            if im.shape!=(512,512):raise ValueError('Disjoint HST crop not fully covered')
            good=np.isfinite(im)
            if dq is not None:good &= np.asarray(dq[bounds])==0
            if err is not None:
                e=np.asarray(err[bounds],float)
                good &= np.isfinite(e)&(e>0)
            stamps,xy,_=extract_stamps(im,good)
            patches[tile]=(im,good,stamps,xy)
    if set(patches)!=set(range(9)):
        raise ValueError('Missing HST tile')
    anchor=split_training_test(patches[4][2],patches[4][3],seed=seed+4)
    if anchor is None:raise ValueError('Insufficient frozen center PRF candidates')
    records=[];tiles=[]
    for i in TILES:
        image,valid,stamps,xy=patches[i]
        rows,status=tile_experiment(image,valid,stamps,xy,anchor['template'],
                                    cluster=cluster,tile_index=i,seed=seed)
        tiles.append(status)
        records.extend(rows)
    sites={(r['tile_index'],r['band'],r['site_index']) for r in records}
    if len(records)!=len(sites)*len(FLUX_RMS_MULTIPLIERS):
        raise ValueError('Every physical injection site must have all 3 flux levels')
    groups=[]
    for band in BANDS:
        for scale in FLUX_RMS_MULTIPLIERS:
            rs=[r for r in records if r['band']==band and
                r['flux_over_tile_robust_pixel_rms']==scale]
            recovered=sum(r['detected'] for r in rs)
            groups.append({'band':band,'flux_over_tile_rms':float(scale),
                'attempted':len(rs),'recovered':recovered,
                'recovery_fraction':recovered/len(rs) if rs else None,
                'median_local_forced_abs_flux_error':float(np.median([
                    r['forced_local_abs_rel_error'] for r in rs])) if rs else None,
                'median_center_forced_abs_flux_error':float(np.median([
                    r['forced_center_abs_rel_error'] for r in rs])) if rs else None,
                'median_matched_position_error_px':float(np.median([
                    r['detector_peak_position_error_px'] for r in rs if r['detected']]))
                    if recovered else None})
    report={
        'version':'0.17.0',
        'kind':'REAL_HST_EXPANDED_PREDECLARED_CROWDING_INJECTIONS',
        'cluster':cluster,'fits_sha256':fits_sha,
        'plate_sha256':sha256_file(plate_path),
        'prf_sha256':sha256_file(prf_path),
        'seed':seed,'bands':{k:list(v) for k,v in BANDS.items()},
        'targets_per_tile_band':PER_TILE_BAND,
        'grid_spacing_px':GRID_SPACING,'minimum_inter_site_separation_px':SITES_MIN_SEPARATION,
        'flux_levels':list(FLUX_RMS_MULTIPLIERS),
        'distinct_selected_injection_positions':len(sites),
        'injection_attempts':len(records),
        'recovered':sum(r['detected'] for r in records),
        'full_test_set_retained':True,
        'all_injection_sites_chosen_before_outcomes':True,
        'all_models_fitted_on_noninjected_training_stars_only':True,
        'source_stamps_from_heldout_real_HST_stars':True,
        'same_injection_positions_shared_by_flux_levels':True,
        'injected_poisson_noise':False,
        'mosaic_noise_physically_calibrated':False,
        'photometric_zeropoint_calibrated':False,
        'stellar_population_completeness_estimated':False,
        'forced_flux_uses_true_injection_position':True,
        'lineage_audit':lineage,
        'noise_metadata_audit':noise_meta,
        'tiles':tiles,'groups':groups,'trials':records,
        'limitations':[
          'Peak-distance bands refer to originally detected 5-sigma peaks; faint and unresolved neighbors are unmeasured',
          'The 3 injected flux levels at one pixel position are repeated conditions, not independent positions',
          'Only one held-out star-shape template per trial; no original-frame detector shot-noise realization',
          'Known truth position in forced PRF flux fit, not blind photometry',
          'No physical source completeness or instrumental magnitude calibration',
          'Mosaic exposure lineage and pixel covariance remain unverified']}
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    (out/'expanded_injections.json').write_text(
        json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    if records:
        with (out/'expanded_injections.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]))
            w.writeheader();w.writerows(records)
    (out/'mosaic_lineage.json').write_text(
        json.dumps({'cluster':cluster,'fits_sha256':fits_sha,
            'lineage_audit':lineage,'noise_metadata_audit':noise_meta},
             indent=2,allow_nan=False),encoding='utf-8')
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description='Actual HST mosaic expanded artificial-star injection')
    for arg in ('fits','plate-report','prf-report','cluster','out'):
        p.add_argument('--'+arg,required=True)
    a=p.parse_args(argv)
    d=run_extended_real_injections(a.fits,a.plate_report,a.prf_report,
                                   a.out,cluster=a.cluster)
    print(json.dumps({'cluster':d['cluster'],'unique_positions':
          d['distinct_selected_injection_positions'],
          'attempted':d['injection_attempts'],'recovered':d['recovered'],
          'bands':[{'band':t['band'],'attempted':t['attempted']}
                   for t in d['groups'] if t['flux_over_tile_rms']==80.],
          'physical_mosaic_noise_verified':False},indent=2))


if __name__=='__main__':
    main()
