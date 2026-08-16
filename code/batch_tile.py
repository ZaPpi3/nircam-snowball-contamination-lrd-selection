'''
General-purpose: given a proposal id + obs_key + list of (id, ra, dec) candidates
in that tile, download F444W uncal files (both nrc long modules) if not already
present, run Detector1Pipeline with SnowblindStep+JumpPlusStep post-hooks, then
check each candidate's position for SATURATED/JUMP_DET DQ flags vs a control.
Idempotent: skips download/pipeline steps if outputs already exist.
'''
import os, sys, glob
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astroquery.mast import Observations
from jwst.pipeline import Detector1Pipeline
from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import ImageModel
from snowblind import SnowblindStep, JumpPlusStep

BASE = os.path.expanduser('~/snowblind_test/batch')
os.makedirs(BASE, exist_ok=True)

STEPS = {
    'jump': {'save_results': True, 'expand_large_events': False,
             'post_hooks': [SnowblindStep, JumpPlusStep]},
    'ramp_fit': {'save_results': True},
    'gain_scale': {'skip': True},
}

def download_tile(proposal_id, obs_key):
    tile_dir = os.path.join(BASE, obs_key, 'uncal')
    os.makedirs(tile_dir, exist_ok=True)
    existing = glob.glob(f'{tile_dir}/**/*_uncal.fits', recursive=True)
    if len(existing) >= 2:
        print(f'  [{obs_key}] uncal already present ({len(existing)} files), skip download')
        return tile_dir
    obs = Observations.query_criteria(obs_collection='JWST', proposal_id=proposal_id,
                                       instrument_name='NIRCAM*', filters='F444W',
                                       obs_id=f'{obs_key}_nircam_clear-f444w')
    if len(obs) == 0:
        print(f'  [{obs_key}] WARNING: no F444W observation found for proposal {proposal_id}')
        return None
    prods = Observations.get_product_list(obs)
    uncal = prods[[str(x).endswith('uncal.fits') for x in prods['productFilename']]]
    print(f'  [{obs_key}] downloading {len(uncal)} uncal files...')
    Observations.download_products(uncal, download_dir=tile_dir)
    return tile_dir

def run_pipeline(obs_key, tile_dir):
    out_dir = os.path.join(BASE, obs_key, 'out')
    os.makedirs(out_dir, exist_ok=True)
    uncal_files = sorted(glob.glob(f'{tile_dir}/**/*_uncal.fits', recursive=True))
    for f in uncal_files:
        expected = os.path.join(out_dir, os.path.basename(f).replace('_uncal.fits', '_rate.fits'))
        if os.path.exists(expected):
            continue
        print(f'  [{obs_key}] processing {os.path.basename(f)}')
        try:
            Detector1Pipeline.call(f, steps=STEPS, output_dir=out_dir, save_results=True)
        except Exception as e:
            print(f'  [{obs_key}] FAILED on {os.path.basename(f)}: {e}')
    return out_dir

def check_candidates(obs_key, out_dir, candidates):
    rate_files = sorted(glob.glob(f'{out_dir}/*_rate.fits'))
    results = {}
    for cand_id, ra, dec in candidates:
        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        hits = []
        for rf in rate_files:
            base = rf.replace('_rate.fits', '')
            jump_f = base + '_jump.fits'
            if not os.path.exists(jump_f):
                continue
            try:
                model = ImageModel(rf)
                wcs_model = AssignWcsStep.call(model, skip=False)
                gwcs = wcs_model.meta.wcs
                x, y = gwcs.world_to_pixel(coord)
                x, y = float(x), float(y)
            except Exception:
                continue
            ny, nx = model.data.shape[-2], model.data.shape[-1]
            if not (np.isfinite(x) and np.isfinite(y) and 15 <= x < nx-215 and 15 <= y < ny-215):
                continue
            xi, yi = int(round(x)), int(round(y))
            with fits.open(jump_f) as hdul:
                groupdq = hdul['GROUPDQ'].data
            flagged_any = np.bitwise_or.reduce(groupdq, axis=(0,1))
            box = flagged_any[yi-15:yi+16, xi-15:xi+16]
            sat = int(((box & 2) > 0).sum())
            jump = int(((box & 4) > 0).sum())
            cbox = flagged_any[yi+185:yi+216, xi+185:xi+216]
            csat = int(((cbox & 2) > 0).sum())
            cjump = int(((cbox & 4) > 0).sum())
            hits.append((os.path.basename(rf), xi, yi, sat, jump, csat, cjump))
        results[cand_id] = hits
        print(f'  [{obs_key}] #{cand_id} (RA {ra:.5f} Dec {dec:.5f}): {len(hits)} exposure(s) on-detector')
        for h in hits:
            print(f'      {h[0]}: SAT={h[3]} JUMP={h[4]} | control SAT={h[5]} JUMP={h[6]}')
    return results

def process_tile(proposal_id, obs_key, candidates):
    print(f'=== TILE {obs_key} ({len(candidates)} candidate(s)) ===')
    tile_dir = download_tile(proposal_id, obs_key)
    if tile_dir is None:
        return {}
    out_dir = run_pipeline(obs_key, tile_dir)
    return check_candidates(obs_key, out_dir, candidates)
