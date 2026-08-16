import sys, os, glob
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import ImageModel

def check(out_dir, ra, dec, label):
    coord = SkyCoord(ra=ra, dec=dec, unit='deg')
    rate_files = sorted(glob.glob(f'{out_dir}/*_rate.fits'))
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
        if not (np.isfinite(x) and np.isfinite(y) and 0 <= x < nx and 0 <= y < ny):
            continue
        xi, yi = int(round(x)), int(round(y))
        with fits.open(jump_f) as hdul:
            groupdq = hdul['GROUPDQ'].data
        flagged_any = np.bitwise_or.reduce(groupdq, axis=(0,1))
        box = flagged_any[max(yi-15,0):yi+16, max(xi-15,0):xi+16]
        sat = int(((box & 2) > 0).sum())
        jump = int(((box & 4) > 0).sum())
        cx, cy = xi + 200, yi + 200
        csat = cjump = None
        if cy+16 < flagged_any.shape[0] and cx+16 < flagged_any.shape[1]:
            cbox = flagged_any[cy-15:cy+16, cx-15:cx+16]
            csat = int(((cbox & 2) > 0).sum())
            cjump = int(((cbox & 4) > 0).sum())
        hits.append((os.path.basename(rf), xi, yi, sat, jump, csat, cjump))
    print(f'=== {label} (RA {ra:.5f} Dec {dec:.5f}): {len(hits)} exposure(s) on-detector ===')
    for h in hits:
        print(f'  {h[0]}: pix({h[1]},{h[2]}) SAT={h[3]} JUMP={h[4]} | control SAT={h[5]} JUMP={h[6]}')
    return hits

if __name__ == '__main__':
    out_dir, ra, dec, label = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    check(out_dir, ra, dec, label)
