import glob, os
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import ImageModel

CAND_RA, CAND_DEC = 189.091878, 62.299907
OUT_DIR = os.path.expanduser('~/snowblind_test/out')

rate_files = sorted(glob.glob(f'{OUT_DIR}/*_rate.fits'))
print('rate files:', len(rate_files))
coord = SkyCoord(ra=CAND_RA, dec=CAND_DEC, unit='deg')

for rf in rate_files:
    base = rf.replace('_rate.fits', '')
    jump_f = base + '_jump.fits'
    if not os.path.exists(jump_f):
        print(os.path.basename(rf), 'no matching jump file, skip')
        continue
    try:
        model = ImageModel(rf)
        wcs_model = AssignWcsStep.call(model, skip=False)
        gwcs = wcs_model.meta.wcs
        x, y = gwcs.world_to_pixel(coord)
        x, y = float(x), float(y)
    except Exception as e:
        print(os.path.basename(rf), 'WCS/assign_wcs failed:', e)
        continue
    ny, nx = model.data.shape[-2], model.data.shape[-1]
    if not (0 <= x < nx and 0 <= y < ny):
        print(f'{os.path.basename(rf)}: off-detector at ({x:.1f},{y:.1f}) of {nx}x{ny}')
        continue
    xi, yi = int(round(x)), int(round(y))
    with fits.open(jump_f) as hdul:
        groupdq = hdul['GROUPDQ'].data
    box = groupdq[..., max(yi-3,0):yi+4, max(xi-3,0):xi+4]
    nonzero = np.count_nonzero(box)
    print(f'{os.path.basename(rf)}: ON-DETECTOR at pixel ({xi},{yi}) of {nx}x{ny} -- '
          f'nonzero GROUPDQ in 7x7 box: {nonzero}/{box.size}')
    if nonzero:
        print('   unique nonzero DQ values:', np.unique(box[box != 0]))
