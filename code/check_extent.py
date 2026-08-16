import os
import numpy as np
from astropy.io import fits

OUT_DIR = os.path.expanduser('~/snowblind_test/out')

# candidate #973 position in exposure 5
f5 = f'{OUT_DIR}/jw01181002001_02101_00005_nrcalong_jump.fits'
f6 = f'{OUT_DIR}/jw01181002001_02101_00006_nrcalong_jump.fits'

for f, (xi, yi) in [(f5, (74, 302)), (f6, (83, 276))]:
    with fits.open(f) as hdul:
        groupdq = hdul['GROUPDQ'].data  # (nints, ngroups, ny, nx)
    # collapse across integrations/groups: any flag set at any read
    flagged_any = np.bitwise_or.reduce(groupdq, axis=(0,1))
    # look at a bigger box (25x25) around the candidate to see extent
    box = flagged_any[max(yi-15,0):yi+16, max(xi-15,0):xi+16]
    sat = (box & 2) > 0
    jump = (box & 4) > 0
    print(f'{os.path.basename(f)} around ({xi},{yi}):')
    print(f'  SATURATED footprint: {sat.sum()} px in 31x31 box, bounding box extent: {np.ptp(np.argwhere(sat), axis=0) if sat.any() else "none"}')
    print(f'  JUMP_DET footprint: {jump.sum()} px in 31x31 box, bounding box extent: {np.ptp(np.argwhere(jump), axis=0) if jump.any() else "none"}')
    # control: same-size box at a random offset position, same exposure
    cx, cy = xi + 200, yi + 200
    if cy+16 < flagged_any.shape[0] and cx+16 < flagged_any.shape[1]:
        cbox = flagged_any[cy-15:cy+16, cx-15:cx+16]
        csat = (cbox & 2) > 0
        cjump = (cbox & 4) > 0
        print(f'  CONTROL (offset +200,+200): SATURATED {csat.sum()} px, JUMP_DET {cjump.sum()} px in 31x31 box')
