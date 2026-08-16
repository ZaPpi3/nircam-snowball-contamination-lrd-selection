"""
Ad hoc cutout generator for the 8 A2744 + 8 GOODS-N candidates flagged by the
Playground LRD-color-systematics check (compact: 0.65<=concentration<=1.0,
break_shape=='smooth_dusty_interloper_like') as a crude LRD-like proxy.
Neither field had ANY visual vetting of its LRD-like-proxy candidates before
this (unlike JADES, where visual vetting is what caught the ~60-80% blocky-
pixel-artifact contamination rate in the equivalent compact+red signature
space -- see NOTES.md 2026-07-16/17). This exists to close that gap, not as
a new general-purpose tool.
"""
import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.visualization import ZScaleInterval, ImageNormalize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CUTOUT_SIZE_ARCSEC = 3.0

A2744_DIR = "F:/Claude/JWST/jwst_data_a2744/mastDownload/JWST"
GOODSN_DIR = "F:/Claude/JWST/jwst_data_goodsn/mastDownload/JWST"
OUT_DIR = "F:/Claude/JWST/candidate_cutouts_lrdproxy"
os.makedirs(OUT_DIR, exist_ok=True)

A2744_CANDIDATES = [
    ("jw02561-o001_t003", 188, 3.566277, -30.420185),
    ("jw02561-o001_t003", 870, 3.519088, -30.371296),
    ("jw02561-o001_t003", 3591, 3.632713, -30.399583),
    ("jw02561-o001_t003", 3822, 3.635645, -30.396787),
    ("jw02561-o001_t003", 4671, 3.634660, -30.376063),
    ("jw02561-o002_t001", 97, 3.654100, -30.453474),
    ("jw02561-o002_t001", 4731, 3.590348, -30.426203),
    ("jw02561-o003_t006", 1753, 3.561357, -30.326939),
]

GOODSN_CANDIDATES = [
    ("jw01181-o002_t006", 973, 189.091878, 62.299907),
    ("jw01181-o002_t006", 3230, 189.056865, 62.267137),
    ("jw01181-o098_t010", 765, 189.356364, 62.257104),
    ("jw01181-o098_t010", 1206, 189.344796, 62.240286),
    ("jw01181-o006_t008", 450, 189.331302, 62.188114),
    ("jw01181-o006_t008", 2779, 189.244822, 62.196313),
    ("jw01181-o001_t001", 268, 189.197422, 62.299916),
    ("jw01181-o001_t001", 1933, 189.204792, 62.245553),
]


def make_cutout(field, data_dir, obs_key, cand_id, ra, dec, bands):
    coord = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')
    fig, axes = plt.subplots(1, len(bands), figsize=(5 * len(bands), 5))
    if len(bands) == 1:
        axes = [axes]
    for ax, (label, filt) in zip(axes, bands):
        path = f"{data_dir}/{obs_key}_nircam_clear-{filt}/{obs_key}_nircam_clear-{filt}_i2d.fits"
        if not os.path.exists(path):
            ax.set_title(f"{label}: missing")
            ax.axis('off')
            continue
        with fits.open(path) as hdul:
            data = hdul['SCI'].data
            header = hdul['SCI'].header
        wcs = WCS(header)
        pixscale = abs(header['CD1_1']) * 3600 if 'CD1_1' in header else abs(header.get('CDELT1', 0.03)) * 3600
        size_px = int(round(CUTOUT_SIZE_ARCSEC / pixscale))
        cutout = Cutout2D(np.nan_to_num(data, nan=0.0), coord, size=size_px, wcs=wcs, mode='partial', fill_value=0.0)
        norm = ImageNormalize(cutout.data, interval=ZScaleInterval())
        ax.imshow(cutout.data, origin='lower', cmap='gray', norm=norm)
        ax.set_title(label)
        ax.axis('off')
        cx, cy = cutout.to_cutout_position(cutout.input_position_original)
        circ = plt.Circle((cx, cy), radius=size_px * 0.06, edgecolor='lime', facecolor='none', linewidth=1.5)
        ax.add_patch(circ)
    fig.suptitle(f"{field} {obs_key}, candidate #{cand_id} - RA {ra:.5f} Dec {dec:.5f}")
    fig.tight_layout()
    outpath = f"{OUT_DIR}/{field.lower()}_{obs_key}_{cand_id}.png"
    fig.savefig(outpath, dpi=130)
    plt.close(fig)
    print(f"Saved {outpath}")


for obs_key, cand_id, ra, dec in A2744_CANDIDATES:
    make_cutout("A2744", A2744_DIR, obs_key, cand_id, ra, dec,
                [("Red (F444W)", "f444w"), ("Mid (F150W)", "f150w"), ("Blue (F115W)", "f115w")])

for obs_key, cand_id, ra, dec in GOODSN_CANDIDATES:
    make_cutout("GOODSN", GOODSN_DIR, obs_key, cand_id, ra, dec,
                [("Red (F444W)", "f444w"), ("Mid (F200W)", "f200w"), ("Blue (F090W)", "f090w")])

print("Done.")
