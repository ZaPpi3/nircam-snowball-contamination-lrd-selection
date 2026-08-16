"""
Second-pass LRD-like proxy search on NEXUS Deep-tier ep01 data (F200W, F444W),
fixing the two real bugs found and documented in
Playground/little-red-dots-nircam-systematics/NOTES.md (2026-08-16) in the
first-pass script (nexus_lrdproxy_search.py):

1. Flat global noise estimate. The first pass used a single sigma-clipped
   std across the whole ~280-megapixel mosaic for significance. Real
   per-pixel error at the 5 spot-checked candidates was 2.6-45x higher than
   that flat estimate. This version uses the real per-pixel ERROR extension
   (`_error.fits.gz`, downloaded alongside DATA) integrated over each
   aperture instead.
2. Coverage-gap bug. The first pass measured F200W flux at each candidate's
   position without checking whether that position actually has real F200W
   coverage; a position landing on a masked/no-data gap silently returned
   near-zero flux, which then satisfied the "faint/red" cut for the wrong
   reason. This version requires a minimum valid-pixel coverage fraction in
   the aperture for both bands before trusting a measurement.

DAOStarFinder detection itself still uses the flat std as a *detection*
threshold (fine for generating a coarse candidate list -- only the
downstream *significance numbers* need the real noise model), but every
reported sigma value in the output CSV is computed from the real error map.
"""
import sys
if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import gzip
import csv
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry

DATA_DIR = "F:/Claude/JWST/nexus_data"
F444W_DATA = f"{DATA_DIR}/nexus_central_deep_ep01_f444w_60mas_i2d_data.fits.gz"
F444W_ERR = f"{DATA_DIR}/nexus_central_deep_ep01_f444w_60mas_i2d_error.fits.gz"
F200W_DATA = f"{DATA_DIR}/nexus_central_deep_ep01_f200w_i2d_data.fits.gz"
F200W_ERR = f"{DATA_DIR}/nexus_central_deep_ep01_f200w_i2d_error.fits.gz"
OUT_CSV = "nexus_ep01_lrdproxy_candidates_v2.csv"

RED_APERTURE_RADIUS_PX = 4.0
CONCENTRATION_OUTER_FACTOR = 2.5
CONCENTRATION_MIN = 0.65
CONCENTRATION_MAX = 1.0
RED_SIGMA_MIN = 10.0
COLOR_SIGMA_MAX = 2.0
MIN_COVERAGE_FRACTION = 0.9  # fraction of aperture pixels that must be real (finite, nonzero) data


def load(path):
    with gzip.open(path, 'rb') as f:
        hdul = fits.open(f)
        data = hdul[0].data.astype(np.float32)
        header = hdul[0].header
    return data, header


def aperture_flux_and_noise(data, err2, valid_mask_u8, wcs, coord, radius_px):
    """Returns flux, real-error-based noise, and coverage fraction for one or more coords.
    valid_mask_u8 must already be a uint8 view/array (not bool), to avoid allocating a
    full-mosaic-sized float32 copy just for the coverage sum -- that alone exceeded
    available memory on the ~1.1-gigapixel F200W array."""
    x, y = wcs.world_to_pixel(coord)
    x, y = np.atleast_1d(x), np.atleast_1d(y)
    n = len(x)
    ny, nx = data.shape
    in_bounds = (x > radius_px) & (x < nx - radius_px) & (y > radius_px) & (y < ny - radius_px)
    flux = np.full(n, np.nan)
    noise = np.full(n, np.nan)
    coverage = np.full(n, 0.0)
    if in_bounds.any():
        xy = np.column_stack([x[in_bounds], y[in_bounds]])
        aps = CircularAperture(xy, r=radius_px)
        flux[in_bounds] = aperture_photometry(data, aps)['aperture_sum']
        var_sum = aperture_photometry(err2, aps)['aperture_sum']
        noise[in_bounds] = np.sqrt(np.clip(var_sum, 0, None))
        cov_sum = aperture_photometry(valid_mask_u8, aps)['aperture_sum']
        coverage[in_bounds] = cov_sum / aps.area
    return flux, noise, coverage


print("Loading F444W data + error...")
data_red, hdr_red = load(F444W_DATA)
err_red_raw, _ = load(F444W_ERR)
wcs_red = WCS(hdr_red)
pixscale_red = proj_plane_pixel_scales(wcs_red)[0]
valid_red = np.isfinite(data_red) & (data_red != 0) & np.isfinite(err_red_raw) & (err_red_raw > 0)
mean_r, median_r, std_r = sigma_clipped_stats(data_red[valid_red], sigma=3.0)
print(f"  shape={data_red.shape} pixscale={pixscale_red*3600*1000:.1f}mas median={median_r:.4g} flat_std={std_r:.4g}")

err2_red = np.nan_to_num(err_red_raw, nan=0.0) ** 2
err2_red[~valid_red] = 0.0

print("Running DAOStarFinder on F444W (coarse candidate generation, flat threshold)...")
data_sub = data_red.copy()
data_sub -= np.float32(median_r)
data_sub[~valid_red] = 0.0
daofind = DAOStarFinder(fwhm=3.0, threshold=4.0 * std_r)
sources = daofind(data_sub)
print(f"  {len(sources) if sources is not None else 0} raw detections")

if sources is None or len(sources) == 0:
    print("No sources found, stopping.")
    sys.exit(0)

xs, ys = np.array(sources['xcentroid']), np.array(sources['ycentroid'])
coords = wcs_red.pixel_to_world(xs, ys)

valid_red_u8 = valid_red.view(np.uint8)
print("Measuring F444W real-error-based significance + concentration + coverage...")
flux_r, noise_r, cov_r = aperture_flux_and_noise(data_sub, err2_red, valid_red_u8, wcs_red, coords, RED_APERTURE_RADIUS_PX)
_, _, cov_r_outer = aperture_flux_and_noise(
    data_sub, err2_red, valid_red_u8, wcs_red, coords, RED_APERTURE_RADIUS_PX * CONCENTRATION_OUTER_FACTOR)
flux_r_outer, _, _ = aperture_flux_and_noise(
    data_sub, err2_red, valid_red_u8, wcs_red, coords, RED_APERTURE_RADIUS_PX * CONCENTRATION_OUTER_FACTOR)
concentration = np.where(flux_r_outer > 0, flux_r / flux_r_outer, np.nan)
red_sigma = np.where(noise_r > 0, flux_r / noise_r, np.nan)

prefilter = (
    (cov_r >= MIN_COVERAGE_FRACTION)
    & (cov_r_outer >= MIN_COVERAGE_FRACTION)
    & (red_sigma > RED_SIGMA_MIN)
    & (concentration >= CONCENTRATION_MIN)
    & (concentration <= CONCENTRATION_MAX)
)
n_pre = int(prefilter.sum())
print(f"  {n_pre} pass real-noise red_sigma>{RED_SIGMA_MIN}, concentration in "
      f"[{CONCENTRATION_MIN},{CONCENTRATION_MAX}], and >={MIN_COVERAGE_FRACTION*100:.0f}% real coverage")

if n_pre == 0:
    print("No candidates survive the red+compact prefilter, stopping.")
    sys.exit(0)

print("Loading F200W data + error for color check on prefiltered candidates only...")
data_blue, hdr_blue = load(F200W_DATA)
err_blue_raw, _ = load(F200W_ERR)
wcs_blue = WCS(hdr_blue)
pixscale_blue = proj_plane_pixel_scales(wcs_blue)[0]
valid_blue = np.isfinite(data_blue) & (data_blue != 0) & np.isfinite(err_blue_raw) & (err_blue_raw > 0)
mean_b, median_b, std_b = sigma_clipped_stats(data_blue[valid_blue], sigma=3.0)
print(f"  shape={data_blue.shape} pixscale={pixscale_blue*3600*1000:.1f}mas median={median_b:.4g} flat_std={std_b:.4g}")

err2_blue = np.nan_to_num(err_blue_raw, nan=0.0) ** 2
err2_blue[~valid_blue] = 0.0
data_blue_sub = data_blue.copy()
data_blue_sub -= np.float32(median_b)
data_blue_sub[~valid_blue] = 0.0

idx_pre = np.where(prefilter)[0]
coords_pre = coords[idx_pre]
blue_aperture_radius_px = RED_APERTURE_RADIUS_PX * (pixscale_red / pixscale_blue)
valid_blue_u8 = valid_blue.view(np.uint8)
flux_b, noise_b, cov_b = aperture_flux_and_noise(
    data_blue_sub, err2_blue, valid_blue_u8, wcs_blue, coords_pre, blue_aperture_radius_px)
blue_sigma = np.where(noise_b > 0, flux_b / noise_b, np.nan)

final_mask = (cov_b >= MIN_COVERAGE_FRACTION) & (blue_sigma < COLOR_SIGMA_MAX)
n_final = int(np.nansum(final_mask))
n_no_coverage = int(((cov_b < MIN_COVERAGE_FRACTION)).sum())
print(f"\n{n_final} candidates pass red+compact+blue-faint (real noise, real coverage) out of {n_pre} prefiltered")
print(f"  ({n_no_coverage} of {n_pre} excluded for insufficient F200W coverage -- the coverage-gap bug this version fixes)\n")

rows = []
for j, i in enumerate(idx_pre):
    if not final_mask[j]:
        continue
    c = coords[i]
    rows.append({
        'id': int(i),
        'ra': c.ra.deg, 'dec': c.dec.deg,
        'red_sigma': float(red_sigma[i]),
        'blue_sigma': float(blue_sigma[j]),
        'concentration': float(concentration[i]),
        'red_coverage': float(cov_r[i]),
        'blue_coverage': float(cov_b[j]),
        'x_f444w': float(xs[i]), 'y_f444w': float(ys[i]),
    })
    print(f"  #{i}: RA/Dec {c.ra.deg:.6f},{c.dec.deg:.6f} red={red_sigma[i]:.1f}sigma "
          f"blue={blue_sigma[j]:.1f}sigma conc={concentration[i]:.2f} "
          f"cov(red/blue)={cov_r[i]:.2f}/{cov_b[j]:.2f}")

with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
    if rows:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
print(f"\nWrote {OUT_CSV}")
