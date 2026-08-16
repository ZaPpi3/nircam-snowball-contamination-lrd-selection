# Snowball Contamination of Compact, Red Source Selections in JWST NIRCam Imaging: A Case Study Motivated by the "Little Red Dot" Population

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Official repository for the paper **"Snowball Contamination of Compact, Red Source Selections in JWST NIRCam Imaging: A Case Study Motivated by the Little Red Dot Population" (2026)** by Paul Jarvis. This work applies a crude proxy for JWST "Little Red Dot" (LRD) selection to NIRCam candidate catalogs already built by this project's dropout-galaxy search, then escalates from visual vetting to a real raw-pipeline check (the field-standard `snowblind` tool) to measure how much of the LRD-like signature is actually a detector artifact rather than a real compact source.

-

## Key Findings

- **The originally-suspected contaminant was not the dominant one.** This study started from a real NIRCam module A/B background offset already found and corrected in this project's own dropout search (`background-systematics-versus-redshift`). All candidates here were already drawn from the current-best, source-masked-background catalogs, so that systematic is controlled for by construction. The dominant contaminant found instead is a different, more mundane one: transient detector artifacts.
- **7 of 22 flagged candidates are confirmed NIRCam "snowball" cosmic-ray artifacts**, not visually suspected but confirmed with real raw-exposure `SATURATED` DQ flags, using the field-standard `snowblind` tool run through the actual `Detector1Pipeline`. The decisive test: a real bright source saturates in every exposure covering its position; a snowball hit is stochastic and does not. The single most extreme candidate in the sample (11238-sigma) is a clean, single-epoch example of exactly this signature.
- **Vetting status, not field identity, explains the result.** The one field whose proxy-flagged candidates came back clean (SMACS J0723.3-7327, 3 of 3) is also the only field that had already been through this project's full multi-round vetting before this study began. The other three fields (A2744, GOODS-N, JADES/GOODS-S) had not, and are roughly 30-45% clean depending on how the open items below are counted.
- **A field-recognized phenomenon, independently measured here.** NIRCam snowball artifacts and the exposure-redundancy mechanism behind them are already documented in the literature (COSMOS-Web's own reduction paper explicitly leaves the resulting false-positive rate unquantified). The 32% confirmed-artifact rate reported here is a genuine, if small-sample, contribution toward that open question, not a claim that any published professional LRD catalog carries a comparable rate.
- **A correction made honestly rather than absorbed silently.** An earlier internal draft of this result miscounted one candidate (GOODS-N #3230) as raw-pipeline-confirmed clean; checking directly against the actual pipeline log output found it was never run at all. The manuscript reports the corrected number (7 of 22 confirmed clean, not 8) and documents the correction in Sec. 4.3.
- **The method generalizes to a new survey.** Section 6 extends the same raw-DQ persistence test to the NEXUS North Ecliptic Pole survey, an independent JWST Treasury program with its own known low-dither-redundancy risk factor. A two-band proxy (NEXUS's public release does not yet include a usable blue band) applied directly to the Stage-3 mosaics, with two real bugs found and fixed along the way (a flat noise estimate that understated true per-pixel noise by 2.6-45x, and a coverage-gap bug that let masked regions register as spuriously faint), finds 2 of 18 (11%) candidates confirmed as transient single-epoch artifacts, 13 clean, 1 inconclusive, and 2 unresolved for lack of raw-exposure coverage.

**Bottom line:** an automated "compact plus smooth red continuum" cut over already-processed NIRCam mosaics is not, by itself, a trustworthy LRD-like selection, even after correcting for a real, previously-known background systematic. The gap is explained by vetting status, not by which side of the AGN-cocoon debate the wider LRD population eventually falls on.

## The Computational Pipeline

1. **LRD-like proxy selection:** concentration index in [0.65, 1.00] plus a `smooth_dusty_interloper_like` mid-band break classification, applied to each field's current-best candidate catalog. *(`code/cutout_lrdproxy_a2744_goodsn.py`)*
2. **Round 2 raw-pipeline verification:** downloads raw `_uncal.fits` exposures per flagged tile and runs the real `Detector1Pipeline` with `snowblind`'s jump-detection post-hook. *(`code/batch_tile.py`, `code/run_all.py`)*
3. **Per-candidate DQ-flag extraction:** resolves which raw exposures actually cover a candidate's sky position via a real distortion-corrected WCS, then extracts `SATURATED`/`JUMP_DET` DQ-flag footprints at the candidate position and a same-exposure control region. *(`code/check_extent.py`, `code/check_flags2.py`, `code/check_candidate.py`)*
4. **Summary figure:** reproduces the classification bar chart from the verified tallies in the manuscript. *(`code/make_summary_figure.py`)*
5. **NEXUS extension (Sec. 6):** two-band proxy selection with real per-pixel error maps and an aperture coverage-fraction check, applied to NEXUS Deep-tier ep01 Stage-3 mosaics. *(`code/nexus_lrdproxy_search_v2.py`)*, then the same raw-pipeline persistence check as step 2-3 above, applied to the 18 surviving candidates across both usable NEXUS episodes. *(`code/run_nexus.py`)*

## Repository Structure

- `main.tex` / `main.pdf` : Manuscript and LaTeX source (plain `article` class, standard packages only).
- `code/` : The scripts above.
- `figures/` : Candidate cutouts (`cutout_973.png`, the decisive artifact case study; `cutout_97.png` and `cutout_34.png`, confirmed-clean candidates) and the classification summary bar chart.
- `requirements.txt` : Pinned dependency versions (the raw-pipeline scripts require a working `jwst`/`snowblind`/CRDS environment; the summary figure script only needs `matplotlib`).

## Reproducing the results

```bash
pip install -r requirements.txt
```

`code/run_all.py` downloads raw JWST exposures from MAST and runs the full `Detector1Pipeline` per tile: expect this to dominate wall-clock time and require a working CRDS setup (`CRDS_PATH`/`CRDS_SERVER_URL`). `code/make_summary_figure.py` only needs the pinned tallies already in the script and `matplotlib`.

```bash
python code/run_all.py                    # raw-pipeline verification, all 22 non-SMACS candidates
python code/check_extent.py                # decisive #973 case study (Sec. 4.1)
python code/make_summary_figure.py          # reproduces figures/classification_summary.png
python code/nexus_lrdproxy_search_v2.py     # NEXUS two-band proxy selection (Sec. 6)
python code/run_nexus.py                    # NEXUS raw-pipeline persistence check (Sec. 6)
```

## License

MIT License. See `LICENSE`.

## Data Availability

JWST data are publicly available via MAST.
