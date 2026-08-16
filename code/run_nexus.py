"""
Raw-exposure DQ-flag persistence check for the 18 NEXUS Deep-tier ep01
candidates from nexus_lrdproxy_search_v2.py (Sec. 6 of the manuscript).

Downloads raw _uncal.fits exposures for the two NEXUS visits with usable
aggregate MAST Level-3 records (o004_t001 = episode 1, o002_t001 = episode 2),
runs Detector1Pipeline + snowblind's SnowblindStep via the same process_tile()
function used for the paper's original four fields, and reports SATURATED/
JUMP_DET DQ flags at each candidate position across all covering exposures
in both episodes.

Requires batch_tile.py (same directory) and a working jwst/snowblind/CRDS
environment; see the top-level README for setup.
"""
import sys
sys.path.insert(0, '.')
from batch_tile import process_tile

# All 18 candidates surviving the real-error, coverage-checked selection in
# nexus_lrdproxy_search_v2.py (id, ra_deg, dec_deg).
candidates = [
    ('3866', 268.3458089730663, 65.15961377761596),
    ('5516', 268.57960436217235, 65.16676644070303),
    ('5605', 268.56313114916065, 65.16714925383185),
    ('6025', 268.5445129953018, 65.16903572139097),
    ('6643', 268.665763, 65.171608),
    ('7996', 268.669969, 65.176732),
    ('9769', 268.305969, 65.183169),
    ('11099', 268.39450240629253, 65.18832747803495),
    ('12003', 268.384905111133, 65.19167788732449),
    ('12608', 268.453401, 65.194034),
    ('14009', 268.2839291824869, 65.19941189320964),
    ('14786', 268.37957585670415, 65.20233221729667),
    ('16564', 268.671507, 65.209519),
    ('20911', 268.59883466528913, 65.2266642163724),
    ('21241', 268.303860, 65.228111),
    ('21519', 268.6090342927112, 65.22935042381854),
    ('24489', 268.5174566046849, 65.24303650929463),
    ('24675', 268.36693726612464, 65.24406629965402),
]

print('=== NEXUS episode 1 (visit o004_t001) ===')
process_tile('9263', 'jw09263-o004_t001', candidates)
print('=== NEXUS episode 2 (visit o002_t001) ===')
process_tile('9263', 'jw09263-o002_t001', candidates)
print('DONE')
