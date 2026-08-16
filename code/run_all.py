import sys
sys.path.insert(0, '/home/jwst/snowblind_test')
from batch_tile import process_tile

TILES = [
    # (proposal_id, obs_key, [(cand_id, ra, dec), ...])
    ('2561', 'jw02561-o001_t003', [
        (188, 3.566277, -30.420185), (870, 3.519088, -30.371296),
        (3591, 3.632713, -30.399583), (3822, 3.635645, -30.396787),
        (4671, 3.634660, -30.376063),
    ]),
    ('2561', 'jw02561-o002_t001', [
        (97, 3.654100, -30.453474), (4731, 3.590348, -30.426203),
    ]),
    ('2561', 'jw02561-o003_t006', [
        (1753, 3.561357, -30.326939),
    ]),
    ('1181', 'jw01181-o098_t010', [
        (765, 189.356364, 62.257104), (1206, 189.344796, 62.240286),
    ]),
    ('1181', 'jw01181-o006_t008', [
        (450, 189.331302, 62.188114), (2779, 189.244822, 62.196313),
    ]),
    ('1181', 'jw01181-o001_t001', [
        (268, 189.197422, 62.299916), (1933, 189.204792, 62.245553),
    ]),
    ('1180', 'jw01180-o026_t028', [
        (812, 53.07945, -27.83363),
    ]),
    ('1180', 'jw01180-o136_t029', [
        (632, 53.14063, -27.88139),
    ]),
    ('1180', 'jw01180-o223_t223', [
        (679, 53.087251, -27.793417), (2300, 53.101792, -27.774275),
    ]),
    ('1180', 'jw01180-o010_t009', [
        (839, 53.20399, -27.77211),
    ]),
    ('1180', 'jw01180-o029_t028', [
        (34, 53.166833, -27.865206),
    ]),
]

all_results = {}
for proposal_id, obs_key, candidates in TILES:
    try:
        res = process_tile(proposal_id, obs_key, candidates)
        all_results[obs_key] = res
    except Exception as e:
        print(f'=== TILE {obs_key} FAILED: {e} ===')
        import traceback; traceback.print_exc()

print()
print('========== FINAL SUMMARY ==========')
for obs_key, res in all_results.items():
    for cand_id, hits in res.items():
        verdict = 'NO ON-DETECTOR EXPOSURE'
        if hits:
            max_sat = max(h[3] for h in hits)
            verdict = f'max_SAT={max_sat} across {len(hits)} exposure(s)'
        print(f'{obs_key} #{cand_id}: {verdict}')
