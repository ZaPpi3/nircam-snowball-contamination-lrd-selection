import matplotlib.pyplot as plt

# Real classification tallies for the 22 LRD-like-proxy candidates, from the
# raw-pipeline verification described in main.tex Sec. 4 (Detector1 + snowblind,
# persistence criterion applied to the actual DQ flags). Source: batch_run3.log's
# FINAL SUMMARY block plus the separate #973 decisive check (check_extent.py).
# #3230 is counted under "visual-pass only" per the 2026-08-13 correction: it was
# never actually run through the raw pipeline (not in run_all.py's candidate list,
# no entry in FINAL SUMMARY), so it cannot be counted as raw-pipeline-confirmed.
categories = [
    "Confirmed artifact\n(partial-exposure\nsaturation)",
    "Plausibly real\n(modest/blended\nsaturation)",
    "Confirmed clean\n(raw pipeline)",
    "Visual pass only\n(#3230, not\nraw-pipeline-tested)",
    "Inconclusive\n(off-detector)",
]
counts = [7, 3, 7, 1, 4]

fig, ax = plt.subplots(figsize=(7, 4.2))
bars = ax.bar(categories, counts, color="0.35", edgecolor="black", linewidth=0.8)
for b, c in zip(bars, counts):
    ax.text(b.get_x() + b.get_width() / 2, c + 0.15, str(c), ha="center", va="bottom", fontsize=10)

ax.set_ylabel("Candidates (of 22 total)")
ax.set_title("Raw-pipeline verification outcome for 22 LRD-like-proxy candidates")
ax.set_ylim(0, 9)
plt.xticks(fontsize=8.5)
fig.tight_layout()
fig.savefig("../figures/classification_summary.png", dpi=300, bbox_inches="tight")
print("Saved classification_summary.png")
