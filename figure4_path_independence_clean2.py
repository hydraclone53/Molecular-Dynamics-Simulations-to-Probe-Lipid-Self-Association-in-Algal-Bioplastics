"""
Path Independence Test 
Compares plateau % self-associating stearic acid between pairs of
simulations that reached the same final temperature and pressure via
different starting paths, using a two-sample t-test for each pair.

Plateau = frames from the 20% mark of the trajectory onward.

CSV files must be named as prefix_rep{number}.csv
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

blue = '#2980b9'  # started at 300K, 1bar
red = '#c0392b'   # started at 373K, 25bar

conditions = [
    # (label, prefix, color of starting configuration)
    ('300K/25bar', '300K25bar_from300K1bar', blue),    # P increase
    ('300K/25bar ', '300K25bar_from373K25bar', red),   # T decrease
    ('373K/1bar', '373K1bar_from300K1bar', blue),       # T increase
    ('373K/1bar ', '373K1bar_from373K25bar', red),      # P decrease
]
replicates = [1, 2, 3]
plateauStartPercent = 0.20

means, ses, labels, colors, allReps = [], [], [], [], []

for label, prefix, color in conditions:
    reps = []
    for rep in replicates:
        df = pd.read_csv(f'{prefix}_rep{rep}.csv')
        tmax = df['time_ns'].max()
        plateau = df[df['time_ns'] >= plateauStartPercent * tmax]['percent_self_associating']
        reps.append(plateau.mean())
    means.append(np.mean(reps))
    ses.append(np.std(reps, ddof=1) / np.sqrt(len(replicates)))
    labels.append(label.strip())
    colors.append(color)
    allReps.append(reps)

    print(label.strip(), "reps =", [round(r, 2) for r in reps],
          "mean =", round(np.mean(reps), 1), "SE =", round(np.std(reps, ddof=1) / np.sqrt(3), 1))

# Two-sample t-tests for the two path-independence comparisons
# bars 0,1 = P increase vs T decrease (both end at 300K/25bar)
# bars 2,3 = T increase vs P decrease (both end at 373K/1bar)
t1, p1 = stats.ttest_ind(allReps[0], allReps[1])
t2, p2 = stats.ttest_ind(allReps[2], allReps[3])
print("300K/25bar pair: t =", round(t1, 2), " p =", p1)
print("373K/1bar pair: t =", round(t2, 2), " p =", p2)

p1Text = "p < 0.0001" if p1 < 0.0001 else "p = " + str(round(p1, 2))
p2Text = "p < 0.0001" if p2 < 0.0001 else "p = " + str(round(p2, 2))

x = np.arange(len(labels))
fig, ax = plt.subplots(figsize=(9.6, 6.5))
ax.bar(x, means, yerr=ses, capsize=6, color=colors, edgecolor='black', linewidth=1.2)

for i, (m, s) in enumerate(zip(means, ses)):
    ax.text(i, m + s + 1, str(round(m, 1)) + "±" + str(round(s, 1)) + "%", ha='center', fontsize=11)

h1 = 72
ax.plot([0, 0, 1, 1], [h1 - 1.5, h1, h1, h1 - 1.5], color='black', linewidth=1.3)
ax.text(0.5, h1 + 1.5, p1Text, ha='center', fontsize=11)

h2 = 65
ax.plot([2, 2, 3, 3], [h2 - 1.5, h2, h2, h2 - 1.5], color='black', linewidth=1.3)
ax.text(2.5, h2 + 1.5, p2Text, ha='center', fontsize=11)

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=11)
ax.set_xlabel('Final Configurations', fontsize=11)
ax.set_ylabel('% Self-Associating Stearic Acid (plateau mean)', fontsize=11)
ax.set_title('Path Independence Test', fontsize=14, fontweight='bold')
ax.set_ylim(0, 85)

ax.legend(handles=[
    Patch(facecolor=red, edgecolor='black', label='373K, 25Bar'),
    Patch(facecolor=blue, edgecolor='black', label='300K, 1Bar'),
], title='Starting Configurations', fontsize=10, title_fontsize=10, loc='upper right')

plt.tight_layout()
plt.savefig('pathIndependenceTest.png', dpi=150)
print("Saved: pathIndependenceTest.png")
