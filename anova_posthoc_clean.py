"""
One-Way ANOVA and Tukey HSD Post-Hoc Test
Tests whether the change in % self-associating stearic acid (delta)
differs significantly across the four MD conditions.

For each replicate:
    start = mean of first 5 frames
    end   = mean of plateau (last 80% of trajectory)
    delta = end - start

CSV files must be named as prefix_rep{number}.csv
Requires: pandas, numpy, scipy, statsmodels
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

conditions = {
    'Pressure increase': '300K25bar_from300K1bar',
    'Temperature increase': '373K1bar_from300K1bar',
    'Pressure decrease': '373K1bar_from373K25bar',
    'Temperature decrease': '300K25bar_from373K25bar',
}
replicates = [1, 2, 3]
plateauStartPercent = 0.20
nStartFrames = 5

def getDelta(path):
    df = pd.read_csv(path)
    pct = df['percent_self_associating'].values
    n = len(pct)
    startMean = np.mean(pct[:nStartFrames])
    plateau = pct[int(n * plateauStartPercent):]
    endMean = np.mean(plateau)
    return endMean - startMean

deltasByCondition = {}
for label, prefix in conditions.items():
    deltas = [getDelta(f'{prefix}_rep{rep}.csv') for rep in replicates]
    deltasByCondition[label] = deltas
    print(label, deltas, "mean =", round(np.mean(deltas), 2))

groups = list(deltasByCondition.values())
fStat, pVal = stats.f_oneway(*groups)
print("\nOne-way ANOVA: F =", round(fStat, 3), " p =", pVal)

allVals = np.concatenate(groups)
allLabels = np.concatenate([[label] * len(replicates) for label in deltasByCondition.keys()])
tukey = pairwise_tukeyhsd(allVals, allLabels, alpha=0.05)
print("\nTukey HSD post-hoc test:")
print(tukey)
