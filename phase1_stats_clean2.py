"""
Convergence Check: First vs Second Half of Plateau

CSV files must be named as prefix_rep{number}.csv
"""

import numpy as np

conditions = {
    '300K/25bar (P increase)': [
        '300K25bar_from300K1bar_rep1.csv',
        '300K25bar_from300K1bar_rep2.csv',
        '300K25bar_from300K1bar_rep3.csv',
    ],
    '373K/1bar (T increase)': [
        '373K1bar_from300K1bar_rep1.csv',
        '373K1bar_from300K1bar_rep2.csv',
        '373K1bar_from300K1bar_rep3.csv',
    ],
    '373K/1bar (P decrease)': [
        '373K1bar_from373K25bar_rep1.csv',
        '373K1bar_from373K25bar_rep2.csv',
        '373K1bar_from373K25bar_rep3.csv',
    ],
    '300K/25bar (T decrease)': [
        '300K25bar_from373K25bar_rep1.csv',
        '300K25bar_from373K25bar_rep2.csv',
        '300K25bar_from373K25bar_rep3.csv',
    ],
}

def loadCsv(path):
    data = np.loadtxt(path, delimiter=',', skiprows=1)
    return data[:, 0], data[:, 1]

def getPlateau(pct, frac=0.2):
    n = len(pct)
    return pct[int(n * frac):]

for condition, paths in conditions.items():
    for i, path in enumerate(paths):
        t, p = loadCsv(path)
        plateau = getPlateau(p)
        n = len(plateau)
        first = np.mean(plateau[:n // 2])
        second = np.mean(plateau[n // 2:])
        diff = abs(second - first)
        status = 'converged' if diff < 5 else 'check convergence'
        print(condition, "rep" + str(i + 1), "first =", round(first, 1),
              "second =", round(second, 1), "diff =", round(diff, 1), status)
