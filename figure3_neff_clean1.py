"""
Effective Sample Size by Condition
===============================================
Calculates the effective sample size (N_eff) for each condition following
the method of Grossfield et al., and plots all four conditions as a
bar chart (mean +/- SE across 3 replicates).

Plateau region = frames from the 20% mark of the trajectory onward
i.e. the last 80%, after discarding the first 20% as equilibration.

CSV files must be named as prefix_rep{number}.csv
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


conditions = [
    #(label, prefix, group)
    ('300K/25bar', '300K25bar_from300K1bar', 'start_300K_1bar'),    
    ('373K/1bar', '373K1bar_from300K1bar', 'start_300K_1bar'),      
    ('373K/1bar ', '373K1bar_from373K25bar', 'start_373K_25bar'),   
    ('300K/25bar ', '300K25bar_from373K25bar', 'start_373K_25bar'), 
]
replicates = [1, 2, 3]
platueStartPercent = 0.20  # plateau begins at the 20% mark, runs to the end

 
"""N_eff via the full autocorrelation sum (Grossfield et al.):
    g = 1 + 2*sum(C(k)) until C(k) first goes negative; N_eff = N/g."""
def get_neff(values):
   
    fluctuations = values - values.mean()
    variance = np.mean(fluctuations ** 2)
    if variance == 0:
        return 1
    max_lag = min(200, len(values) // 4)
    acf_sum = 0
    for k in range(1, max_lag + 1):
        c_k = np.mean(fluctuations[:len(values) - k] * fluctuations[k:]) / variance
        if c_k < 0:
            break
        acf_sum += c_k
    g = 1 + 2 * acf_sum
    return len(values) / g


def main():
    labels, means, ses, colors = [], [], [], []

    for label, prefix, group in conditions:
        neffs = []
        for rep in replicates:
            df = pd.read_csv(f'{prefix}_rep{rep}.csv')
            selfAssociation = df['percent_self_associating'].values
            n = len(selfAssociation)
            plateau = selfAssociation[int(n * platueStartPercent):]
            neffs.append(get_neff(plateau))

        labels.append(label.strip())
        means.append(np.mean(neffs))
        ses.append(np.std(neffs, ddof=1) / np.sqrt(len(replicates)))
        colors.append('#c0392b' if group == 'start_373K_25bar' else '#2980b9')

        print(label.strip(), "N_eff =", [round(x, 1) for x in neffs],
      "mean =", round(np.mean(neffs)), "SE =", round(np.std(neffs, ddof=1) / np.sqrt(3)))

    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.bar(x, means, yerr=ses, capsize=8, color=colors, edgecolor='black',
           linewidth=1.2, width=0.55,
           error_kw=dict(elinewidth=1.5, capthick=1.5))

    for i, (m, s) in enumerate(zip(means, ses)):
        ax.text(x[i], m + s + 5, str(round(m)) + " ± " + str(round(s)),
        ha='center', fontsize=11, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_xlabel('Final Configurations', fontsize=12)
    ax.set_ylabel('Effective Independent Frames ($N_{\\mathrm{eff}}$)', fontsize=12)
    ax.set_title('Effective Sample Size by Condition', fontsize=15, fontweight='bold')
    ax.set_ylim(0, 250)

    ax.legend(handles=[
        Patch(facecolor='#c0392b', edgecolor='black', label='373K, 25Bar'),
        Patch(facecolor='#2980b9', edgecolor='black', label='300K, 1Bar'),
    ], title='Starting Configurations', fontsize=11, title_fontsize=11, loc='upper right')

    plt.tight_layout()
    plt.savefig('effectSampleSize.png', dpi=150)
    print("\nSaved: effectSampleSize.png")


if __name__ == '__main__':
    main()
