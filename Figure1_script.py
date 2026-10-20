"""
Change in STA Self-Association by Condition
============================================
Calculates the change (Delta) in % self-associating molecules from
start to end of each run, with the standard error across replicates.

For each run:
    start = first frame
    end   = mean of plateau (last 80% of trajectory, i.e. 100-500 ns)
    delta = end - start

For each condition:
    mean Delta = mean of the 3 replicate deltas
    SE         = std(3 replicate deltas, ddof=1) / sqrt(3)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

def load_csv(path):
    data = np.loadtxt(path, delimiter=',', skiprows=1)
    return data[:, 0], data[:, 1]

def get_delta(path):
    time, selfAssociation = load_csv(path)
    n = len(selfAssociation)

    startValue = selfAssociation[0]

    plateau = selfAssociation[int(n * 0.2):]
    endMean = np.mean(plateau)

    return endMean - startValue


def analyze_and_plot(conditions, colors, xtick_labels, out_dir='delta_output'):
    os.makedirs(out_dir, exist_ok=True)

    labels = list(conditions.keys())
    all_deltas = []

    print("Condition", "Rep", "Delta%")

    for condition, paths in conditions.items():
        rep_deltas = []

        for i, path in enumerate(paths):
            delta = get_delta(path)
            rep_deltas.append(delta)
            print(condition, "rep" + str(i+1), round(delta, 1))

        all_deltas.append(rep_deltas)

    means_delta = [np.mean(d) for d in all_deltas]
    ses = [np.std(d, ddof=1) / np.sqrt(len(d)) for d in all_deltas]

    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.bar(x, means_delta, yerr=ses, capsize=10, color=colors,
           edgecolor='black', linewidth=1.2, width=0.55,
           error_kw=dict(elinewidth=1.5, capthick=1.5))

    for i, (m, se) in enumerate(zip(means_delta, ses)):
        ax.text(x[i], max(0, m) + se + 0.8, "Δ=" + f"{m:+.1f}" + "%\n±" + str(round(se, 1)) + "%",
            ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax.axhline(0, color='black', linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels, fontsize=11)
    ax.set_xlabel('Final Configurations', fontsize=12)
    ax.set_ylabel('Δ% Self-Associating Stearic Acid Molecules', fontsize=12)
    ax.set_title('Change in Stearic Acid Self-Association by Condition',
                 fontsize=15, fontweight='bold')
    ax.set_ylim(-8, 42)

    ax.legend(handles=[
        Patch(facecolor=red, edgecolor='black', label='373K, 25Bar'),
        Patch(facecolor=blue, edgecolor='black', label='300K, 1Bar'),
    ], title='Starting Configurations', fontsize=11, title_fontsize=11, loc='upper right')

    plt.tight_layout()
    path = os.path.join(out_dir, 'delta_selfassoc.png')
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"\nSaved: {path}")

    return means_delta, ses


conditions = {
    'P increase (300K/1bar→300K/25bar)': [
        '300K25bar_from300K1bar_rep1.csv',
        '300K25bar_from300K1bar_rep2.csv',
        '300K25bar_from300K1bar_rep3.csv',
    ],
    'T increase (300K/1bar→373K/1bar)': [
        '373K1bar_from300K1bar_rep1.csv',
        '373K1bar_from300K1bar_rep2.csv',
        '373K1bar_from300K1bar_rep3.csv',
    ],
    'P decrease (373K/25bar→373K/1bar)': [
        '373K1bar_from373K25bar_rep1.csv',
        '373K1bar_from373K25bar_rep2.csv',
        '373K1bar_from373K25bar_rep3.csv',
    ],
    'T decrease (373K/25bar→300K/25bar)': [
        '300K25bar_from373K25bar_rep1.csv',
        '300K25bar_from373K25bar_rep2.csv',
        '300K25bar_from373K25bar_rep3.csv',
    ],
}

blue = '#2980b9'
red = '#c0392b'
colors = [blue, blue, red, red]
xtickLabels = ['300K/25bar', '373K/1bar', '373K/1bar', '300K/25bar']  # final configurations

if __name__ == '__main__':
    analyze_and_plot(conditions, colors, xtickLabels, out_dir='delta_output')
