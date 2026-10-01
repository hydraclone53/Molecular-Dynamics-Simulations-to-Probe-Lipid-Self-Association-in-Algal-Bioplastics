"""
Change in STA Self-Association with N_eff-Corrected SE
========================================================
Calculates the change (Delta) in % self-associating molecules from
start to end of each run, with N_eff-corrected standard error.

For each run:
    start = mean of first 5 frames
    end   = mean of plateau (last 80% of trajectory)
    delta = end - start
    
    SE_end   = std(plateau) / sqrt(N_eff)
    SE_start = std(first 5 frames) / sqrt(5)
    SE_delta = sqrt(SE_end^2 + SE_start^2)
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

def get_neff(values):
    fluctuations = values - values.mean()
    variance = np.mean(fluctuations ** 2)
    if variance == 0:
        return 1
    max_lag = min(200, len(values) // 4)
    acf_sum = 0
    for k in range(1, max_lag + 1):
        c = np.mean(fluctuations[:len(values) - k] * fluctuations[k:]) / variance
        if c < 0:
            break
        acf_sum += c
    g = 1 + 2 * acf_sum
    return len(values) / g

def get_delta_and_se(path):
    time, selfAssociation = load_csv(path)
    n = len(selfAssociation)

    startFrames = selfAssociation[:5]
    startMean = np.mean(startFrames)
    SE_start = np.std(startFrames, ddof=1) / np.sqrt(5)

    plateau = selfAssociation[int(n * 0.2):]
    endMean = np.mean(plateau)
    neff = get_neff(plateau)
    SE_end = np.std(plateau, ddof=1) / np.sqrt(neff)

    delta = endMean - startMean
    SE_delta = np.sqrt(SE_end ** 2 + SE_start ** 2)

    return delta, SE_delta, neff


def analyze_and_plot(conditions, colors, xtick_labels, out_dir='delta_neff_output'):
    os.makedirs(out_dir, exist_ok=True)

    labels = list(conditions.keys())
    all_deltas = []
    all_ses = []
    all_neffs = []

    print("Condition", "Rep", "Delta%", "SE", "N_eff")

    for condition, paths in conditions.items():
        rep_deltas, rep_ses, rep_neffs = [], [], []

        for i, path in enumerate(paths):
            delta, se_delta, neff = get_delta_and_se(path)
            rep_deltas.append(delta)
            rep_ses.append(se_delta)
            rep_neffs.append(neff)
            print(condition, "rep" + str(i+1), round(delta, 1), round(se_delta, 3), round(neff))

        all_deltas.append(rep_deltas)
        all_ses.append(rep_ses)
        all_neffs.append(rep_neffs)

    means_delta = [np.mean(d) for d in all_deltas]
    ses = [np.mean(s) for s in all_ses]
    neffs_mean = [np.mean(n) for n in all_neffs]

    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.bar(x, means_delta, yerr=ses, capsize=10, color=colors,
           edgecolor='black', linewidth=1.2, width=0.55,
           error_kw=dict(elinewidth=1.5, capthick=1.5))

    for i, (m, se, nf) in enumerate(zip(means_delta, ses, neffs_mean)):
        ax.text(x[i], max(0, m) + se + 0.8, "Δ=" + f"{m:+.1f}" + "%\n±" + str(round(se, 1)) + "%",
            ha='center', va='bottom', fontsize=11, fontweight='bold')
        ax.text(x[i], -2.5, "$N_{eff}$≈" + str(round(nf)),
            ha='center', va='top', fontsize=9, color='dimgray')

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
    path = os.path.join(out_dir, 'delta_selfassoc_neff_corrected.png')
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"\nSaved: {path}")

    return means_delta, ses, neffs_mean


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
    analyze_and_plot(conditions, colors, xtickLabels, out_dir='delta_neff_output')
