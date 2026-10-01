"""
STA COM-Based Self-Association Analysis
Computes percent self-associating stearic acid using center-of-mass
distances between molecules (6 Angstrom COM cutoff).

Usage:
    python sta_com_analysis.py --tpr md_ramp.tpr --xtc md_ramp_center.xtc --out_dir results
"""

import numpy as np
import argparse
import os
from scipy.spatial.distance import cdist

try:
    import MDAnalysis as mda
except ImportError:
    print("MDAnalysis not found. Install with: pip install MDAnalysis --user")
    exit(1)

cutoffAng = 6.0
resname = 'STA'

def computePercentSelfassoc(tprPath, xtcPath, outDir, stride=1):
    os.makedirs(outDir, exist_ok=True)

    print("Loading universe:", tprPath, xtcPath)
    u = mda.Universe(tprPath, xtcPath)

    sta = u.select_atoms('resname ' + resname)
    if len(sta) == 0:
        print("ERROR: no atoms found with resname", resname)
        print("Available residue names:", set(u.residues.resnames))
        exit(1)

    nMolecules = sta.n_residues

    times = []
    pctSelfassoc = []

    for i, ts in enumerate(u.trajectory[::stride]):
        if i % 100 == 0:
            print("Frame", i, "t =", round(ts.time / 1000, 1), "ns")

        coms = np.array([res.atoms.center_of_mass() for res in sta.residues])
        distMatrix = cdist(coms, coms)
        np.fill_diagonal(distMatrix, np.inf)
        hasNeighbor = np.any(distMatrix <= cutoffAng, axis=1)

        pct = (np.sum(hasNeighbor) / nMolecules) * 100
        times.append(ts.time / 1000.0)
        pctSelfassoc.append(pct)

    times = np.array(times)
    pctSelfassoc = np.array(pctSelfassoc)

    csvPath = os.path.join(outDir, 'sta_percent_selfassoc.csv')
    np.savetxt(csvPath, np.column_stack([times, pctSelfassoc]),
               header='time_ns,percent_self_associating', delimiter=',', comments='')
    print("Saved:", csvPath)

    return times, pctSelfassoc


def main():
    parser = argparse.ArgumentParser(description='STA COM self-association analysis')
    parser.add_argument('--tpr', required=True, help='Path to .tpr file')
    parser.add_argument('--xtc', required=True, help='Path to .xtc file')
    parser.add_argument('--out_dir', default='results', help='Output directory')
    parser.add_argument('--stride', type=int, default=1, help='Frame stride')
    args = parser.parse_args()

    computePercentSelfassoc(args.tpr, args.xtc, args.out_dir, stride=args.stride)
    print("Done.")

if __name__ == '__main__':
    main()
