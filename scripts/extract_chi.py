#!/usr/bin/env python3
"""Re-extract chi(k) from Athena project with proper normalization using larch."""

import numpy as np
from larch.io import read_athena
from larch import Group
from larch.xafs import pre_edge, autobk

PRJ_PATH = '/home/nanare/Dokumente/Uni/Master_Wuppertal/ERP_Hausarbeit/data/athena_projects/MEA_athena.prj'
OUT_DIR = '/home/nanare/Dokumente/Uni/Master_Wuppertal/ERP_Hausarbeit/data/chi_extracted'

# Read Athena project
proj = read_athena(PRJ_PATH)

# List all groups in the project
print("Groups in Athena project:")
for name in dir(proj):
    g = getattr(proj, name)
    if isinstance(g, Group) and hasattr(g, 'energy'):
        edge_str = getattr(g, 'bkg_e0', 'unknown')
        print(f"  {name}: E0 = {edge_str}, {len(g.energy)} points, "
              f"E range = {g.energy[0]:.1f} - {g.energy[-1]:.1f} eV")

print("\n--- Processing each spectrum ---\n")

# Map group names to edge labels
edge_map = {}
for name in dir(proj):
    g = getattr(proj, name)
    if isinstance(g, Group) and hasattr(g, 'energy'):
        name_lower = name.lower()
        if 'co' in name_lower and 'pt' not in name_lower:
            edge_map['Co'] = name
        elif 'ni' in name_lower:
            edge_map['Ni'] = name
        elif 'cu' in name_lower:
            edge_map['Cu'] = name

print(f"Identified edges: {edge_map}\n")

for edge, group_name in edge_map.items():
    g = getattr(proj, group_name)

    # Run pre_edge normalization with default (robust) parameters
    pre_edge(g)

    print(f"=== {edge} K-edge ({group_name}) ===")
    print(f"  E0 = {g.e0:.2f} eV")
    print(f"  Edge step = {g.edge_step:.6f}")
    print(f"  Pre-edge range: {getattr(g, 'pre1', '?')} to {getattr(g, 'pre2', '?')}")
    print(f"  Norm range: {getattr(g, 'norm1', '?')} to {getattr(g, 'norm2', '?')}")

    # Run autobk for background removal → chi(k)
    autobk(g, rbkg=1.0)

    k = g.k
    chi = g.chi

    # Check amplitude
    mask = (k >= 3) & (k <= 10)
    chi_k2 = chi * k**2
    print(f"  k range: {k[0]:.2f} - {k[-1]:.2f} Å⁻¹")
    print(f"  max|χ(k)·k²| (k=3-10): {np.max(np.abs(chi_k2[mask])):.4f}")
    print(f"  Points: {len(k)}")

    # Save chi(k) — two columns: k, chi(k)
    outfile = f'{OUT_DIR}/MEA_{edge}_chi_new.dat'
    np.savetxt(outfile, np.column_stack([k, chi]),
               fmt='%.8e', header=f'{edge} K-edge chi(k) extracted with larch autobk')
    print(f"  Saved to: {outfile}\n")

# Also show old data for comparison
print("--- Comparison with old extraction ---")
for edge in ['Co', 'Ni', 'Cu']:
    old = np.loadtxt(f'{OUT_DIR}/MEA_{edge}_chi.dat')
    k_old = old[:, 0]
    chi_old = old[:, 1]
    mask_old = (k_old >= 3) & (k_old <= 10)
    chi_k2_old = chi_old * k_old[mask_old]**2 if len(chi_old) == len(k_old) else chi_old[mask_old]
    # The old files might already be chi*k^n
    print(f"  {edge} old: max|col2| (k=3-10) = {np.max(np.abs(old[mask_old, 1])):.4f}, "
          f"k range = {k_old[0]:.2f}-{k_old[-1]:.2f}")
