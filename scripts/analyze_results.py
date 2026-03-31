#!/usr/bin/env python3
"""Analyse der EvAX-Ergebnisse für MEA CoNiCu."""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.size': 11, 'figure.dpi': 150})

base = Path('/home/nanare/Dokumente/Uni/Master_Wuppertal/ERP_Hausarbeit/MEA_analysis/output_files')
edges = ['Co', 'Ni', 'Cu']
colors = {'Co': '#1f77b4', 'Ni': '#2ca02c', 'Cu': '#d62728'}

# ============================================================
# 1) Convergence plot
# ============================================================
output = np.loadtxt(base / 'output.dat', skiprows=1)
iterations = output[:, 0]
best_residual = output[:, 1]

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(iterations, best_residual, 'k-', lw=2)
ax.set_xlabel('Iteration')
ax.set_ylabel('Best Residual (ξ)')
ax.set_title('EvAX Konvergenz – MEA CoNiCu (3-Kanten simultan)')
ax.grid(True, alpha=0.3)
ax.annotate(f'Final: {best_residual[-1]:.4f}', xy=(iterations[-1], best_residual[-1]),
            xytext=(-80, 20), textcoords='offset points',
            arrowprops=dict(arrowstyle='->', color='red'), color='red', fontweight='bold')
fig.tight_layout()
fig.savefig(base / 'convergence.png')
print(f"Konvergenz: {best_residual[0]:.4f} → {best_residual[-1]:.4f} ({len(iterations)} Iterationen)")

# ============================================================
# 2) EXAFS comparison (exp vs calc) for all edges
# ============================================================
fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
for i, edge in enumerate(edges):
    calc = np.loadtxt(base / edge / 'EXAFS.dat')
    exp_ipol = np.loadtxt(base / edge / 'ipolEXAFS.dat')

    k_calc, chi_calc = calc[:, 0], calc[:, 1]
    k_exp, chi_exp = exp_ipol[:, 0], exp_ipol[:, 1]

    axes[i].plot(k_exp, chi_exp, 'k-', lw=1.5, label='Experiment')
    axes[i].plot(k_calc, chi_calc, '-', color=colors[edge], lw=1.5, label='EvAX Fit', alpha=0.8)
    axes[i].set_ylabel(f'{edge} K-edge\nχ(k)·k² (Å⁻²)')
    axes[i].legend(loc='upper right')
    axes[i].grid(True, alpha=0.3)

    # Residual
    if len(k_calc) == len(k_exp):
        diff = chi_exp - chi_calc
        r_factor = np.sum(diff**2) / np.sum(chi_exp**2)
        axes[i].set_title(f'{edge} K-edge  (R-factor: {r_factor:.4f})')
    else:
        axes[i].set_title(f'{edge} K-edge')

axes[-1].set_xlabel('k (Å⁻¹)')
fig.suptitle('EXAFS Vergleich: Experiment vs. EvAX Fit', fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(base / 'exafs_comparison.png', bbox_inches='tight')
print("EXAFS Vergleich gespeichert.")

# ============================================================
# 3) Fourier Transform comparison
# ============================================================
fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
for i, edge in enumerate(edges):
    ft_calc = np.loadtxt(base / edge / 'ft.dat')
    ft_exp = np.loadtxt(base / edge / 'expFT.dat')

    r_calc, mag_calc = ft_calc[:, 0], ft_calc[:, 3]
    r_exp, mag_exp = ft_exp[:, 0], ft_exp[:, 3]

    axes[i].plot(r_exp, mag_exp, 'k-', lw=1.5, label='Experiment')
    axes[i].plot(r_calc, mag_calc, '-', color=colors[edge], lw=1.5, label='EvAX Fit', alpha=0.8)
    axes[i].set_ylabel(f'{edge} K-edge\n|FT(χ(k)·k²)|')
    axes[i].legend(loc='upper right')
    axes[i].grid(True, alpha=0.3)
    axes[i].set_title(f'{edge} K-edge – Fourier Transform')
    axes[i].set_xlim(0, 6)

axes[-1].set_xlabel('R (Å)')
fig.suptitle('FT Vergleich: Experiment vs. EvAX Fit', fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(base / 'ft_comparison.png', bbox_inches='tight')
print("FT Vergleich gespeichert.")

# ============================================================
# 4) Partial RDF (pair distribution functions)
# ============================================================
fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
pair_labels = {
    'Co': ['Co-Cu', 'Co-Co', 'Co-Ni'],
    'Ni': ['Ni-Cu', 'Ni-Co', 'Ni-Ni'],
    'Cu': ['Cu-Cu', 'Cu-Co', 'Cu-Ni'],
}
pair_colors = ['#d62728', '#1f77b4', '#2ca02c']

for i, edge in enumerate(edges):
    rdf_data = np.loadtxt(base / edge / 'rdf.dat', skiprows=1)
    r = rdf_data[:, 0]

    for j, label in enumerate(pair_labels[edge]):
        axes[i].plot(r, rdf_data[:, j+1], '-', color=pair_colors[j], lw=1.2, label=label)

    axes[i].set_ylabel(f'{edge}-Kante\ng(R)')
    axes[i].legend(loc='upper right')
    axes[i].grid(True, alpha=0.3)
    axes[i].set_title(f'Partielle RDF – {edge} Absorber')
    axes[i].set_xlim(1.5, 6.0)

axes[-1].set_xlabel('R (Å)')
fig.suptitle('Partielle Paarverteilungsfunktionen (RDF)', fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(base / 'rdf_partial.png', bbox_inches='tight')
print("RDF gespeichert.")

# ============================================================
# 5) Koordinationszahlen-Analyse aus RDF
# ============================================================
print("\n" + "="*60)
print("KOORDINATIONSZAHLEN-ANALYSE (1. Schale, R = 2.2-2.8 Å)")
print("="*60)

# Ideal fcc NN distance: a/sqrt(2) = 3.54/1.414 = 2.503 Å
# Ideal coordination: 12 NN, each type ~4 (for equiatomic ternary)
r_min, r_max = 2.2, 2.8
dr = 0.01  # bin width

for edge in edges:
    rdf_data = np.loadtxt(base / edge / 'rdf.dat', skiprows=1)
    r = rdf_data[:, 0]
    mask = (r >= r_min) & (r <= r_max)

    labels = pair_labels[edge]
    print(f"\n{edge}-Kante:")
    total_cn = 0
    for j, label in enumerate(labels):
        # RDF values are already density-like; integrate
        # The RDF from EvAX gives g(R) in histogram form
        # CN ≈ sum(g(R)*dR) but depends on normalization
        cn = np.sum(rdf_data[mask, j+1]) * dr
        total_cn += cn
        print(f"  {label}: CN ≈ {cn:.2f}")
    print(f"  Total 1st shell: {total_cn:.2f} (ideal fcc: 12)")

# ============================================================
# 6) Chemical Short-Range Order (Warren-Cowley α)
# ============================================================
print("\n" + "="*60)
print("CHEMISCHE NAHORDNUNG (Warren-Cowley Parameter α)")
print("="*60)
print("α = 1 - P(ij)/c(j)  |  α=0: zufällig, α<0: Präferenz, α>0: Vermeidung")

# Composition: equiatomic ~1/3 each
c = {'Co': 85/256, 'Ni': 85/256, 'Cu': 86/256}

# Collect partial CN from all edges
cn_matrix = {}
for edge in edges:
    rdf_data = np.loadtxt(base / edge / 'rdf.dat', skiprows=1)
    r = rdf_data[:, 0]
    mask = (r >= r_min) & (r <= r_max)

    labels = pair_labels[edge]
    cns = []
    for j in range(3):
        cn = np.sum(rdf_data[mask, j+1]) * dr
        cns.append(cn)
    cn_matrix[edge] = dict(zip(labels, cns))

# Calculate Warren-Cowley parameters
print(f"\nKompositionen: Co={c['Co']:.3f}, Ni={c['Ni']:.3f}, Cu={c['Cu']:.3f}")
print()

for absorber in edges:
    total = sum(cn_matrix[absorber].values())
    if total < 0.1:
        continue
    print(f"{absorber}-Absorber (Total CN = {total:.1f}):")
    for pair, cn in cn_matrix[absorber].items():
        neighbor = pair.split('-')[1]
        p_ij = cn / total if total > 0 else 0
        alpha = 1 - p_ij / c[neighbor]
        print(f"  {pair}: CN={cn:.2f}, P={p_ij:.3f}, α = {alpha:+.3f}")
    print()

# ============================================================
# 7) Final structure composition check
# ============================================================
print("="*60)
print("FINALE STRUKTUR")
print("="*60)
with open(base / 'final.xyz') as f:
    natom = int(f.readline().strip())
    comment = f.readline().strip()
    elements = []
    for line in f:
        parts = line.split()
        if len(parts) >= 4:
            elements.append(parts[0])
            if len(elements) >= 256:
                break

from collections import Counter
counts = Counter(elements[:256])
print(f"Atome (real): {sum(counts.values())}")
for el, n in sorted(counts.items()):
    print(f"  {el}: {n} ({n/256*100:.1f}%)")

print(f"\nPBC-Kopien: {natom - 256}")
print(f"Gesamt im XYZ: {natom}")

plt.show()
print("\nAnalyse abgeschlossen!")
