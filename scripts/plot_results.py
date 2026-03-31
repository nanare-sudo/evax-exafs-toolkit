#!/usr/bin/env python3
"""Publication-quality plots for MEA CoNiCu EvAX analysis."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# ============================================================
# Publication style setup
# ============================================================
mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times'],
    'mathtext.fontset': 'dejavuserif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'axes.linewidth': 0.8,
    'lines.linewidth': 1.2,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

BASE = Path('/home/nanare/Dokumente/Uni/Master_Wuppertal/ERP_Hausarbeit')
OUT = BASE / 'MEA_analysis' / 'output_files'
FIG = BASE / 'figures'
EDGES = ['Co', 'Ni', 'Cu']

# Consistent color palette (colorblind-friendly)
C_EXP = '#000000'          # black for experiment
C_FIT = {'Co': '#0072B2',  # blue
         'Ni': '#009E73',  # green
         'Cu': '#D55E00'}  # vermillion
C_PAIR = ['#0072B2', '#D55E00', '#009E73']  # for partial RDFs

# S02 correction: with properly normalized chi(k) data (Run 04+), no correction needed
AUTO_S02 = {'Co': 1.0, 'Ni': 1.0, 'Cu': 1.0}


# ============================================================
# Figure 1: Convergence
# ============================================================
def plot_convergence():
    data = np.loadtxt(OUT / 'output.dat', skiprows=1)
    it, res = data[:, 0], data[:, 1]

    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    ax.plot(it, res, 'k-', lw=1.0)
    ax.set_xlabel('Iteration')
    ax.set_ylabel(r'Best residual $\xi$')
    ax.set_xlim(0, it[-1])
    ax.set_ylim(bottom=0.4)

    # Mark start and end
    ax.annotate(f'{res[0]:.3f}', xy=(it[0], res[0]),
                xytext=(40, -5), textcoords='offset points',
                fontsize=8, color='0.4')
    ax.annotate(f'{res[-1]:.3f}', xy=(it[-1], res[-1]),
                xytext=(-45, 10), textcoords='offset points',
                fontsize=8, color='0.4',
                arrowprops=dict(arrowstyle='->', color='0.4', lw=0.6))

    fig.savefig(FIG / 'fig1_convergence.pdf')
    fig.savefig(FIG / 'fig1_convergence.png')
    plt.close(fig)
    print('Fig 1: Convergence saved')


# ============================================================
# Figure 2: EXAFS chi(k)*k^2 comparison (3 panels)
# ============================================================
def plot_exafs():
    fig, axes = plt.subplots(3, 1, figsize=(3.5, 7.5), sharex=True)
    fig.subplots_adjust(hspace=0.08)

    labels = ['(a)', '(b)', '(c)']

    for i, edge in enumerate(EDGES):
        ax = axes[i]
        calc = np.loadtxt(OUT / edge / 'EXAFS.dat')
        exp = np.loadtxt(OUT / edge / 'ipolEXAFS.dat')

        exp_corr = exp[:, 1] * AUTO_S02[edge]
        ax.plot(exp[:, 0], exp_corr, color=C_EXP, lw=1.0, label='Experiment')
        ax.plot(calc[:, 0], calc[:, 1], color=C_FIT[edge], lw=1.0,
                label='EvAX fit', alpha=0.85)

        # Residual below with clear offset
        if len(calc) == len(exp):
            diff = exp_corr - calc[:, 1]
            res_offset = -0.55
            ax.plot(exp[:, 0], diff + res_offset, color='0.5', lw=0.6,
                    label='Residual')
            ax.axhline(res_offset, color='0.8', lw=0.3, ls='-')

        ymax = max(np.max(np.abs(exp_corr)), np.max(np.abs(calc[:, 1])))
        ax.set_ylim(-0.75, ymax * 1.15)
        ax.set_ylabel(r'$\chi(k) \cdot k^2$ ($\mathrm{\AA}^{-2}$)')
        ax.set_xlim(0, 14)

        # Grey shading for excluded k-range
        ylims = ax.get_ylim()
        ax.axvspan(0, 3, alpha=0.08, color='grey', zorder=0)
        ax.axvspan(12, 14, alpha=0.08, color='grey', zorder=0)
        if i == 0:
            ax.text(1.5, ylims[1] * 0.85, 'excluded', ha='center',
                    fontsize=7, color='0.5', style='italic')

        ax.text(0.02, 0.95, f'{labels[i]} {edge} K-edge',
                transform=ax.transAxes, fontsize=10, fontweight='bold',
                va='top')

        if i == 0:
            ax.legend(loc='upper right', frameon=False, handlelength=1.5)

    axes[-1].set_xlabel(r'$k$ ($\mathrm{\AA}^{-1}$)')

    fig.savefig(FIG / 'fig2_exafs_comparison.pdf')
    fig.savefig(FIG / 'fig2_exafs_comparison.png')
    plt.close(fig)
    print('Fig 2: EXAFS comparison saved')


# ============================================================
# Figure 3: Fourier Transform magnitude (3 panels)
# FT computed from chi(k)*k^2 in fit range with Hanning window
# ============================================================
def _compute_ft(k, chi_kw, k_min=3.0, k_max=12.0, r_max=6.0, n_r=600):
    """Compute FT magnitude from chi(k)*k^w using Hanning window."""
    mask = (k >= k_min) & (k <= k_max)
    k_w = k[mask]
    chi_w = chi_kw[mask]

    # Hanning window
    dk = k_w[-1] - k_w[0]
    window = 0.5 * (1 - np.cos(2 * np.pi * (k_w - k_w[0]) / dk))
    chi_windowed = chi_w * window

    r = np.linspace(0, r_max, n_r)
    ft_real = np.zeros_like(r)
    ft_imag = np.zeros_like(r)
    for j, rj in enumerate(r):
        ft_real[j] = np.trapezoid(chi_windowed * np.cos(2 * rj * k_w), k_w)
        ft_imag[j] = np.trapezoid(chi_windowed * np.sin(2 * rj * k_w), k_w)
    ft_mag = np.sqrt(ft_real**2 + ft_imag**2) * np.sqrt(2 / np.pi)
    return r, ft_mag


def plot_ft():
    fig, axes = plt.subplots(3, 1, figsize=(3.5, 6.5), sharex=True)
    fig.subplots_adjust(hspace=0.08)

    labels = ['(a)', '(b)', '(c)']

    for i, edge in enumerate(EDGES):
        ax = axes[i]

        # Load chi(k)*k^2 data and compute FT ourselves (with k-window)
        calc = np.loadtxt(OUT / edge / 'EXAFS.dat')
        exp = np.loadtxt(OUT / edge / 'ipolEXAFS.dat')

        exp_corr = exp[:, 1] * AUTO_S02[edge]

        r_e, ft_e = _compute_ft(exp[:, 0], exp_corr)
        r_c, ft_c = _compute_ft(calc[:, 0], calc[:, 1])

        ax.plot(r_e, ft_e, color=C_EXP, lw=1.0, label='Experiment')
        ax.plot(r_c, ft_c, color=C_FIT[edge], lw=1.0,
                label='EvAX fit', alpha=0.85)

        ax.set_ylabel(r'$|\mathrm{FT}[\chi(k) \cdot k^2]|$')
        ax.set_xlim(0, 6)
        ax.set_ylim(bottom=0)
        ax.text(0.02, 0.92, f'{labels[i]} {edge} K-edge',
                transform=ax.transAxes, fontsize=10, fontweight='bold',
                va='top')

        if i == 0:
            ax.legend(loc='upper right', frameon=False, handlelength=1.5)

    axes[-1].set_xlabel(r'$R$ ($\mathrm{\AA}$)')

    fig.savefig(FIG / 'fig3_ft_comparison.pdf')
    fig.savefig(FIG / 'fig3_ft_comparison.png')
    plt.close(fig)
    print('Fig 3: FT comparison saved')


# ============================================================
# Figure 4: Partial RDFs (3 panels)
# ============================================================
def plot_rdf():
    fig, axes = plt.subplots(3, 1, figsize=(3.5, 6.5), sharex=True)
    fig.subplots_adjust(hspace=0.08)

    labels_panel = ['(a)', '(b)', '(c)']
    pair_labels = {
        'Co': ['Co\u2013Cu', 'Co\u2013Co', 'Co\u2013Ni'],
        'Ni': ['Ni\u2013Cu', 'Ni\u2013Co', 'Ni\u2013Ni'],
        'Cu': ['Cu\u2013Cu', 'Cu\u2013Co', 'Cu\u2013Ni'],
    }

    for i, edge in enumerate(EDGES):
        ax = axes[i]
        rdf = np.loadtxt(OUT / edge / 'rdf.dat', skiprows=1)
        r = rdf[:, 0]

        for j in range(3):
            ax.plot(r, rdf[:, j+1], color=C_PAIR[j], lw=0.9,
                    label=pair_labels[edge][j])

        # Ideal fcc NN distance
        ax.axvline(3.54 / np.sqrt(2), color='0.7', lw=0.5, ls=':', zorder=0)

        ax.set_ylabel(r'$g(R)$')
        ax.set_xlim(1.8, 6.0)
        # Auto-scale y to data range
        rdf_visible = rdf[(r >= 1.8) & (r <= 6.0)]
        ymax = np.max(rdf_visible[:, 1:4]) * 1.15
        ax.set_ylim(0, max(ymax, 5))
        ax.text(0.02, 0.95, f'{labels_panel[i]} {edge} absorber',
                transform=ax.transAxes, fontsize=10, fontweight='bold',
                va='top')
        ax.legend(loc='upper right', frameon=False, handlelength=1.5,
                  fontsize=8, ncol=1)

    axes[-1].set_xlabel(r'$R$ ($\mathrm{\AA}$)')

    fig.savefig(FIG / 'fig4_rdf_partial.pdf')
    fig.savefig(FIG / 'fig4_rdf_partial.png')
    plt.close(fig)
    print('Fig 4: Partial RDF saved')


# ============================================================
# Figure 5: Combined overview (2x2 panels)
# ============================================================
def plot_overview():
    fig = plt.figure(figsize=(7.0, 6.0))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.35)

    # (a) Convergence
    ax_conv = fig.add_subplot(gs[0, 0])
    data = np.loadtxt(OUT / 'output.dat', skiprows=1)
    ax_conv.plot(data[:, 0], data[:, 1], 'k-', lw=0.9)
    ax_conv.set_xlabel('Iteration')
    ax_conv.set_ylabel(r'Best residual $\xi$')
    ax_conv.set_xlim(0, data[-1, 0])
    ax_conv.text(0.03, 0.95, '(a)', transform=ax_conv.transAxes,
                 fontsize=11, fontweight='bold', va='top')

    # (b) EXAFS overlay
    ax_exafs = fig.add_subplot(gs[0, 1])
    offsets = {'Co': 0.4, 'Ni': 0.0, 'Cu': -0.4}
    for edge in EDGES:
        calc = np.loadtxt(OUT / edge / 'EXAFS.dat')
        exp = np.loadtxt(OUT / edge / 'ipolEXAFS.dat')
        off = offsets[edge]
        ax_exafs.plot(exp[:, 0], exp[:, 1] * AUTO_S02[edge] + off, color=C_EXP, lw=0.8)
        ax_exafs.plot(calc[:, 0], calc[:, 1] + off, color=C_FIT[edge],
                      lw=0.8, label=f'{edge} K')
    ax_exafs.set_xlabel(r'$k$ ($\mathrm{\AA}^{-1}$)')
    ax_exafs.set_ylabel(r'$\chi(k) \cdot k^2$ (offset)')
    ax_exafs.set_xlim(0, 14)
    ax_exafs.axvspan(0, 3, alpha=0.08, color='grey', zorder=0)
    ax_exafs.axvspan(12, 14, alpha=0.08, color='grey', zorder=0)
    ax_exafs.legend(loc='upper right', frameon=False, fontsize=8)
    ax_exafs.text(0.03, 0.95, '(b)', transform=ax_exafs.transAxes,
                  fontsize=11, fontweight='bold', va='top')

    # (c) FT overlay
    ax_ft = fig.add_subplot(gs[1, 0])
    ft_offsets = {'Co': 6, 'Ni': 3, 'Cu': 0}
    for edge in EDGES:
        ft_c = np.loadtxt(OUT / edge / 'ft.dat')
        ft_e = np.loadtxt(OUT / edge / 'expFT.dat')
        off = ft_offsets[edge]
        ax_ft.plot(ft_e[:, 0], ft_e[:, 3] * AUTO_S02[edge] + off, color=C_EXP, lw=0.8)
        ax_ft.plot(ft_c[:, 0], ft_c[:, 3] + off, color=C_FIT[edge],
                   lw=0.8, label=f'{edge} K')
    ax_ft.set_xlabel(r'$R$ ($\mathrm{\AA}$)')
    ax_ft.set_ylabel(r'$|\mathrm{FT}|$ (offset)')
    ax_ft.set_xlim(0, 6)
    ax_ft.set_ylim(bottom=0)
    ax_ft.legend(loc='upper right', frameon=False, fontsize=8)
    ax_ft.text(0.03, 0.95, '(c)', transform=ax_ft.transAxes,
               fontsize=11, fontweight='bold', va='top')

    # (d) RDF first shell (all edges combined)
    ax_rdf = fig.add_subplot(gs[1, 1])
    all_pairs = {}
    pair_map = {
        'Co': {0: 'Co-Cu', 1: 'Co-Co', 2: 'Co-Ni'},
        'Ni': {0: 'Ni-Cu', 1: 'Ni-Co', 2: 'Ni-Ni'},
        'Cu': {0: 'Cu-Cu', 1: 'Cu-Co', 2: 'Cu-Ni'},
    }
    # Plot Co absorber RDF as representative
    rdf = np.loadtxt(OUT / 'Co' / 'rdf.dat', skiprows=1)
    r = rdf[:, 0]
    labels_rdf = ['Co\u2013Cu', 'Co\u2013Co', 'Co\u2013Ni']
    for j in range(3):
        ax_rdf.plot(r, rdf[:, j+1], color=C_PAIR[j], lw=0.9,
                    label=labels_rdf[j])
    ax_rdf.axvline(3.54 / np.sqrt(2), color='0.7', lw=0.5, ls=':')
    ax_rdf.set_xlabel(r'$R$ ($\mathrm{\AA}$)')
    ax_rdf.set_ylabel(r'$g(R)$')
    ax_rdf.set_xlim(1.8, 4.5)
    rdf_vis = rdf[(r >= 1.8) & (r <= 4.5)]
    ax_rdf.set_ylim(0, np.max(rdf_vis[:, 1:4]) * 1.15)
    ax_rdf.legend(loc='upper right', frameon=False, fontsize=8)
    ax_rdf.text(0.03, 0.95, '(d)', transform=ax_rdf.transAxes,
                fontsize=11, fontweight='bold', va='top')

    fig.savefig(FIG / 'fig5_overview.pdf')
    fig.savefig(FIG / 'fig5_overview.png')
    plt.close(fig)
    print('Fig 5: Overview saved')


# ============================================================
# Figure 6: Warren-Cowley SRO parameters (bar chart)
# ============================================================
def plot_sro():
    # Calculate from RDFs
    c = {'Co': 85/256, 'Ni': 85/256, 'Cu': 86/256}
    r_min, r_max, dr = 2.2, 2.8, 0.01

    pair_labels_all = {
        'Co': ['Co-Cu', 'Co-Co', 'Co-Ni'],
        'Ni': ['Ni-Cu', 'Ni-Co', 'Ni-Ni'],
        'Cu': ['Cu-Cu', 'Cu-Co', 'Cu-Ni'],
    }

    alphas = {}
    for edge in EDGES:
        rdf = np.loadtxt(OUT / edge / 'rdf.dat', skiprows=1)
        r = rdf[:, 0]
        mask = (r >= r_min) & (r <= r_max)
        cns = [np.sum(rdf[mask, j+1]) * dr for j in range(3)]
        total = sum(cns)
        for j, pair in enumerate(pair_labels_all[edge]):
            neighbor = pair.split('-')[1]
            p_ij = cns[j] / total if total > 0 else 0
            alphas[pair] = 1 - p_ij / c[neighbor]

    # Unique pairs (average symmetric: Co-Ni and Ni-Co)
    unique = ['Co\u2013Co', 'Ni\u2013Ni', 'Cu\u2013Cu',
              'Co\u2013Ni', 'Co\u2013Cu', 'Ni\u2013Cu']
    alpha_vals = [
        alphas['Co-Co'],
        alphas['Ni-Ni'],
        alphas['Cu-Cu'],
        (alphas['Co-Ni'] + alphas['Ni-Co']) / 2,
        (alphas['Co-Cu'] + alphas['Cu-Co']) / 2,
        (alphas['Ni-Cu'] + alphas['Cu-Ni']) / 2,
    ]

    bar_colors = [C_FIT['Co'], C_FIT['Ni'], C_FIT['Cu'],
                  '#666666', '#666666', '#666666']

    fig, ax = plt.subplots(figsize=(4.0, 3.0))
    x = np.arange(len(unique))
    bars = ax.bar(x, alpha_vals, width=0.55, color=bar_colors, edgecolor='k',
                  linewidth=0.5, alpha=0.8)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(unique, fontsize=9, rotation=0)
    ax.set_ylabel(r'Warren-Cowley $\alpha_1$')
    ax.set_ylim(-0.15, 0.17)

    # Add value labels with better spacing
    for bar, val in zip(bars, alpha_vals):
        ypos = val + 0.01 if val >= 0 else val - 0.01
        ax.text(bar.get_x() + bar.get_width()/2, ypos,
                f'{val:+.3f}', ha='center',
                va='bottom' if val >= 0 else 'top',
                fontsize=7.5)

    # Shaded regions for clarity
    ax.axhspan(0, 0.17, alpha=0.03, color='red', zorder=0)
    ax.axhspan(-0.15, 0, alpha=0.03, color='blue', zorder=0)
    ax.text(0.97, 0.95, r'$\alpha > 0$: avoidance', transform=ax.transAxes,
            ha='right', va='top', fontsize=7, color='0.45', style='italic')
    ax.text(0.97, 0.05, r'$\alpha < 0$: preference', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=7, color='0.45', style='italic')

    fig.savefig(FIG / 'fig6_sro_parameters.pdf')
    fig.savefig(FIG / 'fig6_sro_parameters.png')
    plt.close(fig)
    print('Fig 6: SRO parameters saved')


# ============================================================
# Run all
# ============================================================
if __name__ == '__main__':
    print(f'Saving figures to: {FIG}')
    plot_convergence()
    plot_exafs()
    plot_ft()
    plot_rdf()
    plot_overview()
    plot_sro()
    print('\nAll figures saved as PDF + PNG.')
