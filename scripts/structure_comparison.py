#!/usr/bin/env python3
"""
Run EvAX with different starting structures (FCC, BCC, HCP) to determine
which crystal structure best fits the MEA CoNiCu EXAFS data.

Uses optimal parameters from parameter scan + autoES.
Final production runs with full convergence.

Usage:
    python scripts/structure_comparison.py              # run all 3
    python scripts/structure_comparison.py --dry        # show what would run
    python scripts/structure_comparison.py --results    # show results table
    python scripts/structure_comparison.py --plot       # generate all plots
"""

import subprocess
import re
import sys
import time
import shutil
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

BASE = Path('/home/nanare/Dokumente/Uni/Master_Wuppertal/ERP_Hausarbeit')
MDIR = BASE / 'MEA_analysis'
SCAN_DIR = MDIR / 'structure_comparison'
TEMPLATE = MDIR / 'output_files_run05' / 'restart2002'
FIG_DIR = BASE / 'figures'

# Production run settings (~7h max per structure)
SCREENING_ITERATIONS = 100
FROZE_IN = 150

STRUCTURES = {
    'fcc': 'input_files/MEA_CoNiCu.p1',
    'bcc': 'input_files/MEA_CoNiCu_bcc.p1',
    'hcp': 'input_files/MEA_CoNiCu_hcp.p1',
}

STRUCT_LABELS = {'fcc': 'FCC', 'bcc': 'BCC', 'hcp': 'HCP'}
STRUCT_COLORS = {'fcc': '#2ecc71', 'bcc': '#3498db', 'hcp': '#e74c3c'}
EDGES = ['Co', 'Ni', 'Cu']
EDGE_COLORS = {'Co': '#0072B2', 'Ni': '#009E73', 'Cu': '#D55E00'}


# ============================================================
# EvAX run infrastructure
# ============================================================

def setup_workdir(name):
    work_dir = SCAN_DIR / name / 'work'
    work_dir.mkdir(parents=True, exist_ok=True)

    evax_link = work_dir / 'EvAX.exe'
    if not evax_link.exists():
        evax_link.symlink_to(BASE / 'software' / 'EvAX-6.16_linux.exe')

    feff_link = work_dir / 'feff'
    if not feff_link.exists():
        feff_link.symlink_to(BASE / 'feff85exafs' / 'legacy' / 'mod' / 'MONO' / 'feff85L')

    run_feff = work_dir / 'run-FEFF.exe'
    if not run_feff.exists():
        shutil.copy2(MDIR / 'run-FEFF.exe', run_feff)
        run_feff.chmod(0o755)

    inp_link = work_dir / 'input_files'
    if not inp_link.exists():
        inp_link.symlink_to(MDIR / 'input_files')

    (work_dir / 'output_files').mkdir(exist_ok=True)
    return work_dir


def run_evax(name, structure_file):
    result_dir = SCAN_DIR / name

    if (result_dir / 'output.dat').exists():
        lines = (result_dir / 'output.dat').read_text().strip().split('\n')
        if len(lines) > 5:
            return name, 'skipped', 0

    work_dir = setup_workdir(name)

    template = TEMPLATE.read_text(encoding='latin-1')

    def replace_param(text, param, value):
        pattern = rf'(              {param}\s+).*'
        return re.sub(pattern, rf'\g<1>{value}', text)

    text = template
    text = replace_param(text, 'Title', f'structure_{name}')
    text = replace_param(text, 'File_with_structure', structure_file)
    text = replace_param(text, 'S02', '1 1 1')
    text = replace_param(text, 'dE0', '0 0 0')
    text = replace_param(text, 'Stop_after', str(SCREENING_ITERATIONS))
    text = replace_param(text, 'Froze_in', str(FROZE_IN))
    # Optimal params from parameter scan
    text = replace_param(text, 'k_min', '3 3 3')
    text = replace_param(text, 'k_max', '11.5 11.5 11.5')
    text = replace_param(text, 'k_power', '2 2 2')
    text = replace_param(text, 'Space', 'w')
    text = replace_param(text, 'N_legs', '4 4 4')
    text = replace_param(text, 'Maximal_step_length', '0.005')
    text = replace_param(text, 'Maximal_displacement', '0.4')
    text = replace_param(text, 'R_max_for_FEFF', '6 6 6')

    # Reset iteration counter — fresh start
    text = re.sub(r'Total_iterations:\s+\d+', 'Total_iterations:                         0', text)
    text = re.sub(r'Residuals:.*', 'Residuals:', text)
    # Remove old ATOMS section — force fresh start from structure file
    atoms_idx = text.find('\nATOMS\n')
    if atoms_idx > 0:
        text = text[:atoms_idx + 1]

    param_file = work_dir / 'parameters.dat'
    param_file.write_text(text, encoding='latin-1')

    t0 = time.time()
    try:
        result = subprocess.run(
            ['./EvAX.exe', 'parameters.dat', '-autoES'],
            cwd=str(work_dir),
            capture_output=True, timeout=86400,  # 24h max
        )

        result_dir.mkdir(parents=True, exist_ok=True)
        (result_dir / 'evax.log').write_bytes(result.stdout + result.stderr)

        out_dir = work_dir / 'output_files'
        if out_dir.exists():
            # Copy everything from output_files
            for item in out_dir.iterdir():
                dst = result_dir / item.name
                if item.is_dir():
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(item, dst)
                else:
                    shutil.copy2(item, dst)

        shutil.rmtree(work_dir, ignore_errors=True)

        elapsed = time.time() - t0
        return name, 'done', elapsed

    except subprocess.TimeoutExpired:
        return name, 'timeout', time.time() - t0
    except Exception as e:
        return name, f'error: {e}', time.time() - t0


def parse_results(result_dir):
    results = {'total': None, 'Co': None, 'Ni': None, 'Cu': None,
               's02_Co': None, 's02_Ni': None, 's02_Cu': None,
               'dE0_Co': None, 'dE0_Ni': None, 'dE0_Cu': None}

    out_file = result_dir / 'output.dat'
    if out_file.exists():
        lines = out_file.read_text().strip().split('\n')
        if len(lines) > 1:
            vals = lines[-1].split()
            if len(vals) > 1:
                try:
                    results['total'] = float(vals[1])
                except ValueError:
                    pass

    log_file = result_dir / 'evax.log'
    if log_file.exists():
        log = log_file.read_text(encoding='latin-1', errors='replace')
        precise = re.findall(r'Precise/approximate:\s+([\d.]+)', log)
        if len(precise) >= 3:
            results['Co'] = float(precise[-3])
            results['Ni'] = float(precise[-2])
            results['Cu'] = float(precise[-1])
        s02s = re.findall(r'S02\s*=\s*([\d.]+)', log)
        if len(s02s) >= 3:
            results['s02_Co'] = float(s02s[-3])
            results['s02_Ni'] = float(s02s[-2])
            results['s02_Cu'] = float(s02s[-1])
        de0s = re.findall(r'dE0\s*=\s*([+-]?[\d.]+)', log)
        if len(de0s) >= 3:
            results['dE0_Co'] = float(de0s[-3])
            results['dE0_Ni'] = float(de0s[-2])
            results['dE0_Cu'] = float(de0s[-1])

    return results


def print_results():
    if not SCAN_DIR.exists():
        print("No results found!")
        return

    print(f"\n{'='*90}")
    print(f"{'Structure':<10} {'Co':>8} {'Ni':>8} {'Cu':>8} {'TOTAL':>8}  "
          f"{'S02_Co':>7} {'S02_Ni':>7} {'S02_Cu':>7}")
    print(f"{'-'*90}")

    rows = []
    for name in ['fcc', 'bcc', 'hcp']:
        d = SCAN_DIR / name
        if not d.exists():
            continue
        r = parse_results(d)
        rows.append((name, r))
        co = r.get('Co', 0) or 0
        ni = r.get('Ni', 0) or 0
        cu = r.get('Cu', 0) or 0
        tot = r.get('total', 0) or 0
        s_co = r.get('s02_Co', 0) or 0
        s_ni = r.get('s02_Ni', 0) or 0
        s_cu = r.get('s02_Cu', 0) or 0
        print(f"{name.upper():<10} {co:>8.4f} {ni:>8.4f} {cu:>8.4f} {tot:>8.4f}  "
              f"{s_co:>7.3f} {s_ni:>7.3f} {s_cu:>7.3f}")

    if rows:
        best = min(rows, key=lambda x: x[1].get('total', 999) or 999)
        print(f"\n  Best fit: {best[0].upper()} (residual = {best[1]['total']:.4f})")
    print(f"{'='*90}")
    return rows


# ============================================================
# Plotting — individual per structure + combined comparison
# ============================================================

def _compute_ft(k, chi, k_min=3, k_max=11.5, k_power=2, dr=0.02, r_max=6):
    """Compute FT magnitude from chi(k) with Hanning window."""
    mask = (k >= k_min) & (k <= k_max)
    k_w = k[mask]
    chi_w = chi[mask] * k_w**k_power

    # Hanning window
    dk_win = 0.5
    win = np.ones_like(k_w)
    lo = k_w <= k_min + dk_win
    hi = k_w >= k_max - dk_win
    if np.any(lo):
        win[lo] = 0.5 * (1 - np.cos(np.pi * (k_w[lo] - k_min) / dk_win))
    if np.any(hi):
        win[hi] = 0.5 * (1 - np.cos(np.pi * (k_max - k_w[hi]) / dk_win))
    chi_w *= win

    r = np.arange(0, r_max, dr)
    ft_re = np.zeros_like(r)
    ft_im = np.zeros_like(r)
    dk = k_w[1] - k_w[0] if len(k_w) > 1 else 0.05
    for i, ri in enumerate(r):
        ft_re[i] = np.sum(chi_w * np.cos(2 * k_w * ri)) * dk / np.sqrt(np.pi)
        ft_im[i] = np.sum(chi_w * np.sin(2 * k_w * ri)) * dk / np.sqrt(np.pi)

    return r, np.sqrt(ft_re**2 + ft_im**2)


def plot_single_structure(name):
    """Generate individual detail plots for one structure (like fig1-fig6)."""
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    mpl.rcParams.update({
        'font.family': 'serif', 'font.size': 9,
        'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    })

    result_dir = SCAN_DIR / name
    r = parse_results(result_dir)
    label = STRUCT_LABELS[name]

    # --- Convergence plot ---
    out_file = result_dir / 'output.dat'
    if out_file.exists():
        data = np.loadtxt(out_file)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(data[:, 0], data[:, 1], 'k-', lw=1)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Residual')
        ax.set_title(f'{label} — Convergence')
        ax.set_ylim(bottom=0)
        fig.savefig(FIG_DIR / f'struct_{name}_convergence.png')
        fig.savefig(FIG_DIR / f'struct_{name}_convergence.pdf')
        plt.close()

    # --- EXAFS + FT (2×3 grid) ---
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    fig.suptitle(f'{label} Structure — EXAFS Fit', fontsize=13, fontweight='bold')

    for i, edge in enumerate(EDGES):
        color = EDGE_COLORS[edge]

        # Row 0: EXAFS chi(k)*k²
        ax = axes[0, i]
        exafs_file = result_dir / edge / 'EXAFS.dat'
        if exafs_file.exists():
            data = np.loadtxt(exafs_file)
            k = data[:, 0]
            chi_calc = data[:, 1]
            ax.plot(k, chi_calc * k**2, color=color, lw=1.2, label='Calc.')
            if data.shape[1] >= 3:
                chi_exp = data[:, 2]
                ax.plot(k, chi_exp * k**2, 'k--', lw=0.8, label='Exp.')
        ax.axvspan(0, 3, alpha=0.1, color='grey')
        ax.set_xlabel('k ($\\AA^{-1}$)')
        ax.set_ylabel('$\\chi(k) \\cdot k^2$')
        ax.set_title(f'{edge} K-edge')
        ax.legend(fontsize=7)
        ax.set_xlim(0, 14)

        # Row 1: FT
        ax = axes[1, i]
        if exafs_file.exists():
            data = np.loadtxt(exafs_file)
            k = data[:, 0]
            # Calc FT
            r_ft, ft_calc = _compute_ft(k, data[:, 1])
            ax.plot(r_ft, ft_calc, color=color, lw=1.2, label='Calc.')
            if data.shape[1] >= 3:
                r_ft_e, ft_exp = _compute_ft(k, data[:, 2])
                ax.plot(r_ft_e, ft_exp, 'k--', lw=0.8, label='Exp.')
        ax.set_xlabel('R ($\\AA$)')
        ax.set_ylabel('|FT|')
        ax.set_title(f'{edge} K-edge FT')
        ax.legend(fontsize=7)
        ax.set_xlim(0, 6)

    fig.tight_layout()
    fig.savefig(FIG_DIR / f'struct_{name}_exafs_ft.png')
    fig.savefig(FIG_DIR / f'struct_{name}_exafs_ft.pdf')
    plt.close()

    # --- RDF ---
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    fig.suptitle(f'{label} Structure — Radial Distribution Functions', fontsize=12, fontweight='bold')
    pair_colors = {
        'Co': '#0072B2', 'Ni': '#009E73', 'Cu': '#D55E00',
    }
    for i, edge in enumerate(EDGES):
        ax = axes[i]
        for target in EDGES:
            rdf_file = result_dir / edge / f'RDF_{edge}_{target}.dat'
            if not rdf_file.exists():
                rdf_file = result_dir / edge / f'RDF_{target}.dat'
            if rdf_file.exists():
                data = np.loadtxt(rdf_file)
                ax.plot(data[:, 0], data[:, 1], color=pair_colors[target],
                        lw=1.2, label=f'{edge}-{target}')
        ax.set_xlabel('R ($\\AA$)')
        ax.set_ylabel('g(r)')
        ax.set_title(f'{edge} absorber')
        ax.legend(fontsize=7)
        ax.set_xlim(1.5, 5.5)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f'struct_{name}_rdf.png')
    fig.savefig(FIG_DIR / f'struct_{name}_rdf.pdf')
    plt.close()

    print(f"  Individual plots saved for {label}")


def plot_comparison():
    """Publication-quality comparison plot: all 3 structures together."""
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.gridspec import GridSpec

    mpl.rcParams.update({
        'font.family': 'serif', 'font.size': 9,
        'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
        'axes.linewidth': 0.8,
    })

    structs = ['fcc', 'bcc', 'hcp']

    # Load all results
    results = {}
    for s in structs:
        d = SCAN_DIR / s
        if d.exists() and (d / 'output.dat').exists():
            results[s] = parse_results(d)

    if len(results) < 2:
        print("Need at least 2 completed structures for comparison!")
        return

    available = [s for s in structs if s in results]

    # ================================================================
    # Main comparison figure (paper-ready)
    # Layout: 4 rows × 3 columns
    #   Row 0: Residual bar chart (spans all 3 cols)
    #   Row 1: EXAFS chi(k)*k² per edge
    #   Row 2: FT per edge
    #   Row 3: Residual breakdown table + fit parameters
    # ================================================================
    fig = plt.figure(figsize=(7.2, 9.5))  # single-column journal width
    gs = GridSpec(4, 3, figure=fig, height_ratios=[1.0, 1.2, 1.2, 0.8],
                  hspace=0.45, wspace=0.35)

    # --- Row 0: Residual comparison bar chart ---
    ax_res = fig.add_subplot(gs[0, :])
    x = np.arange(len(available))
    w = 0.18
    for i, edge in enumerate(EDGES):
        vals = [results[s].get(edge, 0) or 0 for s in available]
        bars = ax_res.bar(x + (i - 1) * w, vals, w, label=f'{edge} K-edge',
                          color=EDGE_COLORS[edge], edgecolor='black', linewidth=0.4)
    # Total weighted residual as diamond
    totals = [results[s].get('total', 0) or 0 for s in available]
    ax_res.scatter(x, totals, marker='D', s=80, c='black', zorder=5,
                   label='Weighted total')
    ax_res.set_xticks(x)
    ax_res.set_xticklabels([STRUCT_LABELS[s] for s in available],
                           fontsize=11, fontweight='bold')
    ax_res.set_ylabel('R-factor')
    ax_res.set_title('(a) Fit quality by crystal structure', fontsize=10,
                     fontweight='bold', loc='left')
    ax_res.legend(fontsize=7, ncol=4, loc='upper right')

    # Highlight best
    best_idx = np.argmin(totals)
    for i, t in enumerate(totals):
        ax_res.text(i, t + 0.01, f'{t:.3f}', ha='center', va='bottom',
                    fontsize=8, fontweight='bold' if i == best_idx else 'normal',
                    color='green' if i == best_idx else 'black')

    # --- Row 1: EXAFS chi(k)*k² ---
    for col, edge in enumerate(EDGES):
        ax = fig.add_subplot(gs[1, col])
        panel = chr(ord('b') + col)

        for s in available:
            exafs_file = SCAN_DIR / s / edge / 'EXAFS.dat'
            if exafs_file.exists():
                data = np.loadtxt(exafs_file)
                k = data[:, 0]
                ax.plot(k, data[:, 1] * k**2, color=STRUCT_COLORS[s],
                        lw=1.0, label=STRUCT_LABELS[s])

        # Experimental (from first available)
        for s in available:
            exafs_file = SCAN_DIR / s / edge / 'EXAFS.dat'
            if exafs_file.exists():
                data = np.loadtxt(exafs_file)
                if data.shape[1] >= 3:
                    ax.plot(data[:, 0], data[:, 2] * data[:, 0]**2,
                            'k-', lw=0.6, alpha=0.4, label='Exp.')
                break

        ax.axvspan(0, 3, alpha=0.08, color='grey')
        ax.set_xlabel('k ($\\AA^{-1}$)')
        if col == 0:
            ax.set_ylabel('$\\chi(k) \\cdot k^2$')
        ax.set_title(f'({panel}) {edge} K-edge EXAFS', fontsize=9,
                     fontweight='bold', loc='left')
        ax.legend(fontsize=6, loc='upper right')
        ax.set_xlim(2, 12)

    # --- Row 2: Fourier Transform ---
    for col, edge in enumerate(EDGES):
        ax = fig.add_subplot(gs[2, col])
        panel = chr(ord('e') + col)

        for s in available:
            exafs_file = SCAN_DIR / s / edge / 'EXAFS.dat'
            if exafs_file.exists():
                data = np.loadtxt(exafs_file)
                k = data[:, 0]
                r_ft, ft_mag = _compute_ft(k, data[:, 1])
                ax.plot(r_ft, ft_mag, color=STRUCT_COLORS[s],
                        lw=1.0, label=STRUCT_LABELS[s])

        # Experimental FT
        for s in available:
            exafs_file = SCAN_DIR / s / edge / 'EXAFS.dat'
            if exafs_file.exists():
                data = np.loadtxt(exafs_file)
                if data.shape[1] >= 3:
                    r_ft_e, ft_exp = _compute_ft(k, data[:, 2])
                    ax.plot(r_ft_e, ft_exp, 'k-', lw=0.6, alpha=0.4, label='Exp.')
                break

        ax.set_xlabel('R ($\\AA$)')
        if col == 0:
            ax.set_ylabel('|FT($\\chi(k) \\cdot k^2$)|')
        ax.set_title(f'({panel}) {edge} K-edge FT', fontsize=9,
                     fontweight='bold', loc='left')
        ax.legend(fontsize=6, loc='upper right')
        ax.set_xlim(0, 6)

    # --- Row 3: Fit parameters table ---
    ax_tab = fig.add_subplot(gs[3, :])
    ax_tab.axis('off')
    ax_tab.set_title('(h) Fit parameters', fontsize=9, fontweight='bold', loc='left')

    col_labels = ['Structure', 'R-factor\n(total)',
                  'R(Co)', 'R(Ni)', 'R(Cu)',
                  'S$_0^2$(Co)', 'S$_0^2$(Ni)', 'S$_0^2$(Cu)']
    cell_text = []
    cell_colors = []
    min_total = min(results[s].get('total', 999) or 999 for s in available)

    for s in available:
        r = results[s]
        tot = r.get('total', 0) or 0
        row = [
            STRUCT_LABELS[s],
            f"{tot:.4f}",
            f"{r.get('Co', 0) or 0:.4f}",
            f"{r.get('Ni', 0) or 0:.4f}",
            f"{r.get('Cu', 0) or 0:.4f}",
            f"{r.get('s02_Co', 0) or 0:.3f}",
            f"{r.get('s02_Ni', 0) or 0:.3f}",
            f"{r.get('s02_Cu', 0) or 0:.3f}",
        ]
        cell_text.append(row)
        if tot == min_total:
            cell_colors.append(['#d4edda'] * len(col_labels))
        else:
            cell_colors.append(['white'] * len(col_labels))

    table = ax_tab.table(cellText=cell_text, colLabels=col_labels,
                         cellColours=cell_colors,
                         loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.4)

    # Bold header
    for j in range(len(col_labels)):
        table[0, j].set_text_props(fontweight='bold')

    fig.savefig(FIG_DIR / 'fig10_structure_comparison.png')
    fig.savefig(FIG_DIR / 'fig10_structure_comparison.pdf')
    plt.close()
    print("Saved: fig10_structure_comparison")


def plot_all():
    """Generate all plots: individual + comparison."""
    for name in ['fcc', 'bcc', 'hcp']:
        if (SCAN_DIR / name / 'output.dat').exists():
            plot_single_structure(name)
    plot_comparison()


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    if '--results' in sys.argv:
        print_results()
        sys.exit(0)

    if '--plot' in sys.argv:
        plot_all()
        sys.exit(0)

    dry_run = '--dry' in sys.argv
    n_parallel = 3

    if not TEMPLATE.exists():
        print(f"ERROR: Template not found: {TEMPLATE}")
        sys.exit(1)

    SCAN_DIR.mkdir(exist_ok=True)

    jobs = [(name, sf) for name, sf in STRUCTURES.items()]

    print(f"Structure comparison: {len(jobs)} structures (FCC, BCC, HCP)")
    print(f"Screening: {SCREENING_ITERATIONS} iterations, Froze_in: {FROZE_IN}")
    print(f"This is the FINAL production run — expect several hours.\n")

    if dry_run:
        for name, sf in jobs:
            print(f"  {name}: {sf}")
        sys.exit(0)

    t_start = time.time()
    with ProcessPoolExecutor(max_workers=n_parallel) as executor:
        futures = {
            executor.submit(run_evax, name, sf): name
            for name, sf in jobs
        }
        for future in as_completed(futures):
            try:
                name, status, elapsed = future.result()
                r = parse_results(SCAN_DIR / name) if status == 'done' else {}
                res = r.get('total')
                res_str = f"residual={res:.4f}" if res else status
                print(f"  {name}: {res_str} ({elapsed/60:.1f} min)")
            except Exception as e:
                print(f"  {futures[future]}: EXCEPTION {e}")

    total_time = (time.time() - t_start) / 60
    print(f"\nDone! Total time: {total_time:.1f} min ({total_time/60:.1f} h)")

    print_results()
    print("\nGenerating all plots...")
    plot_all()
