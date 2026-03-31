#!/usr/bin/env python3
"""
Automated EvAX parameter scan — sequential execution from MEA_analysis dir.

EvAX must run from MEA_analysis/ (needs feff, run-FEFF.exe, evax_temp/).
Each run gets its own parameter file, and output_files/ is backed up after each run.

Usage:
    python scripts/parameter_scan.py              # run all
    python scripts/parameter_scan.py --dry        # show what would be run
    python scripts/parameter_scan.py --results    # only show results table
    python scripts/parameter_scan.py --plot       # generate comparison plot
"""

import subprocess
import re
import sys
import time
import shutil
import numpy as np
from pathlib import Path
from itertools import product
from concurrent.futures import ProcessPoolExecutor, as_completed

# ============================================================
# Configuration
# ============================================================
BASE = Path('/home/nanare/Dokumente/Uni/Master_Wuppertal/ERP_Hausarbeit')
MDIR = BASE / 'MEA_analysis'
SCAN_DIR = MDIR / 'parameter_scan'
TEMPLATE_RESTART = MDIR / 'output_files_run05' / 'restart2002'

SCREENING_ITERATIONS = 50
FROZE_IN = 50

# ============================================================
# Parameter grid — 2 phases
# ============================================================
# Phase 1: Fit parameters (default optimization settings)
PHASE1_GRID = {
    'k_min':   [2, 3],
    'k_max':   ['11.5 11.5 11.5', '11.5 11.5 14'],
    'k_power': [1, 2, 3],
    'Space':   ['w', 'R', 'k'],
    'N_legs':  [2, 4],
}
PHASE1_DEFAULTS = {
    'Maximal_step_length': 0.003,
    'Maximal_displacement': 0.4,
    'R_max_for_FEFF': '6 6 6',
}

# Phase 2: Optimization parameters (with best-guess fit settings)
PHASE2_FIT = {
    'k_min': 3, 'k_max': '11.5 11.5 14', 'k_power': 2,
    'Space': 'w', 'N_legs': 4,
}
PHASE2_GRID = {
    'Maximal_step_length': [0.001, 0.003, 0.005],
    'Maximal_displacement': [0.2, 0.4],
    'R_max_for_FEFF': ['4 4 4', '6 6 6'],
}


# ============================================================
# Functions
# ============================================================

def build_jobs():
    """Build list of (run_name, params) tuples."""
    jobs = []

    # Phase 1
    for k_min, k_max, k_power, space, n_legs in product(
        PHASE1_GRID['k_min'], PHASE1_GRID['k_max'],
        PHASE1_GRID['k_power'], PHASE1_GRID['Space'], PHASE1_GRID['N_legs'],
    ):
        kmax_short = k_max.replace(' ', '-').replace('.', '')
        name = f'kmin{k_min}_kmax{kmax_short}_kw{k_power}_sp{space}_nl{n_legs}'
        params = {
            'k_min': k_min, 'k_max': k_max, 'k_power': k_power,
            'Space': space, 'N_legs': n_legs, **PHASE1_DEFAULTS,
        }
        jobs.append((name, params))

    # Phase 2
    for step, disp, rfeff in product(
        PHASE2_GRID['Maximal_step_length'],
        PHASE2_GRID['Maximal_displacement'],
        PHASE2_GRID['R_max_for_FEFF'],
    ):
        if step == 0.003 and disp == 0.4 and rfeff == '6 6 6':
            continue  # already in phase 1
        rfeff_short = rfeff.replace(' ', '').replace('.', '')
        name = f'opt_step{step}_disp{disp}_rfeff{rfeff_short}'
        params = {
            **PHASE2_FIT,
            'Maximal_step_length': step,
            'Maximal_displacement': disp,
            'R_max_for_FEFF': rfeff,
        }
        jobs.append((name, params))

    return jobs


def make_param_file(run_name, params):
    """Create parameter file from restart template in MEA_analysis/."""
    template = TEMPLATE_RESTART.read_text(encoding='latin-1')

    def replace_param(text, param, value):
        pattern = rf'(              {param}\s+).*'
        return re.sub(pattern, rf'\g<1>{value}', text)

    text = template
    text = replace_param(text, 'Title', f'scan_{run_name}')
    text = replace_param(text, 'S02', '1 1 1')
    text = replace_param(text, 'dE0', '0 0 0')
    text = replace_param(text, 'Stop_after', str(SCREENING_ITERATIONS))
    text = replace_param(text, 'Froze_in', str(FROZE_IN))
    text = replace_param(text, 'k_min', '{v} {v} {v}'.format(v=params['k_min']))
    text = replace_param(text, 'k_max', params['k_max'])
    text = replace_param(text, 'k_power', '{v} {v} {v}'.format(v=params['k_power']))
    text = replace_param(text, 'Space', params['Space'])
    text = replace_param(text, 'N_legs', '{v} {v} {v}'.format(v=params['N_legs']))
    text = replace_param(text, 'Maximal_step_length', str(params['Maximal_step_length']))
    text = replace_param(text, 'Maximal_displacement', str(params['Maximal_displacement']))
    text = replace_param(text, 'R_max_for_FEFF', params['R_max_for_FEFF'])

    param_file = MDIR / f'parameters_scan_{run_name}.dat'
    param_file.write_text(text, encoding='latin-1')
    return param_file


def setup_workdir(run_name):
    """Create a complete EvAX working directory with all required files."""
    work_dir = SCAN_DIR / run_name / 'work'
    work_dir.mkdir(parents=True, exist_ok=True)

    # Symlink EvAX binary
    evax_link = work_dir / 'EvAX.exe'
    if not evax_link.exists():
        evax_link.symlink_to(BASE / 'software' / 'EvAX-6.16_linux.exe')

    # Symlink feff binary
    feff_link = work_dir / 'feff'
    if not feff_link.exists():
        feff_link.symlink_to(BASE / 'feff85exafs' / 'legacy' / 'mod' / 'MONO' / 'feff85L')

    # Copy run-FEFF.exe (tiny shell script)
    run_feff = work_dir / 'run-FEFF.exe'
    if not run_feff.exists():
        shutil.copy2(MDIR / 'run-FEFF.exe', run_feff)
        run_feff.chmod(0o755)

    # Symlink input_files
    inp_link = work_dir / 'input_files'
    if not inp_link.exists():
        inp_link.symlink_to(MDIR / 'input_files')

    # Create output_files dir
    (work_dir / 'output_files').mkdir(exist_ok=True)

    return work_dir


def run_evax(run_name, params):
    """Run EvAX in isolated working directory. Safe for parallel execution."""
    result_dir = SCAN_DIR / run_name

    # Skip if already completed
    if (result_dir / 'output.dat').exists():
        lines = (result_dir / 'output.dat').read_text().strip().split('\n')
        if len(lines) > 5:
            return run_name, 'skipped', 0

    # Create isolated working directory
    work_dir = setup_workdir(run_name)

    # Create parameter file IN the work dir
    param_file = work_dir / 'parameters.dat'
    template = TEMPLATE_RESTART.read_text(encoding='latin-1')

    def replace_param(text, param, value):
        pattern = rf'(              {param}\s+).*'
        return re.sub(pattern, rf'\g<1>{value}', text)

    text = template
    text = replace_param(text, 'Title', f'scan_{run_name}')
    text = replace_param(text, 'S02', '1 1 1')
    text = replace_param(text, 'dE0', '0 0 0')
    text = replace_param(text, 'Stop_after', str(SCREENING_ITERATIONS))
    text = replace_param(text, 'Froze_in', str(FROZE_IN))
    text = replace_param(text, 'k_min', '{v} {v} {v}'.format(v=params['k_min']))
    text = replace_param(text, 'k_max', params['k_max'])
    text = replace_param(text, 'k_power', '{v} {v} {v}'.format(v=params['k_power']))
    text = replace_param(text, 'Space', params['Space'])
    text = replace_param(text, 'N_legs', '{v} {v} {v}'.format(v=params['N_legs']))
    text = replace_param(text, 'Maximal_step_length', str(params['Maximal_step_length']))
    text = replace_param(text, 'Maximal_displacement', str(params['Maximal_displacement']))
    text = replace_param(text, 'R_max_for_FEFF', params['R_max_for_FEFF'])

    # Reset iteration counter so EvAX starts fresh (not continuing from Run 05)
    text = re.sub(r'Total_iterations:\s+\d+', 'Total_iterations:                         0', text)
    # Reset residuals so EvAX doesn't think it already converged
    text = re.sub(r'Residuals:.*', 'Residuals:', text)

    param_file.write_text(text, encoding='latin-1')

    t0 = time.time()
    try:
        result = subprocess.run(
            ['./EvAX.exe', 'parameters.dat', '-autoES'],
            cwd=str(work_dir),
            capture_output=True, timeout=7200,
        )

        # Save log to result dir
        result_dir.mkdir(parents=True, exist_ok=True)
        (result_dir / 'evax.log').write_bytes(result.stdout + result.stderr)

        # Copy key results out of work_dir
        out_dir = work_dir / 'output_files'
        if out_dir.exists():
            out_dat = out_dir / 'output.dat'
            if out_dat.exists():
                shutil.copy2(out_dat, result_dir / 'output.dat')
            for edge in ['Co', 'Ni', 'Cu']:
                edge_src = out_dir / edge
                if edge_src.exists():
                    edge_dst = result_dir / edge
                    if edge_dst.exists():
                        shutil.rmtree(edge_dst)
                    shutil.copytree(edge_src, edge_dst)

        # Clean up work dir (save disk space)
        shutil.rmtree(work_dir, ignore_errors=True)

        elapsed = time.time() - t0
        return run_name, 'done', elapsed

    except subprocess.TimeoutExpired:
        return run_name, 'timeout', time.time() - t0
    except Exception as e:
        return run_name, f'error: {e}', time.time() - t0


def parse_results(result_dir):
    """Parse EvAX output for residuals."""
    results = {'total': None, 'Co': None, 'Ni': None, 'Cu': None}

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

    return results


def print_results_table():
    """Print sorted comparison table."""
    if not SCAN_DIR.exists():
        print("No scan results found!")
        return

    rows = []
    for d in sorted(SCAN_DIR.iterdir()):
        if not d.is_dir():
            continue
        results = parse_results(d)
        if results['total'] is not None:
            # Parse params from name
            name = d.name
            params = {}
            for p in name.split('_'):
                if p.startswith('kmin'): params['kmin'] = p[4:]
                elif p.startswith('kmax'): params['kmax'] = p[4:]
                elif p.startswith('kw'): params['kw'] = p[2:]
                elif p.startswith('sp'): params['sp'] = p[2:]
                elif p.startswith('nl'): params['nl'] = p[2:]
                elif p.startswith('step'): params['step'] = p[4:]
                elif p.startswith('disp'): params['disp'] = p[4:]
                elif p.startswith('rfeff'): params['rfeff'] = p[5:]
            rows.append({**params, **results, 'name': name})

    if not rows:
        print("No completed runs found!")
        return

    rows.sort(key=lambda r: r.get('total', 999))

    print(f"\n{'='*110}")
    print(f"{'#':>3} {'Run':<45} {'kmin':>4} {'kw':>2} {'Space':>5} {'Nleg':>4} "
          f"{'Co':>7} {'Ni':>7} {'Cu':>7} {'TOTAL':>7}")
    print(f"{'-'*110}")
    for i, r in enumerate(rows):
        marker = '  *** BEST' if i == 0 else ''
        print(f"{i+1:>3} {r['name']:<45} "
              f"{r.get('kmin',''):>4} {r.get('kw',''):>2} "
              f"{r.get('sp',''):>5} {r.get('nl',''):>4} "
              f"{r.get('Co', 0):>7.4f} {r.get('Ni', 0):>7.4f} "
              f"{r.get('Cu', 0):>7.4f} {r.get('total', 0):>7.4f}{marker}")
    print(f"{'='*110}")
    print(f"\nCompleted: {len(rows)} runs")

    # Save CSV
    csv_file = SCAN_DIR / 'results.csv'
    with open(csv_file, 'w') as f:
        f.write('rank,run,k_min,k_power,Space,N_legs,residual_Co,residual_Ni,residual_Cu,residual_total\n')
        for i, r in enumerate(rows):
            f.write(f"{i+1},{r['name']},{r.get('kmin','')},{r.get('kw','')},"
                    f"{r.get('sp','')},{r.get('nl','')},"
                    f"{r.get('Co',0):.6f},{r.get('Ni',0):.6f},{r.get('Cu',0):.6f},"
                    f"{r.get('total',0):.6f}\n")
    print(f"CSV: {csv_file}")
    return rows


def plot_results():
    """Generate comparison plots from scan results."""
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    mpl.rcParams.update({
        'font.family': 'serif', 'font.size': 9,
        'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    })

    rows = []
    for d in sorted(SCAN_DIR.iterdir()):
        if not d.is_dir():
            continue
        results = parse_results(d)
        if results['total'] is not None:
            name = d.name
            params = {}
            for p in name.split('_'):
                if p.startswith('kmin'): params['kmin'] = int(p[4:])
                elif p.startswith('kw'): params['kw'] = int(p[2:])
                elif p.startswith('sp'): params['sp'] = p[2:]
                elif p.startswith('nl'): params['nl'] = int(p[2:])
            rows.append({**params, **results, 'name': name})

    if not rows:
        print("No results to plot!")
        return

    rows.sort(key=lambda r: r['total'])

    # --- Plot 1: Bar chart of all runs sorted by residual ---
    fig, ax = plt.subplots(figsize=(12, 5))
    names = [r['name'][:35] for r in rows]
    totals = [r['total'] for r in rows]
    colors = ['#2ecc71' if t == min(totals) else '#3498db' if t < np.median(totals) else '#e74c3c'
              for t in totals]
    ax.barh(range(len(rows)), totals, color=colors, height=0.7)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(names, fontsize=6)
    ax.set_xlabel('Total Residual')
    ax.set_title('Parameter Scan Results (sorted by residual)')
    ax.invert_yaxis()
    ax.axvline(totals[0], color='green', ls='--', lw=0.5, alpha=0.5)
    fig.savefig(BASE / 'figures' / 'fig8_parameter_scan.png')
    fig.savefig(BASE / 'figures' / 'fig8_parameter_scan.pdf')
    plt.close()

    # --- Plot 2: Heatmaps per parameter ---
    # Only for phase 1 runs (have all fit params)
    p1_rows = [r for r in rows if 'sp' in r and 'kw' in r]
    if p1_rows:
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        fig.subplots_adjust(wspace=0.4)

        # k_power vs Space
        ax = axes[0]
        for sp_i, sp in enumerate(['k', 'R', 'w']):
            for kw in [1, 2, 3]:
                vals = [r['total'] for r in p1_rows if r.get('sp') == sp and r.get('kw') == kw]
                if vals:
                    ax.scatter(kw, sp_i, c=min(vals), cmap='RdYlGn_r',
                              vmin=min(totals), vmax=np.percentile(totals, 75),
                              s=200, edgecolors='k', linewidth=0.5)
                    ax.text(kw, sp_i, f'{min(vals):.3f}', ha='center', va='center', fontsize=7)
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(['k¹', 'k²', 'k³'])
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels(['k-space', 'R-space', 'Wavelet'])
        ax.set_title('k-weight vs Space\n(best residual per combo)')

        # k_min vs N_legs
        ax = axes[1]
        for km_i, km in enumerate([2, 3]):
            for nl in [2, 4]:
                vals = [r['total'] for r in p1_rows if r.get('kmin') == km and r.get('nl') == nl]
                if vals:
                    nl_i = 0 if nl == 2 else 1
                    ax.scatter(km_i, nl_i, c=min(vals), cmap='RdYlGn_r',
                              vmin=min(totals), vmax=np.percentile(totals, 75),
                              s=200, edgecolors='k', linewidth=0.5)
                    ax.text(km_i, nl_i, f'{min(vals):.3f}', ha='center', va='center', fontsize=7)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['k_min=2', 'k_min=3'])
        ax.set_yticks([0, 1])
        ax.set_yticklabels(['N_legs=2', 'N_legs=4'])
        ax.set_title('k_min vs N_legs\n(best residual per combo)')

        # Per-edge residuals for top 10
        ax = axes[2]
        top10 = p1_rows[:min(10, len(p1_rows))]
        x = np.arange(len(top10))
        w = 0.25
        ax.bar(x - w, [r['Co'] for r in top10], w, label='Co', color='#0072B2')
        ax.bar(x, [r['Ni'] for r in top10], w, label='Ni', color='#009E73')
        ax.bar(x + w, [r['Cu'] for r in top10], w, label='Cu', color='#D55E00')
        ax.set_xticks(x)
        ax.set_xticklabels([f"#{i+1}" for i in range(len(top10))], fontsize=7)
        ax.set_ylabel('Residual')
        ax.set_title('Per-edge residuals (Top 10)')
        ax.legend(fontsize=7)

        fig.savefig(BASE / 'figures' / 'fig9_parameter_heatmap.png')
        fig.savefig(BASE / 'figures' / 'fig9_parameter_heatmap.pdf')
        plt.close()

    print("Plots saved: fig8_parameter_scan, fig9_parameter_heatmap")


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    dry_run = '--dry' in sys.argv
    results_only = '--results' in sys.argv
    plot_only = '--plot' in sys.argv

    if results_only:
        print_results_table()
        sys.exit(0)

    if plot_only:
        plot_results()
        sys.exit(0)

    # Parse --parallel N
    n_parallel = 4
    for i, arg in enumerate(sys.argv):
        if arg == '--parallel' and i + 1 < len(sys.argv):
            n_parallel = int(sys.argv[i + 1])

    if not TEMPLATE_RESTART.exists():
        print(f"ERROR: Template not found: {TEMPLATE_RESTART}")
        sys.exit(1)

    SCAN_DIR.mkdir(exist_ok=True)
    jobs = build_jobs()

    est_min = len(jobs) * 15 / n_parallel
    print(f"Parameter scan: {len(jobs)} combinations, {n_parallel} parallel workers")
    print(f"Screening: {SCREENING_ITERATIONS} iterations per run")
    print(f"Estimated time: ~{est_min:.0f} min ({est_min / 60:.1f} h)\n")

    if dry_run:
        for i, (name, params) in enumerate(jobs):
            print(f"  [{i+1:>3}/{len(jobs)}] {name}")
        print(f"\nTotal: {len(jobs)} runs")
        sys.exit(0)

    # Run in parallel with isolated work directories
    completed = 0
    skipped = 0
    failed = 0
    t_start = time.time()

    with ProcessPoolExecutor(max_workers=n_parallel) as executor:
        futures = {
            executor.submit(run_evax, name, params): name
            for name, params in jobs
        }
        for future in as_completed(futures):
            try:
                run_name, status, elapsed = future.result()
                total_done = completed + skipped + failed + 1
                if status == 'skipped':
                    skipped += 1
                    print(f"  [{total_done}/{len(jobs)}] {run_name}: skipped")
                elif status == 'done':
                    completed += 1
                    results = parse_results(SCAN_DIR / run_name)
                    res = results['total']
                    res_str = f"residual={res:.4f}" if res else "done"
                    print(f"  [{total_done}/{len(jobs)}] {run_name}: {res_str} ({elapsed/60:.1f} min)")
                else:
                    failed += 1
                    print(f"  [{total_done}/{len(jobs)}] {run_name}: FAILED ({status})")
            except Exception as e:
                failed += 1
                print(f"  {futures[future]}: EXCEPTION {e}")

    total_time = (time.time() - t_start) / 60
    print(f"\n{'='*60}")
    print(f"Done! {completed} completed, {skipped} skipped, {failed} failed")
    print(f"Total time: {total_time:.1f} min ({total_time/60:.1f} h)")
    print(f"{'='*60}")

    print_results_table()
    print("\nGenerating plots...")
    plot_results()
