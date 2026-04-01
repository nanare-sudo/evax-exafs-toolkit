#!/usr/bin/env python3
"""
EvAX Engine — Kern-Logik fuer alle Analysemodi.

Wird von run.py aufgerufen. Nicht direkt ausfuehren.
"""

import re
import shutil
import subprocess
import sys
import time
import numpy as np
from itertools import product
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    import yaml
except ImportError:
    print("FEHLER: PyYAML ist nicht installiert.")
    print("  pip install pyyaml")
    sys.exit(1)


# ============================================================
# Config loading
# ============================================================

def load_config(config_path):
    """Load and validate a YAML config file."""
    config_path = Path(config_path).resolve()
    base_dir = config_path.parent.parent  # configs/ -> project root

    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    cfg['_config_path'] = config_path
    cfg['_base_dir'] = base_dir
    cfg['_output_dir'] = base_dir / 'output' / cfg['name']
    cfg['_figures_dir'] = base_dir / 'figures'

    return cfg


def _format_param(value, n_edges):
    """Format a parameter value for parameters.dat (space-separated per edge)."""
    if isinstance(value, list):
        return ' '.join(str(v) for v in value)
    else:
        return ' '.join([str(value)] * n_edges)


# ============================================================
# Template management
# ============================================================

def _find_template(cfg):
    """Find a restart template to base parameters.dat on."""
    base = cfg['_base_dir']

    # Look for existing restart files (best starting point)
    candidates = [
        base / 'MEA_analysis' / 'output_files_run05' / 'restart2002',
        base / 'templates' / 'restart_template',
    ]
    # Also search output dir for any restart files
    out_dir = cfg['_output_dir']
    if out_dir.exists():
        for f in sorted(out_dir.glob('restart*')):
            candidates.insert(0, f)

    for c in candidates:
        if c.exists():
            return c

    # No template found — generate a minimal parameters.dat
    return _generate_parameters_dat(cfg)


def _generate_parameters_dat(cfg):
    """Generate a complete parameters.dat from config (no template needed)."""
    base = cfg['_base_dir']
    n = len(cfg['edges'])
    p = cfg['parameters']
    fp = _format_param

    chi_paths = ' '.join(cfg['chi_files'][e] for e in cfg['edges'])

    lines = [
        'BEGIN',
        f"              Title                              {cfg['name']}",
        '###############################################################################',
        '########################## STRUCTURE_MODEL_DEFINITION #########################',
        '###############################################################################',
        '              Number_of_structures                   1',
        '              Scale_coordinates_by                   1',
        f"              Initial_displacement                   {p.get('Initial_displacement', 0.01)}",
        '              Weight                                 1',
        '              Doping_level                           -1',
        f"              File_with_structure                    {cfg['structure_file']}",
        '              Supercell                              1x1x1',
        '              Doping                                 none',
        '              pbc                                    1',
        '              labels_to_element_name                 1',
        '###############################################################################',
        '######################### STRUCTURE_MODEL_OPTIMIZATION ########################',
        '###############################################################################',
        f"              Froze_in                               {p['Froze_in']}",
        f"              Stop_after                             {p['Stop_after']}",
        f"              Number_of_states                       {p.get('Number_of_states', 16)}",
        f"              Maximal_step_length                    {p['Maximal_step_length']}",
        '              Weight_increment                       0.0001',
        '              Acceptance_rate                        0.8',
        '              Energy_weight                          -1',
        '              Move_type                              move_all',
        '              MT_seed                                291 564 837 1110',
        '              adjust_acceptance                      1',
        '              adjust_weights                         0',
        f"              Maximal_displacement                   {p['Maximal_displacement']}",
        '              Potential_type                         box',
        '###############################################################################',
        '############################### EXPERIMENTAL_DATA #############################',
        '###############################################################################',
        f"              Number_of_spectra                      {n}",
        f"              Columns_in_input                       {fp(2, n)}",
        f"              Column_for_k                           {fp(1, n)}",
        f"              Column_for_chi                         {fp(2, n)}",
        f"              S02                                    {fp(1, n)}",
        f"              dE0                                    {fp(0, n)}",
        f"              k_power_in_IO                          {fp(0, n)}",
        f"              Spectrum_weight                        {fp(1, n)}",
        f"              File_with_experimental_EXAFS_signal    {chi_paths}",
        '              adjust_spectra_weights                 1',
        '###############################################################################',
        '############################### THEORETICAL_DATA ##############################',
        '###############################################################################',
        '              Update_basis_every_...                 2',
        '              Maximal_number_of_clusters             10000',
        '              Expansion_accuracy                     -1',
        '              FEFF_path                              feff',
        f"              N_legs                                 {fp(p['N_legs'], n)}",
        f"              R_max_for_FEFF                         {fp(p['R_max_for_FEFF'], n)}",
        f"              FEFF_criteria                          {fp(1, n)}",
        f"              Pathfinder_robustness                  {fp(0.5, n)}",
        f"              Clustering_precision                   {fp(0.001, n)}",
        f"              Edge                                   {' '.join(cfg['edges'])}",
        f"              Edge_type                              {' '.join(cfg['edge_types'])}",
        f"              Polarization                           {' '.join(['none'] * n)}",
        f"              r_scf                                  {fp(4.3356, n)}",
        f"              ca                                     {fp(0.1, n)}",
        f"              n_scf                                  {fp(30, n)}",
        f"              n_mix                                  {fp(1, n)}",
        '###############################################################################',
        '####################### THEORY_AND_EXPERIMENT_COMPARISON ######################',
        '###############################################################################',
        f"              Space                                  {p['Space']}",
        f"              WT_type                                {fp(1, n)}",
        f"              k_min                                  {fp(p['k_min'], n)}",
        f"              k_max                                  {fp(p['k_max'], n)}",
        f"              dk                                     {fp(0.05, n)}",
        f"              k_power                                {fp(p['k_power'], n)}",
        f"              FT_window_type                         {fp(1, n)}",
        f"              R_min                                  {fp(0.5, n)}",
        f"              R_max                                  {fp(5, n)}",
        f"              R_max_for_1st_shell                    {fp(3, n)}",
        f"              R0_for_WT                              {fp(2.75, n)}",
        '###############################################################################',
        '#################################### OUTPUT ###################################',
        '###############################################################################',
        f"              File_for_calculated_EXAFS_signal       {' '.join(f'output_files/{e}/EXAFS.dat' for e in cfg['edges'])}",
        f"              File_for_FT                            {' '.join(f'output_files/{e}/expFT.dat' for e in cfg['edges'])}",
        f"              File_for_BFT                           {' '.join(f'output_files/{e}/expBFT.dat' for e in cfg['edges'])}",
        f"              File_for_WT                            {' '.join(f'output_files/{e}/expWT.dat' for e in cfg['edges'])}",
        f"              File_for_interpolated_EXAFS_signal     {' '.join(f'output_files/{e}/ipolEXAFS.dat' for e in cfg['edges'])}",
        '              Convergence                            output_files/output.dat',
        '              xyz_file                                output_files/final.xyz',
        '              Save_state                             output_files/restart',
        '              Save_state_every                       200',
        '              Total_iterations:                         0',
        '              Residuals:',
    ]

    template_path = base / 'templates' / 'generated_parameters.dat'
    template_path.parent.mkdir(parents=True, exist_ok=True)
    template_path.write_text('\n'.join(lines), encoding='latin-1')
    return template_path


# ============================================================
# Work directory setup
# ============================================================

def _setup_workdir(cfg, work_dir, structure_file=None):
    """Set up an isolated EvAX working directory."""
    base = cfg['_base_dir']
    work_dir.mkdir(parents=True, exist_ok=True)

    # Symlink EvAX binary
    evax_link = work_dir / 'EvAX.exe'
    if not evax_link.exists():
        evax_link.symlink_to(base / cfg['software']['evax'])

    # Symlink feff binary
    feff_link = work_dir / 'feff'
    if not feff_link.exists():
        feff_link.symlink_to(base / cfg['software']['feff'])

    # Copy run-FEFF.exe wrapper
    run_feff = work_dir / 'run-FEFF.exe'
    if not run_feff.exists():
        run_feff.write_text('#!/bin/bash\ncd $1\n./feff > feff.log\n')
        run_feff.chmod(0o755)

    # Symlink input data directories
    # We need chi files and structure files accessible from work_dir
    for subdir in ['data', 'structures']:
        link = work_dir / subdir
        if not link.exists():
            target = base / subdir
            if target.exists():
                link.symlink_to(target)

    # Also link input_files if it exists (legacy compatibility)
    inp_link = work_dir / 'input_files'
    if not inp_link.exists():
        legacy = base / 'MEA_analysis' / 'input_files'
        if legacy.exists():
            inp_link.symlink_to(legacy)

    (work_dir / 'output_files').mkdir(exist_ok=True)
    return work_dir


def _make_parameters_dat(cfg, work_dir, overrides=None, structure_file=None):
    """Create parameters.dat in work_dir from template + config + overrides."""
    template = _find_template(cfg)
    text = template.read_text(encoding='latin-1')
    n = len(cfg['edges'])
    p = cfg['parameters'].copy()

    if overrides:
        p.update(overrides)

    def replace_param(text, param, value):
        pattern = rf'(              {param}\s+).*'
        return re.sub(pattern, rf'\g<1>{value}', text)

    text = replace_param(text, 'Title', cfg['name'])
    if structure_file:
        text = replace_param(text, 'File_with_structure', structure_file)
    else:
        text = replace_param(text, 'File_with_structure', cfg['structure_file'])

    text = replace_param(text, 'S02', _format_param(1, n))
    text = replace_param(text, 'dE0', _format_param(0, n))
    text = replace_param(text, 'Stop_after', str(p['Stop_after']))
    text = replace_param(text, 'Froze_in', str(p['Froze_in']))
    text = replace_param(text, 'k_min', _format_param(p['k_min'], n))
    text = replace_param(text, 'k_max', _format_param(p['k_max'], n))
    text = replace_param(text, 'k_power', _format_param(p['k_power'], n))
    text = replace_param(text, 'Space', str(p['Space']))
    text = replace_param(text, 'N_legs', _format_param(p['N_legs'], n))
    text = replace_param(text, 'Maximal_step_length', str(p['Maximal_step_length']))
    text = replace_param(text, 'Maximal_displacement', str(p['Maximal_displacement']))
    text = replace_param(text, 'R_max_for_FEFF', _format_param(p['R_max_for_FEFF'], n))

    # Reset iteration counter — fresh start
    text = re.sub(r'Total_iterations:\s+\d+',
                  'Total_iterations:                         0', text)
    text = re.sub(r'Residuals:.*', 'Residuals:', text)

    # Remove old ATOMS section — force fresh start
    atoms_idx = text.find('\nATOMS\n')
    if atoms_idx > 0:
        text = text[:atoms_idx + 1]

    param_file = work_dir / 'parameters.dat'
    param_file.write_text(text, encoding='latin-1')
    return param_file


# ============================================================
# EvAX execution
# ============================================================

def _run_evax_in_dir(work_dir, timeout=86400):
    """Run EvAX in the given directory. Returns (log_path, elapsed)."""
    t0 = time.time()
    log_file = work_dir / 'evax_run.log'
    with open(log_file, 'wb') as lf:
        subprocess.run(
            ['./EvAX.exe', 'parameters.dat', '-autoES'],
            cwd=str(work_dir),
            stdout=lf,
            stderr=subprocess.STDOUT,
            timeout=timeout,
        )
    elapsed = time.time() - t0
    return log_file, elapsed


def _collect_results(work_dir, result_dir):
    """Copy output_files from work_dir to result_dir."""
    result_dir.mkdir(parents=True, exist_ok=True)
    out_dir = work_dir / 'output_files'
    if out_dir.exists():
        for item in out_dir.iterdir():
            dst = result_dir / item.name
            if item.is_dir():
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(item, dst)
            else:
                shutil.copy2(item, dst)


def _parse_output(result_dir, edges):
    """Parse EvAX output for residuals and fit parameters."""
    results = {'total': None}
    for e in edges:
        results[e] = None
        results[f's02_{e}'] = None
        results[f'dE0_{e}'] = None

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
        n = len(edges)
        if len(precise) >= n:
            for i, e in enumerate(edges):
                results[e] = float(precise[-(n - i)])
        s02s = re.findall(r'S02\s*=\s*([\d.]+)', log)
        if len(s02s) >= n:
            for i, e in enumerate(edges):
                results[f's02_{e}'] = float(s02s[-(n - i)])
        de0s = re.findall(r'dE0\s*=\s*([+-]?[\d.]+)', log)
        if len(de0s) >= n:
            for i, e in enumerate(edges):
                results[f'dE0_{e}'] = float(de0s[-(n - i)])

    return results


# ============================================================
# Mode: Single production run
# ============================================================

def run_single(cfg, dry_run=False):
    """Run a single EvAX analysis with the config parameters."""
    out_dir = cfg['_output_dir'] / 'production'
    work_dir = out_dir / 'work'

    print(f"Modus: Einzellauf (Produktion)")
    print(f"  Struktur: {cfg['structure_file']}")
    print(f"  Screening: {cfg['parameters']['Stop_after']} Iterationen")
    print(f"  Froze_in:  {cfg['parameters']['Froze_in']}")
    print(f"  Ausgabe:   {out_dir}\n")

    if dry_run:
        print("(Trockenlauf — nichts wird ausgefuehrt)")
        return

    _setup_workdir(cfg, work_dir)
    _make_parameters_dat(cfg, work_dir)

    print("Starte EvAX...")
    log_file, elapsed = _run_evax_in_dir(work_dir)
    shutil.copy2(log_file, out_dir / 'evax.log')
    _collect_results(work_dir, out_dir)
    shutil.rmtree(work_dir, ignore_errors=True)

    r = _parse_output(out_dir, cfg['edges'])
    print(f"\nFertig! ({elapsed / 60:.1f} min)")
    print(f"Residual: {r.get('total', '?')}")
    for e in cfg['edges']:
        print(f"  {e}: {r.get(e, '?')}")


# ============================================================
# Mode: Parameter scan
# ============================================================

def _build_scan_jobs(cfg):
    """Build parameter scan job list from config."""
    scan_cfg = cfg.get('parameter_scan', {})
    p1 = scan_cfg.get('phase1', {})
    p2 = scan_cfg.get('phase2', {})
    n = len(cfg['edges'])
    jobs = []

    # Phase 1: Fit parameters
    if p1:
        grid_keys = sorted(p1.keys())
        grid_vals = [p1[k] for k in grid_keys]
        for combo in product(*grid_vals):
            params = dict(zip(grid_keys, combo))
            # Build name
            parts = []
            for k, v in sorted(params.items()):
                sv = str(v).replace(' ', '').replace('[', '').replace(']', '').replace(',', '-')
                parts.append(f"{k}{sv}")
            name = '_'.join(parts)
            # Use defaults from main config for non-scanned params
            overrides = {
                'Stop_after': scan_cfg.get('screening_iterations', 50),
                'Froze_in': scan_cfg.get('froze_in', 50),
            }
            overrides.update(params)
            jobs.append((f'p1_{name}', overrides))

    # Phase 2: Optimization parameters (with best Phase 1 settings)
    if p2:
        # Use the main config parameters as the "best" fit params
        grid_keys = sorted(p2.keys())
        grid_vals = [p2[k] for k in grid_keys]
        for combo in product(*grid_vals):
            params = dict(zip(grid_keys, combo))
            parts = []
            for k, v in sorted(params.items()):
                sv = str(v).replace(' ', '').replace('[', '').replace(']', '').replace(',', '-')
                parts.append(f"{k}{sv}")
            name = '_'.join(parts)
            overrides = {
                'Stop_after': scan_cfg.get('screening_iterations', 50),
                'Froze_in': scan_cfg.get('froze_in', 50),
            }
            overrides.update(params)
            jobs.append((f'p2_{name}', overrides))

    return jobs


def _run_scan_job(args):
    """Worker function for parallel scan execution."""
    cfg_path, job_name, overrides = args
    cfg = load_config(cfg_path)
    scan_dir = cfg['_output_dir'] / 'parameter_scan' / job_name
    result_dir = scan_dir

    # Skip if already done
    if (result_dir / 'output.dat').exists():
        lines = (result_dir / 'output.dat').read_text().strip().split('\n')
        if len(lines) > 5:
            r = _parse_output(result_dir, cfg['edges'])
            return job_name, 'skipped', 0, r.get('total')

    work_dir = scan_dir / 'work'
    _setup_workdir(cfg, work_dir)
    _make_parameters_dat(cfg, work_dir, overrides=overrides)

    try:
        log_file, elapsed = _run_evax_in_dir(work_dir, timeout=7200)
        shutil.copy2(log_file, result_dir / 'evax.log')
        _collect_results(work_dir, result_dir)
        shutil.rmtree(work_dir, ignore_errors=True)

        r = _parse_output(result_dir, cfg['edges'])
        return job_name, 'done', elapsed, r.get('total')

    except subprocess.TimeoutExpired:
        return job_name, 'timeout', 7200, None
    except Exception as e:
        return job_name, f'error: {e}', 0, None


def run_parameter_scan(cfg, n_parallel=4, dry_run=False):
    """Run automated parameter scan."""
    jobs = _build_scan_jobs(cfg)
    scan_dir = cfg['_output_dir'] / 'parameter_scan'
    scan_dir.mkdir(parents=True, exist_ok=True)

    scan_cfg = cfg.get('parameter_scan', {})
    print(f"Modus: Parameterstudie")
    print(f"  {len(jobs)} Kombinationen, {n_parallel} parallele Instanzen")
    print(f"  Screening: {scan_cfg.get('screening_iterations', 50)} Iterationen")
    print(f"  Ausgabe: {scan_dir}\n")

    if dry_run:
        for i, (name, params) in enumerate(jobs):
            print(f"  [{i+1:>3}/{len(jobs)}] {name}")
        print(f"\n(Trockenlauf — nichts wird ausgefuehrt)")
        return

    t_start = time.time()
    completed = 0
    config_path = cfg['_config_path']

    with ProcessPoolExecutor(max_workers=n_parallel) as executor:
        futures = {
            executor.submit(_run_scan_job, (config_path, name, overrides)): name
            for name, overrides in jobs
        }
        for future in as_completed(futures):
            try:
                name, status, elapsed, residual = future.result()
                completed += 1
                res_str = f"R={residual:.4f}" if residual else status
                print(f"  [{completed}/{len(jobs)}] {name}: {res_str} ({elapsed/60:.1f} min)")
            except Exception as e:
                completed += 1
                print(f"  [{completed}/{len(jobs)}] {futures[future]}: FEHLER {e}")

    total_time = (time.time() - t_start) / 60
    print(f"\nFertig! {total_time:.1f} min ({total_time/60:.1f} h)")
    show_results(cfg, mode='parameter_scan')


# ============================================================
# Mode: Structure comparison
# ============================================================

def _run_structure_job(args):
    """Worker function for parallel structure comparison."""
    cfg_path, struct_name, struct_file = args
    cfg = load_config(cfg_path)
    comp_dir = cfg['_output_dir'] / 'structure_comparison' / struct_name
    work_dir = comp_dir / 'work'

    # Skip if already done
    if (comp_dir / 'output.dat').exists():
        lines = (comp_dir / 'output.dat').read_text().strip().split('\n')
        if len(lines) > 5:
            r = _parse_output(comp_dir, cfg['edges'])
            return struct_name, 'skipped', 0, r.get('total')

    _setup_workdir(cfg, work_dir)

    sc_cfg = cfg.get('structure_compare', {})
    overrides = {
        'Stop_after': sc_cfg.get('screening_iterations', 100),
        'Froze_in': sc_cfg.get('froze_in', 150),
    }
    _make_parameters_dat(cfg, work_dir, overrides=overrides,
                         structure_file=struct_file)

    try:
        log_file, elapsed = _run_evax_in_dir(work_dir, timeout=86400)
        shutil.copy2(log_file, comp_dir / 'evax.log')
        _collect_results(work_dir, comp_dir)
        shutil.rmtree(work_dir, ignore_errors=True)

        r = _parse_output(comp_dir, cfg['edges'])
        return struct_name, 'done', elapsed, r.get('total')

    except subprocess.TimeoutExpired:
        return struct_name, 'timeout', 86400, None
    except Exception as e:
        return struct_name, f'error: {e}', 0, None


def run_structure_comparison(cfg, n_parallel=3, dry_run=False):
    """Run structure comparison (FCC vs BCC vs HCP etc.)."""
    structs = cfg.get('structure_comparison', {})
    if not structs:
        print("FEHLER: Keine Strukturen zum Vergleichen in der Config definiert.")
        print("Fuege 'structure_comparison' in die YAML-Datei ein.")
        sys.exit(1)

    comp_dir = cfg['_output_dir'] / 'structure_comparison'
    comp_dir.mkdir(parents=True, exist_ok=True)

    sc_cfg = cfg.get('structure_compare', {})
    print(f"Modus: Strukturvergleich")
    print(f"  Strukturen: {', '.join(s.get('label', k) for k, s in structs.items())}")
    print(f"  Screening: {sc_cfg.get('screening_iterations', 100)} Iterationen")
    print(f"  Froze_in:  {sc_cfg.get('froze_in', 150)}")
    print(f"  Ausgabe:   {comp_dir}\n")

    jobs = [(name, info['file']) for name, info in structs.items()]

    if dry_run:
        for name, f in jobs:
            label = structs[name].get('label', name)
            print(f"  {label}: {f}")
        print(f"\n(Trockenlauf — nichts wird ausgefuehrt)")
        return

    t_start = time.time()
    config_path = cfg['_config_path']

    with ProcessPoolExecutor(max_workers=n_parallel) as executor:
        futures = {
            executor.submit(_run_structure_job, (config_path, name, sf)): name
            for name, sf in jobs
        }
        for future in as_completed(futures):
            try:
                name, status, elapsed, residual = future.result()
                label = structs[name].get('label', name)
                res_str = f"R={residual:.4f}" if residual else status
                print(f"  {label}: {res_str} ({elapsed/60:.1f} min)")
            except Exception as e:
                print(f"  {futures[future]}: FEHLER {e}")

    total_time = (time.time() - t_start) / 60
    print(f"\nFertig! {total_time:.1f} min ({total_time/60:.1f} h)")
    show_results(cfg, mode='structure_compare')


# ============================================================
# Results display
# ============================================================

def show_results(cfg, mode='single'):
    """Display results table."""
    edges = cfg['edges']

    if mode == 'parameter_scan':
        scan_dir = cfg['_output_dir'] / 'parameter_scan'
        if not scan_dir.exists():
            print("Keine Ergebnisse gefunden.")
            return

        rows = []
        for d in sorted(scan_dir.iterdir()):
            if not d.is_dir() or d.name == 'work':
                continue
            r = _parse_output(d, edges)
            if r['total'] is not None:
                rows.append({'name': d.name, **r})

        if not rows:
            print("Keine abgeschlossenen Runs gefunden.")
            return

        rows.sort(key=lambda x: x.get('total', 999))

        edge_headers = ''.join(f'{e:>8}' for e in edges)
        print(f"\n{'='*90}")
        print(f"{'#':>3} {'Run':<40} {edge_headers} {'TOTAL':>8}")
        print(f"{'-'*90}")
        for i, r in enumerate(rows):
            edge_vals = ''.join(f'{r.get(e, 0) or 0:>8.4f}' for e in edges)
            marker = '  *** BEST' if i == 0 else ''
            print(f"{i+1:>3} {r['name']:<40} {edge_vals} "
                  f"{r.get('total', 0):>8.4f}{marker}")
        print(f"{'='*90}")
        print(f"  {len(rows)} Runs abgeschlossen")

    elif mode == 'structure_compare':
        comp_dir = cfg['_output_dir'] / 'structure_comparison'
        structs = cfg.get('structure_comparison', {})
        if not comp_dir.exists():
            print("Keine Ergebnisse gefunden.")
            return

        edge_headers = ''.join(f'{e:>8}' for e in edges)
        print(f"\n{'='*80}")
        print(f"{'Struktur':<12} {edge_headers} {'TOTAL':>8}")
        print(f"{'-'*80}")

        best_name, best_total = None, 999
        for name in structs:
            d = comp_dir / name
            if not d.exists():
                continue
            r = _parse_output(d, edges)
            label = structs[name].get('label', name.upper())
            edge_vals = ''.join(f'{r.get(e, 0) or 0:>8.4f}' for e in edges)
            tot = r.get('total', 0) or 0
            print(f"{label:<12} {edge_vals} {tot:>8.4f}")
            if tot < best_total and tot > 0:
                best_name, best_total = label, tot

        if best_name:
            print(f"\n  Beste Struktur: {best_name} (R = {best_total:.4f})")
        print(f"{'='*80}")

    else:
        out_dir = cfg['_output_dir'] / 'production'
        if not out_dir.exists():
            print("Keine Ergebnisse gefunden.")
            return
        r = _parse_output(out_dir, edges)
        print(f"\nResidual: {r.get('total', '?')}")
        for e in edges:
            print(f"  {e}: {r.get(e, '?')}")


# ============================================================
# Plotting
# ============================================================

def _compute_ft(k, chi, k_min=3, k_max=11.5, k_power=2, dr=0.02, r_max=6):
    """Compute FT magnitude from chi(k) with Hanning window."""
    mask = (k >= k_min) & (k <= k_max)
    k_w = k[mask]
    chi_w = chi[mask] * k_w**k_power

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
    dk = k_w[1] - k_w[0] if len(k_w) > 1 else 0.05
    ft_re = np.array([np.sum(chi_w * np.cos(2 * k_w * ri)) * dk / np.sqrt(np.pi) for ri in r])
    ft_im = np.array([np.sum(chi_w * np.sin(2 * k_w * ri)) * dk / np.sqrt(np.pi) for ri in r])

    return r, np.sqrt(ft_re**2 + ft_im**2)


def generate_plots(cfg, mode='single'):
    """Generate plots for the specified mode."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        'font.family': 'serif', 'font.size': 9,
        'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    })

    fig_dir = cfg['_figures_dir']
    fig_dir.mkdir(parents=True, exist_ok=True)
    edges = cfg['edges']

    EDGE_COLORS = {}
    palette = ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#56B4E9', '#E69F00']
    for i, e in enumerate(edges):
        EDGE_COLORS[e] = palette[i % len(palette)]

    if mode == 'structure_compare':
        _plot_structure_comparison(cfg, fig_dir, edges, EDGE_COLORS, plt)
    elif mode == 'parameter_scan':
        _plot_parameter_scan(cfg, fig_dir, plt)
    else:
        _plot_single_run(cfg, fig_dir, edges, EDGE_COLORS, plt)


def _plot_single_run(cfg, fig_dir, edges, colors, plt):
    """Plot results from a single production run."""
    out_dir = cfg['_output_dir'] / 'production'
    prefix = cfg['name']

    # Convergence
    out_file = out_dir / 'output.dat'
    if out_file.exists():
        data = np.loadtxt(out_file)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(data[:, 0], data[:, 1], 'k-', lw=1)
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Residual')
        ax.set_title(f'{prefix} — Konvergenz')
        ax.set_ylim(bottom=0)
        fig.savefig(fig_dir / f'{prefix}_convergence.png')
        plt.close()

    # EXAFS + FT
    n = len(edges)
    fig, axes = plt.subplots(2, n, figsize=(4 * n, 7))
    if n == 1:
        axes = axes.reshape(2, 1)
    fig.suptitle(f'{prefix} — EXAFS Fit', fontsize=13, fontweight='bold')

    for i, edge in enumerate(edges):
        exafs_file = out_dir / edge / 'EXAFS.dat'
        if not exafs_file.exists():
            continue
        data = np.loadtxt(exafs_file)
        k = data[:, 0]

        ax = axes[0, i]
        ax.plot(k, data[:, 1] * k**2, color=colors[edge], lw=1.2, label='Calc.')
        if data.shape[1] >= 3:
            ax.plot(k, data[:, 2] * k**2, 'k--', lw=0.8, label='Exp.')
        ax.set_xlabel(r'k ($\AA^{-1}$)')
        ax.set_ylabel(r'$\chi(k) \cdot k^2$')
        ax.set_title(f'{edge}')
        ax.legend(fontsize=7)

        ax = axes[1, i]
        r_ft, ft_calc = _compute_ft(k, data[:, 1])
        ax.plot(r_ft, ft_calc, color=colors[edge], lw=1.2, label='Calc.')
        if data.shape[1] >= 3:
            r_ft_e, ft_exp = _compute_ft(k, data[:, 2])
            ax.plot(r_ft_e, ft_exp, 'k--', lw=0.8, label='Exp.')
        ax.set_xlabel(r'R ($\AA$)')
        ax.set_ylabel('|FT|')
        ax.set_title(f'{edge} FT')
        ax.legend(fontsize=7)
        ax.set_xlim(0, 6)

    fig.tight_layout()
    fig.savefig(fig_dir / f'{prefix}_exafs_ft.png')
    plt.close()
    print(f"Plots gespeichert in {fig_dir}/")


def _plot_structure_comparison(cfg, fig_dir, edges, edge_colors, plt):
    """Plot structure comparison results."""
    from matplotlib.gridspec import GridSpec

    structs = cfg.get('structure_comparison', {})
    comp_dir = cfg['_output_dir'] / 'structure_comparison'
    prefix = cfg['name']

    STRUCT_COLORS = {}
    s_palette = ['#2ecc71', '#3498db', '#e74c3c', '#9b59b6', '#f39c12']
    for i, name in enumerate(structs):
        STRUCT_COLORS[name] = s_palette[i % len(s_palette)]

    # Check which structures have results
    available = []
    results = {}
    for name in structs:
        d = comp_dir / name
        if d.exists() and (d / 'output.dat').exists():
            r = _parse_output(d, edges)
            if r['total'] is not None:
                results[name] = r
                available.append(name)

    if len(available) < 2:
        print(f"Nur {len(available)} Strukturen fertig — mind. 2 noetig fuer Vergleich.")
        return

    # Individual plots per structure
    for name in available:
        label = structs[name].get('label', name.upper())
        d = comp_dir / name

        # EXAFS + FT grid
        n = len(edges)
        fig, axes = plt.subplots(2, n, figsize=(4 * n, 7))
        if n == 1:
            axes = axes.reshape(2, 1)
        fig.suptitle(f'{label} — EXAFS Fit', fontsize=13, fontweight='bold')

        for i, edge in enumerate(edges):
            exafs_file = d / edge / 'EXAFS.dat'
            if not exafs_file.exists():
                continue
            data = np.loadtxt(exafs_file)
            k = data[:, 0]

            ax = axes[0, i]
            ax.plot(k, data[:, 1] * k**2, color=edge_colors[edge], lw=1.2, label='Calc.')
            if data.shape[1] >= 3:
                ax.plot(k, data[:, 2] * k**2, 'k--', lw=0.8, label='Exp.')
            ax.set_xlabel(r'k ($\AA^{-1}$)')
            ax.set_ylabel(r'$\chi(k) \cdot k^2$')
            ax.set_title(f'{edge}')
            ax.legend(fontsize=7)

            ax = axes[1, i]
            r_ft, ft_calc = _compute_ft(k, data[:, 1])
            ax.plot(r_ft, ft_calc, color=edge_colors[edge], lw=1.2, label='Calc.')
            if data.shape[1] >= 3:
                r_ft_e, ft_exp = _compute_ft(k, data[:, 2])
                ax.plot(r_ft_e, ft_exp, 'k--', lw=0.8, label='Exp.')
            ax.set_xlabel(r'R ($\AA$)')
            ax.set_ylabel('|FT|')
            ax.legend(fontsize=7)
            ax.set_xlim(0, 6)

        fig.tight_layout()
        fig.savefig(fig_dir / f'{prefix}_struct_{name}_exafs_ft.png')
        fig.savefig(fig_dir / f'{prefix}_struct_{name}_exafs_ft.pdf')
        plt.close()

    # Combined comparison plot
    n = len(edges)
    fig = plt.figure(figsize=(7.2, 9.5))
    gs = GridSpec(4, n, figure=fig, height_ratios=[1.0, 1.2, 1.2, 0.8],
                  hspace=0.45, wspace=0.35)

    # Row 0: Residual bars
    ax_res = fig.add_subplot(gs[0, :])
    x = np.arange(len(available))
    w = 0.18
    for i, edge in enumerate(edges):
        vals = [results[s].get(edge, 0) or 0 for s in available]
        ax_res.bar(x + (i - len(edges)/2 + 0.5) * w, vals, w,
                   label=f'{edge}', color=edge_colors[edge],
                   edgecolor='black', linewidth=0.4)
    totals = [results[s].get('total', 0) or 0 for s in available]
    ax_res.scatter(x, totals, marker='D', s=80, c='black', zorder=5,
                   label='Total')
    ax_res.set_xticks(x)
    ax_res.set_xticklabels([structs[s].get('label', s.upper()) for s in available],
                           fontsize=11, fontweight='bold')
    ax_res.set_ylabel('R-factor')
    ax_res.set_title('(a) Fit-Qualitaet', fontsize=10, fontweight='bold', loc='left')
    ax_res.legend(fontsize=7, ncol=len(edges)+1, loc='upper right')
    best_idx = np.argmin(totals)
    for i, t in enumerate(totals):
        ax_res.text(i, t + 0.01, f'{t:.3f}', ha='center', va='bottom', fontsize=8,
                    fontweight='bold' if i == best_idx else 'normal',
                    color='green' if i == best_idx else 'black')

    # Row 1+2: EXAFS + FT per edge
    for col, edge in enumerate(edges):
        ax = fig.add_subplot(gs[1, col])
        for s in available:
            exafs_file = comp_dir / s / edge / 'EXAFS.dat'
            if exafs_file.exists():
                data = np.loadtxt(exafs_file)
                k = data[:, 0]
                ax.plot(k, data[:, 1] * k**2, color=STRUCT_COLORS[s],
                        lw=1.0, label=structs[s].get('label', s.upper()))
        # Experimental from first available
        for s in available:
            exafs_file = comp_dir / s / edge / 'EXAFS.dat'
            if exafs_file.exists():
                data = np.loadtxt(exafs_file)
                if data.shape[1] >= 3:
                    ax.plot(data[:, 0], data[:, 2] * data[:, 0]**2,
                            'k-', lw=0.6, alpha=0.4, label='Exp.')
                break
        ax.set_xlabel(r'k ($\AA^{-1}$)')
        if col == 0:
            ax.set_ylabel(r'$\chi(k) \cdot k^2$')
        ax.set_title(f'{edge} EXAFS', fontsize=9, fontweight='bold', loc='left')
        ax.legend(fontsize=6)
        ax.set_xlim(2, 12)

        ax = fig.add_subplot(gs[2, col])
        for s in available:
            exafs_file = comp_dir / s / edge / 'EXAFS.dat'
            if exafs_file.exists():
                data = np.loadtxt(exafs_file)
                k = data[:, 0]
                r_ft, ft_mag = _compute_ft(k, data[:, 1])
                ax.plot(r_ft, ft_mag, color=STRUCT_COLORS[s],
                        lw=1.0, label=structs[s].get('label', s.upper()))
        for s in available:
            exafs_file = comp_dir / s / edge / 'EXAFS.dat'
            if exafs_file.exists():
                data = np.loadtxt(exafs_file)
                if data.shape[1] >= 3:
                    r_ft_e, ft_exp = _compute_ft(k, data[:, 2])
                    ax.plot(r_ft_e, ft_exp, 'k-', lw=0.6, alpha=0.4, label='Exp.')
                break
        ax.set_xlabel(r'R ($\AA$)')
        if col == 0:
            ax.set_ylabel(r'|FT($\chi(k) \cdot k^2$)|')
        ax.set_title(f'{edge} FT', fontsize=9, fontweight='bold', loc='left')
        ax.legend(fontsize=6)
        ax.set_xlim(0, 6)

    # Row 3: Parameter table
    ax_tab = fig.add_subplot(gs[3, :])
    ax_tab.axis('off')
    col_labels = ['Struktur', 'R-factor'] + [f'R({e})' for e in edges]
    cell_text = []
    cell_colors = []
    min_total = min(totals)
    for s in available:
        r = results[s]
        row = [structs[s].get('label', s.upper()), f"{r.get('total', 0) or 0:.4f}"]
        row += [f"{r.get(e, 0) or 0:.4f}" for e in edges]
        cell_text.append(row)
        if (r.get('total', 0) or 0) == min_total:
            cell_colors.append(['#d4edda'] * len(col_labels))
        else:
            cell_colors.append(['white'] * len(col_labels))

    table = ax_tab.table(cellText=cell_text, colLabels=col_labels,
                         cellColours=cell_colors, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.4)
    for j in range(len(col_labels)):
        table[0, j].set_text_props(fontweight='bold')

    fig.savefig(fig_dir / f'{prefix}_structure_comparison.png')
    fig.savefig(fig_dir / f'{prefix}_structure_comparison.pdf')
    plt.close()
    print(f"Vergleichsplot gespeichert: {prefix}_structure_comparison.png")


def _plot_parameter_scan(cfg, fig_dir, plt):
    """Plot parameter scan results."""
    scan_dir = cfg['_output_dir'] / 'parameter_scan'
    prefix = cfg['name']
    edges = cfg['edges']

    rows = []
    for d in sorted(scan_dir.iterdir()):
        if not d.is_dir():
            continue
        r = _parse_output(d, edges)
        if r['total'] is not None:
            rows.append({'name': d.name, **r})

    if not rows:
        print("Keine Ergebnisse zum Plotten.")
        return

    rows.sort(key=lambda x: x['total'])
    totals = [r['total'] for r in rows]

    fig, ax = plt.subplots(figsize=(12, max(5, len(rows) * 0.25)))
    names = [r['name'][:40] for r in rows]
    colors = ['#2ecc71' if t == min(totals) else '#3498db' if t < np.median(totals) else '#e74c3c'
              for t in totals]
    ax.barh(range(len(rows)), totals, color=colors, height=0.7)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(names, fontsize=6)
    ax.set_xlabel('Total Residual')
    ax.set_title(f'{prefix} — Parameterstudie')
    ax.invert_yaxis()
    fig.savefig(fig_dir / f'{prefix}_parameter_scan.png')
    plt.close()

    # CSV
    csv_file = scan_dir / 'results.csv'
    with open(csv_file, 'w') as f:
        f.write('rank,run,' + ','.join(f'residual_{e}' for e in edges) + ',residual_total\n')
        for i, r in enumerate(rows):
            vals = ','.join(f"{r.get(e, 0):.6f}" for e in edges)
            f.write(f"{i+1},{r['name']},{vals},{r['total']:.6f}\n")
    print(f"CSV: {csv_file}")
    print(f"Plot: {prefix}_parameter_scan.png")
