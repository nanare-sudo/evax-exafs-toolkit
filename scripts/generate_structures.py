#!/usr/bin/env python3
"""
Generate POSCAR-format starting structures for EvAX analysis.

Erzeugt Superzellen mit ~256 Atomen fuer verschiedene Kristallstrukturen
(FCC, BCC, HCP, L10) mit zufaelliger Elementverteilung.

Verwendung:
    # FCC-Struktur fuer CoPt (50/50), a=3.80 Angstrom:
    python scripts/generate_structures.py --elements Co Pt --counts 128 128 \
        --structure fcc --lattice 3.80 --output structures/CoPt_fcc.p1

    # BCC-Struktur fuer CoNiCu (33/33/34):
    python scripts/generate_structures.py --elements Co Ni Cu --counts 85 85 86 \
        --structure bcc --lattice 2.87 --output structures/MEA_bcc.p1

    # L10-Struktur fuer CoPt (geordnet, abwechselnde Schichten):
    python scripts/generate_structures.py --elements Co Pt --counts 128 128 \
        --structure l10 --lattice 3.80 --c-over-a 0.98 --output structures/CoPt_l10.p1

    # Alle Strukturen fuer CoPt auf einmal:
    python scripts/generate_structures.py --preset CoPt
"""

import argparse
import numpy as np
from pathlib import Path


def write_poscar(filename, comment, lattice_vectors, frac_coords, labels, elements):
    """Write a POSCAR-format .p1 file."""
    counts = {e: labels.count(e) for e in elements}

    with open(filename, 'w') as f:
        f.write(f'{comment}\n')
        f.write('1.0\n')
        for v in lattice_vectors:
            f.write(f'        {v[0]:.10f}         {v[1]:.10f}         {v[2]:.10f}\n')
        f.write('   ' + '    '.join(elements) + '\n')
        f.write('   ' + '    '.join(str(counts[e]) for e in elements) + '\n')
        f.write('Direct\n')
        for coord, label in zip(frac_coords, labels):
            f.write(f'   {coord[0]:.9f}   {coord[1]:.9f}   {coord[2]:.9f}   {label}\n')

    total = len(labels)
    count_str = ', '.join(f'{counts[e]} {e}' for e in elements)
    print(f'  {filename}: {total} Atome ({count_str})')


def random_labels(elements, counts, seed=42):
    """Create shuffled element labels."""
    labels = []
    for e, n in zip(elements, counts):
        labels.extend([e] * n)
    rng = np.random.default_rng(seed)
    rng.shuffle(labels)
    return labels


def generate_fcc(elements, counts, a, output, seed=42):
    """FCC: 4 atoms/cell, 4x4x4 = 256 atoms."""
    n_total = sum(counts)
    # Find supercell dimensions
    # FCC: 4 atoms per unit cell
    n_cells = n_total // 4
    # Try to make it as cubic as possible
    nx = ny = nz = round(n_cells ** (1/3))
    while nx * ny * nz * 4 != n_total:
        nz += 1
        if nx * ny * nz * 4 > n_total * 2:
            raise ValueError(f'Cannot make FCC supercell with exactly {n_total} atoms. '
                             f'FCC needs multiples of 4. Try {(n_total // 4) * 4}.')

    lattice = [[nx * a, 0, 0], [0, ny * a, 0], [0, 0, nz * a]]
    basis = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]

    coords = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for bx, by, bz in basis:
                    coords.append(((ix + bx) / nx, (iy + by) / ny, (iz + bz) / nz))

    labels = random_labels(elements, counts, seed)
    write_poscar(output, f'{"-".join(elements)} FCC a={a}', lattice, coords, labels, elements)


def generate_bcc(elements, counts, a, output, seed=42):
    """BCC: 2 atoms/cell."""
    n_total = sum(counts)
    n_cells = n_total // 2
    nx = ny = round(n_cells ** (1/3))
    nz = nx
    while nx * ny * nz * 2 < n_total:
        nz += 1
    if nx * ny * nz * 2 != n_total:
        # Try asymmetric
        for nx in range(2, 20):
            for ny in range(2, 20):
                nz = n_total // (2 * nx * ny)
                if nx * ny * nz * 2 == n_total:
                    break
            else:
                continue
            break

    if nx * ny * nz * 2 != n_total:
        raise ValueError(f'Cannot make BCC supercell with exactly {n_total} atoms.')

    lattice = [[nx * a, 0, 0], [0, ny * a, 0], [0, 0, nz * a]]
    basis = [(0, 0, 0), (0.5, 0.5, 0.5)]

    coords = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for bx, by, bz in basis:
                    coords.append(((ix + bx) / nx, (iy + by) / ny, (iz + bz) / nz))

    labels = random_labels(elements, counts, seed + 1)
    write_poscar(output, f'{"-".join(elements)} BCC a={a}', lattice, coords, labels, elements)


def generate_hcp(elements, counts, a, c_over_a, output, seed=42):
    """HCP: 2 atoms/cell, hexagonal lattice."""
    c = a * c_over_a
    n_total = sum(counts)
    n_cells = n_total // 2

    # Find supercell
    nx = ny = round(n_cells ** (1/3))
    nz = nx
    while nx * ny * nz * 2 < n_total:
        nz += 1
    if nx * ny * nz * 2 != n_total:
        for nx in range(2, 20):
            for ny in range(2, 20):
                nz = n_total // (2 * nx * ny)
                if nx * ny * nz * 2 == n_total:
                    break
            else:
                continue
            break

    if nx * ny * nz * 2 != n_total:
        raise ValueError(f'Cannot make HCP supercell with exactly {n_total} atoms.')

    a1 = np.array([a, 0, 0])
    a2 = np.array([-a / 2, a * np.sqrt(3) / 2, 0])
    a3 = np.array([0, 0, c])

    lattice = [(nx * a1).tolist(), (ny * a2).tolist(), (nz * a3).tolist()]
    basis = [(1/3, 2/3, 1/4), (2/3, 1/3, 3/4)]

    coords = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for bx, by, bz in basis:
                    coords.append(((ix + bx) / nx, (iy + by) / ny, (iz + bz) / nz))

    labels = random_labels(elements, counts, seed + 2)
    write_poscar(output, f'{"-".join(elements)} HCP a={a} c={c:.2f}',
                 lattice, coords, labels, elements)


def generate_l10(elements, counts, a, c_over_a, output, seed=42):
    """L10: ordered tetragonal structure (alternating layers).
    Based on FCC but with c/a != 1 and ordered element placement.
    2 atoms per unit cell. Elements alternate in layers along z."""
    if len(elements) != 2:
        raise ValueError('L10 structure requires exactly 2 elements.')

    c = a * c_over_a
    n_total = sum(counts)
    # L10: 2 atoms per cell (one of each element)
    # Basis: A at (0,0,0), B at (0.5,0.5,0.5)
    n_cells = n_total // 2

    nx = ny = round(n_cells ** (1/3))
    nz = nx
    while nx * ny * nz * 2 < n_total:
        nz += 1
    if nx * ny * nz * 2 != n_total:
        for nx in range(2, 20):
            for ny in range(2, 20):
                nz = n_total // (2 * nx * ny)
                if nx * ny * nz * 2 == n_total:
                    break
            else:
                continue
            break

    lattice = [[nx * a, 0, 0], [0, ny * a, 0], [0, 0, nz * c]]

    coords = []
    labels = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                # Site 1: corner → element A
                coords.append(((ix + 0) / nx, (iy + 0) / ny, (iz + 0) / nz))
                labels.append(elements[0])
                # Site 2: face center → element B
                coords.append(((ix + 0.5) / nx, (iy + 0.5) / ny, (iz + 0.5) / nz))
                labels.append(elements[1])

    write_poscar(output, f'{"-".join(elements)} L10 a={a} c/a={c_over_a}',
                 lattice, coords, labels, elements)


# ============================================================
# Presets for common probes
# ============================================================

PRESETS = {
    'MEA_CoNiCu': {
        'elements': ['Co', 'Ni', 'Cu'],
        'counts': [85, 85, 86],
        'structures': {
            'fcc': {'a': 3.54, 'file': 'structures/MEA_CoNiCu.p1'},
            'bcc': {'a': 2.87, 'file': 'structures/MEA_CoNiCu_bcc.p1'},
            'hcp': {'a': 2.51, 'c_over_a': 1.622, 'file': 'structures/MEA_CoNiCu_hcp.p1'},
        },
    },
    'CoPt': {
        'elements': ['Co', 'Pt'],
        'counts': [128, 128],
        'structures': {
            'fcc': {'a': 3.80, 'file': 'structures/CoPt_fcc.p1'},
            'bcc': {'a': 3.01, 'file': 'structures/CoPt_bcc.p1'},
            'hcp': {'a': 2.68, 'c_over_a': 1.622, 'file': 'structures/CoPt_hcp.p1'},
            'l10': {'a': 3.80, 'c_over_a': 0.98, 'file': 'structures/CoPt_l10.p1'},
        },
    },
}


def run_preset(name):
    """Generate all structures for a preset probe."""
    if name not in PRESETS:
        print(f'Unbekanntes Preset: {name}')
        print(f'Verfuegbar: {", ".join(PRESETS.keys())}')
        return

    p = PRESETS[name]
    print(f'Generiere Strukturen fuer {name} ({", ".join(p["elements"])}):')

    for stype, cfg in p['structures'].items():
        Path(cfg['file']).parent.mkdir(parents=True, exist_ok=True)
        if stype == 'fcc':
            generate_fcc(p['elements'], p['counts'], cfg['a'], cfg['file'])
        elif stype == 'bcc':
            generate_bcc(p['elements'], p['counts'], cfg['a'], cfg['file'])
        elif stype == 'hcp':
            generate_hcp(p['elements'], p['counts'], cfg['a'], cfg.get('c_over_a', 1.622), cfg['file'])
        elif stype == 'l10':
            generate_l10(p['elements'], p['counts'], cfg['a'], cfg.get('c_over_a', 0.98), cfg['file'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='POSCAR-Startstrukturen fuer EvAX erzeugen',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Beispiele:
  python scripts/generate_structures.py --preset CoPt
  python scripts/generate_structures.py --preset MEA_CoNiCu
  python scripts/generate_structures.py --elements Co Pt --counts 128 128 \\
      --structure fcc --lattice 3.80 --output structures/CoPt_fcc.p1
        """,
    )
    parser.add_argument('--preset', choices=list(PRESETS.keys()),
                        help='Vordefiniertes Preset (erzeugt alle Strukturen)')
    parser.add_argument('--elements', nargs='+', help='Elementnamen (z.B. Co Pt)')
    parser.add_argument('--counts', nargs='+', type=int, help='Anzahl pro Element')
    parser.add_argument('--structure', choices=['fcc', 'bcc', 'hcp', 'l10'],
                        help='Kristallstruktur')
    parser.add_argument('--lattice', type=float, help='Gitterkonstante a (Angstrom)')
    parser.add_argument('--c-over-a', type=float, default=1.622,
                        help='c/a-Verhaeltnis (fuer HCP/L10, Standard: 1.622)')
    parser.add_argument('--output', help='Ausgabedatei (.p1)')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')

    args = parser.parse_args()

    if args.preset:
        run_preset(args.preset)
    elif args.elements and args.counts and args.structure and args.lattice and args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        if args.structure == 'fcc':
            generate_fcc(args.elements, args.counts, args.lattice, args.output, args.seed)
        elif args.structure == 'bcc':
            generate_bcc(args.elements, args.counts, args.lattice, args.output, args.seed)
        elif args.structure == 'hcp':
            generate_hcp(args.elements, args.counts, args.lattice, args.c_over_a, args.output, args.seed)
        elif args.structure == 'l10':
            generate_l10(args.elements, args.counts, args.lattice, args.c_over_a, args.output, args.seed)
    else:
        parser.print_help()
