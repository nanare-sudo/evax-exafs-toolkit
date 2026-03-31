#!/usr/bin/env python3
"""
Generate BCC and HCP POSCAR-format structure files for CoNiCu MEA
with 256 atoms (85 Co, 85 Ni, 86 Cu) randomly distributed.

Output format matches the EvAX .p1 input format used in this project.
"""

import numpy as np
import os

SEED = 42
OUTPUT_DIR = "/home/nanare/Dokumente/Uni/Master_Wuppertal/ERP_Hausarbeit/MEA_analysis/input_files"

# Composition: 85 Co, 85 Ni, 86 Cu = 256 total
N_CO, N_NI, N_CU = 85, 85, 86
N_TOTAL = N_CO + N_NI + N_CU  # 256


def random_element_assignment(rng):
    """Create a shuffled list of element labels: 85 Co, 85 Ni, 86 Cu."""
    labels = ["Co"] * N_CO + ["Ni"] * N_NI + ["Cu"] * N_CU
    rng.shuffle(labels)
    return labels


def write_poscar(filename, comment, lattice_vectors, frac_coords, labels):
    """Write a POSCAR-format .p1 file matching the EvAX convention."""
    # Count elements in order of appearance (Co, Ni, Cu)
    co_count = labels.count("Co")
    ni_count = labels.count("Ni")
    cu_count = labels.count("Cu")

    with open(filename, "w") as f:
        f.write(f"{comment}\n")
        f.write("1.0\n")
        for v in lattice_vectors:
            f.write(f"        {v[0]:.10f}         {v[1]:.10f}         {v[2]:.10f}\n")
        f.write(f"   Co    Ni    Cu\n")
        f.write(f"   {co_count}    {ni_count}    {cu_count}\n")
        f.write("Direct\n")
        for coord, label in zip(frac_coords, labels):
            f.write(f"   {coord[0]:.9f}   {coord[1]:.9f}   {coord[2]:.9f}   {label}\n")

    print(f"Written {filename} ({len(labels)} atoms: {co_count} Co, {ni_count} Ni, {cu_count} Cu)")


def generate_bcc():
    """
    BCC structure: 4x4x8 supercell = 128 unit cells x 2 atoms/cell = 256 atoms.
    BCC basis: (0,0,0) and (0.5,0.5,0.5) in fractional coords of the unit cell.
    Lattice parameter a = 2.87 A (typical BCC transition metal).
    Supercell dimensions: 4a x 4a x 8a = 11.48 x 11.48 x 22.96 A.
    """
    a_bcc = 2.87  # Angstrom
    nx, ny, nz = 4, 4, 8

    lattice = [
        [nx * a_bcc, 0.0, 0.0],
        [0.0, ny * a_bcc, 0.0],
        [0.0, 0.0, nz * a_bcc],
    ]

    # BCC basis in fractional coordinates of the unit cell
    basis = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]

    frac_coords = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for bx, by, bz in basis:
                    fx = (ix + bx) / nx
                    fy = (iy + by) / ny
                    fz = (iz + bz) / nz
                    frac_coords.append((fx, fy, fz))

    assert len(frac_coords) == N_TOTAL, f"BCC: expected {N_TOTAL}, got {len(frac_coords)}"

    rng = np.random.default_rng(SEED)
    labels = random_element_assignment(rng)

    outfile = os.path.join(OUTPUT_DIR, "MEA_CoNiCu_bcc.p1")
    write_poscar(
        outfile,
        f"CoNiCu MEA bcc a={a_bcc}",
        lattice,
        frac_coords,
        labels,
    )


def generate_hcp():
    """
    HCP structure: 8x4x4 supercell = 128 unit cells x 2 atoms/cell = 256 atoms.
    HCP basis (fractional): (1/3, 2/3, 1/4) and (2/3, 1/3, 3/4).
    Lattice parameters: a = 2.51 A, c = 4.07 A (c/a ~ 1.622).
    Hexagonal lattice vectors:
        a1 = (a, 0, 0)
        a2 = (-a/2, a*sqrt(3)/2, 0)
        a3 = (0, 0, c)
    Supercell: 8*a1, 4*a2, 4*a3.
    """
    a_hcp = 2.51  # Angstrom
    c_hcp = 4.07  # Angstrom
    nx, ny, nz = 8, 4, 4

    # Hexagonal unit cell vectors
    a1 = np.array([a_hcp, 0.0, 0.0])
    a2 = np.array([-a_hcp / 2.0, a_hcp * np.sqrt(3.0) / 2.0, 0.0])
    a3 = np.array([0.0, 0.0, c_hcp])

    # Supercell lattice vectors
    lattice = [
        (nx * a1).tolist(),
        (ny * a2).tolist(),
        (nz * a3).tolist(),
    ]

    # HCP basis in fractional coordinates of the unit cell
    basis = [(1.0 / 3.0, 2.0 / 3.0, 1.0 / 4.0), (2.0 / 3.0, 1.0 / 3.0, 3.0 / 4.0)]

    frac_coords = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for bx, by, bz in basis:
                    fx = (ix + bx) / nx
                    fy = (iy + by) / ny
                    fz = (iz + bz) / nz
                    frac_coords.append((fx, fy, fz))

    assert len(frac_coords) == N_TOTAL, f"HCP: expected {N_TOTAL}, got {len(frac_coords)}"

    rng = np.random.default_rng(SEED + 1)  # Different seed for different random distribution
    labels = random_element_assignment(rng)

    outfile = os.path.join(OUTPUT_DIR, "MEA_CoNiCu_hcp.p1")
    write_poscar(
        outfile,
        f"CoNiCu MEA hcp a={a_hcp} c={c_hcp}",
        lattice,
        frac_coords,
        labels,
    )


if __name__ == "__main__":
    generate_bcc()
    generate_hcp()
    print("Done.")
