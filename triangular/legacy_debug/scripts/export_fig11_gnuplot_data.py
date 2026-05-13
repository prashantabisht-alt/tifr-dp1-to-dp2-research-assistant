"""
Export Fig. 11 probability-distribution data for gnuplot.

This script keeps the scientific calculation in Python, but writes plain text
tables that can be plotted by fig11_corrected.gnu:

    outputs/fig11_theory_corrected_xyz.txt
    outputs/fig11_kmc_xyz.txt
    outputs/fig11_cross_section.txt

The heatmap tables use isotropic triangular display coordinates:

    X = 2 a n1 + a n2
    Y = sqrt(3) a n2

while the probability values are still the lattice probabilities P[n2, n1].
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# Parameters matching kmc_triangular_jmvr.f90 and Fig. 11 reproduction.
gamma = 0.01
epsilon = 0.15
L = 30
t_final = 50.0

# Theory convention from Dipanjan / draft equations.
a_theory = 1.0
b_theory = 1.0

# Display convention for true 60-degree triangular geometry.
a_disp = 1.0
b_disp = np.sqrt(3.0)

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)


def load_fortran_kmc(path: Path) -> tuple[np.ndarray, int]:
    """Load columns n1 n2 count into P[n2, n1]."""
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run ./run_kmc.sh first.")

    counts = np.zeros((L, L), dtype=np.int64)
    n_walkers = None
    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                if "N=" in line:
                    try:
                        n_walkers = int(line.split("N=")[1].strip())
                    except ValueError:
                        n_walkers = None
                continue

            parts = line.split()
            if len(parts) != 3:
                continue
            n1, n2, count = int(parts[0]), int(parts[1]), int(parts[2])
            counts[n2, n1] = count

    if n_walkers is None:
        n_walkers = int(counts.sum())
    return counts, n_walkers


def build_Mk(k1: float, k2: float, fix_c3: bool) -> np.ndarray:
    """6x6 continuous-time generator M(k)."""
    bulk = (1.0 / 3.0) * (
        np.cos(2.0 * a_theory * k1)
        + np.cos(a_theory * k1 + b_theory * k2)
        + np.cos(a_theory * k1 - b_theory * k2)
    )

    c1 = bulk - 1.0 - gamma + 2.0j * epsilon * np.sin(2.0 * a_theory * k1)
    c2 = bulk - 1.0 - gamma + 2.0j * epsilon * np.sin(a_theory * k1 + b_theory * k2)
    sign3 = -1.0 if fix_c3 else 1.0
    c3 = bulk - 1.0 - gamma + sign3 * 2.0j * epsilon * np.sin(
        a_theory * k1 - b_theory * k2
    )

    diag = [c1, c2, c3, np.conj(c1), np.conj(c2), np.conj(c3)]
    M = np.zeros((6, 6), dtype=complex)
    for src in range(6):
        M[src, src] = diag[src]
        M[(src + 1) % 6, src] += gamma / 2.0
        M[(src - 1) % 6, src] += gamma / 2.0
    return M


def theory_P_on_lattice(fix_c3: bool) -> np.ndarray:
    """Return exact P[n2, n1] at t_final."""
    eig_data: list[tuple[float, float, complex]] = []

    for m1 in range(L):
        for m2 in range(L):
            kx = np.pi * m1 / (a_theory * L)
            ky = np.pi * (2.0 * m2 - m1) / (b_theory * L)
            M = build_Mk(kx, ky, fix_c3=fix_c3)
            evals, evecs = np.linalg.eig(M)
            initial_director_distribution = (1.0 / 6.0) * np.ones(6, dtype=complex)
            prefactor = np.linalg.solve(evecs, initial_director_distribution)
            ptilde_vec = evecs @ (prefactor * np.exp(evals * t_final))
            eig_data.append((kx, ky, np.sum(ptilde_vec)))

    P = np.zeros((L, L), dtype=float)
    for n1 in range(L):
        for n2 in range(L):
            x_cart = 2.0 * a_theory * n1 + a_theory * n2
            y_cart = b_theory * n2
            total = 0.0
            for kx, ky, ptilde in eig_data:
                total += np.real(ptilde * np.exp(-1j * (kx * x_cart + ky * y_cart)))
            P[n2, n1] = total / (L * L)
    return P


def write_xyz_grid(path: Path, P_centered: np.ndarray) -> None:
    """
    Write X Y P table with blank lines between rows.

    The blank lines let gnuplot pm3d understand this as a 2D grid.
    P_centered is stored as P[n2, n1].
    """
    with path.open("w") as f:
        f.write("# X  Y  P\n")
        f.write("# X=2*n1+n2, Y=sqrt(3)*n2, P=P[n2,n1]\n")
        for n2 in range(L):
            for n1 in range(L):
                x = 2.0 * a_disp * n1 + a_disp * n2
                y = b_disp * n2
                f.write(f"{x:.12g} {y:.12g} {P_centered[n2, n1]:.16e}\n")
            f.write("\n")


def write_cross_section(
    path: Path,
    P_kmc_centered: np.ndarray,
    P_bug_centered: np.ndarray,
    P_fix_centered: np.ndarray,
) -> None:
    """Write n1, KMC, buggy theory, corrected theory through n2=L/2."""
    n2_slice = L // 2
    with path.open("w") as f:
        f.write("# n1  P_kmc  P_buggy_theory  P_corrected_theory\n")
        for n1 in range(L):
            f.write(
                f"{n1:d} "
                f"{P_kmc_centered[n2_slice, n1]:.16e} "
                f"{P_bug_centered[n2_slice, n1]:.16e} "
                f"{P_fix_centered[n2_slice, n1]:.16e}\n"
            )


def main() -> None:
    counts, n_walkers = load_fortran_kmc(HERE / "kmc_triangular_counts.txt")
    P_kmc = counts / n_walkers

    print(f"Loaded KMC counts: N = {n_walkers:,}")
    print("Computing corrected theory...")
    P_fix = theory_P_on_lattice(fix_c3=True)
    print("Computing buggy theory...")
    P_bug = theory_P_on_lattice(fix_c3=False)

    shift = (L // 2, L // 2)
    P_kmc_c = np.roll(P_kmc, shift=shift, axis=(0, 1))
    P_fix_c = np.roll(P_fix, shift=shift, axis=(0, 1))
    P_bug_c = np.roll(P_bug, shift=shift, axis=(0, 1))

    write_xyz_grid(OUT / "fig11_theory_corrected_xyz.txt", P_fix_c)
    write_xyz_grid(OUT / "fig11_kmc_xyz.txt", P_kmc_c)
    write_cross_section(OUT / "fig11_cross_section.txt", P_kmc_c, P_bug_c, P_fix_c)

    rms_bug = float(np.sqrt(np.mean((P_kmc - P_bug) ** 2)))
    rms_fix = float(np.sqrt(np.mean((P_kmc - P_fix) ** 2)))
    mc_noise = float(np.sqrt(P_kmc.max() / n_walkers))

    print("Saved:")
    print(f"  {OUT / 'fig11_theory_corrected_xyz.txt'}")
    print(f"  {OUT / 'fig11_kmc_xyz.txt'}")
    print(f"  {OUT / 'fig11_cross_section.txt'}")
    print(f"RMS buggy  vs KMC = {rms_bug:.3e}")
    print(f"RMS fixed  vs KMC = {rms_fix:.3e}")
    print(f"MC noise floor    = {mc_noise:.3e}")


if __name__ == "__main__":
    main()
