"""
Export corrected Fig. 11 data in the same visual coordinate style as the draft.

The draft's Fig. 11 is not a smooth triangular-coordinate pm3d surface. It is
a colored point plot in Dipanjan's rectangular Cartesian embedding:

    x = (2 n1 + n2) mod 2L
    y = n2

Only the parity sublattice x+y even is physical for a walker started at (0,0).
This representation makes the triangular lattice visible as a dot pattern and
keeps the domain rectangular, closer to the original draft figure.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from export_fig11_gnuplot_data import (
    HERE,
    L,
    OUT,
    load_fortran_kmc,
    theory_P_on_lattice,
)


def lattice_to_draft_xy(P_lattice: np.ndarray) -> np.ndarray:
    """
    Convert P[n2, n1] to Pxy[x, y] in Dipanjan's rectangular embedding.

    Off-parity sites are stored as NaN so the gnuplot point plot draws only the
    physical triangular-lattice sites.
    """
    Pxy = np.full((2 * L, L), np.nan, dtype=float)
    for n2 in range(L):
        for n1 in range(L):
            x = (2 * n1 + n2) % (2 * L)
            y = n2
            Pxy[x, y] = P_lattice[n2, n1]
    return Pxy


def center_rectangular_domain(Pxy: np.ndarray) -> np.ndarray:
    """
    Put the initial site near the center of the rectangular display.

    The unshifted walker starts at (x,y)=(0,0). The draft-style centered view
    uses a rectangular translation by (L, L/2).
    """
    return np.roll(Pxy, shift=(L, L // 2), axis=(0, 1))


def write_points(path: Path, Pxy_centered: np.ndarray) -> None:
    """Write x y P for physical sites only."""
    with path.open("w") as f:
        f.write("# x  y  P\n")
        f.write("# Dipanjan rectangular embedding; physical parity sites only\n")
        for y in range(L):
            for x in range(2 * L):
                val = Pxy_centered[x, y]
                if np.isfinite(val):
                    f.write(f"{x:d} {y:d} {val:.16e}\n")
            f.write("\n")


def write_cross_section(
    path: Path,
    P_kmc_xy: np.ndarray,
    P_bug_xy: np.ndarray,
    P_fix_xy: np.ndarray,
) -> None:
    """Write x, KMC, buggy, corrected along y=L/2."""
    y = L // 2
    with path.open("w") as f:
        f.write("# x  P_kmc  P_buggy_theory  P_corrected_theory\n")
        f.write(f"# rectangular cross-section y={y}\n")
        for x in range(2 * L):
            vals = (P_kmc_xy[x, y], P_bug_xy[x, y], P_fix_xy[x, y])
            if all(np.isfinite(v) for v in vals):
                f.write(f"{x:d} {vals[0]:.16e} {vals[1]:.16e} {vals[2]:.16e}\n")


def main() -> None:
    counts, n_walkers = load_fortran_kmc(HERE / "kmc_triangular_counts.txt")
    P_kmc = counts / n_walkers

    print(f"Loaded KMC counts: N = {n_walkers:,}")
    print("Computing corrected theory...")
    P_fix = theory_P_on_lattice(fix_c3=True)
    print("Computing buggy theory...")
    P_bug = theory_P_on_lattice(fix_c3=False)

    P_kmc_xy = center_rectangular_domain(lattice_to_draft_xy(P_kmc))
    P_fix_xy = center_rectangular_domain(lattice_to_draft_xy(P_fix))
    P_bug_xy = center_rectangular_domain(lattice_to_draft_xy(P_bug))

    write_points(OUT / "fig11_draft_theory_points.txt", P_fix_xy)
    write_points(OUT / "fig11_draft_kmc_points.txt", P_kmc_xy)
    write_cross_section(OUT / "fig11_draft_cross_section.txt", P_kmc_xy, P_bug_xy, P_fix_xy)

    print("Saved:")
    print(f"  {OUT / 'fig11_draft_theory_points.txt'}")
    print(f"  {OUT / 'fig11_draft_kmc_points.txt'}")
    print(f"  {OUT / 'fig11_draft_cross_section.txt'}")


if __name__ == "__main__":
    main()
