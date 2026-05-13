"""
Export the final Fig. 11 data for gnuplot.

This is the PI-friendly workflow:

    Python computes exact theory / reads KMC / writes plain .txt files.
    gnuplot makes the final figure from those .txt files.

The heatmap coordinates match fig11_final_hex.py: periodic-image triangular
Cartesian coordinates tiled across the same rectangular viewing window. The
cross-section still uses minimum-image coordinates, so the one-dimensional
slice runs cleanly through the walker.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from triangular_jmvr_corrected import (
    L_DEFAULT,
    T_DEFAULT,
    build_Mk_corrected,
    build_Mk_dipanjan,
    load_kmc_counts,
    theory_probability_on_lattice,
)


HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)

gamma = 0.01
epsilon = 0.15
L = L_DEFAULT
t_final = T_DEFAULT


def min_image_coordinates(L: int) -> tuple[np.ndarray, np.ndarray]:
    """Return X[n2,n1], Y[n2,n1] for closest periodic image to origin."""
    X = np.zeros((L, L), dtype=float)
    Y = np.zeros((L, L), dtype=float)
    sqrt3 = np.sqrt(3.0)

    for n2 in range(L):
        for n1 in range(L):
            best = (0.0, 0.0, np.inf)
            for dn1 in (-1, 0, 1):
                for dn2 in (-1, 0, 1):
                    nn1 = n1 + dn1 * L
                    nn2 = n2 + dn2 * L
                    x = 2.0 * nn1 + nn2
                    y = sqrt3 * nn2
                    r2 = x * x + y * y
                    if r2 < best[2]:
                        best = (x, y, r2)
            X[n2, n1] = best[0]
            Y[n2, n1] = best[1]
    return X, Y


def periodic_tiled_lattice(
    L: int,
    x_range: tuple[float, float] = (-28.0, 28.0),
    y_range: tuple[float, float] = (-25.0, 25.0),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return all periodic-image sites inside the final Fig. 11 window."""
    sqrt3 = np.sqrt(3.0)
    T1x, T1y = 2.0 * L, 0.0
    T2x, T2y = float(L), sqrt3 * L
    n_images = 3

    xs, ys, n1s, n2s = [], [], [], []
    for n1 in range(L):
        for n2 in range(L):
            x_base = 2.0 * n1 + n2
            y_base = sqrt3 * n2
            for i1 in range(-n_images, n_images + 1):
                for i2 in range(-n_images, n_images + 1):
                    x = x_base + i1 * T1x + i2 * T2x
                    y = y_base + i1 * T1y + i2 * T2y
                    if x_range[0] <= x <= x_range[1] and y_range[0] <= y <= y_range[1]:
                        xs.append(x)
                        ys.append(y)
                        n1s.append(n1)
                        n2s.append(n2)
    return (
        np.array(xs, dtype=float),
        np.array(ys, dtype=float),
        np.array(n1s, dtype=int),
        np.array(n2s, dtype=int),
    )


def write_points(path: Path, x: np.ndarray, y: np.ndarray, p: np.ndarray) -> None:
    with path.open("w") as f:
        f.write("# x  y  P\n")
        for row in zip(x.ravel(), y.ravel(), p.ravel()):
            f.write(f"{row[0]:.12g} {row[1]:.12g} {row[2]:.16e}\n")


def horizontal_cross_section(
    P: np.ndarray, X: np.ndarray, Y: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    mask = np.isclose(Y, 0.0)
    x = X[mask]
    p = P[mask]
    order = np.argsort(x)
    return x[order], p[order]


def write_cross_section(
    path: Path,
    X: np.ndarray,
    Y: np.ndarray,
    P_kmc: np.ndarray,
    P_bug: np.ndarray,
    P_fix: np.ndarray,
) -> None:
    x, kmc = horizontal_cross_section(P_kmc, X, Y)
    _, bug = horizontal_cross_section(P_bug, X, Y)
    _, fix = horizontal_cross_section(P_fix, X, Y)
    with path.open("w") as f:
        f.write("# x  P_kmc  P_buggy_theory  P_corrected_theory\n")
        for row in zip(x, kmc, bug, fix):
            f.write(
                f"{row[0]:.12g} {row[1]:.16e} {row[2]:.16e} {row[3]:.16e}\n"
            )


def main() -> None:
    counts, n_walkers = load_kmc_counts(HERE / "kmc_triangular_counts.txt", L=L)
    P_kmc = counts / n_walkers

    print(f"Loaded KMC counts: N = {n_walkers:,}")
    print("Computing corrected theory...")
    P_fix = theory_probability_on_lattice(
        build_Mk_corrected, gamma, epsilon, L=L, t_final=t_final
    )
    print("Computing buggy theory...")
    P_bug = theory_probability_on_lattice(
        build_Mk_dipanjan, gamma, epsilon, L=L, t_final=t_final
    )

    X_tile, Y_tile, n1_tile, n2_tile = periodic_tiled_lattice(L)
    write_points(
        OUT / "fig11_final_hex_theory_points.txt",
        X_tile,
        Y_tile,
        P_fix[n2_tile, n1_tile],
    )
    write_points(
        OUT / "fig11_final_hex_kmc_points.txt",
        X_tile,
        Y_tile,
        P_kmc[n2_tile, n1_tile],
    )

    X, Y = min_image_coordinates(L)
    write_cross_section(
        OUT / "fig11_final_hex_cross_section.txt", X, Y, P_kmc, P_bug, P_fix
    )

    rms_bug = float(np.sqrt(np.mean((P_kmc - P_bug) ** 2)))
    rms_fix = float(np.sqrt(np.mean((P_kmc - P_fix) ** 2)))
    mc_noise = float(np.sqrt(P_kmc.max() / n_walkers))

    print("Saved:")
    print(f"  {OUT / 'fig11_final_hex_theory_points.txt'}")
    print(f"  {OUT / 'fig11_final_hex_kmc_points.txt'}")
    print(f"  {OUT / 'fig11_final_hex_cross_section.txt'}")
    print(f"RMS buggy  vs KMC = {rms_bug:.3e}")
    print(f"RMS fixed  vs KMC = {rms_fix:.3e}")
    print(f"MC noise floor    = {mc_noise:.3e}")


if __name__ == "__main__":
    main()
