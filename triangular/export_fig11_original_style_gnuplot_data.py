"""
Export original-Confinement-style Fig. 11 data for gnuplot.

This is the visual-reproduction figure, not the finite-L=30 verification
figure.  The old draft Fig. 11 shows an annular/hexagonal active front with
low probability at the center and edges.  To reproduce that visual convention
honestly, we use a larger L=60 torus at t=50 so the front has not visibly
wrapped around the finite box.

Outputs plain text:
    outputs/fig11_original_style_theory_points.txt
    outputs/fig11_original_style_kmc_points.txt     (if L60 KMC exists)
    outputs/fig11_original_style_cross_section.txt

The cross-section is an off-center fixed-y slice, matching the kind of slice
printed by Dipanjan's notebook for the old Fig. 11.  We include both old
failure modes: the draft-like old c3-sign curve and the full historical buggy
theory (rectangular notebook k-grid plus old c3 sign).  The corrected curve
should sit on top of KMC.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from triangular_jmvr_corrected import (
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
L = 60
t_final = 50.0


def theory_probability_rectangular_notebook(builder) -> np.ndarray:
    """Exact-looking P[n2,n1] using Dipanjan's rectangular 2L x L k-grid.

    This intentionally reproduces the old notebook convention for comparison.
    It is not the corrected triangular-torus Fourier grid.
    """
    eig_data: list[tuple[float, float, complex]] = []
    initial = (1.0 / 6.0) * np.ones(6, dtype=complex)

    for nx in range(2 * L):
        for ny in range(L):
            k1 = 2.0 * np.pi * nx / (2.0 * L)
            k2 = 2.0 * np.pi * ny / L
            M = builder(gamma, epsilon, k1, k2, 1.0, 1.0)
            evals, evecs = np.linalg.eig(M)
            prefactor = np.linalg.solve(evecs, initial)
            ptilde = np.sum(evecs @ (prefactor * np.exp(evals * t_final)))
            eig_data.append((k1, k2, ptilde))

    P = np.zeros((L, L), dtype=float)
    for n1 in range(L):
        for n2 in range(L):
            x = 2.0 * n1 + n2
            y = n2
            total = 0.0
            for k1, k2, ptilde in eig_data:
                total += np.real(ptilde * np.exp(-1j * (k1 * x + k2 * y)))
            P[n2, n1] = total / (2 * L * L)
    return P


def display_points(
    L: int,
    x_range: tuple[float, float] = (-1.0, 60.0),
    y_range: tuple[float, float] = (-1.0, 52.0),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return centered axial triangular display points for the draft-style panel.

    Important: our simulation coordinates are axial triangular coordinates
    (n1,n2), not odd-row square coordinates.  The correct Cartesian
    plotting map is

        x = n1 + n2/2,
        y = sqrt(3) n2 / 2,

    up to an overall display shift.  A previous draft-style plot used
    ``x = n1 + (n2 mod 2)/2``; that is an odd-row plotting convention and it
    visibly tilted the annulus.  This function uses the true axial map and
    includes nearby periodic images so the rectangular plotting window is
    filled with valid low-probability lattice sites.
    """
    sqrt3 = np.sqrt(3.0)
    center_x = L / 2.0
    center_y = 0.5 * sqrt3 * (L / 2.0)

    xs, ys, n1s, n2s = [], [], [], []
    for n1 in range(L):
        for n2 in range(L):
            for i1 in (-1, 0, 1):
                for i2 in (-1, 0, 1):
                    u = n1 + i1 * L
                    v = n2 + i2 * L
                    x = center_x + u + 0.5 * v
                    y = center_y + 0.5 * sqrt3 * v
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


def write_cross_section(
    path: Path,
    x: np.ndarray,
    y: np.ndarray,
    n1: np.ndarray,
    n2: np.ndarray,
    P_full_bug: np.ndarray,
    P_sign_bug: np.ndarray,
    P_fix: np.ndarray,
    P_kmc: np.ndarray | None,
) -> None:
    # Dipanjan's notebook printed a fixed off-center y-slice, not the center
    # row of the annulus.  This row reproduces the old unequal two-peak curve.
    target_y = 0.5 * np.sqrt(3.0) * 20.0
    mask = np.isclose(y, target_y)
    order = np.argsort(x[mask])
    xs = x[mask][order]
    n1s = n1[mask][order]
    n2s = n2[mask][order]
    with path.open("w") as f:
        if P_kmc is None:
            f.write("# x  P_full_original_buggy  P_old_c3_sign  P_corrected_theory\n")
            for xi, i1, i2 in zip(xs, n1s, n2s):
                f.write(
                    f"{xi:.12g} {P_full_bug[i2, i1]:.16e} "
                    f"{P_sign_bug[i2, i1]:.16e} {P_fix[i2, i1]:.16e}\n"
                )
        else:
            f.write(
                "# x  P_kmc  P_full_original_buggy  P_old_c3_sign  "
                "P_corrected_theory\n"
            )
            for xi, i1, i2 in zip(xs, n1s, n2s):
                f.write(
                    f"{xi:.12g} {P_kmc[i2, i1]:.16e} "
                    f"{P_full_bug[i2, i1]:.16e} {P_sign_bug[i2, i1]:.16e} "
                    f"{P_fix[i2, i1]:.16e}\n"
                )


def main() -> None:
    print(f"Computing L={L} corrected theory...")
    P_fix = theory_probability_on_lattice(
        build_Mk_corrected, gamma, epsilon, L=L, t_final=t_final
    )
    print(f"Computing L={L} old c3-sign theory on corrected triangular grid...")
    P_sign_bug = theory_probability_on_lattice(
        build_Mk_dipanjan, gamma, epsilon, L=L, t_final=t_final
    )
    print(f"Computing L={L} full historical buggy theory...")
    P_full_bug = theory_probability_rectangular_notebook(build_Mk_dipanjan)

    X, Y, n1_display, n2_display = display_points(L)
    write_points(
        OUT / "fig11_original_style_theory_points.txt",
        X,
        Y,
        P_fix[n2_display, n1_display],
    )

    kmc_path = HERE / "kmc_triangular_counts_L60.txt"
    P_kmc_disp = None
    if kmc_path.exists():
        counts, n_walkers = load_kmc_counts(kmc_path, L=L)
        P_kmc = counts / n_walkers
        P_kmc_disp = P_kmc
        write_points(
            OUT / "fig11_original_style_kmc_points.txt",
            X,
            Y,
            P_kmc[n2_display, n1_display],
        )
        rms_fix = float(np.sqrt(np.mean((P_kmc - P_fix) ** 2)))
        rms_sign_bug = float(np.sqrt(np.mean((P_kmc - P_sign_bug) ** 2)))
        rms_full_bug = float(np.sqrt(np.mean((P_kmc - P_full_bug) ** 2)))
        mc_noise = float(np.sqrt(P_kmc.max() / n_walkers))
        print(f"Loaded L={L} KMC counts: N = {n_walkers:,}")
        print(f"RMS full original buggy vs KMC = {rms_full_bug:.3e}")
        print(f"RMS old c3 sign vs KMC         = {rms_sign_bug:.3e}")
        print(f"RMS corrected vs KMC           = {rms_fix:.3e}")
        print(f"MC noise floor                 = {mc_noise:.3e}")
    else:
        print("No L60 KMC file found yet; gnuplot will use theory in both heatmap slots.")

    write_cross_section(
        OUT / "fig11_original_style_cross_section.txt",
        X,
        Y,
        n1_display,
        n2_display,
        P_full_bug,
        P_sign_bug,
        P_fix,
        P_kmc_disp,
    )
    print("Saved:")
    print(f"  {OUT / 'fig11_original_style_theory_points.txt'}")
    print(f"  {OUT / 'fig11_original_style_cross_section.txt'}")
    if P_kmc_disp is not None:
        print(f"  {OUT / 'fig11_original_style_kmc_points.txt'}")


if __name__ == "__main__":
    main()
