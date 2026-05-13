"""Two-bug forensic comparison against the 100M-walker Fortran KMC.

Dipanjan's notebook uses a rectangular ``2L x L`` Fourier grid. The correct
triangular torus uses the sheared ``L x L`` reciprocal grid. Independently, the
notebook/draft has the old wrong sign in the c3 chirality term.

This script compares all four combinations:

  1. rectangular notebook grid + old c3 sign
  2. rectangular notebook grid + corrected c3 sign
  3. sheared triangular-torus grid + old c3 sign
  4. sheared triangular-torus grid + corrected c3 sign

Only the fourth case reaches the Monte Carlo noise floor.
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

gamma   = 0.01
epsilon = 0.15
a, b    = 1.0, 1.0
L       = 30
t_final = 50.0


def theory_P_rectangular_notebook(builder):
    """Exact P[n2,n1] using Dipanjan's rectangular 2L x L k-grid."""
    eig_data = []
    initial = (1.0 / 6.0) * np.ones(6, dtype=complex)
    for nx in range(2 * L):
        for ny in range(L):
            k1 = 2.0 * np.pi * nx / (2.0 * a * L)
            k2 = 2.0 * np.pi * ny / (b * L)
            M = builder(gamma, epsilon, k1, k2, a, b)
            evals, evecs = np.linalg.eig(M)
            prefact = np.linalg.solve(evecs, initial)
            ptilde = np.sum(evecs @ (prefact * np.exp(evals * t_final)))
            eig_data.append((k1, k2, ptilde))

    P = np.zeros((L, L))
    for n1 in range(L):
        for n2 in range(L):
            x_cart = 2.0 * a * n1 + a * n2
            y_cart = b * n2
            s = 0.0
            for k1, k2, pt in eig_data:
                s += np.real(pt * np.exp(-1j * (k1 * x_cart + k2 * y_cart)))
            P[n2, n1] = s / (2 * L * L)
    return P


def theory_P_sheared(builder):
    """Exact P[n2,n1] using the correct triangular-torus reciprocal grid."""
    return theory_probability_on_lattice(
        builder,
        gamma,
        epsilon,
        L=L,
        t_final=t_final,
        a=a,
        b=b,
    )


def main() -> None:
    counts, n_walkers = load_kmc_counts(Path("kmc_triangular_counts.txt"), L=L)
    P_kmc = counts / n_walkers
    mc_noise = float(np.sqrt(P_kmc.max() / n_walkers))

    cases = (
        ("rectangular notebook grid + BUGGY c3", theory_P_rectangular_notebook, build_Mk_dipanjan),
        ("rectangular notebook grid + FIXED c3", theory_P_rectangular_notebook, build_Mk_corrected),
        ("sheared triangular grid   + BUGGY c3", theory_P_sheared, build_Mk_dipanjan),
        ("sheared triangular grid   + FIXED c3", theory_P_sheared, build_Mk_corrected),
    )

    print(f"KMC: N = {n_walkers:,}, MC noise floor = {mc_noise:.3e}")
    print(f"{'case':<44}  {'RMS vs KMC':>12}   ratio/noise")
    print("-" * 74)
    for label, method, builder in cases:
        P = method(builder)
        rms = float(np.sqrt(np.mean((P_kmc - P) ** 2)))
        print(f"  {label:<42}  {rms:12.3e}   {rms/mc_noise:.2f}")


if __name__ == "__main__":
    main()
