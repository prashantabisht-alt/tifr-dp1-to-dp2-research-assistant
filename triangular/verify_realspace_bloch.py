"""Finite-torus real-space generator versus Bloch matrix check.

This is the most direct check that the sheared k-grid is the correct grid for
the triangular L x L torus used by the KMC.

For a small L, build the full real-space generator on states (n1,n2,d) with
periodic lattice-coordinate boundaries. Then compare its eigenvalues with the
union of 6x6 Bloch spectra over:

  - the corrected sheared triangular-torus grid,
  - Dipanjan's rectangular notebook grid,
  - the old-sign matrix on the sheared grid.

Expected result:
  - real-space spectrum == sheared corrected Bloch spectrum,
  - rectangular notebook grid contains extra off-torus modes,
  - old c3 sign does not match the real-space generator.
"""
from __future__ import annotations

import numpy as np

from triangular_jmvr_corrected import build_Mk_corrected, build_Mk_dipanjan


gamma = 0.01
epsilon = 0.15
a = 1.0
b = 1.0

NN = (
    (1, 0),
    (0, 1),
    (-1, 1),
    (-1, 0),
    (0, -1),
    (1, -1),
)


def state_index(n1: int, n2: int, d: int, L: int) -> int:
    return ((n2 * L + n1) * 6 + d)


def real_space_generator(L: int) -> np.ndarray:
    """Build the full 6L^2 x 6L^2 real-space generator."""
    G = np.zeros((6 * L * L, 6 * L * L), dtype=complex)
    for n1 in range(L):
        for n2 in range(L):
            for director in range(6):
                src = state_index(n1, n2, director, L)
                G[src, src] -= 1.0 + gamma

                # Director switching.
                G[state_index(n1, n2, (director + 1) % 6, L), src] += gamma / 2.0
                G[state_index(n1, n2, (director - 1) % 6, L), src] += gamma / 2.0

                # Translation to all six nearest neighbours.
                for hop_dir, (dn1, dn2) in enumerate(NN):
                    if hop_dir == director:
                        rate = 1.0 / 6.0 + epsilon
                    elif hop_dir == (director + 3) % 6:
                        rate = 1.0 / 6.0 - epsilon
                    else:
                        rate = 1.0 / 6.0
                    dst = state_index((n1 + dn1) % L, (n2 + dn2) % L, director, L)
                    G[dst, src] += rate
    return G


def bloch_spectra_sheared(builder, L: int) -> np.ndarray:
    """Union of Bloch spectra over the correct triangular-torus grid."""
    vals = []
    for m1 in range(L):
        for m2 in range(L):
            k1 = np.pi * m1 / (a * L)
            k2 = np.pi * (2.0 * m2 - m1) / (b * L)
            vals.extend(np.linalg.eigvals(builder(gamma, epsilon, k1, k2, a, b)))
    return np.array(vals)


def bloch_spectra_rectangular_notebook(builder, L: int) -> np.ndarray:
    """Union of Bloch spectra over Dipanjan's rectangular 2L x L notebook grid."""
    vals = []
    for nx in range(2 * L):
        for ny in range(L):
            k1 = 2.0 * np.pi * nx / (2.0 * a * L)
            k2 = 2.0 * np.pi * ny / (b * L)
            vals.extend(np.linalg.eigvals(builder(gamma, epsilon, k1, k2, a, b)))
    return np.array(vals)


def directed_max_nearest_distance(a_vals: np.ndarray, b_vals: np.ndarray) -> float:
    """For every a, find the nearest b; return the largest such distance."""
    return float(max(min(abs(a - b) for b in b_vals) for a in a_vals))


def report(L: int) -> None:
    real_vals = np.linalg.eigvals(real_space_generator(L))
    shear_fix = bloch_spectra_sheared(build_Mk_corrected, L)
    shear_old = bloch_spectra_sheared(build_Mk_dipanjan, L)
    rect_fix = bloch_spectra_rectangular_notebook(build_Mk_corrected, L)

    print(f"L = {L}")
    print(f"  real-space dimension        = {len(real_vals)}")
    print(f"  sheared Bloch modes         = {len(shear_fix)}")
    print(f"  rectangular notebook modes  = {len(rect_fix)}")

    print("  real <-> sheared corrected:")
    print(
        "    real to shear = "
        f"{directed_max_nearest_distance(real_vals, shear_fix):.3e}"
    )
    print(
        "    shear to real = "
        f"{directed_max_nearest_distance(shear_fix, real_vals):.3e}"
    )

    print("  real <-> rectangular corrected:")
    print(
        "    real to rect  = "
        f"{directed_max_nearest_distance(real_vals, rect_fix):.3e}"
    )
    print(
        "    rect to real  = "
        f"{directed_max_nearest_distance(rect_fix, real_vals):.3e}"
        "  (extra off-torus modes)"
    )

    print("  real <-> sheared old-sign:")
    print(
        "    real to old   = "
        f"{directed_max_nearest_distance(real_vals, shear_old):.3e}"
    )
    print(
        "    old to real   = "
        f"{directed_max_nearest_distance(shear_old, real_vals):.3e}"
    )
    print()


def main() -> None:
    for L in (3, 4, 5):
        report(L)


if __name__ == "__main__":
    main()
