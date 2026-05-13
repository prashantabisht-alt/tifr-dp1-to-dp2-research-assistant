"""
Diagnostic that killed the "PBC bug" hypothesis.

Kabir said Dipanjan probably implemented PBC wrong. My first guess was
that the rectangular k-grid in rtp_tl_2.nb was wrong for the triangular
torus, and a sheared k-grid would give a different (correct) answer.

This script computes P(n1, n2, t) two ways and shows they agree to
machine precision — proving PBC is NOT the bug. The bug must therefore
be inside the Bloch matrix M(k) itself, which is where we then found
the c_3 sign error.

Run:
    python3 forensic_pbc_is_fine.py
Expected output:
    max |P_rect - P_shear|  ≈  1e-16
"""
from __future__ import annotations
import numpy as np

# Parameters from rtp_tl_2.nb
gamma   = 0.01
epsilon = 0.15
a, b    = 1.0, 1.0
L       = 30
t_final = 50.0


def build_Mk(k1, k2):
    """6x6 Bloch matrix (using Dipanjan's BUGGY c_3 with +sin — same in both
    grids; the point is to test PBC, not the bug)."""
    bulk = (1.0 / 3) * (
        np.cos(2 * a * k1) + np.cos(a * k1 + b * k2) + np.cos(a * k1 - b * k2)
    )
    c1 = bulk - 1 - gamma + 2j * epsilon * np.sin(2 * a * k1)
    c2 = bulk - 1 - gamma + 2j * epsilon * np.sin(a * k1 + b * k2)
    c3 = bulk - 1 - gamma + 2j * epsilon * np.sin(a * k1 - b * k2)
    diag = [c1, c2, c3, np.conj(c1), np.conj(c2), np.conj(c3)]
    M = np.zeros((6, 6), dtype=complex)
    for i in range(6):
        M[i, i] = diag[i]
        M[i, (i + 1) % 6] += gamma / 2
        M[i, (i - 1) % 6] += gamma / 2
    return M


def inverse_fourier(grid_kx, grid_ky):
    """Given a list of (kx, ky) pairs, evolve and inverse-FT to P(n1, n2)."""
    P = np.zeros((L, L))
    eig_data = []
    for kx, ky in zip(grid_kx, grid_ky):
        M = build_Mk(kx, ky)
        evals, evecs = np.linalg.eig(M)
        rhtvct = (1.0 / 6) * np.ones(6, dtype=complex)
        prefact = np.linalg.solve(evecs, rhtvct)
        ptilde_vec = evecs @ (prefact * np.exp(evals * t_final))
        pt = np.sum(ptilde_vec)
        eig_data.append((kx, ky, pt))
    for n1 in range(L):
        for n2 in range(L):
            x_cart = 2 * a * n1 + a * n2
            y_cart = b * n2
            s = 0.0
            for kx, ky, pt in eig_data:
                s += np.real(pt * np.exp(-1j * (kx * x_cart + ky * y_cart)))
            P[n2, n1] = s / (L * L)
    return P


# Way A — Dipanjan's RECTANGULAR k-grid (what's in rtp_tl_2.nb)
kx_rect, ky_rect = [], []
for m1 in range(L):
    for m2 in range(L):
        kx_rect.append(2 * np.pi * m1 / (2 * a * L))
        ky_rect.append(2 * np.pi * m2 / (b * L))

# Way B — SHEARED triangular k-grid (my "fix" hypothesis)
kx_shear, ky_shear = [], []
for m1 in range(L):
    for m2 in range(L):
        kx_shear.append(np.pi * m1 / (a * L))
        ky_shear.append(np.pi * (2 * m2 - m1) / (b * L))

print("Computing P with Dipanjan's rectangular k-grid...")
P_rect = inverse_fourier(kx_rect, ky_rect)

print("Computing P with sheared triangular k-grid...")
P_shear = inverse_fourier(kx_shear, ky_shear)

diff = np.max(np.abs(P_rect - P_shear))
print(f"\n  max |P_rect - P_shear|  =  {diff:.3e}")
print(f"  max |P|                 =  {P_rect.max():.3e}")
print(f"  relative                =  {diff / P_rect.max():.3e}")

if diff < 1e-12:
    print("\n  => RECT and SHEAR grids give bit-identical P.")
    print("     PBC is NOT the bug. The two grids enumerate the same")
    print("     k-vectors modulo the reciprocal lattice; relabelling")
    print("     cannot change the inverse Fourier sum.")
    print("     The bug must be inside the Bloch matrix M(k).")
else:
    print("\n  => Grids disagree — investigate.")
