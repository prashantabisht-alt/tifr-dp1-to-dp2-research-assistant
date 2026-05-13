"""
LEGACY forensic test: Dipanjan's rectangular k-grid vs a sheared
triangular-torus k-grid.

This predates the c3 sign-error diagnosis. It intentionally reproduces
Dipanjan's old Bloch matrix and should not be used for the current corrected
JMVR theory. The current source of truth is triangular_jmvr_corrected.py.

Reproduces Dipanjan's Mathematica calculation exactly, then runs the
corrected version. Compares P(x, 4, t=10) for x = 0, 2, ..., 2L-2.

If the rectangular-vs-sheared bug analysis is right:
  - Both versions should give similar but not identical results.
  - Discrepancies are finite-size effects of the wrong PBC topology.
  - They should converge as L -> infinity.

Parameters match Dipanjan exactly.
"""
import numpy as np

# Parameters (matching Dipanjan)
gamma   = 0.01
epsilon = 0.15
a, b    = 1.0, 1.0
L       = 30
x0, y0  = 0, 0
t       = 10.0


def build_Mk(k1, k2):
    """6x6 Bloch matrix — matches Dipanjan's rtp_tl_2.nb structure exactly."""
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


def compute_P(target_x, target_y, k_points, N_normalize, t_eval=t):
    """
    Inverse Fourier reconstruction of P(target_x, target_y, t) over
    the given k-grid, divided by N_normalize.
    """
    pxy = 0.0
    for kx, ky in k_points:
        M = build_Mk(kx, ky)
        evals, evecs = np.linalg.eig(M)
        # Initial condition in Fourier: (1/6) exp(i k . r0) for all 6 directors.
        rhtvct = (1.0 / 6) * np.exp(1j * (kx * x0 + ky * y0)) * np.ones(6, dtype=complex)
        # Project onto eigenbasis: prefact = V^{-1} rhtvct
        prefact = np.linalg.solve(evecs, rhtvct)
        # Time-evolved P~_m(k, t)
        ptilde_vec = evecs @ (prefact * np.exp(evals * t_eval))
        # Sum over directors to get total P~(k, t)
        pt = np.sum(ptilde_vec)
        # Add to inverse Fourier
        pxy += np.real(pt * np.exp(-1j * (kx * target_x + ky * target_y)))
    return pxy / N_normalize


# ──────────────────────────────────────────────────────────────────────────
# Two k-grids
# ──────────────────────────────────────────────────────────────────────────

# Dipanjan's rectangular grid (his Mathematica notebook):
#   k_x = 2π nx / (2 a L), nx ∈ [0, 2L)
#   k_y = 2π ny / (b L),   ny ∈ [0, L)
# Total: 2 L² k-points; normalize by 2 L².
rect_grid = [
    (2 * np.pi * nx / (2 * a * L), 2 * np.pi * ny / (b * L))
    for nx in range(2 * L)
    for ny in range(L)
]
rect_N = 2 * L ** 2

# CORRECTED grid for true skew-torus PBC.
# Reciprocal lattice of primitive vectors a1=(2a,0), a2=(a,b):
#   b1 = (π/a, -π/b)
#   b2 = (0,   2π/b)
# Allowed k on L×L torus: k = (m1/L) b1 + (m2/L) b2, m1,m2 ∈ [0, L)
# In Cartesian: k_x = π m1 / (a L),  k_y = π (2 m2 - m1) / (b L)
# Total: L² k-points; normalize by L².
shear_grid = [
    (np.pi * m1 / (a * L), np.pi * (2 * m2 - m1) / (b * L))
    for m1 in range(L)
    for m2 in range(L)
]
shear_N = L ** 2


# ──────────────────────────────────────────────────────────────────────────
# Run the comparison
# ──────────────────────────────────────────────────────────────────────────

print("=" * 72)
print("Forensic comparison: Dipanjan's rectangular k-grid vs corrected sheared")
print("=" * 72)
print(f"L = {L}, γ = {gamma}, ε = {epsilon}, a = b = {a}, t = {t}")
print(f"Initial condition at ({x0}, {y0}); evaluated along row y = 4")
print()
print(f"Rectangular grid (Dipanjan):  {len(rect_grid):>6d} k-points,  divide by {rect_N}")
print(f"Sheared grid (corrected):     {len(shear_grid):>6d} k-points,  divide by {shear_N}")
print()
print(f"{'x':>4}  {'P_rect (Dipanjan)':>22}  {'P_shear (corrected)':>22}  {'|Δ|':>11}  {'|Δ|/P_shear':>12}")
print("-" * 80)

for x in range(0, 2 * L, 2):
    Pr = compute_P(x, 4, rect_grid, rect_N)
    Ps = compute_P(x, 4, shear_grid, shear_N)
    d = abs(Pr - Ps)
    rel = d / abs(Ps) if abs(Ps) > 1e-30 else float("inf")
    print(f"{x:>4}  {Pr:>22.6e}  {Ps:>22.6e}  {d:>11.2e}  {rel:>12.2e}")
