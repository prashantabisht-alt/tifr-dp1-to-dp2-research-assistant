"""
Forensic test #2: the BUG is a sign error in c_3 (chirality term).

Compares:
  (a) Dipanjan's c_3 as written:           + 2 i epsilon sin(a k1 - b k2)
  (b) Corrected c_3 (state 2 forward = (-a,b)): - 2 i epsilon sin(a k1 - b k2)

Plots P(x, y, t=50) over a 2D region for both versions. The corrected
version should be 6-fold symmetric; the buggy version should not.

Parameters match Fig 11 of the Confinement_enhanced_clustering draft:
  gamma = 0.01, epsilon = 0.15, t = 50, L = 30 grid for k-sum.
"""
import numpy as np
import matplotlib.pyplot as plt

# Parameters from the draft's Fig 11
gamma   = 0.01
epsilon = 0.15
a, b    = 1.0, 1.0
L       = 30
x0, y0  = 0, 0
t       = 50.0


def build_Mk(k1, k2, fix_c3_sign=False):
    """
    6x6 Bloch matrix.

    fix_c3_sign=False -> Dipanjan's original (buggy: c_3 has +sin)
    fix_c3_sign=True  -> Corrected (c_3 has -sin, as the master eqn demands)
    """
    bulk = (1.0 / 3) * (
        np.cos(2 * a * k1) + np.cos(a * k1 + b * k2) + np.cos(a * k1 - b * k2)
    )
    c1 = bulk - 1 - gamma + 2j * epsilon * np.sin(2 * a * k1)
    c2 = bulk - 1 - gamma + 2j * epsilon * np.sin(a * k1 + b * k2)
    # The bug: Dipanjan writes +sin, but state 2 forward=(-a,b) so it should be -sin
    chir3_sign = -1.0 if fix_c3_sign else +1.0
    c3 = bulk - 1 - gamma + chir3_sign * 2j * epsilon * np.sin(a * k1 - b * k2)
    diag = [c1, c2, c3, np.conj(c1), np.conj(c2), np.conj(c3)]
    M = np.zeros((6, 6), dtype=complex)
    for i in range(6):
        M[i, i] = diag[i]
        M[i, (i + 1) % 6] += gamma / 2
        M[i, (i - 1) % 6] += gamma / 2
    return M


def Pxy_grid(target_xs, target_ys, fix_c3_sign):
    """Compute P(x, y, t) on a 2D grid of (target_x, target_y) sites."""
    # Use sheared k-grid (matches Dipanjan's rect, just less wasteful)
    k_points = []
    for nx in range(2 * L):
        for ny in range(L):
            k_points.append((2 * np.pi * nx / (2 * a * L),
                              2 * np.pi * ny / (b * L)))
    N = 2 * L ** 2

    # Precompute eigendecomposition at each k
    eig_data = []
    for kx, ky in k_points:
        M = build_Mk(kx, ky, fix_c3_sign=fix_c3_sign)
        evals, evecs = np.linalg.eig(M)
        rhtvct = (1.0 / 6) * np.exp(1j * (kx * x0 + ky * y0)) * np.ones(6, dtype=complex)
        prefact = np.linalg.solve(evecs, rhtvct)
        ptilde_vec = evecs @ (prefact * np.exp(evals * t))
        pt = np.sum(ptilde_vec)
        eig_data.append((kx, ky, pt))

    Pxy = np.zeros((len(target_ys), len(target_xs)))
    for ix, x in enumerate(target_xs):
        for iy, y in enumerate(target_ys):
            s = 0.0
            for kx, ky, pt in eig_data:
                s += np.real(pt * np.exp(-1j * (kx * x + ky * y)))
            Pxy[iy, ix] = s / N
    return Pxy


# Use centered grid like Fig 11: x ∈ [0, 60], y ∈ [0, 30]
# Step by 2 in x (only physical sites on this sublattice convention)
# and by 1 in y. Re-center around origin (initial condition at 0).
target_xs = np.arange(0, 60, 2)   # 30 even x values
target_ys = np.arange(0, 30, 1)   # 30 y values

print("Computing buggy version...")
P_bug = Pxy_grid(target_xs, target_ys, fix_c3_sign=False)
print(f"  range: [{P_bug.min():.3e}, {P_bug.max():.3e}]")

print("Computing corrected version...")
P_fix = Pxy_grid(target_xs, target_ys, fix_c3_sign=True)
print(f"  range: [{P_fix.min():.3e}, {P_fix.max():.3e}]")

# Plot side-by-side
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

im0 = axes[0].imshow(P_bug, origin="lower",
                     extent=[target_xs.min(), target_xs.max(),
                             target_ys.min(), target_ys.max()],
                     aspect="auto", cmap="hot")
axes[0].set_title(rf"(a) Dipanjan's $c_3$ (BUGGY: $+\sin$)  $\gamma$={gamma}, $\epsilon$={epsilon}, t={int(t)}")
axes[0].set_xlabel("x"); axes[0].set_ylabel("y")
plt.colorbar(im0, ax=axes[0])

im1 = axes[1].imshow(P_fix, origin="lower",
                     extent=[target_xs.min(), target_xs.max(),
                             target_ys.min(), target_ys.max()],
                     aspect="auto", cmap="hot")
axes[1].set_title(rf"(b) Corrected $c_3$ (FIXED: $-\sin$)")
axes[1].set_xlabel("x"); axes[1].set_ylabel("y")
plt.colorbar(im1, ax=axes[1])

diff = P_fix - P_bug
vmax = np.max(np.abs(diff))
im2 = axes[2].imshow(diff, origin="lower",
                     extent=[target_xs.min(), target_xs.max(),
                             target_ys.min(), target_ys.max()],
                     aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
axes[2].set_title("(c) Difference  (fixed - buggy)")
axes[2].set_xlabel("x"); axes[2].set_ylabel("y")
plt.colorbar(im2, ax=axes[2])

plt.tight_layout()
plt.savefig("signbug_demo.png", dpi=180, bbox_inches="tight")
print("\nSaved signbug_demo.png")

# Also dump a numerical summary
print("\nNumerical summary at the row y = 4 (matching Dipanjan's original output):")
iy_4 = np.argmin(np.abs(target_ys - 4))
print(f"  {'x':>4}  {'buggy P':>15}  {'fixed P':>15}  {'|diff|/P':>10}")
for ix, x in enumerate(target_xs):
    pb = P_bug[iy_4, ix]
    pf = P_fix[iy_4, ix]
    rel = abs(pb - pf) / max(abs(pf), 1e-30)
    print(f"  {x:>4}  {pb:>15.4e}  {pf:>15.4e}  {rel:>10.2%}")
