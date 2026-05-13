"""
Fig-11 replacement, plotted in TRUE 60°-isotropic triangular Cartesian
coordinates (b = √3 a).

Important fact: the lattice-coordinate P[n_1, n_2] computed by both the
CTRW master-equation theory AND the kinetic Monte Carlo walker is
INVARIANT under the choice of a, b. The KMC histogram is in lattice
coordinates; the Fourier theory's k-grid and inverse-FT both use lattice
indices. Only the Cartesian visualization depends on a, b.

So this script reuses the same 100 M Fortran KMC output
(kmc_triangular_counts.txt) and the same theory calculation, but displays
the result with the isotropic Cartesian transform that makes the
walker's spread look properly hexagonal.

Uses scatter with hexagonal markers (matplotlib marker='h') so each
lattice site shows up as a clean hexagon in the visualization.

Output: fig11_corrected_isotropic.png — publication-quality
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib.pyplot as plt

# Match the model parameters used everywhere else
gamma   = 0.01
epsilon = 0.15
L       = 30
t_final = 50.0

# Lattice parameters for the THEORY computation (must match Fortran code)
a_theory, b_theory = 1.0, 1.0

# Lattice parameters for VISUALIZATION (isotropic — true triangular)
a_disp, b_disp = 1.0, np.sqrt(3)


# ──────────────────────────────────────────────────────────────────────────
# Load Fortran KMC output
# ──────────────────────────────────────────────────────────────────────────
def load_fortran_kmc(path="kmc_triangular_counts.txt"):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found — run kmc_triangular first")
    counts = np.zeros((L, L), dtype=np.int64)
    N = None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if "N=" in line:
                    try:
                        N = int(line.split("N=")[1].strip())
                    except Exception:
                        pass
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            n1, n2, c = int(parts[0]), int(parts[1]), int(parts[2])
            counts[n2, n1] = c
    if N is None:
        N = int(counts.sum())
    return counts, N


counts, N_walkers = load_fortran_kmc()
P_kmc = counts / N_walkers
print(f"Loaded {N_walkers:,}-walker KMC.   P range: [{P_kmc.min():.3e}, {P_kmc.max():.3e}]")


# ──────────────────────────────────────────────────────────────────────────
# Theory (same as before)
# ──────────────────────────────────────────────────────────────────────────
def build_Mk(k1, k2, fix_c3):
    bulk = (1.0 / 3) * (
        np.cos(2 * a_theory * k1)
        + np.cos(a_theory * k1 + b_theory * k2)
        + np.cos(a_theory * k1 - b_theory * k2)
    )
    c1 = bulk - 1 - gamma + 2j * epsilon * np.sin(2 * a_theory * k1)
    c2 = bulk - 1 - gamma + 2j * epsilon * np.sin(a_theory * k1 + b_theory * k2)
    sign3 = -1.0 if fix_c3 else +1.0
    c3 = bulk - 1 - gamma + sign3 * 2j * epsilon * np.sin(a_theory * k1 - b_theory * k2)
    diag = [c1, c2, c3, np.conj(c1), np.conj(c2), np.conj(c3)]
    M = np.zeros((6, 6), dtype=complex)
    for i in range(6):
        M[i, i] = diag[i]
        M[i, (i + 1) % 6] += gamma / 2
        M[i, (i - 1) % 6] += gamma / 2
    return M


def theory_P_on_lattice(fix_c3):
    """Returns P[n2, n1] for n1, n2 in [0, L) — lattice coordinate output."""
    P = np.zeros((L, L))
    eig_data = []
    for m1 in range(L):
        for m2 in range(L):
            kx = np.pi * m1 / (a_theory * L)
            ky = np.pi * (2 * m2 - m1) / (b_theory * L)
            M = build_Mk(kx, ky, fix_c3=fix_c3)
            evals, evecs = np.linalg.eig(M)
            rhtvct = (1.0 / 6) * np.ones(6, dtype=complex)
            prefact = np.linalg.solve(evecs, rhtvct)
            ptilde_vec = evecs @ (prefact * np.exp(evals * t_final))
            pt = np.sum(ptilde_vec)
            eig_data.append((kx, ky, pt))
    for n1 in range(L):
        for n2 in range(L):
            x_cart = 2 * a_theory * n1 + a_theory * n2
            y_cart = b_theory * n2
            s = 0.0
            for kx, ky, pt in eig_data:
                s += np.real(pt * np.exp(-1j * (kx * x_cart + ky * y_cart)))
            P[n2, n1] = s / (L * L)
    return P


print("Computing theory (corrected)...")
P_fix = theory_P_on_lattice(fix_c3=True)
print("Computing theory (buggy)...")
P_bug = theory_P_on_lattice(fix_c3=False)

# Center the walker in the display
shift = (L // 2, L // 2)
P_kmc_c = np.roll(P_kmc, shift=shift, axis=(0, 1))
P_fix_c = np.roll(P_fix, shift=shift, axis=(0, 1))
P_bug_c = np.roll(P_bug, shift=shift, axis=(0, 1))

# Compare metrics
rms_bug = float(np.sqrt(np.mean((P_kmc - P_bug) ** 2)))
rms_fix = float(np.sqrt(np.mean((P_kmc - P_fix) ** 2)))
mc_noise = float(np.sqrt(P_kmc.max() / N_walkers))
print(f"\n  BUGGY-KMC RMS = {rms_bug:.3e}")
print(f"  FIXED-KMC RMS = {rms_fix:.3e}")
print(f"  MC noise floor = {mc_noise:.3e}")


# ──────────────────────────────────────────────────────────────────────────
# Build the figure with isotropic Cartesian display
# ──────────────────────────────────────────────────────────────────────────
# Cartesian positions of each lattice site (n1, n2) under display lattice:
#   X = 2 a_disp n1 + a_disp n2
#   Y = b_disp n2     with b_disp = √3 a_disp = √3
N1, N2 = np.meshgrid(np.arange(L), np.arange(L), indexing="ij")
X_cart = 2 * a_disp * N1 + a_disp * N2          # shape (L, L)  [n1, n2]
Y_cart = b_disp * N2                            # shape (L, L)

# Flatten for scatter
X_flat = X_cart.flatten()
Y_flat = Y_cart.flatten()


def plot_panel(ax, P_centered, title, vmax, cbar_label=None):
    """Hexagonal-marker scatter heatmap of P_centered[n2, n1]."""
    # We had P_centered stored as P[n2, n1]; flatten matching the meshgrid (n1, n2) order
    P_for_scatter = P_centered.T.flatten()           # T to get [n1, n2] ordering
    sc = ax.scatter(X_flat, Y_flat, c=P_for_scatter, s=120, marker="h",
                    cmap="hot", vmin=0, vmax=vmax, edgecolor="none")
    ax.set_aspect("equal")
    ax.set_xlim(X_flat.min() - 1, X_flat.max() + 1)
    ax.set_ylim(Y_flat.min() - 1, Y_flat.max() + 1)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(title, fontsize=12)
    cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    if cbar_label:
        cb.set_label(cbar_label)


# ──────────────────────────────────────────────────────────────────────────
# The figure — 3 panels matching draft Fig 11 layout
# ──────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 11))
gs = fig.add_gridspec(2, 2, height_ratios=[1.4, 1])

ax_theory = fig.add_subplot(gs[0, 0])
ax_kmc    = fig.add_subplot(gs[0, 1])
ax_cross  = fig.add_subplot(gs[1, :])

vmax = max(P_fix_c.max(), P_kmc_c.max())

plot_panel(ax_theory, P_fix_c, r"(a) Theory  (corrected $c_3$)", vmax, "Probability")
plot_panel(ax_kmc,    P_kmc_c, f"(b) Kinetic Monte Carlo  (N = {N_walkers:,})", vmax, "Probability")

# Cross-section through the center of the walker
# Walker is at lattice (n1, n2) = (L/2, L/2) after rolling
# In Cartesian: x_center = 2*L/2 + L/2 = 1.5*L, y_center = b_disp*L/2
# For the cross-section, pick a horizontal slice at the row containing the walker
# i.e. n2 = L/2 (in centered indexing)
n2_slice = L // 2
n1_axis  = np.arange(L)

ax_cross.plot(n1_axis, P_kmc_c[n2_slice, :], "ko", ms=6,
              label=f"KMC (N = {N_walkers:,})", zorder=3)
ax_cross.plot(n1_axis, P_bug_c[n2_slice, :], "r-", lw=2,
              label=f"Buggy theory  (RMS = {rms_bug:.2e})", zorder=2)
ax_cross.plot(n1_axis, P_fix_c[n2_slice, :], "b-", lw=2,
              label=f"Fixed theory  (RMS = {rms_fix:.2e})", zorder=2)
ax_cross.set_xlabel(r"$n_1$  (walker centered at $n_1 = L/2 = $" + f"{L//2})", fontsize=11)
ax_cross.set_ylabel(r"$P(n_1, n_2 = L/2, t)$", fontsize=11)
ax_cross.set_title(rf"(c) Cross-section through the walker,  "
                   rf"$\gamma = {gamma}$, $\epsilon = {epsilon}$, $t = {t_final}$",
                   fontsize=12)
ax_cross.legend(loc="upper right", fontsize=10, framealpha=0.95)
ax_cross.grid(alpha=0.3, which="both")
ax_cross.set_xlim(0, L)

fig.suptitle(
    "Fig 11 — corrected, isotropic triangular display ($b = \\sqrt{3}\\,a$).  "
    "Theory matches KMC to MC noise; sign error in $c_3$ identified.",
    fontsize=13, y=0.995
)
plt.tight_layout()
plt.savefig("fig11_corrected_isotropic.png", dpi=180, bbox_inches="tight")
print("\nSaved fig11_corrected_isotropic.png  — draft-matching hexagonal layout.")
