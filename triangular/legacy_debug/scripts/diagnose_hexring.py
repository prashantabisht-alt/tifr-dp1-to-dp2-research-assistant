"""
Diagnose why our heatmap doesn't show the dramatic hexagonal ring
that Fig 11 of the Confinement draft shows, at SAME parameters
(gamma=0.01, eps=0.15, t=50).

Hypotheses:
  H1: Our data IS hexagonal-ring-shaped, but rendering masks it
      (imshow in sheared lattice coords squashes hex into rhombus).
  H2: Their figure uses a larger L so the ring is fully inside and
      shows hexagonal corners; our L=30 is too small.
  H3: We're plotting raw lattice P[n2,n1] which lives on a sheared
      grid; a proper hex-grid Cartesian rendering will show 6-fold.

This script:
  1. Loads existing 100M-walker Fortran KMC counts
  2. Plots in 3 ways:
     (a) raw P[n2,n1] with imshow         -> looks like a blob/rhombus
     (b) Cartesian (x,y) scatter with hex markers  -> should show hexagonal ring
     (c) angular profile P(r,theta) integrated over r in the ring band
         -> if 6-fold modulation exists, this is the clearest test
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib.pyplot as plt

L       = 30
a, b    = 1.0, np.sqrt(3)         # ISOTROPIC display
gamma   = 0.01
epsilon = 0.15
t_final = 50.0

# ── Load 100M-walker Fortran KMC ───────────────────────────────────────
counts = np.zeros((L, L), dtype=np.int64)
N = None
path = "kmc_triangular_counts.txt"
with open(path) as f:
    for line in f:
        if line.startswith("#"):
            if "N=" in line:
                N = int(line.split("N=")[1].strip())
            continue
        parts = line.split()
        if len(parts) != 3:
            continue
        n1, n2, c = int(parts[0]), int(parts[1]), int(parts[2])
        counts[n2, n1] = c
if N is None:
    N = int(counts.sum())
P = counts / N
print(f"Loaded {N:,} walkers. P range: [{P.min():.3e}, {P.max():.3e}]")
print(f"P at start (0,0): {P[0,0]:.3e}")
print(f"P at far (15,15): {P[15,15]:.3e}")

# ── Cartesian positions of every lattice site, WITHOUT roll ──────────
# Walker started at (n1, n2) = (0, 0) which is x=0, y=0
# Apply MINIMUM-IMAGE convention: each lattice site (n1, n2) has periodic
# images in 9 surrounding torus cells; pick the one closest to (0,0).
N1, N2 = np.meshgrid(np.arange(L), np.arange(L), indexing="ij")

def min_image_cart(n1, n2):
    """Return Cartesian (x, y) of lattice site (n1, n2) using the
    periodic image closest to the origin."""
    best_x, best_y, best_d2 = None, None, np.inf
    for dn1 in (-1, 0, 1):
        for dn2 in (-1, 0, 1):
            nn1 = n1 + dn1 * L
            nn2 = n2 + dn2 * L
            x = 2 * a * nn1 + a * nn2
            y = b * nn2
            d2 = x * x + y * y
            if d2 < best_d2:
                best_d2 = d2
                best_x, best_y = x, y
    return best_x, best_y

X = np.zeros_like(P)
Y = np.zeros_like(P)
for i in range(L):
    for j in range(L):
        X[i, j], Y[i, j] = min_image_cart(j, i)   # P[n2=i, n1=j]

R = np.sqrt(X**2 + Y**2)
THETA = np.arctan2(Y, X)

# ── 1) Angular profile in a ring band ─────────────────────────────────
# Ring should be at r ~ epsilon * t in Cartesian-distance units.
# eps*t = 7.5 in lattice; in Cartesian, drift in director m happens
# at speed eps along direction d_hat_m with |d_hat_m| = lattice spacing.
# Since hop lengths are 1 in NN units (Cartesian distance between
# adjacent sites with a=1, b=sqrt(3) is 1), the drift distance is
# eps*t = 7.5 in Cartesian.
r_ring = epsilon * t_final
print(f"\nExpected ring radius (Cartesian): {r_ring:.2f}")
print(f"Max in-torus Cartesian distance: {R.max():.2f}")

# Histogram P over angle, weighted, for points with r in [r_ring-3, r_ring+3]
mask_ring = (R >= r_ring - 3.0) & (R <= r_ring + 3.0)
print(f"Number of sites in ring band: {mask_ring.sum()}")

theta_bins = np.linspace(-np.pi, np.pi, 25)
theta_centers = 0.5 * (theta_bins[1:] + theta_bins[:-1])
angular_P = np.zeros(len(theta_centers))
for k in range(len(theta_centers)):
    bin_mask = mask_ring & (THETA >= theta_bins[k]) & (THETA < theta_bins[k + 1])
    if bin_mask.sum() > 0:
        angular_P[k] = P[bin_mask].mean()

# ── 2) Radial profile (1D cross-section after isotropy averaging) ────
r_bins = np.linspace(0, R.max(), 30)
r_centers = 0.5 * (r_bins[1:] + r_bins[:-1])
radial_P = np.zeros(len(r_centers))
for k in range(len(r_centers)):
    rb_mask = (R >= r_bins[k]) & (R < r_bins[k + 1])
    if rb_mask.sum() > 0:
        radial_P[k] = P[rb_mask].mean()

# ── 3) The three plots ────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel A: Cartesian scatter with hex markers (closest to draft Fig 11)
sc = axes[0].scatter(X.flatten(), Y.flatten(), c=P.flatten(),
                      s=150, marker="h", cmap="hot",
                      vmin=0, vmax=P.max(), edgecolor="none")
axes[0].set_aspect("equal")
axes[0].set_title(
    f"(a) Hex-marker Cartesian scatter — walker at origin\n"
    f"L={L}, γ={gamma}, ε={epsilon}, t={t_final}", fontsize=11
)
axes[0].set_xlabel("x"); axes[0].set_ylabel("y")
axes[0].axhline(0, color="cyan", ls="--", lw=0.5, alpha=0.5)
axes[0].axvline(0, color="cyan", ls="--", lw=0.5, alpha=0.5)
# Mark expected ring radius
theta_circle = np.linspace(0, 2*np.pi, 100)
axes[0].plot(r_ring * np.cos(theta_circle), r_ring * np.sin(theta_circle),
             "c--", lw=1.5, alpha=0.7, label=f"expected ring r={r_ring:.1f}")
axes[0].legend(loc="upper left", fontsize=9)
plt.colorbar(sc, ax=axes[0], fraction=0.046)

# Panel B: Radial profile P(r)
axes[1].plot(r_centers, radial_P, "bo-", lw=2, ms=6)
axes[1].axvline(r_ring, color="red", ls="--", label=f"ε·t = {r_ring:.1f}")
axes[1].set_xlabel("Cartesian radius r")
axes[1].set_ylabel(r"$\langle P \rangle_\theta(r)$")
axes[1].set_title("(b) Radial profile (angle-averaged)\nshould peak at r = ε·t if ring exists",
                   fontsize=11)
axes[1].legend(fontsize=10)
axes[1].grid(alpha=0.3)

# Panel C: Angular profile in ring band
axes[2].plot(np.degrees(theta_centers), angular_P, "ro-", lw=2, ms=6)
# Mark 6-fold director angles (in 60° increments)
for k in range(6):
    axes[2].axvline(60 * k - 180 + 30, color="gray", ls=":", alpha=0.5)
axes[2].set_xlabel("angle θ (degrees)")
axes[2].set_ylabel(r"$P(r \approx \epsilon t, \theta)$")
axes[2].set_title("(c) Angular profile in ring band\n6 peaks → hexagonal; flat → circular",
                   fontsize=11)
axes[2].grid(alpha=0.3)

plt.tight_layout()
plt.savefig("diagnose_hexring.png", dpi=140, bbox_inches="tight")
print("\nSaved diagnose_hexring.png")
print("\nINTERPRETATION:")
print("  - Panel (a): if you see a hexagonal ring with depleted center,")
print("    our data has the same structure as their Fig 11.")
print("  - Panel (b): radial profile should peak at r ≈ ε·t = 7.5.")
print("  - Panel (c): if there are 6 visible peaks (one per director),")
print("    the hexagonal modulation is real and present in our data.")
