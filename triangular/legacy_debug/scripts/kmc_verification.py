"""
KMC verification of the sign-error bug.

Builds a continuous-time kinetic Monte Carlo walker that obeys master
equations B2-B7 of the Confinement draft (no Bloch matrix shortcuts).
Compares the histogram of walker positions at t=50 against:
  (a) the BUGGY theory  (c_3 has  +sin(a k1 - b k2))
  (b) the FIXED theory  (c_3 has  -sin(a k1 - b k2))

If our diagnosis is correct, KMC ~ FIXED, KMC ≠ BUGGY.

Parameters match the Fig 11 of the draft:
    gamma=0.01, epsilon=0.15, t=50, L=30 torus.

The KMC uses lattice coordinates (n1, n2) with PBC mod L. This makes
the triangular-torus topology trivial. Cartesian is only used for
plotting/comparison.

Authors:
  KMC and analysis: Prashant Bisht, with Claude
"""
import numpy as np
import matplotlib.pyplot as plt
import time

# ──────────────────────────────────────────────────────────────────────────
# Parameters (match Fig 11 of the Confinement draft)
# ──────────────────────────────────────────────────────────────────────────
gamma   = 0.01
epsilon = 0.15
a, b    = 1.0, 1.0
L       = 30
t_final = 50.0


# ──────────────────────────────────────────────────────────────────────────
# KMC walker
# ──────────────────────────────────────────────────────────────────────────

# 6 NN displacements in lattice coordinates (n1, n2):
#   director 0 = +a1 direction = (+1, 0)
#   director 1 = +a2          = (0, +1)
#   director 2 = -a1 + a2     = (-1, +1)
#   director 3 = -a1          = (-1, 0)
#   director 4 = -a2          = (0, -1)
#   director 5 = +a1 - a2     = (+1, -1)
NN_LATTICE = np.array([
    (+1, 0),
    (0, +1),
    (-1, +1),
    (-1, 0),
    (0, -1),
    (+1, -1),
], dtype=np.int32)

# Hop rate matrices per director: hop_rates[m, d] = rate to hop in dir d when in state m
# Forward (d = m): 1/6 + ε
# Backward (d = (m+3) mod 6): 1/6 - ε
# Sideways: 1/6
HOP_RATES = np.full((6, 6), 1.0 / 6, dtype=np.float64)
for m in range(6):
    HOP_RATES[m, m] = 1.0 / 6 + epsilon
    HOP_RATES[m, (m + 3) % 6] = 1.0 / 6 - epsilon

# Cumulative hop probabilities for fast direction sampling
HOP_CUM = np.cumsum(HOP_RATES, axis=1)  # each row sums to 1


def run_walker(rng):
    """One walker, continuous-time KMC, until t_final. Returns final (n1, n2)."""
    n1, n2 = 0, 0
    m = int(rng.integers(6))
    t = 0.0
    rate_total = 1.0 + gamma  # translation total = 1, rotation total = γ
    while True:
        dt = rng.exponential(1.0 / rate_total)
        if t + dt >= t_final:
            return n1, n2
        t += dt
        # Pick event: rotation or translation?
        if rng.random() < gamma / rate_total:
            # Rotation: ±1 mod 6 with equal prob (achiral rotation)
            m = (m + 1) % 6 if rng.random() < 0.5 else (m - 1) % 6
        else:
            # Translation: pick direction by cumulative
            u = rng.random()
            # Find first d where HOP_CUM[m, d] >= u
            d = int(np.searchsorted(HOP_CUM[m], u))
            n1 = (n1 + NN_LATTICE[d, 0]) % L
            n2 = (n2 + NN_LATTICE[d, 1]) % L


# ──────────────────────────────────────────────────────────────────────────
# Theory (Fourier / Bloch matrix)
# ──────────────────────────────────────────────────────────────────────────

def build_Mk(k1, k2, fix_c3=False):
    bulk = (1.0 / 3) * (
        np.cos(2 * a * k1) + np.cos(a * k1 + b * k2) + np.cos(a * k1 - b * k2)
    )
    c1 = bulk - 1 - gamma + 2j * epsilon * np.sin(2 * a * k1)
    c2 = bulk - 1 - gamma + 2j * epsilon * np.sin(a * k1 + b * k2)
    sign3 = -1.0 if fix_c3 else +1.0
    c3 = bulk - 1 - gamma + sign3 * 2j * epsilon * np.sin(a * k1 - b * k2)
    diag = [c1, c2, c3, np.conj(c1), np.conj(c2), np.conj(c3)]
    M = np.zeros((6, 6), dtype=complex)
    for i in range(6):
        M[i, i] = diag[i]
        M[i, (i + 1) % 6] += gamma / 2
        M[i, (i - 1) % 6] += gamma / 2
    return M


def theory_P_on_lattice(fix_c3):
    """Compute P(n1, n2, t_final) on the L×L lattice using the sheared BZ."""
    # Sheared reciprocal lattice: k = (m1/L) b1 + (m2/L) b2
    #   b1 = (π/a, -π/b)
    #   b2 = (0, 2π/b)
    # In Cartesian:
    #   k_x = π m1 / (a L)
    #   k_y = π (2 m2 - m1) / (b L)
    P = np.zeros((L, L))
    eig_data = []
    for m1 in range(L):
        for m2 in range(L):
            kx = np.pi * m1 / (a * L)
            ky = np.pi * (2 * m2 - m1) / (b * L)
            M = build_Mk(kx, ky, fix_c3=fix_c3)
            evals, evecs = np.linalg.eig(M)
            # Initial condition: P_m(r=0, 0) = 1/6 each → tilde P_m(k, 0) = 1/6 each
            rhtvct = (1.0 / 6) * np.ones(6, dtype=complex)
            prefact = np.linalg.solve(evecs, rhtvct)
            ptilde_vec = evecs @ (prefact * np.exp(evals * t_final))
            pt = np.sum(ptilde_vec)
            eig_data.append((kx, ky, pt))
    # Inverse Fourier at Cartesian site (x = 2 n1 + n2, y = n2)
    for n1 in range(L):
        for n2 in range(L):
            x_cart = 2 * a * n1 + a * n2
            y_cart = b * n2
            s = 0.0
            for kx, ky, pt in eig_data:
                s += np.real(pt * np.exp(-1j * (kx * x_cart + ky * y_cart)))
            P[n2, n1] = s / (L * L)
    return P


# ──────────────────────────────────────────────────────────────────────────
# Run everything
# ──────────────────────────────────────────────────────────────────────────

N_walkers = 200_000

print("=" * 64)
print(f"KMC verification of the c_3 sign bug")
print(f"  Parameters: γ = {gamma}, ε = {epsilon}, t = {t_final}, L = {L}")
print(f"  KMC walkers: {N_walkers:,}")
print("=" * 64)

t0 = time.time()
rng = np.random.default_rng(20260511)
counts = np.zeros((L, L), dtype=np.int64)
for i in range(N_walkers):
    n1, n2 = run_walker(rng)
    counts[n2, n1] += 1
    if (i + 1) % 50000 == 0:
        print(f"  KMC: {i+1:,}/{N_walkers:,}  ({time.time()-t0:.1f}s elapsed)")

P_kmc = counts / N_walkers
print(f"\nKMC done in {time.time()-t0:.1f}s")
print(f"  KMC P range: [{P_kmc.min():.3e}, {P_kmc.max():.3e}]")
print(f"  Sum P_kmc = {P_kmc.sum():.6f}  (should be 1)")

print("\nComputing buggy theory (Dipanjan's c_3 with +sin)...")
t0 = time.time()
P_bug = theory_P_on_lattice(fix_c3=False)
print(f"  done in {time.time()-t0:.1f}s, range [{P_bug.min():.3e}, {P_bug.max():.3e}]")

print("Computing fixed theory (corrected c_3 with -sin)...")
t0 = time.time()
P_fix = theory_P_on_lattice(fix_c3=True)
print(f"  done in {time.time()-t0:.1f}s, range [{P_fix.min():.3e}, {P_fix.max():.3e}]")

# Compare
rms_bug = np.sqrt(np.mean((P_kmc - P_bug) ** 2))
rms_fix = np.sqrt(np.mean((P_kmc - P_fix) ** 2))
max_bug = np.max(np.abs(P_kmc - P_bug))
max_fix = np.max(np.abs(P_kmc - P_fix))
mc_noise_floor = np.sqrt(P_kmc.max() / N_walkers)

print("\n" + "=" * 64)
print("THE TEST")
print("=" * 64)
print(f"  KMC vs BUGGY theory:  RMS = {rms_bug:.3e},  max = {max_bug:.3e}")
print(f"  KMC vs FIXED theory:  RMS = {rms_fix:.3e},  max = {max_fix:.3e}")
print(f"  MC noise floor (~√(P_peak/N)): {mc_noise_floor:.3e}")
print()
if rms_fix < rms_bug * 0.5:
    print(f"  ✓ FIXED is closer to KMC than BUGGY by {rms_bug/rms_fix:.1f}× — sign fix CONFIRMED.")
else:
    print(f"  ✗ FIXED is not significantly better than BUGGY — re-examine the diagnosis.")

# Plot
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
ax_kmc, ax_bug, ax_fix, ax_diff_bug, ax_diff_fix, ax_ratio = axes.flatten()

vmax = max(P_kmc.max(), P_bug.max(), P_fix.max())
for ax, P, ttl in [(ax_kmc, P_kmc, "KMC (ground truth)"),
                    (ax_bug, P_bug, "Buggy theory (+sin)"),
                    (ax_fix, P_fix, "Fixed theory (−sin)")]:
    im = ax.imshow(P, origin="lower", cmap="hot", vmin=0, vmax=vmax)
    ax.set_title(ttl); ax.set_xlabel("n1"); ax.set_ylabel("n2")
    plt.colorbar(im, ax=ax)

# Difference plots
diff_bug = P_bug - P_kmc
diff_fix = P_fix - P_kmc
vlim = max(np.max(np.abs(diff_bug)), np.max(np.abs(diff_fix)))

for ax, D, ttl in [(ax_diff_bug, diff_bug, f"BUGGY − KMC (RMS={rms_bug:.2e})"),
                    (ax_diff_fix, diff_fix, f"FIXED − KMC (RMS={rms_fix:.2e})")]:
    im = ax.imshow(D, origin="lower", cmap="RdBu_r", vmin=-vlim, vmax=vlim)
    ax.set_title(ttl); ax.set_xlabel("n1"); ax.set_ylabel("n2")
    plt.colorbar(im, ax=ax)

# Center cross-section
ax_ratio.plot(np.arange(L), P_kmc[L // 2, :], "ko-", ms=3, label="KMC", lw=0.8)
ax_ratio.plot(np.arange(L), P_bug[L // 2, :], "r-", label="buggy")
ax_ratio.plot(np.arange(L), P_fix[L // 2, :], "b-", label="fixed")
ax_ratio.set_xlabel("n1"); ax_ratio.set_ylabel("P")
ax_ratio.set_title(f"Cross-section at n2 = {L//2}")
ax_ratio.legend(); ax_ratio.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("kmc_verification.png", dpi=180, bbox_inches="tight")
print("\nSaved kmc_verification.png")
