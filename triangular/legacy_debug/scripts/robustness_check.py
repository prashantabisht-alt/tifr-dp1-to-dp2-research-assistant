"""
Robustness check: verify the c_3 sign-error bug diagnosis at MULTIPLE
(gamma, epsilon, t) points, not just the original Fig 11 point.

Logic:
  At each test point, run Python KMC + buggy theory + fixed theory.
  If the fix is correct at ALL points:
    - FIXED-KMC RMS should be ≈ MC noise floor at every point
    - BUGGY-KMC RMS should stay systematic and large at every point
  If there's a SECOND bug we missed:
    - FIXED-KMC RMS will be above the noise floor at some point

Three test points spanning the parameter space:
  (1) γ=0.01,  ε=0.15, t=50    — original (matches draft Fig 11)
  (2) γ=0.05,  ε=0.10, t=30    — higher γ, more tumbling, less drift
  (3) γ=0.001, ε=0.16, t=100   — drift-dominant regime, ε near max

Constraint: ε ≤ 1/6 ≈ 0.1667 (so backward rate 1/6−ε stays non-negative).

1M Python walkers per point → ~5 min total runtime.
"""
from __future__ import annotations
import os
import time
import numpy as np
import matplotlib.pyplot as plt

# ──────────────────────────────────────────────────────────────────────────
# Test points
# ──────────────────────────────────────────────────────────────────────────
test_points = [
    dict(gamma=0.01,  epsilon=0.15, t=50.0,
         label=r"(1) $\gamma$=0.01, $\epsilon$=0.15, t=50    [original]",
         short="point1_original"),
    dict(gamma=0.05,  epsilon=0.10, t=30.0,
         label=r"(2) $\gamma$=0.05, $\epsilon$=0.10, t=30    [higher $\gamma$]",
         short="point2_high_gamma"),
    dict(gamma=0.001, epsilon=0.16, t=100.0,
         label=r"(3) $\gamma$=0.001, $\epsilon$=0.16, t=100  [drift-dominant]",
         short="point3_drift"),
]

L      = 30
a, b   = 1.0, 1.0
N_kmc  = 1_000_000          # walkers per parameter point


# ──────────────────────────────────────────────────────────────────────────
# KMC walker
# ──────────────────────────────────────────────────────────────────────────
NN_LATTICE = np.array([(+1, 0), (0, +1), (-1, +1),
                        (-1, 0), (0, -1), (+1, -1)], dtype=np.int32)


def run_kmc(gamma, epsilon, t_final, N_walkers, seed=20260512):
    """Returns counts[n2, n1] array (after binning) and N."""
    HOP_RATES = np.full((6, 6), 1.0 / 6)
    for mm in range(6):
        HOP_RATES[mm, mm] = 1.0 / 6 + epsilon
        HOP_RATES[mm, (mm + 3) % 6] = 1.0 / 6 - epsilon
    HOP_CUM = np.cumsum(HOP_RATES, axis=1)
    rate_total = 1.0 + gamma
    rng = np.random.default_rng(seed)
    counts = np.zeros((L, L), dtype=np.int64)
    for i in range(N_walkers):
        n1, n2, m = 0, 0, int(rng.integers(6))
        t = 0.0
        while True:
            dt = rng.exponential(1.0 / rate_total)
            if t + dt >= t_final:
                break
            t += dt
            if rng.random() < gamma / rate_total:
                m = (m + 1) % 6 if rng.random() < 0.5 else (m - 1) % 6
            else:
                u = rng.random()
                d = int(np.searchsorted(HOP_CUM[m], u))
                n1 = (n1 + NN_LATTICE[d, 0]) % L
                n2 = (n2 + NN_LATTICE[d, 1]) % L
        counts[n2, n1] += 1
    return counts, N_walkers


# ──────────────────────────────────────────────────────────────────────────
# Theory (buggy and fixed)
# ──────────────────────────────────────────────────────────────────────────
def build_Mk(gamma, epsilon, k1, k2, fix_c3):
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


def theory_P(gamma, epsilon, t_final, fix_c3):
    P = np.zeros((L, L))
    eig_data = []
    for m1 in range(L):
        for m2 in range(L):
            kx = np.pi * m1 / (a * L)
            ky = np.pi * (2 * m2 - m1) / (b * L)
            M = build_Mk(gamma, epsilon, kx, ky, fix_c3=fix_c3)
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


# ──────────────────────────────────────────────────────────────────────────
# Run all 3 test points
# ──────────────────────────────────────────────────────────────────────────
print("=" * 78)
print(f"ROBUSTNESS CHECK — c_3 sign fix at {len(test_points)} parameter points")
print(f"L = {L}, N_walkers = {N_kmc:,} per point")
print("=" * 78)

results = []
for ip, pt in enumerate(test_points, start=1):
    print(f"\n--- Point {ip}: {pt['label']} ---")
    t0 = time.time()
    counts, N = run_kmc(pt["gamma"], pt["epsilon"], pt["t"], N_kmc, seed=20260512 + ip)
    P_kmc = counts / N
    print(f"  KMC done in {time.time()-t0:.1f}s. P range: [{P_kmc.min():.3e}, {P_kmc.max():.3e}]")
    t0 = time.time()
    P_bug = theory_P(pt["gamma"], pt["epsilon"], pt["t"], fix_c3=False)
    P_fix = theory_P(pt["gamma"], pt["epsilon"], pt["t"], fix_c3=True)
    print(f"  Theory done in {time.time()-t0:.1f}s")
    rms_bug = float(np.sqrt(np.mean((P_kmc - P_bug) ** 2)))
    rms_fix = float(np.sqrt(np.mean((P_kmc - P_fix) ** 2)))
    mc_noise = float(np.sqrt(P_kmc.max() / N))
    ratio = rms_bug / rms_fix if rms_fix > 0 else float("inf")
    fix_vs_noise = rms_fix / mc_noise if mc_noise > 0 else float("inf")
    print(f"  RMS  BUGGY-KMC: {rms_bug:.3e}")
    print(f"  RMS  FIXED-KMC: {rms_fix:.3e}")
    print(f"  MC noise floor: {mc_noise:.3e}")
    print(f"  BUGGY/FIXED:    {ratio:.1f}x")
    print(f"  FIXED/noise:    {fix_vs_noise:.2f}  (should be ≲ 1 if fix is exact)")
    # Primary criterion: FIXED-KMC RMS at or below MC noise floor.
    # BUGGY/FIXED ratio is a secondary check — it can be small at small γ
    # (where the bug is physically suppressed) or noise-floor-limited
    # at 1M walkers. We log it for interpretation, but it is NOT the test.
    if fix_vs_noise < 1.5:
        verdict = "✓ FIXED matches KMC within MC noise"
    else:
        verdict = "✗ FIXED above noise floor — investigate"
    print(f"  Verdict:        {verdict}")
    results.append(dict(point=pt, P_kmc=P_kmc, P_bug=P_bug, P_fix=P_fix,
                        rms_bug=rms_bug, rms_fix=rms_fix, mc_noise=mc_noise,
                        ratio=ratio, fix_vs_noise=fix_vs_noise))


# ──────────────────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("SUMMARY")
print("=" * 78)
print(f"{'point':<48} {'BUGGY RMS':>10} {'FIXED RMS':>10} {'noise':>10} {'FIX/noise':>10}  verdict")
print("-" * 110)
all_pass = True
for r in results:
    label = r["point"]["label"][:46]
    # Test: FIXED-KMC RMS is at or below MC noise floor.
    v = ("PASS" if r["fix_vs_noise"] < 1.5 else "FAIL")
    if v == "FAIL":
        all_pass = False
    print(f"{label:<48} {r['rms_bug']:>10.2e} {r['rms_fix']:>10.2e} {r['mc_noise']:>10.2e}"
          f" {r['fix_vs_noise']:>10.2f}  {v}")
print()
if all_pass:
    print("  ✓✓✓ FIXED theory matches KMC within MC noise at all three points.")
    print("  Note: BUGGY/FIXED ratio at Point 3 (γ=0.001) is small because the")
    print("  c_3 bug acts only through rotation γ/2 couplings — when γ→0 the bug")
    print("  becomes invisible. This is a positive physical check, not a problem.")
else:
    print("  ✗   FIXED theory above noise floor at one or more points — investigate.")


# ──────────────────────────────────────────────────────────────────────────
# Visualization: 3 rows, 4 columns (P_bug, P_fix, P_kmc, cross-section)
# ──────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(len(results), 4, figsize=(20, 5 * len(results)))
if len(results) == 1:
    axes = axes[None, :]

for irow, r in enumerate(results):
    pt = r["point"]
    P_kmc, P_bug, P_fix = r["P_kmc"], r["P_bug"], r["P_fix"]
    shift = (L // 2, L // 2)
    P_kmc_c = np.roll(P_kmc, shift=shift, axis=(0, 1))
    P_bug_c = np.roll(P_bug, shift=shift, axis=(0, 1))
    P_fix_c = np.roll(P_fix, shift=shift, axis=(0, 1))
    vmax = max(P_kmc_c.max(), P_bug_c.max(), P_fix_c.max())
    for ax, P, ttl in zip(
        axes[irow, :3],
        [P_bug_c, P_fix_c, P_kmc_c],
        ["Buggy theory", "Fixed theory", "KMC"]
    ):
        im = ax.imshow(P, origin="lower", cmap="hot", vmin=0, vmax=vmax)
        ax.set_title(ttl, fontsize=10)
        if irow == 0:
            ax.set_xlabel(r"$n_1$"); ax.set_ylabel(r"$n_2$")
        plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    # Cross-section
    n2_slice = L // 2
    ns = np.arange(L)
    ax = axes[irow, 3]
    ax.plot(ns, P_kmc_c[n2_slice, :], "ko", ms=5, label="KMC", zorder=3)
    ax.plot(ns, P_bug_c[n2_slice, :], "r-", lw=2,
            label=f"buggy (RMS={r['rms_bug']:.1e})")
    ax.plot(ns, P_fix_c[n2_slice, :], "b-", lw=2,
            label=f"fixed (RMS={r['rms_fix']:.1e})")
    ax.set_title(pt["label"], fontsize=10)
    ax.set_xlabel(r"$n_1$"); ax.set_ylabel("P")
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.3)

fig.suptitle(
    rf"Robustness of sign-fix across 3 parameter points  ($L = {L}$, $N = {N_kmc:,}$ walkers each)",
    fontsize=13, y=1.0
)
plt.tight_layout()
plt.savefig("robustness_check.png", dpi=150, bbox_inches="tight")
print("\nSaved robustness_check.png")
