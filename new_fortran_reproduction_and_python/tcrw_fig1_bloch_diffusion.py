"""
Exact D(ω, D_r) via small-k Bloch expansion + Fortran MC cross-check
====================================================================

For an unbounded TCRW the long-time diffusion coefficient D(ω, D_r) is
determined by the small-k expansion of the slowest Bloch eigenvalue of
the 4×4 Fourier transition matrix P(k) (paper Methods Eq 1):

    λ₁(k) = 1 - D(ω, D_r) |k|² + O(|k|⁴)
    ⇒  D = lim_{k→0} (1 - Re λ₁(k)) / |k|²

Since `build_Pk` is verified element-wise = authors' real-space operator
in Fourier space (fig4b crosscheck, Hausdorff ≤ 1.2e-14), the D value
extracted this way is the **exact** paper-target D for this model — no
Monte Carlo, no fitting, no statistical noise.

Use this script to cross-check the Fortran Fig 1(d) MC:
    Fortran outputs `tcrw_fig1d_D_dr2.txt` and `tcrw_fig1d_D_dr3.txt`
    (one row per ω with the fitted D from MSD-vs-t slope).
The MC values must converge to the exact Bloch values as T → ∞.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Locate fig4b_paper.py and import build_Pk
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
for cand in (HERE,
             "/sessions/elegant-wizardly-einstein/mnt/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python",
             "/Users/prashantbisht/Documents/Claude/Projects/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python"):
    if os.path.isfile(os.path.join(cand, "tcrw_fig4b_paper.py")):
        sys.path.insert(0, cand)
        ROOT = cand
        break

from tcrw_fig4b_paper import build_Pk     # 4×4 Bloch matrix


# ---------------------------------------------------------------------------
# 1. Exact diffusion coefficient at a single (ω, D_r)
# ---------------------------------------------------------------------------
def D_exact(omega: float, D_r: float, k_eps: float = 1e-2,
            return_diag: bool = False):
    """
    Extract D from the slowest Bloch eigenvalue at small k.

    Computes the eigenvalue closest to 1 at k = (k_eps, 0) and at
    k = (0, k_eps); averages the two for robustness (they should agree
    by 4-fold rotational symmetry of the model).

    Parameters
    ----------
    omega, D_r : floats in [0, 1]
    k_eps      : small momentum at which to evaluate.  Default 1e-2 is
                 numerically cleanest: shift (1 - λ₁) is of order D·1e-4,
                 well above np.linalg.eig noise (~1e-15) yet still small
                 enough for the leading-order |k|² to dominate.
    return_diag: if True, also return D from a second k (k_eps/2) — use
                 to verify Richardson convergence.

    Returns
    -------
    D_val (float)   : best estimate of D(ω, D_r)
    diag (dict)     : if return_diag=True, dict with D_x, D_y, D_x_half, D_y_half
    """
    def _D_at(kx, ky):
        P = build_Pk(omega, D_r, kx, ky)
        ev = np.linalg.eigvals(P)
        idx = np.argmin(np.abs(ev - 1.0))
        return (1.0 - ev[idx].real) / (kx * kx + ky * ky)

    D_x = _D_at(k_eps, 0.0)
    D_y = _D_at(0.0,   k_eps)
    D_val = 0.5 * (D_x + D_y)

    if return_diag:
        D_x_half = _D_at(k_eps / 2, 0.0)
        D_y_half = _D_at(0.0, k_eps / 2)
        return D_val, dict(D_x=D_x, D_y=D_y,
                           D_x_half=D_x_half, D_y_half=D_y_half,
                           D_half=0.5 * (D_x_half + D_y_half))
    return D_val


# ---------------------------------------------------------------------------
# 2. Sanity: convergence vs k_eps + 4-fold symmetry
# ---------------------------------------------------------------------------
def validate_small_k(omegas=(0.0, 0.3, 0.5, 0.7, 1.0), D_r=1e-2):
    """
    Make sure D extracted at k_eps = {1e-1, 1e-2, 1e-3} agrees, and that
    D_x ≡ D_y (the model has 4-fold spatial symmetry, so D is isotropic).
    """
    print("Convergence test (D_r = {:.0e}):".format(D_r))
    print(f"  {'ω':>5}  {'D(k=1e-1)':>11}  {'D(k=1e-2)':>11}  {'D(k=1e-3)':>11}  "
          f"{'rel diff k1↔k2':>15}  {'D_x − D_y':>11}")
    for w in omegas:
        Ds = []
        for k in (1e-1, 1e-2, 1e-3):
            P = build_Pk(w, D_r, k, 0.0)
            ev = np.linalg.eigvals(P); idx = np.argmin(np.abs(ev - 1.0))
            Dx = (1.0 - ev[idx].real) / (k * k)
            P = build_Pk(w, D_r, 0.0, k)
            ev = np.linalg.eigvals(P); idx = np.argmin(np.abs(ev - 1.0))
            Dy = (1.0 - ev[idx].real) / (k * k)
            Ds.append(0.5 * (Dx + Dy))
        rel = abs(Ds[1] - Ds[2]) / max(Ds[1], 1e-30)

        # symmetry at k_eps = 1e-2
        P = build_Pk(w, D_r, 1e-2, 0.0); ev = np.linalg.eigvals(P)
        Dx_check = (1 - ev[np.argmin(np.abs(ev - 1))].real) / 1e-4
        P = build_Pk(w, D_r, 0.0, 1e-2); ev = np.linalg.eigvals(P)
        Dy_check = (1 - ev[np.argmin(np.abs(ev - 1))].real) / 1e-4
        sym = Dx_check - Dy_check

        print(f"  {w:>5.2f}  {Ds[0]:>11.6f}  {Ds[1]:>11.6f}  {Ds[2]:>11.6f}  "
              f"{rel:>15.2e}  {sym:>11.2e}")


# ---------------------------------------------------------------------------
# 3. Load Fortran fig1d summary file (handles its quirky scientific notation)
# ---------------------------------------------------------------------------
def load_fig1d_summary(path: str):
    """columns: omega   D_fit   D_err   slope   slope_err   D_end"""
    rows = []
    D_r = None
    with open(path) as f:
        for line in f:
            if "D_r =" in line and D_r is None:
                # parse '# TCRW Fig 1(d) | D_r =  1.000E-03 ...'
                tok = line.split("D_r =")[1].split()[0]
                D_r = float(tok)
                continue
            if line.lstrip().startswith("#") or not line.strip():
                continue
            rows.append([float(x) for x in line.split()])
    arr = np.array(rows, dtype=float)
    return D_r, arr


# ---------------------------------------------------------------------------
# 4. Cross-check + plot
# ---------------------------------------------------------------------------
def crosscheck(summary_files, savepath="tcrw_fig1d_crosscheck.png"):
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = ["tab:purple", "tab:pink", "tab:orange"]
    omega_grid = np.linspace(0.0, 1.0, 51)

    print("\nCross-check Fortran MC vs exact Bloch:")
    for color, path in zip(colors, summary_files):
        D_r, arr = load_fig1d_summary(path)
        if arr.size == 0:
            continue
        omegas = arr[:, 0]
        D_mc   = arr[:, 1]
        D_err  = arr[:, 2]

        # exact Bloch on a denser grid for the line
        D_line = np.array([D_exact(w, D_r) for w in omega_grid])
        # exact at the SAME ω as the MC, for a numerical comparison
        D_at_mc = np.array([D_exact(float(w), float(D_r)) for w in omegas])

        ax.errorbar(omegas, D_mc, yerr=D_err, fmt="o",
                    color=color, ms=6, capsize=3, label=f"MC  D_r={D_r:.0e}")
        ax.plot(omega_grid, D_line, color=color, lw=1.6, alpha=0.7,
                label=f"Bloch exact  D_r={D_r:.0e}")

        rel = np.max(np.abs(D_mc - D_at_mc) / np.maximum(D_at_mc, 1e-12))
        med = np.median(np.abs(D_mc - D_at_mc) / np.maximum(D_at_mc, 1e-12))
        print(f"  D_r = {D_r:.0e}:  max rel diff = {rel:.2e},  median = {med:.2e}")

    ax.set_xlabel(r"$\omega$"); ax.set_ylabel(r"$D(\omega, D_r)$")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(0, None)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.legend(fontsize=9, ncol=2)
    ax.set_title("Fig 1(d) — Fortran MC (markers) vs exact Bloch (lines)")
    ax.grid(alpha=0.3)
    plt.tight_layout(); plt.savefig(savepath, dpi=200, bbox_inches="tight")
    print(f"  saved {savepath}")


# ---------------------------------------------------------------------------
# 5. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 64)
    print(" TCRW Fig 1(d) cross-check: MC vs exact Bloch diffusion coefficient")
    print("=" * 64)

    print("\n[1] Convergence + symmetry of small-k extrapolation")
    validate_small_k(omegas=(0.0, 0.2, 0.5, 0.8, 1.0), D_r=1e-2)
    validate_small_k(omegas=(0.0, 0.2, 0.5, 0.8, 1.0), D_r=1e-3)

    print("\n[2] Cross-check Fortran MC against exact Bloch")
    files = []
    for tag in ("dr3", "dr2"):
        f = os.path.join(ROOT, f"tcrw_fig1d_D_{tag}.txt")
        if os.path.exists(f):
            files.append(f)
        else:
            print(f"  [skip] {f} not found")
    if files:
        crosscheck(files, os.path.join(ROOT, "tcrw_fig1d_crosscheck.png"))

    print("\nDone.")
