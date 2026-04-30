"""
TCRW Fig 2 — using authors' TRW.py directly (no re-implementation)
===================================================================

Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020.

This script reproduces the data panels of Fig 2 by calling the authors'
own functions verbatim from `TRW._original_code_by_paperauthors.py`:

    TRW.build_sparse_transition_matrix(L, ω, D_r, defects)
    TRW.solve_steady_state_sparse(W, L, defects)
    TRW.calculate_J1_J2_with_boundaries(π, W, L, D_r, ω, defects)

These are the exact functions that produced the paper's "black line"
exact calculations (the MC dots come from a separate single-walker
loop in the same notebook).  No re-implementation by us — this is
literally the authors' code in action.

Use this script when:
  • You want to demonstrate the paper's exact calculation as-is.
  • You need a sanity reference for our self-contained `tcrw_fig2_pymc.py`
    (which has been bit-verified against this).

Use `tcrw_fig2_pymc.py` when:
  • You want a self-contained file with no TRW.py dependency.
  • You're going to fork it for DP2 (jerky walker) and modify the rules.

Both produce identical numbers to floating-point precision.

Convention
----------
Authors' L = max-index ⇒ L+1 sites per axis.  Paper L = 10 (10 sites)
corresponds to passing L = 9 here.  Defects are (i, j) tuples on the
0..L grid.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import sys
import importlib.util
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
for cand in (HERE,
             "/sessions/elegant-wizardly-einstein/mnt/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python",
             "/Users/prashantbisht/Documents/Claude/Projects/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python"):
    p = os.path.join(cand, "TRW._original_code_by_paperauthors.py")
    if os.path.isfile(p):
        ROOT = cand
        TRW_PATH = p
        break

# Import authors' TRW module by file path (filename has dots, can't import directly)
spec = importlib.util.spec_from_file_location("TRW", TRW_PATH)
TRW = importlib.util.module_from_spec(spec)
spec.loader.exec_module(TRW)


# ---------------------------------------------------------------------------
# 1. Pull P(X, Y), J_Dr, J_ω at given (ω, D_r, L, defects) using AUTHORS' code
# ---------------------------------------------------------------------------
def fig2_data(omega: float, D_r: float, L: int, defects=None):
    """Return P(X,Y), J_Dr, J_omega from authors' TRW.py directly."""
    W = TRW.build_sparse_transition_matrix(L, omega, D_r, defects=defects)
    pi, _ = TRW.solve_steady_state_sparse(W, L, defects=defects)
    J_Dr, J_omega = TRW.calculate_J1_J2_with_boundaries(
        pi, W, L, D_r, omega, defects=defects
    )

    # Reduce π to P(X, Y) by summing over the 4 directors at each site
    P_xy = np.zeros((L + 1, L + 1))
    for i in range(L + 1):
        for j in range(L + 1):
            for d_name in TRW.directions:
                P_xy[i, j] += pi[TRW.index(i, j, d_name, L)]

    return P_xy, J_Dr, J_omega


# ---------------------------------------------------------------------------
# 2. Plot one row of Fig 2 (P heatmap + 3 current quivers)
# ---------------------------------------------------------------------------
def _plot_row(axs_row, omega, D_r, L, defects, title_prefix,
              quiver_scale: float, mag_threshold: float):
    """
    Plot one row.  `quiver_scale` is shared across all panels of all rows
    so arrow length is comparable globally — paper-style.
    `mag_threshold` masks out arrows below this magnitude (suppresses
    bulk-noise residuals at ω = 0.5 and bulk sites at ω ∈ {0, 1}).
    """
    P_xy, J_Dr, J_omega = fig2_data(omega, D_r, L, defects=defects)
    J_tot = J_Dr + J_omega

    # Panel: P(X, Y)
    ax = axs_row[0]
    valid = P_xy > 0
    vmin = max(P_xy[valid].min(), 1e-10) if valid.any() else 1e-10
    vmax = P_xy.max()
    im = ax.imshow(
        P_xy.T, origin="lower", cmap="hot",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        interpolation="nearest",
    )
    ax.set_xlabel("X"); ax.set_ylabel("Y")
    ax.set_title(rf"{title_prefix}  $P(X,Y)$,  $\omega = {omega}$", fontsize=11)
    if defects:
        for (di, dj) in defects:
            ax.scatter(di, dj, c="red", marker="x", s=40, zorder=5)
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)

    # Quiver panels — shared scale, magnitude threshold to hide bulk noise
    Xq, Yq = np.meshgrid(np.arange(L + 1), np.arange(L + 1), indexing="ij")
    for ax, J, lbl in zip(axs_row[1:],
                           (J_tot, J_omega, J_Dr),
                           (r"$\vec{J}$", r"$\vec{J}_\omega$", r"$\vec{J}_{D_r}$")):
        mag = np.sqrt(J[..., 0] ** 2 + J[..., 1] ** 2)
        # mask sites below threshold (zero out so quiver draws nothing)
        keep = mag > mag_threshold
        Jx_plot = np.where(keep, J[..., 0], 0.0)
        Jy_plot = np.where(keep, J[..., 1], 0.0)
        mag_plot = np.where(keep, mag, np.nan)
        ax.quiver(Xq, Yq, Jx_plot, Jy_plot, mag_plot,
                   cmap="copper_r", scale=quiver_scale, width=0.005,
                   angles="xy", scale_units="xy", pivot="mid")
        ax.set_xlabel("X"); ax.set_ylabel("Y")
        ax.set_title(lbl, fontsize=11)
        ax.set_xlim(-0.5, L + 0.5); ax.set_ylim(-0.5, L + 0.5)
        ax.set_aspect("equal")
        if defects:
            for (di, dj) in defects:
                ax.scatter(di, dj, c="red", marker="x", s=40, zorder=5)


def make_fig2(L: int = 9, D_r: float = 1e-3,
              omega_list=(1.0, 0.5, 0.0),
              defects_list=(None, None,
                            [(3, 3), (3, 4), (3, 5), (3, 6),
                             (4, 3), (5, 3), (6, 3)]),
              savepath: str = "tcrw_fig2_authors.png"):
    nrows = len(omega_list)

    # PRE-PASS: compute all panels' currents to set a global quiver scale.
    # This is what the paper does (implicitly) — fixed scale across panels
    # so bulk arrows (tiny) are below visible length while edge arrows
    # (large) show clearly.
    all_J = []
    for omega, defects in zip(omega_list, defects_list):
        _, J_Dr, J_omega = fig2_data(omega, D_r, L, defects=defects)
        all_J.extend([J_Dr, J_omega, J_Dr + J_omega])
    global_max = max(np.sqrt(J[..., 0]**2 + J[..., 1]**2).max() for J in all_J)
    # Scale: longest arrow ≈ 0.8 grid cells.
    quiver_scale = global_max / 0.8 if global_max > 0 else 1.0
    # Threshold: hide arrows below 1% of global max  (suppresses ω=0.5
    # noise + bulk residuals at ω=1, 0).
    mag_threshold = 0.01 * global_max
    print(f"  global |J| max = {global_max:.3e}")
    print(f"  quiver scale   = {quiver_scale:.3e}  (longest arrow ~ 0.8 cell)")
    print(f"  mag threshold  = {mag_threshold:.3e}  (1% of max)")

    fig, axs = plt.subplots(nrows, 4, figsize=(15, 3.7 * nrows),
                             constrained_layout=True)
    if nrows == 1:
        axs = axs[None, :]

    titles = [f"({chr(ord('a') + 5*r)})" for r in range(nrows)]
    for r, (omega, defects, title_prefix) in enumerate(
            zip(omega_list, defects_list, titles)):
        _plot_row(axs[r], omega, D_r, L, defects, title_prefix,
                  quiver_scale=quiver_scale,
                  mag_threshold=mag_threshold)

    fig.suptitle(
        rf"Fig 2 (using authors' TRW.py) — $D_r = {D_r}$, paper L = {L+1}",
        fontsize=13,
    )
    plt.savefig(savepath, dpi=200, bbox_inches="tight")
    print(f"saved {savepath}")


# ---------------------------------------------------------------------------
# 3. Sanity check: confirm we're calling the right authors' code
# ---------------------------------------------------------------------------
def _sanity_print():
    print(f"using authors' TRW from: {TRW_PATH}")
    print(f"  has build_sparse_transition_matrix: {hasattr(TRW, 'build_sparse_transition_matrix')}")
    print(f"  has solve_steady_state_sparse:      {hasattr(TRW, 'solve_steady_state_sparse')}")
    print(f"  has calculate_J1_J2_with_boundaries:{hasattr(TRW, 'calculate_J1_J2_with_boundaries')}")
    print(f"  has directions list:                {hasattr(TRW, 'directions')}  → {TRW.directions}")


# ---------------------------------------------------------------------------
# 4. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 64)
    print(" TCRW Fig 2 — exact, using authors' TRW.py directly")
    print("=" * 64)

    _sanity_print()
    print()

    make_fig2(
        L=9, D_r=1e-3,
        omega_list=(1.0, 0.5, 0.0),
        defects_list=(None, None,
                      [(3, 3), (3, 4), (3, 5), (3, 6),
                       (4, 3), (5, 3), (6, 3)]),
        savepath=os.path.join(ROOT, "tcrw_fig2_authors.png"),
    )

    print("\nDone.")
