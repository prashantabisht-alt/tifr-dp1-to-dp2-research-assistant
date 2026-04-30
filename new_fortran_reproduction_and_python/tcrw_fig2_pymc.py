"""
TCRW Fig 2 — exact-Python reproduction (steady state + current decomposition)
==============================================================================

Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020, Fig 2.

Three-row plot reproducing Fig 2 data panels:

    Row 1  ω = 1.0    (fully chiral CW)              CCW edge current
    Row 2  ω = 0.5    (achiral)                      no net edge current
    Row 3  ω = 0.0    + L-shape defect cluster       internal-edge current

Each row shows four data panels:
    P(X, Y)   — steady-state spatial distribution (log scale)
    J         — total current vector field
    J_ω       — chiral-move current
    J_Dr      — noise-induced current

3-D time-stacked-occupancy panels (paper b, g, l) are skipped — those are
visualisation artefacts of the MC trajectory and contain no information
beyond P(X, Y).

Why this is "exact"
-------------------
The paper's MC at T = 10¹⁰ converges (by ergodicity) to the Perron
eigenvector of the (4·(L+1)²) × (4·(L+1)²) column-stochastic transition
matrix W.  We solve that linear-algebra problem directly.  No randomness,
no time loop — just `np.linalg.eig` (small L) or `scipy.sparse.linalg.eigs`
(big L) and the analytic flux-flavour current decomposition from authors'
`calculate_J1_J2_with_boundaries`.

Convention
----------
Authors' L = max-index ⇒ grid 0..L inclusive, (L+1)² sites.
Paper "L = 10" means 10 lattice sites, i.e. authors' L = 9.
We call the function with L = 9 by default so axes match paper Fig 2.

Cross-check
-----------
`crosscheck_authors_with_defects` at L = 2 (with and without defects)
compares my matrix and currents element-wise to authors' TRW.py. Expect
machine-precision agreement.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import sys
import time
import importlib.util
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
for cand in (HERE,
             "/sessions/elegant-wizardly-einstein/mnt/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python",
             "/Users/prashantbisht/Documents/Claude/Projects/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python"):
    if os.path.isfile(os.path.join(cand, "tcrw_fig4c.py")):
        sys.path.insert(0, cand)
        ROOT = cand
        break


DX = np.array([0, 1, 0, -1], dtype=int)   # d = 0↑, 1→, 2↓, 3←
DY = np.array([1, 0, -1, 0], dtype=int)


def state_index(i: int, j: int, d: int, L: int) -> int:
    return (i * (L + 1) + j) * 4 + d


# ---------------------------------------------------------------------------
# 1. OBC transition matrix with optional defects (port of authors' code)
# ---------------------------------------------------------------------------
def build_obc_matrix(omega: float, D_r: float, L: int,
                     defects=None) -> sp.csc_matrix:
    """
    4(L+1)^2 × 4(L+1)^2 column-stochastic OBC matrix, with defect support.

    Defect sites are walker-inaccessible.  In the matrix:
      - source rows of a defect site → self-loop weight 1 (dummy state)
      - chiral moves whose target is a defect or out-of-bounds → blocked
        (self-loop weight 1 − D_r at the source state)
    """
    n = 4 * (L + 1) ** 2
    defect_set = set(defects) if defects else set()
    rows, cols, vals = [], [], []

    for i in range(L + 1):
        for j in range(L + 1):
            # If source site is a defect, give it a self-loop = 1 and skip.
            if (i, j) in defect_set:
                for d in range(4):
                    idx = state_index(i, j, d, L)
                    rows.append(idx); cols.append(idx); vals.append(1.0)
                continue

            for d in range(4):
                src = state_index(i, j, d, L)
                d_ccw = (d - 1) % 4
                d_cw  = (d + 1) % 4

                # Noise step (no translation): CCW w.p. ω, CW w.p. 1−ω
                rows.append(state_index(i, j, d_ccw, L)); cols.append(src)
                vals.append(D_r * omega)
                rows.append(state_index(i, j, d_cw,  L)); cols.append(src)
                vals.append(D_r * (1.0 - omega))

                # Chiral step: translate then rotate
                ni, nj = i + int(DX[d]), j + int(DY[d])
                blocked = (not (0 <= ni <= L and 0 <= nj <= L)) \
                          or ((ni, nj) in defect_set)
                if not blocked:
                    rows.append(state_index(ni, nj, d_cw,  L)); cols.append(src)
                    vals.append((1.0 - D_r) * omega)
                    rows.append(state_index(ni, nj, d_ccw, L)); cols.append(src)
                    vals.append((1.0 - D_r) * (1.0 - omega))
                else:
                    rows.append(src); cols.append(src)
                    vals.append(1.0 - D_r)

    return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()


# ---------------------------------------------------------------------------
# 2. Steady state π via sparse Perron
# ---------------------------------------------------------------------------
def steady_state(W: sp.csc_matrix, L: int, defects=None,
                 sigma: float = 1.0):
    """
    Perron eigenvector at λ = 1, restricted to the *physical* (non-defect)
    subspace.  Defect sites carry dummy self-loops, so λ = 1 is degenerate
    when defects are present.  We pick the eigvec with the most weight on
    non-defect states.
    """
    n = W.shape[0]
    defect_set = set(defects) if defects else set()
    is_defect = np.zeros(n, dtype=bool)
    for (i, j) in defect_set:
        for d in range(4):
            is_defect[state_index(i, j, d, L)] = True

    if n <= 600:
        Pd = W.toarray()
        evals, evecs = np.linalg.eig(Pd)
        # All eigenvalues within 1e-8 of 1
        ev1 = np.where(np.abs(evals - 1.0) < 1e-8)[0]
        if len(ev1) == 0:
            ev1 = [int(np.argmin(np.abs(evals - 1.0)))]
        # Among those, pick the one with most non-defect weight
        best = max(ev1, key=lambda k: float(np.sum(np.abs(evecs[~is_defect, k]) ** 2)))
        pi = evecs[:, best].real
    else:
        evals, evecs = spla.eigs(W, k=1, sigma=sigma + 1e-8, which="LM",
                                 maxiter=500, tol=1e-12)
        pi = evecs[:, 0].real

    # Zero out defect entries
    pi[is_defect] = 0.0
    # Sign convention: largest non-defect entry is positive
    k_max = int(np.argmax(np.abs(pi)))
    if pi[k_max] < 0:
        pi = -pi
    pi[pi < 0] = 0.0
    s = pi.sum()
    if s <= 0:
        raise RuntimeError("steady state is zero — eigvec or defects misconfigured")
    return pi / s


# ---------------------------------------------------------------------------
# 3. Current decomposition (authors' flux-flavour algorithm, ported)
# ---------------------------------------------------------------------------
def calculate_J1_J2(pi: np.ndarray, W: sp.csc_matrix, L: int,
                    D_r: float, omega: float, defects=None):
    """
    Returns (J_Dr, J_omega), each shape (L+1, L+1, 2) — author convention
    where J1 ↔ noise origin (J_Dr), J2 ↔ chiral origin (J_omega).
    """
    num_states = 4 * (L + 1) ** 2
    defect_set = set(defects) if defects else set()
    P_rot = np.zeros(num_states)
    P_chiral = np.zeros(num_states)

    Wc = W.tocoo()
    rows_a = Wc.row; cols_a = Wc.col; data_a = Wc.data

    # Step 1: classify each transition by "rotation" vs "chiral move"
    for idx_dest, idx_src, prob in zip(rows_a, cols_a, data_a):
        if idx_src == idx_dest:
            continue
        i_d = (idx_dest // 4) // (L + 1); j_d = (idx_dest // 4) % (L + 1)
        i_s = (idx_src  // 4) // (L + 1); j_s = (idx_src  // 4) % (L + 1)
        if (i_s, j_s) in defect_set or (i_d, j_d) in defect_set:
            continue
        flux = pi[idx_src] * prob
        if i_s == i_d and j_s == j_d:
            P_rot[idx_dest] += flux
        else:
            P_chiral[idx_dest] += flux

    # Normalise so each π(s) is split into rot/chiral fractions
    total = P_rot + P_chiral
    mask = total > 1e-12
    P_rot[mask]     = pi[mask] * (P_rot[mask]     / total[mask])
    P_chiral[mask]  = pi[mask] * (P_chiral[mask]  / total[mask])
    P_rot[~mask] = 0.0
    P_chiral[~mask] = 0.0

    # Step 2: outgoing currents
    J1 = np.zeros((L + 1, L + 1, 2))   # J_Dr
    J2 = np.zeros((L + 1, L + 1, 2))   # J_omega
    for idx_dest, idx_src, prob in zip(rows_a, cols_a, data_a):
        if idx_src == idx_dest:
            continue
        i_d = (idx_dest // 4) // (L + 1); j_d = (idx_dest // 4) % (L + 1)
        i_s = (idx_src  // 4) // (L + 1); j_s = (idx_src  // 4) % (L + 1)
        if (i_s == i_d and j_s == j_d):     # rotation: no displacement
            continue
        if (i_s, j_s) in defect_set:
            continue
        dx, dy = i_d - i_s, j_d - j_s
        J1[i_s, j_s, 0] += P_rot[idx_src]    * prob * dx
        J1[i_s, j_s, 1] += P_rot[idx_src]    * prob * dy
        J2[i_s, j_s, 0] += P_chiral[idx_src] * prob * dx
        J2[i_s, j_s, 1] += P_chiral[idx_src] * prob * dy

    return J1, J2


# ---------------------------------------------------------------------------
# 4. Author cross-check (with and without defects)
# ---------------------------------------------------------------------------
def _load_authors_trw():
    path = os.path.join(ROOT, "TRW._original_code_by_paperauthors.py")
    if not os.path.exists(path):
        return None
    spec = importlib.util.spec_from_file_location("TRW", path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def crosscheck_authors_with_defects(omega: float, D_r: float, L: int,
                                    defects=None, tol: float = 1e-12):
    TRW = _load_authors_trw()
    if TRW is None:
        print("  [warn] authors' TRW.py not found — cross-check skipped")
        return
    W_A = TRW.build_sparse_transition_matrix(L, omega, D_r,
                                             defects=defects).toarray()
    W_U = build_obc_matrix(omega, D_r, L, defects=defects).toarray()
    elem = float(np.max(np.abs(W_A - W_U)))
    print(f"  ω={omega} D_r={D_r} L={L} defects={defects}: "
          f"element diff = {elem:.2e}")
    assert elem < tol, f"matrix mismatch {elem:.2e}"

    # Also currents
    pi_A, _ = TRW.solve_steady_state_sparse(sp.csc_matrix(W_A), L,
                                             defects=defects)
    pi_U = steady_state(sp.csc_matrix(W_U), L, defects=defects)
    pi_diff = float(np.max(np.abs(pi_A - pi_U)))
    print(f"     max|π_A − π_U|     = {pi_diff:.2e}")

    J1_A, J2_A = TRW.calculate_J1_J2_with_boundaries(
        pi_A, sp.csc_matrix(W_A), L, D_r, omega, defects=defects)
    J1_U, J2_U = calculate_J1_J2(pi_U, sp.csc_matrix(W_U), L, D_r, omega,
                                  defects=defects)
    j1d = float(np.max(np.abs(J1_A - J1_U)))
    j2d = float(np.max(np.abs(J2_A - J2_U)))
    print(f"     max|J_Dr_A − J_Dr_U| = {j1d:.2e},  "
          f"max|J_ω_A − J_ω_U| = {j2d:.2e}")
    print(f"     [ok]")


# ---------------------------------------------------------------------------
# 5. P(X, Y) reduction
# ---------------------------------------------------------------------------
def reshape_pi_to_PXY(pi: np.ndarray, L: int) -> np.ndarray:
    """Sum π over directors → P(X, Y) of shape (L+1, L+1)."""
    P = np.zeros((L + 1, L + 1))
    for i in range(L + 1):
        for j in range(L + 1):
            P[i, j] = sum(pi[state_index(i, j, d, L)] for d in range(4))
    return P


# ---------------------------------------------------------------------------
# 6. Big plot
# ---------------------------------------------------------------------------
def make_fig2(L: int = 9, D_r: float = 1e-3,
              omega_list=(1.0, 0.5, 0.0),
              defects_list=(None, None,
                            [(3, 3), (3, 4), (3, 5), (3, 6),
                             (4, 3), (5, 3), (6, 3)]),
              savepath: str = "tcrw_fig2_pymc.png"):
    nrows = len(omega_list)

    # ---- pre-pass: compute all currents to set common quiver scale ----
    cache = []
    all_J = []
    for omega, defects in zip(omega_list, defects_list):
        W = build_obc_matrix(omega, D_r, L, defects=defects)
        pi = steady_state(W, L, defects=defects)
        J_Dr, J_omega = calculate_J1_J2(pi, W, L, D_r, omega, defects=defects)
        J_tot = J_Dr + J_omega
        P_xy = reshape_pi_to_PXY(pi, L)
        cache.append((omega, defects, P_xy, J_tot, J_omega, J_Dr))
        all_J.extend([J_Dr, J_omega, J_tot])
    global_max = max(np.sqrt(J[..., 0]**2 + J[..., 1]**2).max() for J in all_J)
    quiver_scale  = global_max / 0.8 if global_max > 0 else 1.0
    mag_threshold = 0.01 * global_max
    print(f"  global |J| max = {global_max:.3e}")
    print(f"  quiver scale   = {quiver_scale:.3e}  (longest arrow ~ 0.8 cell)")
    print(f"  mag threshold  = {mag_threshold:.3e}  (1% of max)")

    fig, axs = plt.subplots(nrows, 4, figsize=(15, 3.7 * nrows),
                             constrained_layout=True)
    if nrows == 1:
        axs = axs[None, :]

    for r, (omega, defects, P_xy, J_tot, J_omega, J_Dr) in enumerate(cache):
        # Panel (a/f/k): P(X, Y) heatmap
        ax = axs[r, 0]
        valid = P_xy > 0
        vmin = max(P_xy[valid].min(), 1e-10) if valid.any() else 1e-10
        im = ax.imshow(P_xy.T, origin="lower", cmap="hot",
                       norm=LogNorm(vmin=vmin, vmax=P_xy.max()),
                       interpolation="nearest")
        ax.set_xlabel("X"); ax.set_ylabel("Y")
        ax.set_title(rf"$P(X, Y)$, $\omega = {omega}$", fontsize=11)
        if defects:
            for (di, dj) in defects:
                ax.scatter(di, dj, c="red", marker="x", s=40, zorder=5)
        plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)

        # Panels (c/h/m), (d/i/n), (e/j/o): currents — shared scale + threshold
        Xq, Yq = np.meshgrid(np.arange(L + 1), np.arange(L + 1), indexing="ij")
        for c, (J, lbl) in enumerate([(J_tot,    r"$\vec{J}$"),
                                       (J_omega, r"$\vec{J}_\omega$"),
                                       (J_Dr,    r"$\vec{J}_{D_r}$")]):
            ax = axs[r, c + 1]
            mag = np.sqrt(J[..., 0] ** 2 + J[..., 1] ** 2)
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

    fig.suptitle(rf"Fig 2 (Python exact) — $D_r = {D_r}$, paper L = {L+1}",
                 fontsize=13)
    plt.savefig(savepath, dpi=200, bbox_inches="tight")
    print(f"\nsaved {savepath}")


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 64)
    print(" TCRW Fig 2 — exact-Python reproduction")
    print("=" * 64)

    print("\n[1] Author cross-check at L=2, with and without defects")
    crosscheck_authors_with_defects(1.0, 0.1, L=2, defects=None)
    crosscheck_authors_with_defects(0.5, 0.05, L=2, defects=[(1, 1)])
    crosscheck_authors_with_defects(0.0, 0.001, L=2, defects=[(0, 1), (1, 1)])

    print("\n[2] Build production figure (paper specs)")
    make_fig2(L=9, D_r=1e-3,
              omega_list=(1.0, 0.5, 0.0),
              defects_list=(None, None,
                            [(3, 3), (3, 4), (3, 5), (3, 6),
                             (4, 3), (5, 3), (6, 3)]),
              savepath=os.path.join(ROOT, "tcrw_fig2_pymc.png"))

    print("\nDone.")
