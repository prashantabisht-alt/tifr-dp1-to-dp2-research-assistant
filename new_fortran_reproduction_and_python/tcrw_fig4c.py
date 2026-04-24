"""
TCRW Fig 4(c) — OBC spectrum on L = 2 lattice, 3-D surfaces over (D_r, ω)
========================================================================

Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020, Fig 4(c).

Physics
-------
On an OBC grid of linear size L (authors' convention: sites 0..L inclusive,
so (L+1)^2 spatial sites), the one-walker Markov chain is a
4 (L+1)^2 × 4 (L+1)^2 column-stochastic matrix P.  For L = 2 we have
9 sites × 4 directors = 36 states.

Fig 4(c) plots the real and imaginary parts of the 36 eigenvalues of P
as 3-D surfaces over the parameter plane (D_r, ω) ∈ [0, 1]^2, each
point coloured by the edge weight of the corresponding right
eigenvector.  The "edge weight" of a mode is

    w_edge = Σ_{s on boundary site} |v(s)|^2 / Σ_{s} |v(s)|^2

summed over every state s = (i, j, d) whose lattice site (i, j) lies on
the outer frame (i ∈ {0, L} or j ∈ {0, L}, any director d).  Paper
caption (c) uses a single "edge" colorbar, 0 → 1.

Convention
----------
Direction indexing matches the authors' TRW.py exactly:
    d = 0 → ↑  (ΔX, ΔY) = ( 0, +1)
    d = 1 → →  (ΔX, ΔY) = (+1,  0)
    d = 2 → ↓  (ΔX, ΔY) = ( 0, -1)
    d = 3 → ←  (ΔX, ΔY) = (-1,  0)

State index: s = (i * (L + 1) + j) * 4 + d, i.e. site-major / dir-minor,
identical to the authors' `index(i, j, d, L)`.  This choice makes the
author cross-check element-wise (not just eigenvalue-equivalent).

Rotation rule:
    CCW: d → (d - 1) mod 4 = (d + 3) mod 4
    CW : d → (d + 1) mod 4

Transition probabilities (paper / authors):
    Noise  (prob D_r       ):  CCW w.p. ω,     CW w.p. 1-ω       (no move)
    Chiral (prob 1 - D_r   ):  CW  w.p. ω,     CCW w.p. 1-ω      (translate THEN rotate)
    Blocked chiral: self-loop with weight (1 - D_r); no move, no rotation.

Output
------
Single PNG with two side-by-side 3-D `plot_surface` plots (Re, Im).
Colormap: `turbo` (closest matplotlib match to the paper's blue → green →
red scale).  Every one of the 36 eigenvalues is rendered as its own
semi-transparent sheet; the sheets are sorted independently for the Re
panel (by Re(λ)) and the Im panel (by Im(λ)).

Cross-check
-----------
`crosscheck_authors` densifies the authors' sparse transition matrix at
a given (ω, D_r, L) and compares, element-by-element, to our
`build_obc_matrix`.  Expect 0 difference (not just eigenvalue-equivalent).

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import sys
import importlib.util
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3-D projection)

# ---------------------------------------------------------------------------
# Direction vectors — authors' convention
# ---------------------------------------------------------------------------
DX = np.array([0, 1, 0, -1], dtype=int)   # dx for d = 0↑, 1→, 2↓, 3←
DY = np.array([1, 0, -1, 0], dtype=int)


def state_index(i: int, j: int, d: int, L: int) -> int:
    """Authors' state indexing: s = (i * (L + 1) + j) * 4 + d."""
    return (i * (L + 1) + j) * 4 + d


# ---------------------------------------------------------------------------
# 1. OBC transition matrix (self-contained)
# ---------------------------------------------------------------------------
def build_obc_matrix(omega: float, D_r: float, L: int) -> sp.csc_matrix:
    """
    4 (L+1)^2 × 4 (L+1)^2 column-stochastic OBC transition matrix.

    Element-wise identical to authors' TRW.build_sparse_transition_matrix.
    """
    n = 4 * (L + 1) ** 2
    rows, cols, vals = [], [], []

    for i in range(L + 1):
        for j in range(L + 1):
            for d in range(4):
                src = state_index(i, j, d, L)
                d_ccw = (d - 1) % 4
                d_cw  = (d + 1) % 4

                # ---- noise step (stay, rotate) ----
                # CCW with prob D_r·ω
                rows.append(state_index(i, j, d_ccw, L))
                cols.append(src)
                vals.append(D_r * omega)
                # CW  with prob D_r·(1-ω)
                rows.append(state_index(i, j, d_cw, L))
                cols.append(src)
                vals.append(D_r * (1.0 - omega))

                # ---- chiral step (translate then rotate) ----
                ni, nj = i + int(DX[d]), j + int(DY[d])
                if 0 <= ni <= L and 0 <= nj <= L:
                    # CW with prob (1 - D_r)·ω
                    rows.append(state_index(ni, nj, d_cw, L))
                    cols.append(src)
                    vals.append((1.0 - D_r) * omega)
                    # CCW with prob (1 - D_r)·(1 - ω)
                    rows.append(state_index(ni, nj, d_ccw, L))
                    cols.append(src)
                    vals.append((1.0 - D_r) * (1.0 - omega))
                else:
                    # Blocked chiral: self-loop, no move, no rotation
                    rows.append(src)
                    cols.append(src)
                    vals.append(1.0 - D_r)

    return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()


# ---------------------------------------------------------------------------
# 2. Edge mask — paper Fig 4(a) definition
# ---------------------------------------------------------------------------
# At each boundary site, a state s = (x, y, d) is an "edge state" if d is
# either "to-edge" (director points into the adjacent wall) or
# "CCW-along-edge" (director points in the CCW-circulation direction for
# that edge, using ω=1 as the reference chirality).  All other directors
# at a boundary site are bulk (e.g. CW-along-edge, bulk-directed).
# Interior sites contribute zero edge states.
#
# Direction labels:
#   d = 0 ↑   d = 1 →   d = 2 ↓   d = 3 ←
#
# For each edge, the two edge-director indices (to-edge, CCW-along):
#   left   (x = 0): to-edge = ← (d=3), CCW-along = ↓ (d=2)
#   right  (x = L): to-edge = → (d=1), CCW-along = ↑ (d=0)
#   bottom (y = 0): to-edge = ↓ (d=2), CCW-along = → (d=1)
#   top    (y = L): to-edge = ↑ (d=0), CCW-along = ← (d=3)
#
# At corners, two edges apply; the edge-state set is the union.
# ---------------------------------------------------------------------------

_EDGE_DIRS_CCW = {
    # edge name: (to-edge director, CCW-along director)
    "left":   (3, 2),
    "right":  (1, 0),
    "bottom": (2, 1),
    "top":    (0, 3),
}


def edge_mask(L: int) -> np.ndarray:
    """
    Boolean mask of length 4 (L+1)^2 identifying paper-defined edge states
    (to-edge ∪ CCW-along-edge, using ω=1 as the reference chirality).

    At L = 2 (3x3 grid) this yields 20 edge states out of 36:
      4 corners × 3 edge states each   = 12
      4 mid-edges × 2 edge states each = 8
      1 interior site × 0              = 0
    """
    n = 4 * (L + 1) ** 2
    mask = np.zeros(n, dtype=bool)
    for i in range(L + 1):
        for j in range(L + 1):
            applies = []
            if i == 0: applies.append("left")
            if i == L: applies.append("right")
            if j == 0: applies.append("bottom")
            if j == L: applies.append("top")
            for edge in applies:
                d_to, d_along = _EDGE_DIRS_CCW[edge]
                mask[state_index(i, j, d_to,    L)] = True
                mask[state_index(i, j, d_along, L)] = True
    return mask


# ---------------------------------------------------------------------------
# 3. OBC spectrum with edge weights
# ---------------------------------------------------------------------------
def obc_spectrum_with_weights(omega: float, D_r: float, L: int):
    """Dense eig; returns (eigvals, edge_weights)."""
    P = build_obc_matrix(omega, D_r, L).toarray()
    evals, evecs = np.linalg.eig(P)
    abs2 = np.abs(evecs) ** 2
    total = abs2.sum(axis=0)
    total = np.where(total > 1e-30, total, 1.0)
    edge = edge_mask(L)
    w_edge = abs2[edge, :].sum(axis=0) / total
    return evals, w_edge


# ---------------------------------------------------------------------------
# 4. Author cross-check — element-wise
# ---------------------------------------------------------------------------
def _load_authors_trw():
    """Import the authors' TRW.py (filename has dots, so importlib)."""
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(here, "TRW._original_code_by_paperauthors.py"),
        os.path.join(here, "new_fortran_reproduction_and_python",
                     "TRW._original_code_by_paperauthors.py"),
        os.path.join(os.path.dirname(here),
                     "new_fortran_reproduction_and_python",
                     "TRW._original_code_by_paperauthors.py"),
    ]
    for path in candidates:
        if os.path.exists(path):
            spec = importlib.util.spec_from_file_location("TRW_authors", path)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod
    return None


def crosscheck_authors(omega: float, D_r: float, L: int = 2,
                       tol: float = 1e-12, verbose: bool = True):
    """
    Compare our OBC matrix to authors' sparse matrix element-wise.

    Returns (max_elem_diff, eigval_hausdorff).  Raises AssertionError if
    either exceeds `tol`.
    """
    TRW = _load_authors_trw()
    if TRW is None:
        print("  [warn] authors' TRW module not found; skipping cross-check")
        return None, None

    W_A = TRW.build_sparse_transition_matrix(L, omega, D_r).toarray()
    W_U = build_obc_matrix(omega, D_r, L).toarray()
    max_elem = float(np.max(np.abs(W_A - W_U)))

    ev_A = np.linalg.eigvals(W_A)
    ev_U = np.linalg.eigvals(W_U)
    d_UA = float(np.min(np.abs(ev_U[:, None] - ev_A[None, :]), axis=1).max())
    d_AU = float(np.min(np.abs(ev_A[:, None] - ev_U[None, :]), axis=1).max())
    hausdorff = max(d_UA, d_AU)

    if verbose:
        print(f"  crosscheck  ω={omega:.3g}  D_r={D_r:.3g}  L={L}")
        print(f"    element-wise max |W_A - W_U| = {max_elem:.2e}")
        print(f"    eigenvalue Hausdorff         = {hausdorff:.2e}")

    assert max_elem < tol, f"element mismatch {max_elem:.2e} > tol {tol:.2e}"
    assert hausdorff < tol, f"eigval mismatch {hausdorff:.2e} > tol {tol:.2e}"
    print(f"    [ok] matches authors (elem diff = {max_elem:.2e})")
    return max_elem, hausdorff


# ---------------------------------------------------------------------------
# 5. Scan over (ω, D_r) grid
# ---------------------------------------------------------------------------
def scan_grid(L: int = 2, N_w: int = 40, N_D: int = 40):
    """
    Loop over (ω, D_r) ∈ [0, 1]^2; diagonalize at each grid point.

    Returns
    -------
    w_vals       : (N_w,) — ω axis
    D_vals       : (N_D,) — D_r axis
    evals_grid   : (N_w, N_D, n_states)  complex (unsorted, native eig order)
    weights_grid : (N_w, N_D, n_states)  real    (edge weights, aligned with evals_grid)
    """
    w_vals = np.linspace(0.0, 1.0, N_w)
    D_vals = np.linspace(0.0, 1.0, N_D)
    n_states = 4 * (L + 1) ** 2

    evals_grid = np.zeros((N_w, N_D, n_states), dtype=complex)
    weights_grid = np.zeros((N_w, N_D, n_states))

    edge = edge_mask(L)

    print(f"  scanning {N_w} × {N_D} grid   (matrix {n_states}×{n_states})")
    for i, w in enumerate(w_vals):
        for j, D in enumerate(D_vals):
            P = build_obc_matrix(w, D, L).toarray()
            evals, evecs = np.linalg.eig(P)
            abs2 = np.abs(evecs) ** 2
            total = abs2.sum(axis=0)
            total = np.where(total > 1e-30, total, 1.0)
            evals_grid[i, j, :] = evals
            weights_grid[i, j, :] = abs2[edge, :].sum(axis=0) / total
        if (i + 1) % max(1, N_w // 10) == 0:
            print(f"    row {i + 1}/{N_w}")
    return w_vals, D_vals, evals_grid, weights_grid


# ---------------------------------------------------------------------------
# 6. Figure
# ---------------------------------------------------------------------------
def make_figure(N_w: int = 40, N_D: int = 40, L: int = 2,
                savepath: str = "tcrw_fig4c_paper.png"):
    """Build and save the Fig 4(c) figure."""
    w_vals, D_vals, evals, weights = scan_grid(L=L, N_w=N_w, N_D=N_D)

    # Sort each (i, j) slice independently for the two panels.
    # (Each panel's sheets are visually smoother if ordered by the plotted axis.)
    sort_re = np.argsort(evals.real, axis=2, kind="stable")
    sort_im = np.argsort(evals.imag, axis=2, kind="stable")

    evals_re = np.take_along_axis(evals,   sort_re, axis=2).real
    weights_re = np.take_along_axis(weights, sort_re, axis=2)
    evals_im = np.take_along_axis(evals,   sort_im, axis=2).imag
    weights_im = np.take_along_axis(weights, sort_im, axis=2)

    D_mesh, W_mesh = np.meshgrid(D_vals, w_vals, indexing="xy")

    fig = plt.figure(figsize=(14, 6))
    ax_re = fig.add_subplot(1, 2, 1, projection="3d")
    ax_im = fig.add_subplot(1, 2, 2, projection="3d")

    cmap = plt.cm.turbo
    n_states = evals.shape[2]

    # Plot each of the n_states sheets.  alpha=0.5 keeps overlapping sheets
    # legible; shade=False preserves facecolors; rstride=ccount=1 uses every
    # grid cell for colour (no subsampling).
    for k in range(n_states):
        C_re = cmap(weights_re[:, :, k])
        C_im = cmap(weights_im[:, :, k])
        ax_re.plot_surface(
            D_mesh, W_mesh, evals_re[:, :, k],
            facecolors=C_re, shade=False, alpha=0.5,
            linewidth=0, antialiased=True,
            rcount=N_w, ccount=N_D,
        )
        ax_im.plot_surface(
            D_mesh, W_mesh, evals_im[:, :, k],
            facecolors=C_im, shade=False, alpha=0.5,
            linewidth=0, antialiased=True,
            rcount=N_w, ccount=N_D,
        )

    for ax, zlabel in [(ax_re, r"$\mathrm{Re}(\lambda)$"),
                       (ax_im, r"$\mathrm{Im}(\lambda)$")]:
        ax.set_xlabel(r"$D_r$",    labelpad=8)
        ax.set_ylabel(r"$\omega$", labelpad=8)
        ax.set_zlabel(zlabel,      labelpad=6)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_zlim(-1.2, 1.2)
        ax.set_xticks([0.0, 0.5, 1.0])
        ax.set_yticks([0.0, 0.5, 1.0])
        ax.set_zticks([-1, 0, 1])
        # Paper view: looking from ~+x +y, slightly down.  azim=-65, elev=18.
        ax.view_init(elev=18, azim=-65)

    # Panel titles like the paper ("Re(λ)" and "Im(λ)" floating at top)
    ax_re.set_title(r"$\mathrm{Re}(\lambda)$", fontsize=13, pad=2)
    ax_im.set_title(r"$\mathrm{Im}(\lambda)$", fontsize=13, pad=2)

    # Shared horizontal colorbar above the plots, centred
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.0, 1.0))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=[ax_re, ax_im], orientation="horizontal",
                        shrink=0.22, aspect=25, pad=-0.04,
                        location="top")
    cbar.set_label("edge", labelpad=2)
    cbar.set_ticks([0.0, 0.5, 1.0])

    plt.savefig(savepath, dpi=220, bbox_inches="tight")
    print(f"  saved {savepath}")
    return fig


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 62)
    print(" TCRW Fig 4(c) — OBC L=2, 3-D spectrum surfaces over (D_r, ω)")
    print("=" * 62)

    print("\n[1] Cross-check our OBC matrix vs authors' TRW.py")
    crosscheck_authors(0.3, 0.1, L=2)
    crosscheck_authors(0.5, 0.5, L=2)
    crosscheck_authors(1.0, 1e-3, L=2)

    print("\n[2] Building figure at 40×40 grid")
    make_figure(N_w=40, N_D=40, L=2, savepath="tcrw_fig4c_paper.png")

    print("\nDone.")
