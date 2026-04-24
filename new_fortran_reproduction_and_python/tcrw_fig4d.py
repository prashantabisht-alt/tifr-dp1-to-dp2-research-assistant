"""
TCRW Fig 4(d) — OBC spectrum Re(λ) vs D_r at ω = 1, L = 10
===========================================================

Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020, Fig 4(d).

Physics
-------
For a fully chiral OBC walker (ω = 1) we scan the noise probability D_r
and plot every eigenvalue's Re(λ) as a function of D_r, colouring by the
edge weight of its right eigenvector.

At L = 10 (authors' convention: grid 0..L inclusive → 11×11 playground),
the transition matrix is 4·(L+1)² = 484×484.  Dense diagonalisation at
each D_r; ~50 ms per point, trivial.

Limits (useful as sanity checks)
--------------------------------
D_r → 0 (deterministic chiral): only chiral steps, no noise.  Steady state
     traps the walker on edge-CCW cycles; eigenvalues collapse to
     {+1, +i, -1, -i} (the fourth roots of unity from period-4 cycles),
     each with high multiplicity.  Re(λ) → {-1, 0, +1}.

D_r = 1 (pure noise): matrix is block-diagonal per site, each block = 4×4
     all-CCW permutation at ω = 1.  Eigenvalues = {1, i, -1, -i} with
     multiplicity (L+1)² = 121 each; Re(λ) = {-1, 0 (doubly), +1}.

So the spectrum is pinned at Re = {-1, 0, +1} at *both* endpoints and fans
out in between — this is the "coalescence" mentioned in the caption.

Convention
----------
Direction / indexing identical to tcrw_fig4c.py and authors' TRW.py:
    d = 0 ↑,  1 →,  2 ↓,  3 ←
    s = (i * (L+1) + j) * 4 + d   (authors' indexing)

Edge-mask partition (paper Fig 4(a), ω = 1 reference chirality):
    edge = to-edge ∪ CCW-along-edge   (everything else = bulk)

Plot
----
2-D scatter of (D_r, Re(λ)) with marker colour = edge weight under the
`turbo` colormap (matches the paper's blue → green → red scale).  At
N_D = 200 × 484 points the density reads as continuous coloured curves.

Cross-check
-----------
`crosscheck_authors` at L = 2 — the small matrix is fast and the
element-wise comparison is the strongest possible check.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import sys
import time
import importlib.util
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Direction vectors — authors' convention
# ---------------------------------------------------------------------------
DX = np.array([0, 1, 0, -1], dtype=int)   # dx for d = 0↑, 1→, 2↓, 3←
DY = np.array([1, 0, -1, 0], dtype=int)


def state_index(i: int, j: int, d: int, L: int) -> int:
    return (i * (L + 1) + j) * 4 + d


# ---------------------------------------------------------------------------
# 1. OBC transition matrix (self-contained, any L)
# ---------------------------------------------------------------------------
def build_obc_matrix(omega: float, D_r: float, L: int) -> sp.csc_matrix:
    """4(L+1)² × 4(L+1)² column-stochastic OBC matrix, matches authors'."""
    n = 4 * (L + 1) ** 2
    rows, cols, vals = [], [], []
    for i in range(L + 1):
        for j in range(L + 1):
            for d in range(4):
                src = state_index(i, j, d, L)
                d_ccw = (d - 1) % 4
                d_cw  = (d + 1) % 4
                # Noise step
                rows.append(state_index(i, j, d_ccw, L)); cols.append(src)
                vals.append(D_r * omega)
                rows.append(state_index(i, j, d_cw,  L)); cols.append(src)
                vals.append(D_r * (1.0 - omega))
                # Chiral step
                ni, nj = i + int(DX[d]), j + int(DY[d])
                if 0 <= ni <= L and 0 <= nj <= L:
                    rows.append(state_index(ni, nj, d_cw,  L)); cols.append(src)
                    vals.append((1.0 - D_r) * omega)
                    rows.append(state_index(ni, nj, d_ccw, L)); cols.append(src)
                    vals.append((1.0 - D_r) * (1.0 - omega))
                else:
                    rows.append(src); cols.append(src)
                    vals.append(1.0 - D_r)
    return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()


# ---------------------------------------------------------------------------
# 2. Edge mask — paper Fig 4(a) definition
# ---------------------------------------------------------------------------
_EDGE_DIRS_CCW = {
    "left":   (3, 2),   # to-edge=←, CCW-along=↓
    "right":  (1, 0),   # to-edge=→, CCW-along=↑
    "bottom": (2, 1),   # to-edge=↓, CCW-along=→
    "top":    (0, 3),   # to-edge=↑, CCW-along=←
}


def edge_mask(L: int) -> np.ndarray:
    n = 4 * (L + 1) ** 2
    mask = np.zeros(n, dtype=bool)
    for i in range(L + 1):
        for j in range(L + 1):
            applies = []
            if i == 0: applies.append("left")
            if i == L: applies.append("right")
            if j == 0: applies.append("bottom")
            if j == L: applies.append("top")
            for e in applies:
                d_to, d_along = _EDGE_DIRS_CCW[e]
                mask[state_index(i, j, d_to,    L)] = True
                mask[state_index(i, j, d_along, L)] = True
    return mask


# ---------------------------------------------------------------------------
# 3. OBC spectrum with edge weights
# ---------------------------------------------------------------------------
def obc_spectrum_with_weights(omega: float, D_r: float, L: int):
    P = build_obc_matrix(omega, D_r, L).toarray()
    evals, evecs = np.linalg.eig(P)
    abs2 = np.abs(evecs) ** 2
    total = abs2.sum(axis=0)
    total = np.where(total > 1e-30, total, 1.0)
    edge = edge_mask(L)
    w_edge = abs2[edge, :].sum(axis=0) / total
    return evals, w_edge


# ---------------------------------------------------------------------------
# 4. Author cross-check (at L = 2, any ω, any D_r)
# ---------------------------------------------------------------------------
def _load_authors_trw():
    here = os.path.dirname(os.path.abspath(__file__))
    for path in [
        os.path.join(here, "TRW._original_code_by_paperauthors.py"),
        os.path.join(here, "new_fortran_reproduction_and_python",
                     "TRW._original_code_by_paperauthors.py"),
    ]:
        if os.path.exists(path):
            spec = importlib.util.spec_from_file_location("TRW_authors", path)
            mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
            return mod
    return None


def crosscheck_authors(omega: float, D_r: float, L: int,
                       tol: float = 1e-12, verbose: bool = True):
    TRW = _load_authors_trw()
    if TRW is None:
        print("  [warn] authors' TRW module not found; skipping cross-check")
        return None
    W_A = TRW.build_sparse_transition_matrix(L, omega, D_r).toarray()
    W_U = build_obc_matrix(omega, D_r, L).toarray()
    max_elem = float(np.max(np.abs(W_A - W_U)))
    ev_A = np.linalg.eigvals(W_A); ev_U = np.linalg.eigvals(W_U)
    d_UA = float(np.min(np.abs(ev_U[:, None] - ev_A[None, :]), axis=1).max())
    d_AU = float(np.min(np.abs(ev_A[:, None] - ev_U[None, :]), axis=1).max())
    hausdorff = max(d_UA, d_AU)
    if verbose:
        print(f"  crosscheck  ω={omega:.3g}  D_r={D_r:.3g}  L={L}")
        print(f"    element-wise max |W_A - W_U| = {max_elem:.2e}")
        print(f"    eigenvalue Hausdorff         = {hausdorff:.2e}")
    assert max_elem  < tol, f"element mismatch {max_elem:.2e} > tol"
    assert hausdorff < tol, f"eigval mismatch {hausdorff:.2e} > tol"
    print(f"    [ok] matches authors (elem diff = {max_elem:.2e})")
    return max_elem, hausdorff


# ---------------------------------------------------------------------------
# 5. Scan over D_r at fixed ω
# ---------------------------------------------------------------------------
def scan_D_r(omega: float = 1.0, L: int = 10, N_D: int = 200,
             D_min: float = 0.0, D_max: float = 1.0):
    """
    Diagonalize at each D_r; collect (D_r, eigvals, weights).

    Returns
    -------
    D_vals : (N_D,)
    evals  : (N_D, n_states)  complex
    wts    : (N_D, n_states)  real in [0, 1]
    """
    D_vals = np.linspace(D_min, D_max, N_D)
    n_states = 4 * (L + 1) ** 2
    evals = np.zeros((N_D, n_states), dtype=complex)
    wts   = np.zeros((N_D, n_states))
    edge = edge_mask(L)

    print(f"  scanning {N_D} D_r points   (matrix {n_states}×{n_states}, ω={omega}, L={L})")
    t0 = time.time()
    for i, D in enumerate(D_vals):
        P = build_obc_matrix(omega, D, L).toarray()
        e, v = np.linalg.eig(P)
        a2 = np.abs(v) ** 2
        tot = a2.sum(axis=0)
        tot = np.where(tot > 1e-30, tot, 1.0)
        evals[i, :] = e
        wts[i, :]   = a2[edge, :].sum(axis=0) / tot
        if (i + 1) % max(1, N_D // 10) == 0:
            print(f"    {i + 1}/{N_D}   cpu {time.time() - t0:.1f}s")
    return D_vals, evals, wts


# ---------------------------------------------------------------------------
# 6. Figure
# ---------------------------------------------------------------------------
def make_figure(omega: float = 1.0, L: int = 10, N_D: int = 200,
                savepath: str = "tcrw_fig4d_paper.png"):
    D_vals, evals, wts = scan_D_r(omega=omega, L=L, N_D=N_D)

    # Flatten for scatter:  x = D_r (repeated), y = Re(λ), c = edge weight
    N_D_, n_states = evals.shape
    D_rep = np.broadcast_to(D_vals[:, None], (N_D_, n_states)).ravel()
    y_re  = evals.real.ravel()
    c_all = wts.ravel()

    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    cmap = plt.cm.turbo
    sc = ax.scatter(D_rep, y_re, c=c_all, cmap=cmap, vmin=0.0, vmax=1.0,
                    s=3, alpha=0.75, edgecolors="none")

    ax.set_xlabel(r"$D_r$")
    ax.set_ylabel(r"$\mathrm{Re}(\lambda)$")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-1.05, 1.05)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_yticks([-1, 0, 1])
    ax.axhline(0, color="0.6", lw=0.5, ls="--")
    ax.set_title(rf"(d)  OBC Re($\lambda$) vs $D_r$   at  $\omega={omega}$,  $L={L}$",
                 fontsize=11)

    cbar = fig.colorbar(sc, ax=ax, shrink=0.85, aspect=25)
    cbar.set_label("edge")
    cbar.set_ticks([0.0, 0.5, 1.0])

    plt.tight_layout()
    plt.savefig(savepath, dpi=220, bbox_inches="tight")
    print(f"  saved {savepath}")
    return fig, D_vals, evals, wts


# ---------------------------------------------------------------------------
# 7. Physics spot-checks (hand-calculated limits)
# ---------------------------------------------------------------------------
def spot_checks(L: int = 10, verbose: bool = True):
    """
    At D_r = 1 (pure noise, ω = 1) every spatial site decouples and its
    4×4 block is the CCW permutation with eigenvalues {1, i, -1, -i}.
    So the full 484 eigenvalues are {±1, ±i} each with multiplicity 121.
    Re(λ) should be exactly {-1, 0, +1} with counts (121, 242, 121).
    """
    e, _ = obc_spectrum_with_weights(1.0, 1.0, L)
    re = np.sort(e.real)
    n_m1 = int(np.sum(np.abs(re + 1) < 1e-10))
    n_0  = int(np.sum(np.abs(re    ) < 1e-10))
    n_p1 = int(np.sum(np.abs(re - 1) < 1e-10))
    ns   = (L + 1) ** 2
    if verbose:
        print(f"  D_r = 1 spot-check  (expect (121, 242, 121) for L=10)")
        print(f"    #{{Re(λ) = -1}} = {n_m1}   expected {ns}")
        print(f"    #{{Re(λ) =  0}} = {n_0}    expected {2 * ns}")
        print(f"    #{{Re(λ) = +1}} = {n_p1}   expected {ns}")
    assert n_m1 == ns and n_0 == 2 * ns and n_p1 == ns

    # λ = 1 must be present at every D_r (column-stochasticity)
    for D in (0.0, 0.1, 0.5, 0.9, 1.0):
        e, _ = obc_spectrum_with_weights(1.0, D, L)
        dist = np.min(np.abs(e - 1.0))
        assert dist < 1e-10, f"λ=1 missing at D_r={D}: min |λ-1| = {dist:.2e}"
    print(f"  λ = 1 present at every checked D_r")


# ---------------------------------------------------------------------------
# 8. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 64)
    print(" TCRW Fig 4(d) — OBC Re(λ) vs D_r at ω = 1, L = 10")
    print("=" * 64)

    print("\n[1] Author cross-check at L = 2 (fast)")
    crosscheck_authors(1.0, 0.3, L=2)
    crosscheck_authors(1.0, 0.05, L=2)
    crosscheck_authors(1.0, 0.8, L=2)

    print("\n[2] Physics spot-checks at L = 10")
    spot_checks(L=10)

    print("\n[3] Production figure (N_D = 200)")
    make_figure(omega=1.0, L=10, N_D=200, savepath="tcrw_fig4d_paper.png")

    print("\nDone.")
