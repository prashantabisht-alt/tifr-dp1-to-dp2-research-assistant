"""
TCRW Fig 4(f) — OBC complex-plane spectrum at ω = 1, L = 10
============================================================

Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020, Fig 4(f).

Physics
-------
For a fully chiral OBC walker (ω = 1) and L = 10 (4·11² = 484 states),
plot every eigenvalue in the complex plane (Re λ, Im λ) at three values
of D_r ∈ {0.65, 0.5, 0.35} (paper columns).  Two coloring schemes (paper
rows), each revealing a different facet of Fig 4(a)'s partition of
boundary states:

    "in edge" (top row)    — fraction of |v|^2 on TO-EDGE states only
                             (director points into the adjacent wall)
    "CCW"     (bottom row) — fraction of |v|^2 on CCW-ALONG-EDGE states
                             only (director goes in CCW direction, i.e.
                             CCW edge-current carriers)

These two sub-categories are disjoint; a state at a boundary site is
classified as exactly one of {to-edge, CCW-along, CW-along or
bulk-directed}.  Interior sites contribute zero in either coloring.

Physics expectation (sanity before plot)
----------------------------------------
Unit-circle eigenvalues (|λ| ≈ 1) are the long-lived topological edge
currents.  They live on CCW-along-edge states at ω = 1:
    "in edge" weight → 0   (dark blue  on top row)
    "CCW"      weight → 1   (bright on bottom row)

Short-lived bulk-ish modes near the origin:
    "in edge" → non-zero   (lighter on top row)
    "CCW"      → small      (dark on bottom row)

Convention
----------
Direction / indexing identical to tcrw_fig4c/d/e.py and authors' TRW.py.
State index  s = (i · (L+1) + j) · 4 + d.

Cross-check
-----------
`crosscheck_authors` at L = 2: element-wise zero vs authors' sparse.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import time
import importlib.util
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Direction vectors — authors' convention
# ---------------------------------------------------------------------------
DX = np.array([0, 1, 0, -1], dtype=int)
DY = np.array([1, 0, -1, 0], dtype=int)


def state_index(i: int, j: int, d: int, L: int) -> int:
    return (i * (L + 1) + j) * 4 + d


# ---------------------------------------------------------------------------
# 1. OBC transition matrix
# ---------------------------------------------------------------------------
def build_obc_matrix(omega: float, D_r: float, L: int) -> sp.csc_matrix:
    n = 4 * (L + 1) ** 2
    rows, cols, vals = [], [], []
    for i in range(L + 1):
        for j in range(L + 1):
            for d in range(4):
                src = state_index(i, j, d, L)
                d_ccw = (d - 1) % 4
                d_cw  = (d + 1) % 4
                rows.append(state_index(i, j, d_ccw, L)); cols.append(src)
                vals.append(D_r * omega)
                rows.append(state_index(i, j, d_cw,  L)); cols.append(src)
                vals.append(D_r * (1.0 - omega))
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
# 2. Two disjoint edge-type masks — paper Fig 4(a)
# ---------------------------------------------------------------------------
# At each boundary site:
#   left   (x=0): to-edge = ← (d=3),  CCW-along = ↓ (d=2)
#   right  (x=L): to-edge = → (d=1),  CCW-along = ↑ (d=0)
#   bottom (y=0): to-edge = ↓ (d=2),  CCW-along = → (d=1)
#   top    (y=L): to-edge = ↑ (d=0),  CCW-along = ← (d=3)
# At corners two edges apply; the mask is the union.
_EDGE_DIRS_CCW = {
    "left":   (3, 2),
    "right":  (1, 0),
    "bottom": (2, 1),
    "top":    (0, 3),
}


def to_edge_mask(L: int) -> np.ndarray:
    """States where director points INTO a wall."""
    n = 4 * (L + 1) ** 2
    m = np.zeros(n, dtype=bool)
    for i in range(L + 1):
        for j in range(L + 1):
            applies = []
            if i == 0: applies.append("left")
            if i == L: applies.append("right")
            if j == 0: applies.append("bottom")
            if j == L: applies.append("top")
            for e in applies:
                d_to, _ = _EDGE_DIRS_CCW[e]
                m[state_index(i, j, d_to, L)] = True
    return m


def ccw_along_mask(L: int) -> np.ndarray:
    """
    States where director points in the CCW circulation direction on that
    edge, EXCLUDING states that are already classified as to-edge from
    another wall (this happens at corners, where one director is
    simultaneously "to-edge from bottom" and "CCW-along from left", for
    example).  Paper Fig 4(a) resolves the overlap by labelling those
    states as to-edge (blue) rather than CCW-along (red).
    """
    n = 4 * (L + 1) ** 2
    m = np.zeros(n, dtype=bool)
    te = to_edge_mask(L)
    for i in range(L + 1):
        for j in range(L + 1):
            applies = []
            if i == 0: applies.append("left")
            if i == L: applies.append("right")
            if j == 0: applies.append("bottom")
            if j == L: applies.append("top")
            for e in applies:
                _, d_along = _EDGE_DIRS_CCW[e]
                s = state_index(i, j, d_along, L)
                if not te[s]:
                    m[s] = True
    return m


# ---------------------------------------------------------------------------
# 3. Spectrum + two weight fields
# ---------------------------------------------------------------------------
def obc_spectrum_with_weights(omega: float, D_r: float, L: int):
    P = build_obc_matrix(omega, D_r, L).toarray()
    evals, evecs = np.linalg.eig(P)
    abs2 = np.abs(evecs) ** 2
    total = abs2.sum(axis=0)
    total = np.where(total > 1e-30, total, 1.0)
    m_to  = to_edge_mask(L)
    m_ccw = ccw_along_mask(L)
    w_to  = abs2[m_to,  :].sum(axis=0) / total
    w_ccw = abs2[m_ccw, :].sum(axis=0) / total
    return evals, w_to, w_ccw


# ---------------------------------------------------------------------------
# 4. Author cross-check
# ---------------------------------------------------------------------------
def _load_authors_trw():
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "TRW._original_code_by_paperauthors.py")
    if os.path.exists(path):
        spec = importlib.util.spec_from_file_location("TRW_authors", path)
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
        return mod
    return None


def crosscheck_authors(omega: float, D_r: float, L: int, tol: float = 1e-12):
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
    print(f"  crosscheck ω={omega:.3g} D_r={D_r:.3g} L={L}")
    print(f"    element diff = {max_elem:.2e},  Hausdorff = {hausdorff:.2e}")
    assert max_elem  < tol
    assert hausdorff < tol
    print(f"    [ok] matches authors (elem diff = {max_elem:.2e})")


# ---------------------------------------------------------------------------
# 5. Physics spot-checks
# ---------------------------------------------------------------------------
def spot_checks(L: int = 10):
    """
    1. λ = 1 present at each D_r in paper list.
    2. |λ| ≤ 1 for every eigenvalue (stochasticity).
    3. Masks are disjoint: to_edge ∩ ccw_along = ∅, so w_to + w_ccw ≤ 1.
    4. At D_r = 1 the spectrum is {1, i, -1, -i} each ×(L+1)^2; hand-checkable.
    """
    m_to  = to_edge_mask(L)
    m_ccw = ccw_along_mask(L)
    overlap = int(np.sum(m_to & m_ccw))
    print(f"  to_edge ∩ ccw_along = {overlap} states  (must be 0)")
    assert overlap == 0, "masks overlap!"
    n_to  = int(np.sum(m_to));  n_ccw = int(np.sum(m_ccw))
    # Expected counts at L:
    #   to-edge: each corner contributes 2, each mid-edge 1.  Corners = 4,
    #            mid-edges = 4(L-1), so |to-edge| = 8 + 4(L-1) = 4L + 4.
    #   ccw-along (disjoint): each corner contributes 1 (the non-overlap
    #            direction), each mid-edge 1.  |CCW-along| = 4 + 4(L-1) = 4L.
    print(f"  |to_edge| = {n_to}  (expect {4*L + 4} for L={L})")
    print(f"  |ccw_along| = {n_ccw}  (expect {4*L} for L={L})")
    assert n_to  == 4*L + 4
    assert n_ccw == 4*L

    for D_r in (0.65, 0.5, 0.35):
        e, w_to, w_ccw = obc_spectrum_with_weights(1.0, D_r, L)
        assert np.min(np.abs(e - 1.0)) < 1e-10
        assert np.max(np.abs(e)) <= 1.0 + 1e-12
        assert np.all(w_to + w_ccw <= 1.0 + 1e-12)
    print(f"  λ = 1 present; |λ| ≤ 1; w_to + w_ccw ≤ 1 at all three D_r")


# ---------------------------------------------------------------------------
# 6. Figure
# ---------------------------------------------------------------------------
def make_figure(L: int = 10, D_r_list=(0.65, 0.5, 0.35),
                savepath: str = "tcrw_fig4f_paper.png"):
    omega = 1.0
    # Compute once per D_r
    data = []
    print(f"  computing 2x{len(D_r_list)} panels  (ω={omega}, L={L})")
    t0 = time.time()
    for D in D_r_list:
        e, w_to, w_ccw = obc_spectrum_with_weights(omega, D, L)
        data.append((D, e, w_to, w_ccw))
        print(f"    D_r = {D}   cpu {time.time()-t0:.1f}s")

    ncols = len(D_r_list)
    fig, axs = plt.subplots(2, ncols, figsize=(4.2 * ncols, 7),
                             sharex=True, sharey=True)

    cmap_in  = plt.cm.cividis   # blue → yellow ("in edge")
    cmap_ccw = plt.cm.copper    # black → orange ("CCW")

    for col, (D, e, w_to, w_ccw) in enumerate(data):
        ax_t = axs[0, col]; ax_b = axs[1, col]
        sc_t = ax_t.scatter(e.real, e.imag, c=w_to,  cmap=cmap_in,
                            vmin=0, vmax=1, s=18, edgecolors="none", alpha=0.9)
        sc_b = ax_b.scatter(e.real, e.imag, c=w_ccw, cmap=cmap_ccw,
                            vmin=0, vmax=1, s=18, edgecolors="none", alpha=0.9)
        # Column title: D_r value
        ax_t.set_title(rf"$D_r = {D}$", fontsize=11)

    for ax in axs.ravel():
        ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05)
        ax.set_aspect("equal")
        ax.axhline(0, color="0.7", lw=0.4, ls="--")
        ax.axvline(0, color="0.7", lw=0.4, ls="--")

    # Only outer axes show labels
    for col in range(ncols):
        axs[1, col].set_xlabel(r"$\mathrm{Re}(\lambda)$")
    axs[0, 0].set_ylabel(r"$\mathrm{Im}(\lambda)$")
    axs[1, 0].set_ylabel(r"$\mathrm{Im}(\lambda)$")

    # Two colorbars: one per row, on the right
    cax_t = fig.add_axes([0.92, 0.55, 0.013, 0.33])
    cax_b = fig.add_axes([0.92, 0.14, 0.013, 0.33])
    sm_t = plt.cm.ScalarMappable(cmap=cmap_in,  norm=plt.Normalize(0, 1)); sm_t.set_array([])
    sm_b = plt.cm.ScalarMappable(cmap=cmap_ccw, norm=plt.Normalize(0, 1)); sm_b.set_array([])
    cb_t = fig.colorbar(sm_t, cax=cax_t); cb_t.set_label("in edge")
    cb_b = fig.colorbar(sm_b, cax=cax_b); cb_b.set_label("CCW")
    cb_t.set_ticks([0, 1]); cb_b.set_ticks([0, 1])

    fig.suptitle(r"(f)  OBC spectrum at $\omega = 1$, $L = 10$", fontsize=12, y=0.95)
    plt.subplots_adjust(left=0.08, right=0.90, top=0.88, bottom=0.10, wspace=0.05, hspace=0.18)
    plt.savefig(savepath, dpi=220, bbox_inches="tight")
    print(f"  saved {savepath}")


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 62)
    print(" TCRW Fig 4(f) — OBC complex-plane at ω = 1, L = 10")
    print(" paper D_r = 0.65, 0.5, 0.35 ; colorings: in-edge + CCW")
    print("=" * 62)

    print("\n[1] Author cross-check at L = 2")
    crosscheck_authors(1.0, 0.65, L=2)
    crosscheck_authors(1.0, 0.50, L=2)
    crosscheck_authors(1.0, 0.35, L=2)

    print("\n[2] Physics spot-checks at L = 10")
    spot_checks(L=10)

    print("\n[3] Production figure")
    make_figure(L=10, D_r_list=(0.65, 0.5, 0.35), savepath="tcrw_fig4f_paper.png")

    print("\nDone.")
