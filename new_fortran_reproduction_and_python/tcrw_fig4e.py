"""
TCRW Fig 4(e) — OBC spectrum Re(λ) vs ω at D_r = 0.1, L = 10
=============================================================

Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020, Fig 4(e).

Physics
-------
Fix D_r = 10⁻¹ = 0.1 and scan ω ∈ [0, 1].  Every eigenvalue's Re(λ)
is plotted as a function of ω with colour = edge weight of its right
eigenvector.  Gap closes at ω = 0.5 (the paper's topological transition).

Matrix is 4(L+1)² = 484×484 for L = 10.  Dense diagonalisation; trivial.

Symmetry
--------
Under ω → 1 − ω the chiral rule swaps CW ↔ CCW; this is equivalent to
a global mirror reflection of the director.  Consequence: the eigenvalue
*set* at ω and 1 − ω agree as *complex-conjugate pairs*.  Re(λ) is
therefore symmetric about ω = 0.5, and Im(λ) is antisymmetric.  The plot
must be mirror-symmetric under ω ↔ 1 − ω.  Spot-check asserts this
numerically.

Edge-mask convention
--------------------
Same paper definition (Fig 4(a), ω = 1 reference): edge = to-edge ∪
CCW-along-edge.  This partition is held fixed across the ω-scan — so
at ω < 0.5 the walker is biased toward CW-along-edge (NOT in the edge
set), which is why the edge weight on bulk-ish modes dips lower on the
ω < 0.5 side in the paper figure.

Plot
----
Scatter (ω, Re(λ)) with colour = edge weight under `turbo`.  N_w = 300
gives curve-like rendering.

Cross-check
-----------
`crosscheck_authors` at L = 2 against authors' sparse matrix — expect
element-wise 0.

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
# 1. OBC transition matrix (self-contained, any L)
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
# 2. Edge mask — paper Fig 4(a) definition, ω = 1 reference chirality
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
    return evals, abs2[edge, :].sum(axis=0) / total


# ---------------------------------------------------------------------------
# 4. Author cross-check
# ---------------------------------------------------------------------------
def _load_authors_trw():
    here = os.path.dirname(os.path.abspath(__file__))
    for path in [
        os.path.join(here, "TRW._original_code_by_paperauthors.py"),
    ]:
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
# 5. Scan over ω at fixed D_r
# ---------------------------------------------------------------------------
def scan_omega(D_r: float = 0.1, L: int = 10, N_w: int = 300,
               w_min: float = 0.0, w_max: float = 1.0):
    w_vals = np.linspace(w_min, w_max, N_w)
    n_states = 4 * (L + 1) ** 2
    evals = np.zeros((N_w, n_states), dtype=complex)
    wts   = np.zeros((N_w, n_states))
    edge = edge_mask(L)
    print(f"  scanning {N_w} ω points (matrix {n_states}×{n_states}, D_r={D_r}, L={L})")
    t0 = time.time()
    for i, w in enumerate(w_vals):
        P = build_obc_matrix(w, D_r, L).toarray()
        e, v = np.linalg.eig(P)
        a2 = np.abs(v) ** 2
        tot = a2.sum(axis=0); tot = np.where(tot > 1e-30, tot, 1.0)
        evals[i, :] = e
        wts[i, :]   = a2[edge, :].sum(axis=0) / tot
        if (i + 1) % max(1, N_w // 10) == 0:
            print(f"    {i + 1}/{N_w}   cpu {time.time() - t0:.1f}s")
    return w_vals, evals, wts


# ---------------------------------------------------------------------------
# 6. Physics spot-checks
# ---------------------------------------------------------------------------
def spot_checks(L: int = 10, D_r: float = 0.1):
    """
    1. λ = 1 present at every ω (column-stochasticity).
    2. ω ↔ 1 − ω: eigenvalue sets are complex-conjugate pairs.
       Re(λ) distribution at ω should match Re(λ) distribution at 1 − ω.
    """
    print("  [1] λ = 1 at every ω")
    for w in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0):
        e, _ = obc_spectrum_with_weights(w, D_r, L)
        dist = np.min(np.abs(e - 1.0))
        assert dist < 1e-10, f"λ=1 missing at ω={w}: {dist:.2e}"
    print("     ok")

    print("  [2] ω ↔ 1 − ω symmetry — Re(λ) sets match")
    for w in (0.2, 0.3, 0.4):
        e1, _ = obc_spectrum_with_weights(w,         D_r, L)
        e2, _ = obc_spectrum_with_weights(1.0 - w,   D_r, L)
        # The spectrum at 1 − ω should equal the complex conjugate of ω.
        # So sorted Re sets should match; sorted Im sets should be negatives.
        re1_sorted = np.sort(e1.real); re2_sorted = np.sort(e2.real)
        d_re = float(np.max(np.abs(re1_sorted - re2_sorted)))
        im1_sorted = np.sort(e1.imag); im2_sorted_neg = np.sort(-e2.imag)
        d_im = float(np.max(np.abs(im1_sorted - im2_sorted_neg)))
        print(f"     ω = {w} vs 1 − ω = {1-w}:  "
              f"max|ΔRe| = {d_re:.2e}   max|ΔIm+Im'| = {d_im:.2e}")
        assert d_re < 1e-10, f"Re symmetry fails at ω={w}"
        assert d_im < 1e-10, f"Im anti-symmetry fails at ω={w}"
    print("     ok")


# ---------------------------------------------------------------------------
# 7. Figure
# ---------------------------------------------------------------------------
def make_figure(D_r: float = 0.1, L: int = 10, N_w: int = 300,
                savepath: str = "tcrw_fig4e_paper.png"):
    w_vals, evals, wts = scan_omega(D_r=D_r, L=L, N_w=N_w)

    N_w_, n_states = evals.shape
    w_rep = np.broadcast_to(w_vals[:, None], (N_w_, n_states)).ravel()
    y_re  = evals.real.ravel()
    c_all = wts.ravel()

    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    cmap = plt.cm.turbo
    sc = ax.scatter(w_rep, y_re, c=c_all, cmap=cmap, vmin=0.0, vmax=1.0,
                    s=3, alpha=0.75, edgecolors="none")

    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"$\mathrm{Re}(\lambda)$")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-1.05, 1.05)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_yticks([-1, 0, 1])
    ax.axhline(0, color="0.6", lw=0.5, ls="--")
    ax.set_title(rf"(e)  OBC Re($\lambda$) vs $\omega$   at  $D_r={D_r}$,  $L={L}$",
                 fontsize=11)

    cbar = fig.colorbar(sc, ax=ax, shrink=0.85, aspect=25)
    cbar.set_label("edge"); cbar.set_ticks([0.0, 0.5, 1.0])

    plt.tight_layout()
    plt.savefig(savepath, dpi=220, bbox_inches="tight")
    print(f"  saved {savepath}")


# ---------------------------------------------------------------------------
# 8. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 64)
    print(" TCRW Fig 4(e) — OBC Re(λ) vs ω at D_r = 0.1, L = 10")
    print("=" * 64)

    print("\n[1] Author cross-check at L = 2")
    crosscheck_authors(0.2, 0.1, L=2)
    crosscheck_authors(0.5, 0.1, L=2)
    crosscheck_authors(0.9, 0.1, L=2)

    print("\n[2] Physics spot-checks at L = 10")
    spot_checks(L=10, D_r=0.1)

    print("\n[3] Production figure (N_w = 300)")
    make_figure(D_r=0.1, L=10, N_w=300, savepath="tcrw_fig4e_paper.png")

    print("\nDone.")
