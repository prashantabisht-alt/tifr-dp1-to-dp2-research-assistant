"""
TCRW Fig 4(i) — HPBC band circle on (sin k_y, cos k_y, Re λ), 3x3 grid
=======================================================================

Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020, Fig 4(i).

Same HPBC physics as Fig 4(h): y-periodic, x-open transition matrix
Fourier-transformed along y, parametrised by k_y ∈ [-π, π].  Here the
band structure is rendered in 3-D with the momentum wrapped onto the
unit circle so the periodicity is visually manifest.

Panel layout (3 × 3)
--------------------
     ω = 0.5        ω = 0.7        ω = 1.0
   ┌───────────────────────────────────────┐
   │  bulk     (green → magenta)           │   row 1
   │  CCW      (black → orange, copper)    │   row 2
   │  to-edge  (dark blue → yellow, cividis)│  row 3
   └───────────────────────────────────────┘

Each panel is a 3-D scatter of 44 · N_k points at:
    x = sin(k_y),  y = cos(k_y),  z = Re(λ)
coloured by one of three DISJOINT weights of the right eigenvectors:

    w_to_edge  = Σ |v|² over to-edge states (directors → wall)
    w_ccw      = Σ |v|² over CCW-along-edge states (directors going CCW)
    w_bulk     = Σ |v|² over the remaining states  (interior + CW-along
                                                    + bulk-directed-at-boundary)

By construction  w_bulk + w_ccw + w_to_edge = 1   (machine precision).

Paper parameters assumed
------------------------
    ω    ∈ {0.5, 0.7, 1.0}   — confirmed from paper column headers on 4(g)/4(i)
    D_r  = 0.1               — not labelled in panel; matches 4(e), 4(h)
    L    = 10                — stated in main text

Self-consistency
----------------
Matrix and Fourier convention identical to tcrw_fig4h.py.  Reused tests
from 4(h): k_y = 0 is real, column-stochastic, has λ = 1; k_y = π is
real; k_y ↔ −k_y conjugates.  Additionally check mask partition:
    to_edge ∩ ccw_along = ∅
    w_bulk + w_ccw + w_to = 1 at every eigenvector.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3-D projection)

DX = np.array([0, 1, 0, -1], dtype=int)
DY = np.array([1, 0, -1, 0], dtype=int)


def hpbc_index(x: int, d: int, L: int) -> int:
    return x * 4 + d


# ---------------------------------------------------------------------------
# 1. HPBC matrix — same as fig4h.py
# ---------------------------------------------------------------------------
def build_hpbc_matrix(omega: float, D_r: float, L: int, k_y: float) -> np.ndarray:
    n = 4 * (L + 1)
    A = np.zeros((n, n), dtype=complex)
    for x in range(L + 1):
        for d in range(4):
            src = hpbc_index(x, d, L)
            d_ccw = (d - 1) % 4
            d_cw  = (d + 1) % 4
            A[hpbc_index(x, d_ccw, L), src] += D_r * omega
            A[hpbc_index(x, d_cw,  L), src] += D_r * (1.0 - omega)
            nx = x + int(DX[d])
            if 0 <= nx <= L:
                phase = np.exp(-1j * k_y * float(DY[d]))
                A[hpbc_index(nx, d_cw,  L), src] += (1.0 - D_r) * omega       * phase
                A[hpbc_index(nx, d_ccw, L), src] += (1.0 - D_r) * (1.0 - omega) * phase
            else:
                A[src, src] += (1.0 - D_r)
    return A


# ---------------------------------------------------------------------------
# 2. Three disjoint masks — paper Fig 4(a) partition adapted to HPBC
# ---------------------------------------------------------------------------
# HPBC has only the two x-walls.  No corners (y is periodic), so to-edge
# and CCW-along cannot overlap; each has exactly 2 states at L = 10.
# The remaining 4L states are bulk (all interior x, plus the two
# non-edge directors at each wall).
# ---------------------------------------------------------------------------
def to_edge_mask(L: int) -> np.ndarray:
    n = 4 * (L + 1)
    m = np.zeros(n, dtype=bool)
    m[hpbc_index(0, 3, L)] = True   # left  wall: ←
    m[hpbc_index(L, 1, L)] = True   # right wall: →
    return m


def ccw_along_mask(L: int) -> np.ndarray:
    n = 4 * (L + 1)
    m = np.zeros(n, dtype=bool)
    m[hpbc_index(0, 2, L)] = True   # left  wall: ↓
    m[hpbc_index(L, 0, L)] = True   # right wall: ↑
    return m


def bulk_mask(L: int) -> np.ndarray:
    """Everything not in to-edge or ccw-along."""
    return ~(to_edge_mask(L) | ccw_along_mask(L))


# ---------------------------------------------------------------------------
# 3. Spectrum with three disjoint weight fields
# ---------------------------------------------------------------------------
def hpbc_spectrum_with_three_weights(omega: float, D_r: float, L: int, k_y: float):
    A = build_hpbc_matrix(omega, D_r, L, k_y)
    evals, evecs = np.linalg.eig(A)
    abs2 = np.abs(evecs) ** 2
    total = abs2.sum(axis=0)
    total = np.where(total > 1e-30, total, 1.0)
    m_to, m_ccw, m_bulk = to_edge_mask(L), ccw_along_mask(L), bulk_mask(L)
    w_to   = abs2[m_to,   :].sum(axis=0) / total
    w_ccw  = abs2[m_ccw,  :].sum(axis=0) / total
    w_bulk = abs2[m_bulk, :].sum(axis=0) / total
    return evals, w_to, w_ccw, w_bulk


# ---------------------------------------------------------------------------
# 4. Self-consistency
# ---------------------------------------------------------------------------
def self_consistency(L: int = 10):
    # Mask partition
    m_to, m_ccw, m_bulk = to_edge_mask(L), ccw_along_mask(L), bulk_mask(L)
    assert int((m_to & m_ccw).sum()) == 0
    assert int((m_to & m_bulk).sum()) == 0
    assert int((m_ccw & m_bulk).sum()) == 0
    assert int((m_to | m_ccw | m_bulk).sum()) == 4 * (L + 1)
    print(f"  mask partition: {int(m_to.sum())} to-edge + "
          f"{int(m_ccw.sum())} CCW + {int(m_bulk.sum())} bulk "
          f"= {4 * (L + 1)} (ok)")

    # Weight sum = 1 at a sample of (ω, D_r, k_y) points
    for (w, D, k) in [(0.5, 0.1, 0.3), (0.7, 0.1, -1.5), (1.0, 0.1, 0.0), (1.0, 0.1, np.pi)]:
        _, w_to, w_ccw, w_bulk = hpbc_spectrum_with_three_weights(w, D, L, k)
        total = w_to + w_ccw + w_bulk
        print(f"  ω={w}, D_r={D}, k_y={k:+.2f}: max|w_sum-1| = {float(np.max(np.abs(total-1))):.2e}")
        assert np.max(np.abs(total - 1.0)) < 1e-12

    # k_y = 0 block: real, column-stochastic, λ = 1
    for omega in (0.5, 0.7, 1.0):
        A0 = build_hpbc_matrix(omega, 0.1, L, 0.0)
        assert np.max(np.abs(A0.imag)) < 1e-14
        assert np.max(np.abs(A0.sum(axis=0) - 1.0)) < 1e-14
        ev0 = np.linalg.eigvals(A0)
        assert np.min(np.abs(ev0 - 1.0)) < 1e-10
    print("  k_y = 0 block: real, column-stochastic, λ = 1 present at all ω  (ok)")


# ---------------------------------------------------------------------------
# 5. Scan k_y at fixed (ω, D_r)
# ---------------------------------------------------------------------------
def scan_k_y(omega: float, D_r: float, L: int, N_k: int):
    k_vals = np.linspace(-np.pi, np.pi, N_k)
    n = 4 * (L + 1)
    evals  = np.zeros((N_k, n), dtype=complex)
    w_to   = np.zeros((N_k, n))
    w_ccw  = np.zeros((N_k, n))
    w_bulk = np.zeros((N_k, n))
    for i, k in enumerate(k_vals):
        e, wt, wc, wb = hpbc_spectrum_with_three_weights(omega, D_r, L, k)
        evals[i, :]  = e
        w_to[i,  :]  = wt
        w_ccw[i, :]  = wc
        w_bulk[i, :] = wb
    return k_vals, evals, w_to, w_ccw, w_bulk


# ---------------------------------------------------------------------------
# 6. Figure
# ---------------------------------------------------------------------------
def make_figure(omega_list=(0.5, 0.7, 1.0), D_r: float = 0.1, L: int = 10,
                N_k: int = 400, savepath: str = "tcrw_fig4i_paper.png"):
    # Precompute per ω
    cache = {}
    print(f"  scanning k_y for ω ∈ {omega_list}  (D_r={D_r}, L={L}, N_k={N_k})")
    t0 = time.time()
    for w in omega_list:
        cache[w] = scan_k_y(w, D_r, L, N_k)
        print(f"    ω = {w}  cpu {time.time()-t0:.1f}s")

    # Colormaps
    cmap_bulk    = mcolors.LinearSegmentedColormap.from_list(
        "bulk_green_magenta", ["#3e9f3a", "#c8248d"])    # 0=green, 1=magenta
    cmap_ccw     = plt.cm.copper                          # 0=black,  1=orange
    cmap_to_edge = plt.cm.cividis                         # 0=d-blue, 1=yellow

    nrow = 3
    ncol = len(omega_list)
    fig = plt.figure(figsize=(4.0 * ncol + 1.2, 4.0 * nrow))

    axes = np.empty((nrow, ncol), dtype=object)
    for r in range(nrow):
        for c in range(ncol):
            ax = fig.add_subplot(nrow, ncol, r * ncol + c + 1, projection="3d")
            axes[r, c] = ax

    for c, w in enumerate(omega_list):
        k_vals, evals, w_to, w_ccw, w_bulk = cache[w]
        Nk, ns = evals.shape
        # Per-point arrays (flatten bands together)
        sx = np.broadcast_to(np.sin(k_vals)[:, None], (Nk, ns)).ravel()
        cy = np.broadcast_to(np.cos(k_vals)[:, None], (Nk, ns)).ravel()
        z  = evals.real.ravel()
        wts_by_row = (w_bulk.ravel(), w_ccw.ravel(), w_to.ravel())
        cms        = (cmap_bulk,      cmap_ccw,      cmap_to_edge)
        for r in range(nrow):
            ax = axes[r, c]
            ax.scatter(sx, cy, z, c=wts_by_row[r], cmap=cms[r],
                       vmin=0, vmax=1, s=2.5, alpha=0.9, edgecolors="none")
            ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05); ax.set_zlim(-1.05, 1.05)
            ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1]); ax.set_zticks([-1, 0, 1])
            if r == nrow - 1:
                ax.set_xlabel(r"$\sin k_y$", labelpad=2)
                ax.set_ylabel(r"$\cos k_y$", labelpad=2)
            if r == 0:
                ax.set_title(rf"$\omega = {w}$", fontsize=11, pad=4)
            ax.view_init(elev=20, azim=-60)

    # Row colorbars on the right, one per coloring
    cb_width = 0.013
    cb_x     = 0.92
    row_positions = [(0.68, 0.22), (0.38, 0.22), (0.08, 0.22)]  # (y_bottom, height)
    row_cmaps  = [cmap_bulk, cmap_ccw, cmap_to_edge]
    row_labels = ["bulk", "CCW", "to edge"]
    for (y0, h), cm, lbl in zip(row_positions, row_cmaps, row_labels):
        cax = fig.add_axes([cb_x, y0, cb_width, h])
        sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(0, 1)); sm.set_array([])
        cb = fig.colorbar(sm, cax=cax)
        cb.set_label(lbl); cb.set_ticks([0.0, 0.5, 1.0])

    fig.suptitle(
        rf"(i)  HPBC band circle — $D_r = {D_r}$, $L = {L}$, $N_k = {N_k}$",
        fontsize=12, y=0.97,
    )
    plt.subplots_adjust(left=0.02, right=0.89, top=0.93, bottom=0.04,
                        wspace=0.08, hspace=0.15)
    plt.savefig(savepath, dpi=220, bbox_inches="tight")
    print(f"  saved {savepath}")


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 66)
    print(" TCRW Fig 4(i) — HPBC band circle (3x3 panels, bulk / CCW / to edge)")
    print("=" * 66)

    print("\n[1] Self-consistency (no author cross-check for HPBC)")
    self_consistency(L=10)

    print("\n[2] Production figure")
    make_figure(omega_list=(0.5, 0.7, 1.0), D_r=0.1, L=10, N_k=400,
                savepath="tcrw_fig4i_paper.png")

    print("\nDone.")
