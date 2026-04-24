"""
TCRW Fig 4(h) — HPBC spectrum vs k_y  at  ω = 1, D_r = 0.1, L = 10
===================================================================

Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020, Fig 4(h).
Caption: "Spectrum of a model with hybrid boundary conditions (periodic
along the y-direction and open along the x-direction).  Note the bands
localized on the edge.  There are two bands, each localized at one of
the edges and with opposite winding."

Paper parameters assumed (not labelled on panel):
    ω    = 1      — fully chiral; matches 4(d) / 4(f) style
    D_r  = 0.1    — matches 4(e)'s explicit "D_r = 10⁻¹"
    L    = 10     — stated in main text "system of size L = 10"

Physics
-------
With y periodic, Fourier-transform the 2-D transition matrix along y:
    v(x, d, k_y) = Σ_y exp(-i k_y y) ψ(x, y, d)
Each k_y block is a 4(L+1) × 4(L+1) complex matrix A(k_y) with

    noise steps            : Δy = 0                  → no k_y phase
    chiral step (valid)    : Δy = DY[d]              → factor exp(-i k_y DY[d])
    chiral blocked (x OBC) : no displacement         → no phase (self-loop)

At k_y = 0 every phase = 1, and A reduces to the y-summed real-space
matrix (still column-stochastic and real, with λ = 1).  At k_y = π every
phase is ±1, so A is real but not column-stochastic.

Self-consistency tests (the only sanity we have — authors ship no HPBC)
----------------------------------------------------------------------
1. k_y = 0:  max |Im A| = 0   and   max |col_sum - 1| = 0   and   λ = 1 in spectrum.
2. k_y = π:  max |Im A| = 0.
3. k_y ↔ -k_y:  spectrum conjugates (Re(λ) even in k_y, Im(λ) odd).
4. Matrix size sanity: 4(L+1) = 44 at L = 10.
5. At each k_y, max |λ| ≤ 1 (stochastic spectrum).

Edge-mask convention (single colorbar "edge" like 4(c/d/e))
-----------------------------------------------------------
HPBC has only the two x-walls (left at x = 0, right at x = L):
    left  (x = 0): to-edge = ← (d = 3),  CCW-along = ↓ (d = 2)
    right (x = L): to-edge = → (d = 1),  CCW-along = ↑ (d = 0)
The edge-state set is their union (4 edge states out of 4(L+1)=44 at L=10).

Plot
----
2-row × 1-col figure:  top = Re(λ) vs k_y,  bottom = Im(λ) vs k_y.
Every one of the 44 bands plotted; colour = edge weight under `turbo`.
The paper text says only two of them are edge-localised; those two
bands should appear red-colored while the other 42 are blue-ish.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import time
import importlib.util
import numpy as np
import matplotlib.pyplot as plt

DX = np.array([0, 1, 0, -1], dtype=int)
DY = np.array([1, 0, -1, 0], dtype=int)


def hpbc_index(x: int, d: int, L: int) -> int:
    """State index at fixed k_y: s = x * 4 + d,  x ∈ 0..L,  d ∈ 0..3."""
    return x * 4 + d


# ---------------------------------------------------------------------------
# 1. HPBC transition matrix at fixed k_y
# ---------------------------------------------------------------------------
def build_hpbc_matrix(omega: float, D_r: float, L: int, k_y: float) -> np.ndarray:
    """
    Dense 4(L+1) × 4(L+1) complex matrix A(k_y) for the y-Fourier block.
    Convention:  v(x, d, k_y) = Σ_y e^{-i k_y y} ψ(x, y, d)
    → chiral-step entries pick up e^{-i k_y DY[d]}.
    """
    n = 4 * (L + 1)
    A = np.zeros((n, n), dtype=complex)
    for x in range(L + 1):
        for d in range(4):
            src = hpbc_index(x, d, L)
            d_ccw = (d - 1) % 4
            d_cw  = (d + 1) % 4

            # Noise step: stay at x, rotate.  No k_y phase.
            A[hpbc_index(x, d_ccw, L), src] += D_r * omega
            A[hpbc_index(x, d_cw,  L), src] += D_r * (1.0 - omega)

            # Chiral step: translate (DX[d], DY[d]), rotate.  OBC along x.
            nx = x + int(DX[d])
            if 0 <= nx <= L:
                phase = np.exp(-1j * k_y * float(DY[d]))
                A[hpbc_index(nx, d_cw,  L), src] += (1.0 - D_r) * omega       * phase
                A[hpbc_index(nx, d_ccw, L), src] += (1.0 - D_r) * (1.0 - omega) * phase
            else:
                A[src, src] += (1.0 - D_r)        # blocked chiral self-loop
    return A


# ---------------------------------------------------------------------------
# 2. Edge mask (HPBC: only x-walls contribute)
# ---------------------------------------------------------------------------
def edge_mask(L: int) -> np.ndarray:
    n = 4 * (L + 1)
    m = np.zeros(n, dtype=bool)
    # left wall x = 0 : to-edge ← (d=3), CCW-along ↓ (d=2)
    for d in (3, 2):
        m[hpbc_index(0, d, L)] = True
    # right wall x = L: to-edge → (d=1), CCW-along ↑ (d=0)
    for d in (1, 0):
        m[hpbc_index(L, d, L)] = True
    return m


# ---------------------------------------------------------------------------
# 3. Spectrum + edge weights at fixed k_y
# ---------------------------------------------------------------------------
def hpbc_spectrum(omega: float, D_r: float, L: int, k_y: float):
    A = build_hpbc_matrix(omega, D_r, L, k_y)
    evals, evecs = np.linalg.eig(A)
    abs2 = np.abs(evecs) ** 2
    total = abs2.sum(axis=0)
    total = np.where(total > 1e-30, total, 1.0)
    edge = edge_mask(L)
    return evals, abs2[edge, :].sum(axis=0) / total


# ---------------------------------------------------------------------------
# 4. Author cross-check at k_y = 0 vs OBC-in-x + PBC-in-y at k_y = 0
# ---------------------------------------------------------------------------
# At k_y = 0, the HPBC matrix equals the y-summed version of the full 2-D
# matrix.  We compare against the authors' full 2-D (both-direction-OBC)
# matrix TRACED over y — but authors only ship 2-D OBC, so there's no
# direct author cross-check.  The physically-useful check is that at
# k_y = 0 the HPBC matrix is real, column-stochastic, and has λ = 1.
# A deeper check exists against an independent PBC torus build — not
# included here (the fig4b torus code is a 2-D PBC anyway, different
# problem).
# ---------------------------------------------------------------------------
def self_consistency(omega: float = 1.0, D_r: float = 0.1, L: int = 10):
    # 1. k_y = 0: real, column-stochastic, λ = 1 present
    A0 = build_hpbc_matrix(omega, D_r, L, 0.0)
    max_im = float(np.max(np.abs(A0.imag)))
    col = A0.sum(axis=0)
    max_col = float(np.max(np.abs(col - 1.0)))
    ev0 = np.linalg.eigvals(A0)
    dist1 = float(np.min(np.abs(ev0 - 1.0)))
    print(f"  k_y = 0:   max |Im A| = {max_im:.2e}")
    print(f"             max |col_sum - 1| = {max_col:.2e}")
    print(f"             min |λ - 1| = {dist1:.2e}")
    assert max_im < 1e-14
    assert max_col < 1e-14
    assert dist1 < 1e-10

    # 2. k_y = π: real (all phases are ±1)
    Ap = build_hpbc_matrix(omega, D_r, L, np.pi)
    max_im_pi = float(np.max(np.abs(Ap.imag)))
    print(f"  k_y = π:   max |Im A| = {max_im_pi:.2e}")
    assert max_im_pi < 1e-14

    # 3. Conjugate symmetry k_y ↔ -k_y
    for kv in (0.3, 1.1, 2.0):
        e_p = np.linalg.eigvals(build_hpbc_matrix(omega, D_r, L, +kv))
        e_m = np.linalg.eigvals(build_hpbc_matrix(omega, D_r, L, -kv))
        re_diff = float(np.max(np.abs(np.sort(e_p.real) - np.sort(e_m.real))))
        im_diff = float(np.max(np.abs(np.sort(e_p.imag) - np.sort(-e_m.imag))))
        print(f"  k_y = ±{kv}:  max|ΔRe| = {re_diff:.2e}, max|ΔIm+Im'| = {im_diff:.2e}")
        assert re_diff < 1e-10
        assert im_diff < 1e-10

    # 4. Matrix size
    n_expect = 4 * (L + 1)
    assert A0.shape == (n_expect, n_expect)
    print(f"  matrix size: {A0.shape}  (expected {n_expect}x{n_expect})")

    # 5. |λ| ≤ 1 at a few k_y
    for kv in (-np.pi, -1.0, 0.0, 1.0, np.pi):
        ev = np.linalg.eigvals(build_hpbc_matrix(omega, D_r, L, kv))
        assert np.max(np.abs(ev)) <= 1.0 + 1e-12
    print(f"  |λ| ≤ 1 at k_y ∈ {{-π, -1, 0, 1, π}}")


# ---------------------------------------------------------------------------
# 5. Scan k_y
# ---------------------------------------------------------------------------
def scan_k_y(omega: float = 1.0, D_r: float = 0.1, L: int = 10, N_k: int = 200):
    k_vals = np.linspace(-np.pi, np.pi, N_k)
    n = 4 * (L + 1)
    evals = np.zeros((N_k, n), dtype=complex)
    wts   = np.zeros((N_k, n))
    edge = edge_mask(L)

    print(f"  scanning {N_k} k_y points  (matrix {n}x{n}, ω={omega}, D_r={D_r}, L={L})")
    t0 = time.time()
    for i, k in enumerate(k_vals):
        A = build_hpbc_matrix(omega, D_r, L, k)
        e, v = np.linalg.eig(A)
        a2 = np.abs(v) ** 2
        tot = a2.sum(axis=0); tot = np.where(tot > 1e-30, tot, 1.0)
        evals[i, :] = e
        wts[i, :]   = a2[edge, :].sum(axis=0) / tot
        if (i + 1) % max(1, N_k // 10) == 0:
            print(f"    {i+1}/{N_k}  cpu {time.time()-t0:.1f}s")
    return k_vals, evals, wts


# ---------------------------------------------------------------------------
# 6. Figure
# ---------------------------------------------------------------------------
def make_figure(omega: float = 1.0, D_r: float = 0.1, L: int = 10,
                N_k: int = 400, savepath: str = "tcrw_fig4h_paper.png"):
    k_vals, evals, wts = scan_k_y(omega, D_r, L, N_k)

    N_k_, n_states = evals.shape
    k_rep = np.broadcast_to(k_vals[:, None], (N_k_, n_states)).ravel()
    re_flat = evals.real.ravel()
    im_flat = evals.imag.ravel()
    w_flat  = wts.ravel()

    fig, (ax_re, ax_im) = plt.subplots(2, 1, figsize=(6.5, 7),
                                        sharex=True, gridspec_kw=dict(hspace=0.15))
    cmap = plt.cm.turbo

    sc_r = ax_re.scatter(k_rep, re_flat, c=w_flat, cmap=cmap,
                         vmin=0, vmax=1, s=4, alpha=0.85, edgecolors="none")
    ax_im.scatter(k_rep, im_flat, c=w_flat, cmap=cmap,
                   vmin=0, vmax=1, s=4, alpha=0.85, edgecolors="none")

    for ax, ylab in [(ax_re, r"$\mathrm{Re}(\lambda)$"),
                     (ax_im, r"$\mathrm{Im}(\lambda)$")]:
        ax.set_ylabel(ylab)
        ax.set_ylim(-1.05, 1.05)
        ax.set_xlim(-np.pi, np.pi)
        ax.set_xticks([-np.pi, 0, np.pi])
        ax.set_xticklabels([r"$-\pi$", "0", r"$\pi$"])
        ax.set_yticks([-1, 0, 1])
        ax.axhline(0, color="0.6", lw=0.4, ls="--")
    ax_im.set_xlabel(r"$k_y$")
    ax_re.set_title(rf"(h)  HPBC spectrum at $\omega={omega}$, $D_r={D_r}$, $L={L}$",
                    fontsize=11)

    cbar = fig.colorbar(sc_r, ax=[ax_re, ax_im], shrink=0.7, aspect=25, pad=0.03)
    cbar.set_label("edge")
    cbar.set_ticks([0.0, 0.5, 1.0])

    plt.savefig(savepath, dpi=220, bbox_inches="tight")
    print(f"  saved {savepath}")


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 62)
    print(" TCRW Fig 4(h) — HPBC Re/Im(λ) vs k_y at ω = 1, D_r = 0.1, L = 10")
    print("=" * 62)

    print("\n[1] Self-consistency (no author cross-check available)")
    self_consistency(omega=1.0, D_r=0.1, L=10)

    print("\n[2] Production figure")
    make_figure(omega=1.0, D_r=0.1, L=10, N_k=400, savepath="tcrw_fig4h_paper.png")

    print("\nDone.")
