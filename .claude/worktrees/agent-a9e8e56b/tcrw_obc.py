"""
TCRW Phase 2: OBC — Edge localization and probability distributions
=====================================================================

Reproduces:
  Fig 2(a): P(X,Y) heatmap for omega=1 (chiral), OBC, L=10, D_r=1e-3
  Fig 2(f): P(X,Y) heatmap for omega=0 (achiral), OBC, L=10, D_r=1e-3
  Fig 3(a): P_edge/P_bulk vs D_r for omega=1, various L
  Fig 3(f): P_edge/P_bulk vs omega for D_r=1e-3

Two approaches:
  1. Exact steady state: build the full (4*L^2) x (4*L^2) column-stochastic
     transition matrix P, find its right eigenvector with eigenvalue 1
     (stationary distribution satisfying P @ pi = pi).
  2. Long-trajectory MC: run T~10^8-10^10 steps, histogram visits.

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import time

# Direction vectors: d=0(↑), 1(→), 2(↓), 3(←)
DX = np.array([0, 1, 0, -1])
DY = np.array([1, 0, -1, 0])


# ============================================================
# Exact transition matrix approach
# ============================================================

def state_index(x, y, d, L):
    """Map (x, y, d) -> linear index.  d ∈ {0,1,2,3}, x,y ∈ {0,..,L-1}."""
    return d * L * L + y * L + x


def build_transition_matrix(omega, D_r, L):
    """
    Build the (4*L^2) x (4*L^2) transition matrix P for OBC.

    P[j, i] = probability of going from state i to state j.
    (Column-stochastic: columns sum to 1.)

    States: index = d*L^2 + y*L + x, where d=0(↑),1(→),2(↓),3(←).

    Rules:
      Noise step (prob D_r): stay at (x,y), rotate d.
        CCW (d-1 mod 4) with prob omega.
        CW  (d+1 mod 4) with prob (1-omega).

      Chiral step (prob 1-D_r): translate in direction d, then rotate.
        New position: (x + DX[d], y + DY[d]).
        If new position is outside [0,L-1]^2: BLOCKED (stay, no rotation).
        If valid: move there, then rotate with OPPOSITE chirality:
          CW  (d+1 mod 4) with prob omega.
          CCW (d-1 mod 4) with prob (1-omega).
    """
    N = 4 * L * L
    rows = []
    cols = []
    vals = []

    for x in range(L):
        for y in range(L):
            for d in range(4):
                i = state_index(x, y, d, L)

                # --- Noise step (prob D_r) ---
                # CCW: d -> (d-1) % 4, prob omega * D_r
                d_ccw = (d - 1) % 4
                j_ccw = state_index(x, y, d_ccw, L)
                rows.append(j_ccw)
                cols.append(i)
                vals.append(omega * D_r)

                # CW: d -> (d+1) % 4, prob (1-omega) * D_r
                d_cw = (d + 1) % 4
                j_cw = state_index(x, y, d_cw, L)
                rows.append(j_cw)
                cols.append(i)
                vals.append((1 - omega) * D_r)

                # --- Chiral step (prob 1-D_r) ---
                nx = x + DX[d]
                ny = y + DY[d]

                if 0 <= nx < L and 0 <= ny < L:
                    # Valid move: translate, then rotate (opposite chirality)
                    # CW: d -> (d+1) % 4 with prob omega * (1-D_r)
                    d_cw_c = (d + 1) % 4
                    j1 = state_index(nx, ny, d_cw_c, L)
                    rows.append(j1)
                    cols.append(i)
                    vals.append(omega * (1 - D_r))

                    # CCW: d -> (d-1) % 4 with prob (1-omega) * (1-D_r)
                    d_ccw_c = (d - 1) % 4
                    j2 = state_index(nx, ny, d_ccw_c, L)
                    rows.append(j2)
                    cols.append(i)
                    vals.append((1 - omega) * (1 - D_r))
                else:
                    # Blocked: stay at (x,y,d), no rotation
                    j_stay = state_index(x, y, d, L)
                    rows.append(j_stay)
                    cols.append(i)
                    vals.append(1 - D_r)

    P = sp.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsc()
    return P


def exact_steady_state(omega, D_r, L):
    """
    Find the steady-state distribution pi of the transition matrix P.

    P is column-stochastic: sum_j P[j,i] = 1 for all i.
    Stationary distribution satisfies pi^T P = pi^T,
    i.e., P^T pi = pi  (LEFT eigenvector of P = RIGHT eigenvector of P^T).

    Strategy:
      - L <= 30 (~3600 states): dense eigensolve (reliable)
      - L > 30: sparse eigensolve with shift-invert near eigenvalue 1

    Returns P(x,y) = sum_d pi(x,y,d), normalized so sum = 1.
    """
    P = build_transition_matrix(omega, D_r, L)
    N = P.shape[0]

    if L <= 30:
        # Dense eigensolve — fully reliable
        P_dense = P.toarray()
        # Stationary distribution: right eigenvector of P with eigenvalue 1
        # P is column-stochastic: P[j,i] = prob(i→j), cols sum to 1
        # We need P @ pi = pi, i.e., RIGHT eigenvector of P
        eigenvalues, eigenvectors = np.linalg.eig(P_dense)
        # Find the eigenvalue closest to 1
        idx = np.argmin(np.abs(eigenvalues - 1.0))
        pi = eigenvectors[:, idx].real
    else:
        # Sparse: shift-invert around sigma=1 for robustness
        # Find right eigenvector of P with eigenvalue 1
        try:
            eigenvalues, eigenvectors = spla.eigs(P, k=1, sigma=1.0,
                                                   which='LM', tol=1e-10)
            pi = eigenvectors[:, 0].real
        except Exception:
            # Fallback: power iteration (P @ pi converges to steady state)
            pi = np.ones(N) / N
            for _ in range(50000):
                pi_new = P @ pi
                pi_new /= np.sum(np.abs(pi_new))
                if np.max(np.abs(pi_new - pi)) < 1e-14:
                    break
                pi = pi_new
            pi = pi_new

    # Ensure positive and normalize
    if np.sum(pi) < 0:
        pi = -pi
    pi = np.abs(pi)  # remove numerical noise
    pi = pi / np.sum(pi)

    # Sum over director states to get P(x,y)
    Pxy = np.zeros((L, L))
    for d in range(4):
        offset = d * L * L
        Pxy += pi[offset:offset + L * L].reshape(L, L)

    return Pxy, pi


def compute_edge_bulk_ratio(Pxy, L, edge_width=1, per_site=True):
    """
    Compute P_edge / P_bulk from the spatial probability distribution.

    Edge sites: within edge_width of the boundary.
    Bulk sites: everything else.

    If per_site=True (default), returns the ratio of AVERAGE per-site
    probabilities: (P_edge/n_edge) / (P_bulk/n_bulk).
    This is independent of system size L (as in the paper).

    If per_site=False, returns the ratio of TOTAL probabilities.
    """
    mask_edge = np.zeros((L, L), dtype=bool)
    mask_edge[:edge_width, :] = True   # left edge (x=0..ew-1)
    mask_edge[-edge_width:, :] = True  # right edge
    mask_edge[:, :edge_width] = True   # bottom edge (y=0..ew-1)
    mask_edge[:, -edge_width:] = True  # top edge

    n_edge = mask_edge.sum()
    n_bulk = (~mask_edge).sum()

    P_edge = np.sum(Pxy[mask_edge])
    P_bulk = np.sum(Pxy[~mask_edge])

    if P_bulk == 0 or n_bulk == 0:
        return np.inf

    if per_site:
        return (P_edge / n_edge) / (P_bulk / n_bulk)
    else:
        return P_edge / P_bulk


# ============================================================
# Long-trajectory MC approach (for verification / large L)
# ============================================================

def mc_steady_state(omega, D_r, L, T_steps=10**8, N_traj=1, seed=42):
    """
    MC steady state using the vectorized OBC simulator from tcrw_core.
    Uses N_traj parallel walkers for faster convergence.
    Returns visit histogram P(x,y) normalized to sum=1.
    """
    from tcrw_core import simulate_tcrw_obc
    res = simulate_tcrw_obc(omega, D_r, L, T_steps, N_traj=N_traj,
                             seed=seed, track_visits=True)
    Pxy = res['visits'].astype(np.float64)
    Pxy /= Pxy.sum()
    return Pxy


# ============================================================
# Figures
# ============================================================

def fig2_heatmaps():
    """
    P(X,Y) for omega=1 (CW chiral), omega=0.5 (achiral), omega=0 (CCW chiral).
    L=10, D_r=1e-3, OBC.

    Use exact transition matrix (L=10 → 400×400 matrix, trivial).
    Corresponds to Figs 2(a) and 2(f) in the paper.
    """
    L = 10
    D_r = 1e-3
    omegas = [1.0, 0.5, 0.0]
    labels = ['(a) $\\omega=1$ (CW chiral)',
              '$\\omega=0.5$ (achiral)',
              '(f) $\\omega=0$ (CCW chiral)']

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    for ax, omega, label in zip(axes, omegas, labels):
        print(f"  Computing exact steady state: omega={omega}, D_r={D_r}, L={L}...")
        t0 = time.time()
        Pxy, _ = exact_steady_state(omega, D_r, L)
        print(f"    Done in {time.time()-t0:.1f}s")

        ratio = compute_edge_bulk_ratio(Pxy, L, per_site=True)
        print(f"    P_edge/P_bulk (per-site) = {ratio:.2f}")

        # Heatmap — use log scale like the paper
        im = ax.imshow(Pxy.T, origin='lower', cmap='hot',
                        norm=LogNorm(vmin=Pxy[Pxy > 0].min(), vmax=Pxy.max()),
                        extent=[-0.5, L-0.5, -0.5, L-0.5])
        ax.set_xlabel('X', fontsize=12)
        ax.set_ylabel('Y', fontsize=12)
        ax.set_title(f'{label}\n$D_r = {D_r}$, $L = {L}$', fontsize=12)
        fig.colorbar(im, ax=ax, label='$P(X,Y)$', shrink=0.85)

        # Annotate P_edge/P_bulk
        ax.text(0.03, 0.03, f'$P_{{edge}}/P_{{bulk}}={ratio:.1f}$',
                transform=ax.transAxes, fontsize=9, color='white',
                bbox=dict(boxstyle='round', facecolor='black', alpha=0.6),
                verticalalignment='bottom')

    plt.suptitle('Fig 2: Steady-state probability $P(X,Y)$ (OBC)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('tcrw_fig2_Pxy_heatmaps.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig2_Pxy_heatmaps.png")


def fig3a_Pedge_vs_Dr():
    """
    Fig 3(a): P_edge/P_bulk vs D_r for omega=1, various L.
    Uses exact transition matrix for all sizes.

    Paper uses: L = 4, 9, 19, 49, 99, 199, 499.
    Uses per-site normalization: (P_edge/n_edge) / (P_bulk/n_bulk),
    which makes the ratio independent of system size L.
    """
    omega = 1.0
    # D_r grid: logarithmic from 1e-4 to 1
    D_r_values = np.logspace(-4, 0, 25)

    # L values — go as far as feasible
    L_values = [4, 9, 19, 49, 99]
    markers = ['o', 's', '^', 'D', 'v']
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(L_values)))

    fig, ax = plt.subplots(figsize=(7, 5))

    for L, marker, color in zip(L_values, markers, colors):
        ratios = []
        print(f"\n  L = {L} (matrix size {4*L*L}×{4*L*L}):")
        for D_r in D_r_values:
            try:
                Pxy, _ = exact_steady_state(omega, D_r, L)
                r = compute_edge_bulk_ratio(Pxy, L)
                ratios.append(r)
                print(f"    D_r={D_r:.4f}: P_edge/P_bulk = {r:.4f}")
            except Exception as e:
                print(f"    D_r={D_r:.4f}: FAILED ({e})")
                ratios.append(np.nan)

        ax.loglog(D_r_values, ratios, marker=marker, color=color, ms=5, lw=1.2,
                  label=f'$L = {L}$', markerfacecolor='none')

    ax.set_xlabel('$D_r$', fontsize=13)
    ax.set_ylabel('$P_{\\rm edge} / P_{\\rm bulk}$', fontsize=13)
    ax.set_title('Fig 3(a): Edge localization vs $D_r$ ($\\omega = 1$)', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(5e-5, 2)

    plt.tight_layout()
    plt.savefig('tcrw_fig3a_Pedge_vs_Dr.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("\nSaved: tcrw_fig3a_Pedge_vs_Dr.png")


def fig3f_Pedge_vs_omega():
    """
    Fig 3(f): P_edge/P_bulk vs omega for D_r = 1e-3.
    Should be approximately independent of omega!

    Use exact steady state with moderate L.
    """
    D_r = 1e-3
    omega_values = np.linspace(0.0, 1.0, 21)
    L_values = [10, 19, 49]
    markers = ['o', 's', '^']
    colors = ['tab:blue', 'tab:orange', 'tab:green']

    fig, ax = plt.subplots(figsize=(7, 5))

    for L, marker, color in zip(L_values, markers, colors):
        ratios = []
        print(f"\n  L = {L}:")
        for omega in omega_values:
            Pxy, _ = exact_steady_state(omega, D_r, L)
            r = compute_edge_bulk_ratio(Pxy, L)
            ratios.append(r)
            print(f"    omega={omega:.2f}: P_edge/P_bulk = {r:.2f}")

        ax.plot(omega_values, ratios, marker=marker, color=color, ms=6, lw=1.5,
                label=f'$L = {L}$', markerfacecolor='none')

    ax.set_xlabel('$\\omega$', fontsize=13)
    ax.set_ylabel('$P_{\\rm edge} / P_{\\rm bulk}$', fontsize=13)
    ax.set_title('Fig 3(f): Edge localization vs $\\omega$ ($D_r = 10^{-3}$)', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(-0.02, 1.02)

    plt.tight_layout()
    plt.savefig('tcrw_fig3f_Pedge_vs_omega.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("\nSaved: tcrw_fig3f_Pedge_vs_omega.png")


def fig2_mc_verification():
    """
    Verify exact steady state against MC simulation for L=10.
    Overlay MC histogram and exact P(x,y).
    """
    L = 10
    D_r = 1e-3
    omega = 1.0
    T_mc = 10**8  # 10^8 steps — should be enough for L=10

    print(f"  Exact steady state (omega={omega}, D_r={D_r}, L={L})...")
    Pxy_exact, _ = exact_steady_state(omega, D_r, L)

    # Use multiple shorter trajectories for faster convergence via vectorization
    N_mc = 100  # parallel walkers
    T_per = T_mc // N_mc
    print(f"  MC simulation: {N_mc} walkers × {T_per:.0e} steps...")
    t0 = time.time()
    Pxy_mc = mc_steady_state(omega, D_r, L, T_steps=T_per, N_traj=N_mc, seed=42)
    print(f"    MC done in {time.time()-t0:.1f}s")

    # Compare: scatter plot of exact vs MC for each site
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # (a) Exact heatmap
    im0 = axes[0].imshow(Pxy_exact.T, origin='lower', cmap='hot',
                          norm=LogNorm(), extent=[-0.5, L-0.5, -0.5, L-0.5])
    axes[0].set_title('Exact steady state', fontsize=12)
    axes[0].set_xlabel('X'); axes[0].set_ylabel('Y')
    fig.colorbar(im0, ax=axes[0], shrink=0.85)

    # (b) MC heatmap
    im1 = axes[1].imshow(Pxy_mc.T, origin='lower', cmap='hot',
                          norm=LogNorm(), extent=[-0.5, L-0.5, -0.5, L-0.5])
    axes[1].set_title(f'MC ($T = {T_mc:.0e}$)', fontsize=12)
    axes[1].set_xlabel('X'); axes[1].set_ylabel('Y')
    fig.colorbar(im1, ax=axes[1], shrink=0.85)

    # (c) Scatter: exact vs MC
    pe = Pxy_exact.flatten()
    pm = Pxy_mc.flatten()
    axes[2].scatter(pe, pm, s=15, alpha=0.7, edgecolors='k', linewidths=0.5)
    pmin, pmax = min(pe.min(), pm.min()), max(pe.max(), pm.max())
    axes[2].plot([pmin, pmax], [pmin, pmax], 'r--', lw=1)
    axes[2].set_xlabel('$P(x,y)$ exact', fontsize=12)
    axes[2].set_ylabel('$P(x,y)$ MC', fontsize=12)
    axes[2].set_title('Exact vs MC verification', fontsize=12)
    axes[2].set_xscale('log'); axes[2].set_yscale('log')

    rel_err = np.max(np.abs(pe - pm) / (pe + 1e-20))
    axes[2].text(0.05, 0.92, f'max rel. error: {rel_err:.2%}',
                 transform=axes[2].transAxes, fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.suptitle(f'Phase 2 verification: $\\omega={omega}$, $D_r={D_r}$, $L={L}$',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('tcrw_fig2_mc_verification.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig2_mc_verification.png")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("=" * 60)
    print("TCRW Phase 2: OBC — Edge localization")
    print("=" * 60)

    print("\n--- Fig 2: P(X,Y) heatmaps (exact) ---")
    fig2_heatmaps()

    print("\n--- Fig 3(a): P_edge/P_bulk vs D_r ---")
    fig3a_Pedge_vs_Dr()

    print("\n--- Fig 3(f): P_edge/P_bulk vs omega ---")
    fig3f_Pedge_vs_omega()

    print("\n--- MC verification ---")
    fig2_mc_verification()

    print("\n" + "=" * 60)
    print("Phase 2 complete.")
    print("=" * 60)
