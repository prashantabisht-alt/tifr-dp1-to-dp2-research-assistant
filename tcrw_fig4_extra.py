"""
TCRW Extra figures: Advanced spectrum analysis
==============================================

Reproduces specialized spectrum figures:
  Fig 4(d): Re(λ) of OBC spectrum vs D_r for ω=1, with edge localization coloring
  Fig 4(e): Re(λ) of OBC spectrum vs ω for fixed D_r=0.1, with edge localization coloring
  Fig 4(i): PBC band structure in (cos kx, cos ky) plane
  Fig 8(d)-(l): 3×3 grid of OBC complex-plane spectra at 9 parameter points

Each eigenvalue is colored by its edge localization weight:
  - Edge weight = fraction of |eigenvector|² on boundary sites
  - Yellow (weight ~ 0) = bulk-localized
  - Red (weight ~ 1) = edge-localized

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
import matplotlib.cm as cm
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tcrw_obc import build_transition_matrix, state_index
from tcrw_spectrum import build_Pk, pbc_eigenvalues


# ============================================================
# Helper functions
# ============================================================

def compute_edge_weight(eigenvector, L):
    """
    Compute edge localization weight for a single eigenvector.

    Edge sites: (x,y) with x==0 or x==L-1 or y==0 or y==L-1.
    Edge weight = sum of |v|^2 over edge sites / total sum of |v|^2.

    Parameters
    ----------
    eigenvector : ndarray of shape (4*L*L,)
        A single eigenvector (right eigenvector).
    L : int
        System size.

    Returns
    -------
    edge_weight : float in [0, 1]
        Fraction of probability on edge sites.
    """
    edge_mask = np.zeros(4 * L * L, dtype=bool)
    for d in range(4):
        offset = d * L * L
        for x in range(L):
            for y in range(L):
                if x == 0 or x == L - 1 or y == 0 or y == L - 1:
                    edge_mask[offset + y * L + x] = True

    prob = np.abs(eigenvector) ** 2
    total = prob.sum()
    if total > 0:
        edge_weight = prob[edge_mask].sum() / total
    else:
        edge_weight = 0.0
    return edge_weight


def obc_spectrum_with_weights(omega, D_r, L):
    """
    Compute full OBC eigenvalue spectrum with edge weights for each eigenvector.

    Returns
    -------
    eigenvalues : ndarray, shape (4*L*L,)
        All eigenvalues of the transition matrix.
    edge_weights : ndarray, shape (4*L*L,)
        Edge localization weight for each eigenvalue.
    """
    P = build_transition_matrix(omega, D_r, L)
    P_dense = P.toarray()

    eigenvalues, eigenvectors = np.linalg.eig(P_dense)

    edge_weights = np.zeros(len(eigenvalues))
    for n in range(len(eigenvalues)):
        v = eigenvectors[:, n]
        edge_weights[n] = compute_edge_weight(v, L)

    return eigenvalues, edge_weights


# ============================================================
# Fig 4(d): Re(λ) vs D_r
# ============================================================

def fig4d_spectrum_vs_Dr():
    """
    Fig 4(d): Re(λ) of OBC spectrum vs D_r for ω=1.

    For each D_r in [0.01, 0.5] (~25 values), compute all eigenvalues
    of the (4L²)×(4L²) OBC matrix, plot Re(λ) as scatter points colored
    by edge weight.
    """
    L = 10
    omega = 1.0
    D_r_values = np.linspace(0.01, 0.5, 25)

    print("Fig 4(d): Re(λ) vs D_r (ω=1)")
    print(f"  L={L}, ω={omega}")
    print(f"  D_r sweep: {D_r_values.min():.3f} to {D_r_values.max():.3f} ({len(D_r_values)} points)")

    fig, ax = plt.subplots(figsize=(10, 6))

    # Collect all eigenvalues and weights for scatter plot
    all_re_lambda = []
    all_edge_weights = []
    all_D_r = []

    for D_r in D_r_values:
        print(f"  D_r = {D_r:.4f}...", end=" ", flush=True)
        try:
            evals, weights = obc_spectrum_with_weights(omega, D_r, L)
            print(f"✓ ({len(evals)} states)")
            all_re_lambda.extend(evals.real)
            all_edge_weights.extend(weights)
            all_D_r.extend([D_r] * len(evals))
        except Exception as e:
            print(f"✗ ({e})")

    # Create scatter plot with colormap
    all_re_lambda = np.array(all_re_lambda)
    all_edge_weights = np.array(all_edge_weights)
    all_D_r = np.array(all_D_r)

    sc = ax.scatter(all_D_r, all_re_lambda, c=all_edge_weights, cmap='YlOrRd',
                    s=15, alpha=0.6, edgecolors='none', vmin=0, vmax=1)

    ax.set_xlabel(r'$D_r$', fontsize=13)
    ax.set_ylabel(r'Re($\lambda$)', fontsize=13)
    ax.set_title(f'Fig 4(d): OBC spectrum Re($\\lambda$) vs $D_r$ ($\\omega={omega}$, $L={L}$)',
                 fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(D_r_values.min() - 0.02, D_r_values.max() + 0.02)

    cbar = plt.colorbar(sc, ax=ax, label='Edge weight', shrink=0.85)

    plt.tight_layout()
    plt.savefig('tcrw_fig4d_spectrum_vs_Dr.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig4d_spectrum_vs_Dr.png\n")


# ============================================================
# Fig 4(e): Re(λ) vs ω
# ============================================================

def fig4e_spectrum_vs_omega():
    """
    Fig 4(e): Re(λ) of OBC spectrum vs ω for fixed D_r=0.1.

    For each ω in [0, 1] (~25 values), compute all eigenvalues,
    plot Re(λ) as scatter points colored by edge weight.
    """
    L = 10
    D_r = 0.1
    omega_values = np.linspace(0.0, 1.0, 25)

    print("Fig 4(e): Re(λ) vs ω (D_r=0.1)")
    print(f"  L={L}, D_r={D_r}")
    print(f"  ω sweep: {omega_values.min():.3f} to {omega_values.max():.3f} ({len(omega_values)} points)")

    fig, ax = plt.subplots(figsize=(10, 6))

    all_re_lambda = []
    all_edge_weights = []
    all_omega = []

    for omega in omega_values:
        print(f"  ω = {omega:.4f}...", end=" ", flush=True)
        try:
            evals, weights = obc_spectrum_with_weights(omega, D_r, L)
            print(f"✓ ({len(evals)} states)")
            all_re_lambda.extend(evals.real)
            all_edge_weights.extend(weights)
            all_omega.extend([omega] * len(evals))
        except Exception as e:
            print(f"✗ ({e})")

    all_re_lambda = np.array(all_re_lambda)
    all_edge_weights = np.array(all_edge_weights)
    all_omega = np.array(all_omega)

    sc = ax.scatter(all_omega, all_re_lambda, c=all_edge_weights, cmap='YlOrRd',
                    s=15, alpha=0.6, edgecolors='none', vmin=0, vmax=1)

    ax.set_xlabel(r'$\omega$', fontsize=13)
    ax.set_ylabel(r'Re($\lambda$)', fontsize=13)
    ax.set_title(f'Fig 4(e): OBC spectrum Re($\\lambda$) vs $\\omega$ ($D_r={D_r}$, $L={L}$)',
                 fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.03, 1.03)
    ax.axvline(x=0.5, color='red', ls='--', lw=1, alpha=0.5, label='$\\omega=0.5$ (critical)')
    ax.legend(fontsize=10)

    cbar = plt.colorbar(sc, ax=ax, label='Edge weight', shrink=0.85)

    plt.tight_layout()
    plt.savefig('tcrw_fig4e_spectrum_vs_omega.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig4e_spectrum_vs_omega.png\n")


# ============================================================
# Fig 4(i): PBC band structure on (cos kx, cos ky) plane
# ============================================================

def fig4i_band_circle():
    """
    Fig 4(i): PBC band structure plotted in (cos kx, cos ky) space.

    Sweep k-space on a 60×60 grid in [-π, π)².
    For each k-point, compute 4 eigenvalues of P(k).
    Plot in (cos kx, cos ky) plane, colored by Re(λ) or |λ|.

    The (cos kx, cos ky) representation maps the Brillouin zone to a
    filled circle (the "extended zone" diagram).
    """
    L = 10  # Not used for PBC, but we keep naming consistent
    omega = 1.0
    D_r = 0.1
    Nk = 60

    print("Fig 4(i): PBC band structure in (cos kx, cos ky) plane")
    print(f"  ω={omega}, D_r={D_r}")
    print(f"  k-grid: {Nk}×{Nk}")

    # Generate k-grid
    kx_arr = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    ky_arr = np.linspace(-np.pi, np.pi, Nk, endpoint=False)

    # Collect eigenvalues
    cos_kx_arr = []
    cos_ky_arr = []
    re_lambda_arr = []

    for kx in kx_arr:
        for ky in ky_arr:
            evals = pbc_eigenvalues(omega, D_r, kx, ky)
            cos_kx = np.cos(kx)
            cos_ky = np.cos(ky)

            for lam in evals:
                cos_kx_arr.append(cos_kx)
                cos_ky_arr.append(cos_ky)
                re_lambda_arr.append(lam.real)

    cos_kx_arr = np.array(cos_kx_arr)
    cos_ky_arr = np.array(cos_ky_arr)
    re_lambda_arr = np.array(re_lambda_arr)

    # Create figure
    fig, ax = plt.subplots(figsize=(9, 8.5))

    sc = ax.scatter(cos_kx_arr, cos_ky_arr, c=re_lambda_arr, cmap='RdYlBu_r',
                    s=10, alpha=0.6, edgecolors='none', vmin=-1.0, vmax=1.0)

    ax.set_xlabel(r'$\cos k_x$', fontsize=13)
    ax.set_ylabel(r'$\cos k_y$', fontsize=13)
    ax.set_title(f'Fig 4(i): PBC spectrum ($\\cos k_x$, $\\cos k_y$) space\n'
                 f'$\\omega={omega}$, $D_r={D_r}$, {Nk}×{Nk} k-grid',
                 fontsize=13)

    ax.set_aspect('equal')
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)

    # Draw unit circle for reference
    theta = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(theta), np.sin(theta), 'k--', lw=0.5, alpha=0.3, label='unit circle')
    ax.legend(fontsize=10)

    cbar = plt.colorbar(sc, ax=ax, label=r'Re($\lambda$)', shrink=0.85)

    plt.tight_layout()
    plt.savefig('tcrw_fig4i_band_circle.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig4i_band_circle.png\n")


# ============================================================
# Fig 8(d)-(l): 3×3 grid of OBC spectra
# ============================================================

def fig8_spectrum_grid():
    """
    Fig 8(d)-(l): 3×3 grid of OBC complex-plane spectra.

    9 parameter points: ω ∈ {0.2, 0.5, 0.8} × D_r ∈ {0.05, 0.3, 0.7}
    For each: compute OBC spectrum at L=10, plot in complex plane
    with colorbar for edge weight (yellow→orange→red).
    """
    L = 10
    omega_vals = [0.2, 0.5, 0.8]
    D_r_vals = [0.05, 0.3, 0.7]

    print("Fig 8: 3×3 grid of OBC complex-plane spectra")
    print(f"  L={L}")
    print(f"  ω ∈ {omega_vals}, D_r ∈ {D_r_vals}")

    fig, axes = plt.subplots(3, 3, figsize=(14, 13))

    all_weights_for_colorbar = []

    for i, D_r in enumerate(D_r_vals):
        for j, omega in enumerate(omega_vals):
            ax = axes[i, j]
            print(f"  ω={omega}, D_r={D_r}...", end=" ", flush=True)

            try:
                evals, weights = obc_spectrum_with_weights(omega, D_r, L)
                print(f"✓ ({len(evals)} states)")

                all_weights_for_colorbar.extend(weights)

                # Scatter plot with edge weight color (YlOrRd: yellow=bulk, red=edge)
                sc = ax.scatter(evals.real, evals.imag, c=weights, cmap='YlOrRd',
                               s=20, alpha=0.7, edgecolors='none', vmin=0, vmax=1)

                ax.set_xlabel(r'Re($\lambda$)', fontsize=10)
                ax.set_ylabel(r'Im($\lambda$)', fontsize=10)
                ax.set_title(f'$\\omega={omega}$, $D_r={D_r}$', fontsize=11, fontweight='bold')
                ax.set_aspect('equal')
                ax.set_xlim(-1.15, 1.15)
                ax.set_ylim(-0.65, 0.65)
                ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
                ax.axvline(0, color='gray', lw=0.5, alpha=0.3)

                # Draw unit circle reference
                theta = np.linspace(0, 2 * np.pi, 200)
                ax.plot(np.cos(theta), np.sin(theta), 'k--', lw=0.5, alpha=0.15)

                # Annotate phase: topological if D_r<0.5 AND omega!=0.5
                if abs(omega - 0.5) < 0.01:
                    phase = 'critical'
                elif D_r < 0.5:
                    phase = 'topological'
                else:
                    phase = 'trivial'
                ax.text(0.02, 0.98, phase, transform=ax.transAxes, fontsize=8,
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            except Exception as e:
                print(f"✗ ({e})")
                ax.text(0.5, 0.5, f'ERROR: {str(e)[:20]}', transform=ax.transAxes,
                       ha='center', va='center', fontsize=10, color='red')

    # Add single colorbar spanning all subplots
    fig.colorbar(sc, ax=axes.ravel().tolist(), label='Edge weight',
                pad=0.01, shrink=0.8, aspect=30)

    fig.suptitle(f'Fig 8: OBC spectrum grid ($L={L}$, complex plane)',
                fontsize=14, y=0.995)

    plt.tight_layout()
    plt.savefig('tcrw_fig8_spectrum_grid.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig8_spectrum_grid.png\n")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("=" * 70)
    print("TCRW Extra Figures: Advanced Spectrum Analysis")
    print("=" * 70)

    print("\n--- Fig 4(d): Re(λ) vs D_r ---")
    fig4d_spectrum_vs_Dr()

    print("--- Fig 4(e): Re(λ) vs ω ---")
    fig4e_spectrum_vs_omega()

    print("--- Fig 4(i): PBC band structure in (cos kx, cos ky) plane ---")
    fig4i_band_circle()

    print("--- Fig 8: 3×3 grid of OBC complex-plane spectra ---")
    fig8_spectrum_grid()

    print("=" * 70)
    print("All figures generated successfully!")
    print("=" * 70)
