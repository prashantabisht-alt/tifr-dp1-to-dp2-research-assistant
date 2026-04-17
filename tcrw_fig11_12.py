"""
TCRW Figures 11 & 12: Boundary effects and nested spectra
===========================================================

Reproduces two publication-quality figures showing:

Figure 11 (4 panels):
  (a) OBC lattice L=10 visualization with boundary types
  (b) OBC spectrum colored by external edge localization weight
  (c) PBC + internal defect lattice visualization
  (d) Spectrum of PBC+defect colored by internal boundary weight

Figure 12 (4 panels):
  (a) OBC lattice L=15 with two separated 3×3 holes
  (b) Spectrum colored by EXTERNAL boundary weight
  (c) Internal boundary sites highlighted on the lattice
  (d) Spectrum colored by INTERNAL boundary weight (nested spectrum)

Parameters: omega=1, D_r=0.1
Uses exact eigensolve with dense matrices (feasible for L≤15).

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tcrw_geometry import (RectangleMask, PBCWithDefects, RectangleWithHoles,
                           paper_defect_internal_block, paper_defect_two_holes,
                           paper_pbc_with_internal_defect)
from tcrw_obc import build_transition_matrix_generic


# ============================================================
# Edge/boundary weight calculation
# ============================================================

def compute_edge_weights(P_dense, mask, boundary_type='all'):
    """
    For each eigenvector of P, compute fraction of |v|^2 on boundary sites.

    Parameters
    ----------
    P_dense : ndarray, shape (n_states, n_states)
        Dense transition matrix (column-stochastic)
    mask : LatticeMask
        Lattice geometry defining which sites are boundaries
    boundary_type : str
        'all' = external + internal
        'external' = only outer edge boundaries
        'internal' = only defect/hole boundaries

    Returns
    -------
    eigenvalues : ndarray, shape (n_states,)
        All eigenvalues of P
    weights : ndarray, shape (n_states,)
        For each eigenvector, fraction of |v|^2 on specified boundary sites
    """
    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eig(P_dense)

    n_sites = mask.n_sites

    # Identify boundary site state indices
    bnd_sites = mask.get_boundary_sites(boundary_type)
    bnd_state_indices = set()

    # Each state index i is expanded to 4 * n_sites total states
    # State index = d * n_sites + site_index
    for x, y in bnd_sites:
        si = mask.site_to_index(x, y)
        if si is not None:
            for d in range(4):
                bnd_state_indices.add(d * n_sites + si)

    # Compute weights
    weights = np.zeros(len(eigenvalues))
    for i in range(len(eigenvalues)):
        v = eigenvectors[:, i]
        prob = np.abs(v) ** 2
        prob_sum = np.sum(prob)
        if prob_sum > 1e-30:
            prob = prob / prob_sum
            weights[i] = sum(prob[j] for j in bnd_state_indices)
        else:
            weights[i] = 0.0

    return eigenvalues, weights


# ============================================================
# Figure 11: Boundary effects (4 panels)
# ============================================================

def fig11_boundary_effects():
    """
    Figure 11: Effect of boundaries on spectrum.

    (a) OBC lattice L=10, colored by boundary type
    (b) OBC spectrum, eigenvalues colored by external boundary weight
    (c) PBC + internal defect lattice, colored by boundary type
    (d) PBC+defect spectrum, colored by internal boundary weight
    """
    print("\n" + "=" * 70)
    print("FIGURE 11: Boundary effects on spectrum")
    print("=" * 70)

    omega = 1.0
    D_r = 0.1
    L_obc = 10
    L_pbc = 10

    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    # ========== Panel (a): OBC lattice visualization ==========
    ax_a = fig.add_subplot(gs[0, 0])
    print(f"\nPanel (a): OBC lattice L={L_obc}")
    mask_obc = RectangleMask(L_obc)
    mask_obc.visualize(ax=ax_a, show_types=True)
    ax_a.set_title(f'(a) OBC lattice $L={L_obc}$\n({mask_obc.n_sites} sites, {mask_obc.n_states} states)',
                   fontsize=12, fontweight='bold')

    # ========== Panel (b): OBC spectrum (external boundary weight) ==========
    ax_b = fig.add_subplot(gs[0, 1])
    print(f"\nPanel (b): OBC spectrum (external boundary weight)")
    print(f"  Building transition matrix (omega={omega}, D_r={D_r})...")
    t0 = time.time()
    P_obc = build_transition_matrix_generic(omega, D_r, mask_obc)
    P_obc_dense = P_obc.toarray()
    print(f"  Matrix built in {time.time()-t0:.2f}s, shape {P_obc_dense.shape}")

    print(f"  Computing eigenvalues and edge weights...")
    t0 = time.time()
    evals_obc, weights_obc = compute_edge_weights(P_obc_dense, mask_obc, 'external')
    print(f"  Done in {time.time()-t0:.2f}s")

    # Scatter plot: Re(λ) vs Im(λ), colored by external boundary weight
    norm = Normalize(vmin=0, vmax=1)
    cmap = cm.get_cmap('RdYlBu_r')
    scatter = ax_b.scatter(evals_obc.real, evals_obc.imag, c=weights_obc,
                          cmap=cmap, norm=norm, s=40, alpha=0.7, edgecolors='k', linewidths=0.5)
    cbar_b = plt.colorbar(scatter, ax=ax_b)
    cbar_b.set_label('Edge weight (external)', fontsize=10)

    ax_b.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_b.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_b.set_xlabel('Re($\lambda$)', fontsize=11)
    ax_b.set_ylabel('Im($\lambda$)', fontsize=11)
    ax_b.set_title(f'(b) OBC spectrum ($\omega={omega}$, $D_r={D_r}$)\ncolored by external boundary weight',
                   fontsize=12, fontweight='bold')
    ax_b.grid(True, alpha=0.3)

    # ========== Panel (c): PBC + defect lattice visualization ==========
    ax_c = fig.add_subplot(gs[1, 0])
    print(f"\nPanel (c): PBC + internal defect lattice L={L_pbc}")
    mask_pbc = paper_pbc_with_internal_defect(L_pbc)
    mask_pbc.visualize(ax=ax_c, show_types=True)
    ax_c.set_title(f'(c) PBC + internal defect $L={L_pbc}$\n({mask_pbc.n_sites} sites, {mask_pbc.n_states} states)',
                   fontsize=12, fontweight='bold')

    # ========== Panel (d): PBC+defect spectrum (internal boundary weight) ==========
    ax_d = fig.add_subplot(gs[1, 1])
    print(f"\nPanel (d): PBC+defect spectrum (internal boundary weight)")
    print(f"  Building transition matrix...")
    t0 = time.time()
    P_pbc = build_transition_matrix_generic(omega, D_r, mask_pbc)
    P_pbc_dense = P_pbc.toarray()
    print(f"  Matrix built in {time.time()-t0:.2f}s, shape {P_pbc_dense.shape}")

    print(f"  Computing eigenvalues and internal boundary weights...")
    t0 = time.time()
    evals_pbc, weights_pbc = compute_edge_weights(P_pbc_dense, mask_pbc, 'internal')
    print(f"  Done in {time.time()-t0:.2f}s")

    # Scatter plot: colored by internal boundary weight
    scatter_d = ax_d.scatter(evals_pbc.real, evals_pbc.imag, c=weights_pbc,
                            cmap=cmap, norm=norm, s=40, alpha=0.7, edgecolors='k', linewidths=0.5)
    cbar_d = plt.colorbar(scatter_d, ax=ax_d)
    cbar_d.set_label('Boundary weight (internal)', fontsize=10)

    ax_d.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_d.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_d.set_xlabel('Re($\lambda$)', fontsize=11)
    ax_d.set_ylabel('Im($\lambda$)', fontsize=11)
    ax_d.set_title(f'(d) PBC+defect spectrum ($\omega={omega}$, $D_r={D_r}$)\ncolored by internal boundary weight',
                   fontsize=12, fontweight='bold')
    ax_d.grid(True, alpha=0.3)

    plt.suptitle('Figure 11: Boundary effects on chiral active particle spectrum',
                 fontsize=14, fontweight='bold', y=0.995)

    print(f"\nSaving figure 11...")
    plt.savefig('tcrw_fig11_boundary_effects.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig11_boundary_effects.png")


# ============================================================
# Figure 12: Disconnected boundaries and nested spectrum
# ============================================================

def fig12_nested_spectrum():
    """
    Figure 12: Two disconnected internal boundaries.

    (a) OBC lattice L=15 with two 3×3 holes
    (b) Spectrum colored by EXTERNAL boundary weight
    (c) Lattice with internal boundary sites highlighted
    (d) Spectrum colored by INTERNAL boundary weight (nested ovals)
    """
    print("\n" + "=" * 70)
    print("FIGURE 12: Disconnected boundaries and nested spectrum")
    print("=" * 70)

    omega = 1.0
    D_r = 0.1
    L = 15

    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    # ========== Panel (a): OBC lattice with two holes ==========
    ax_a = fig.add_subplot(gs[0, 0])
    print(f"\nPanel (a): OBC lattice L={L} with two 3×3 holes")
    mask_holes = paper_defect_two_holes(L)
    mask_holes.visualize(ax=ax_a, show_types=True)
    ax_a.set_title(f'(a) OBC with two 3×3 holes, $L={L}$\n({mask_holes.n_sites} sites, {mask_holes.n_states} states)',
                   fontsize=12, fontweight='bold')

    # ========== Panel (b): Spectrum colored by EXTERNAL boundary weight ==========
    ax_b = fig.add_subplot(gs[0, 1])
    print(f"\nPanel (b): Spectrum (external boundary weight)")
    print(f"  Building transition matrix (omega={omega}, D_r={D_r})...")
    t0 = time.time()
    P = build_transition_matrix_generic(omega, D_r, mask_holes)
    P_dense = P.toarray()
    print(f"  Matrix built in {time.time()-t0:.2f}s, shape {P_dense.shape}")

    print(f"  Computing eigenvalues and external boundary weights...")
    t0 = time.time()
    evals, weights_ext = compute_edge_weights(P_dense, mask_holes, 'external')
    print(f"  Done in {time.time()-t0:.2f}s")

    # Scatter colored by external boundary weight
    norm = Normalize(vmin=0, vmax=1)
    cmap = cm.get_cmap('RdYlBu_r')
    scatter_b = ax_b.scatter(evals.real, evals.imag, c=weights_ext,
                            cmap=cmap, norm=norm, s=40, alpha=0.7, edgecolors='k', linewidths=0.5)
    cbar_b = plt.colorbar(scatter_b, ax=ax_b)
    cbar_b.set_label('Boundary weight (external)', fontsize=10)

    ax_b.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_b.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_b.set_xlabel('Re($\lambda$)', fontsize=11)
    ax_b.set_ylabel('Im($\lambda$)', fontsize=11)
    ax_b.set_title(f'(b) Spectrum colored by external boundary weight\n($\omega={omega}$, $D_r={D_r}$)',
                   fontsize=12, fontweight='bold')
    ax_b.grid(True, alpha=0.3)

    # ========== Panel (c): Lattice with internal boundaries highlighted ==========
    ax_c = fig.add_subplot(gs[1, 0])
    print(f"\nPanel (c): Internal boundary sites highlighted")

    # Draw the lattice with custom coloring
    colors_custom = {'bulk': '#d4e6f1', 'external': '#f5b041',
                     'internal': '#e74c3c', 'invalid': '#2c3e50'}

    for x in range(L):
        for y in range(L):
            if mask_holes.is_valid(x, y):
                bt = mask_holes.boundary_type(x, y)
                rect = mpatches.Rectangle((x - 0.5, y - 0.5), 1, 1,
                                         facecolor=colors_custom[bt],
                                         edgecolor='gray', lw=0.5)
                ax_c.add_patch(rect)
            else:
                rect = mpatches.Rectangle((x - 0.5, y - 0.5), 1, 1,
                                         facecolor=colors_custom['invalid'],
                                         edgecolor='gray', lw=0.5, alpha=0.3)
                ax_c.add_patch(rect)

    ax_c.set_xlim(-0.6, L - 0.4)
    ax_c.set_ylim(-0.6, L - 0.4)
    ax_c.set_aspect('equal')
    ax_c.set_xlabel('x', fontsize=11)
    ax_c.set_ylabel('y', fontsize=11)
    ax_c.set_title(f'(c) Internal boundary sites (defect interfaces)',
                   fontsize=12, fontweight='bold')

    legend_elements = [
        mpatches.Patch(facecolor=colors_custom['bulk'], label='Bulk'),
        mpatches.Patch(facecolor=colors_custom['external'], label='External boundary'),
        mpatches.Patch(facecolor=colors_custom['internal'], label='Internal boundary'),
        mpatches.Patch(facecolor=colors_custom['invalid'], alpha=0.3, label='Blocked/removed'),
    ]
    ax_c.legend(handles=legend_elements, loc='upper right', fontsize=8)

    # ========== Panel (d): Spectrum colored by INTERNAL boundary weight ==========
    ax_d = fig.add_subplot(gs[1, 1])
    print(f"\nPanel (d): Spectrum (internal boundary weight - nested)")
    print(f"  Computing internal boundary weights...")
    t0 = time.time()
    _, weights_int = compute_edge_weights(P_dense, mask_holes, 'internal')
    print(f"  Done in {time.time()-t0:.2f}s")

    # Scatter colored by internal boundary weight
    scatter_d = ax_d.scatter(evals.real, evals.imag, c=weights_int,
                            cmap=cmap, norm=norm, s=40, alpha=0.7, edgecolors='k', linewidths=0.5)
    cbar_d = plt.colorbar(scatter_d, ax=ax_d)
    cbar_d.set_label('Boundary weight (internal)', fontsize=10)

    ax_d.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_d.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_d.set_xlabel('Re($\lambda$)', fontsize=11)
    ax_d.set_ylabel('Im($\lambda$)', fontsize=11)
    ax_d.set_title(f'(d) NESTED spectrum colored by internal boundary weight\n(different defects show distinct oval signatures)',
                   fontsize=12, fontweight='bold')
    ax_d.grid(True, alpha=0.3)

    plt.suptitle('Figure 12: Disconnected internal boundaries and nested spectral signatures',
                 fontsize=14, fontweight='bold', y=0.995)

    print(f"\nSaving figure 12...")
    plt.savefig('tcrw_fig12_nested_spectrum.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig12_nested_spectrum.png")


# ============================================================
# Figures 11 & 12 — paper-matched D_r sweep versions
# ============================================================

def fig11_dr_sweep(D_r_values=(0.5, 0.4, 0.25, 0.2, 0.15), omega=1.0, L=10,
                   outfile='tcrw_fig11_boundary_Dr_sweep.png'):
    """
    Fig 11(b), (d) with paper's D_r sweep ω=1, D_r ∈ {0.5, 0.4, 0.25, 0.2, 0.15}.
    Row 1: OBC spectra colored by external-edge weight.
    Row 2: PBC+defect spectra colored by internal-boundary weight.
    """
    print("\n" + "=" * 70)
    print("FIGURE 11 (D_r sweep): paper-matched params")
    print("=" * 70)

    mask_obc = RectangleMask(L)
    mask_pbc = paper_pbc_with_internal_defect(L)

    n = len(D_r_values)
    fig, axes = plt.subplots(2, n, figsize=(3.3 * n, 7))
    cmap = cm.get_cmap('RdYlBu_r')
    norm = Normalize(vmin=0, vmax=1)

    last_sc = None
    for col, D_r in enumerate(D_r_values):
        print(f"  [{col+1}/{n}] D_r={D_r}")
        P_obc = build_transition_matrix_generic(omega, D_r, mask_obc).toarray()
        ev_o, w_o = compute_edge_weights(P_obc, mask_obc, 'external')
        P_pbc = build_transition_matrix_generic(omega, D_r, mask_pbc).toarray()
        ev_p, w_p = compute_edge_weights(P_pbc, mask_pbc, 'internal')

        for row, (ev, w, title) in enumerate([(ev_o, w_o, 'OBC ext'),
                                                (ev_p, w_p, 'PBC+def int')]):
            ax = axes[row, col]
            last_sc = ax.scatter(ev.real, ev.imag, c=w, cmap=cmap, norm=norm,
                                  s=22, alpha=0.75, edgecolors='k', linewidths=0.3)
            ax.axhline(0, color='gray', lw=0.3, alpha=0.4)
            ax.axvline(0, color='gray', lw=0.3, alpha=0.4)
            ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1); ax.set_aspect('equal')
            if row == 0:
                ax.set_title(f'$D_r={D_r}$', fontsize=11)
            if row == 1:
                ax.set_xlabel(r'$\mathrm{Re}(\lambda)$', fontsize=10)
            if col == 0:
                ax.set_ylabel(r'$\mathrm{Im}(\lambda)$' + f'\n{title}', fontsize=10)

    fig.colorbar(last_sc, ax=axes.ravel().tolist(), shrink=0.6, pad=0.02,
                  label='Boundary weight')
    plt.suptitle(f'Fig 11 (paper params): $\\omega={omega}$, $L={L}$, $D_r$ sweep',
                  fontsize=13, y=0.995)
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outfile}")


def fig12_dr_sweep(D_r_values=(0.5, 0.4, 0.25, 0.2, 0.15), omega=1.0, L=15,
                   outfile='tcrw_fig12_nested_Dr_sweep.png'):
    """
    Fig 12(b), (d) with paper's D_r sweep ω=1, D_r ∈ {0.5, 0.4, 0.25, 0.2, 0.15}.
    Row 1: OBC+holes colored by external boundary.
    Row 2: same geometry colored by internal boundary -> nested ovals.
    """
    print("\n" + "=" * 70)
    print("FIGURE 12 (D_r sweep): paper-matched params")
    print("=" * 70)

    mask = RectangleWithHoles(L, [(3, 3, 3, 3), (L-6, L-6, 3, 3)])

    n = len(D_r_values)
    fig, axes = plt.subplots(2, n, figsize=(3.3 * n, 7))
    cmap = cm.get_cmap('RdYlBu_r')
    norm = Normalize(vmin=0, vmax=1)

    last_sc = None
    for col, D_r in enumerate(D_r_values):
        print(f"  [{col+1}/{n}] D_r={D_r}")
        P = build_transition_matrix_generic(omega, D_r, mask).toarray()
        ev_e, w_e = compute_edge_weights(P, mask, 'external')
        ev_i, w_i = compute_edge_weights(P, mask, 'internal')

        for row, (ev, w, title) in enumerate([(ev_e, w_e, 'ext boundary'),
                                                (ev_i, w_i, 'int boundary')]):
            ax = axes[row, col]
            last_sc = ax.scatter(ev.real, ev.imag, c=w, cmap=cmap, norm=norm,
                                  s=22, alpha=0.75, edgecolors='k', linewidths=0.3)
            ax.axhline(0, color='gray', lw=0.3, alpha=0.4)
            ax.axvline(0, color='gray', lw=0.3, alpha=0.4)
            ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1); ax.set_aspect('equal')
            if row == 0:
                ax.set_title(f'$D_r={D_r}$', fontsize=11)
            if row == 1:
                ax.set_xlabel(r'$\mathrm{Re}(\lambda)$', fontsize=10)
            if col == 0:
                ax.set_ylabel(r'$\mathrm{Im}(\lambda)$' + f'\n{title}', fontsize=10)

    fig.colorbar(last_sc, ax=axes.ravel().tolist(), shrink=0.6, pad=0.02,
                  label='Boundary weight')
    plt.suptitle(f'Fig 12 (paper params): $\\omega={omega}$, $L={L}$, $D_r$ sweep',
                  fontsize=13, y=0.995)
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outfile}")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("\n" + "=" * 70)
    print("TCRW FIGURES 11 & 12: Boundary effects and nested spectra")
    print("=" * 70)

    fig11_boundary_effects()
    fig12_nested_spectrum()
    # Paper-matched sweep versions (new):
    fig11_dr_sweep()
    fig12_dr_sweep()

    print("\n" + "=" * 70)
    print("All figures completed successfully!")
    print("=" * 70)
    print("\nGenerated files:")
    print("  - tcrw_fig11_boundary_effects.png")
    print("  - tcrw_fig12_nested_spectrum.png")
