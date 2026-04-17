"""
TCRW Phase 4: Transition matrix spectrum — the heart of the topology
=====================================================================

Two regimes:

4A — PBC (Fourier space):
  The infinite-lattice transition matrix block-diagonalises in momentum.
  At each k = (k_x, k_y), it reduces to a 4×4 matrix P(k) acting on the
  four director states d ∈ {↑, →, ↓, ←}.

  P(k)[j,i] = prob of i → j, with spatial translation encoded as phases.

  Parameters:
    C1 = (1−ω)(1−D_r)   chiral, CCW rotation after translation
    C2 = ω(1−D_r)        chiral, CW  rotation after translation
    R1 = ω D_r            noise,  CCW rotation (no translation)
    R2 = (1−ω) D_r        noise,  CW  rotation (no translation)

  P(k) =
    | 0                R1+C1·e^{ikx}   0                R2+C2·e^{-ikx} |
    | R2+C2·e^{iky}   0                R1+C1·e^{-iky}   0              |
    | 0                R2+C2·e^{ikx}   0                R1+C1·e^{-ikx} |
    | R1+C1·e^{iky}   0                R2+C2·e^{-iky}   0              |

  Key properties:
    - Bipartite structure ⇒ eigenvalues come in ±λ pairs (from P² block-diag)
    - One pair is always purely imaginary (±ib): Re = 0 to machine precision
    - The other pair can be real (±a) or complex, depending on parameters
    - Gap ≡ min_k |Re(λ)| of the non-imaginary pair
    - Gap closes exactly at ω = 0.5 (topological phase transition)

4B — OBC (real space):
  Full (4L²)×(4L²) transition matrix. Eigenvalues become complex.
  Imaginary parts encode edge currents. Bulk-boundary correspondence:
  edge modes appear inside the bulk gap when D_r < 0.5.

Reproduces:
  Fig 4(b): Re(λ) band structure along high-symmetry path
  Fig 4(c): OBC spectrum in complex plane (small L)
  Fig 4(f)-(g): OBC spectrum for L=10, colored by edge-localization weight

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tcrw_obc import build_transition_matrix, state_index


# ============================================================
# 4A: Fourier-space transition matrix
# ============================================================

def build_Pk(omega, D_r, kx, ky):
    """
    Build the 4×4 transition matrix P(k) in Fourier space.

    Rows/columns indexed by director: 0=↑, 1=→, 2=↓, 3=←.
    P(k)[j,i] = amplitude for state i → state j at wavevector (kx, ky).

    Returns complex 4×4 numpy array.
    """
    C1 = (1 - omega) * (1 - D_r)   # chiral CCW
    C2 = omega * (1 - D_r)          # chiral CW
    R1 = omega * D_r                 # noise CCW
    R2 = (1 - omega) * D_r          # noise CW

    ekx = np.exp(1j * kx)
    eky = np.exp(1j * ky)

    P = np.array([
        [0,                R1 + C1 * ekx,    0,                R2 + C2 / ekx],
        [R2 + C2 * eky,    0,                R1 + C1 / eky,    0],
        [0,                R2 + C2 * ekx,    0,                R1 + C1 / ekx],
        [R1 + C1 * eky,    0,                R2 + C2 / eky,    0],
    ], dtype=complex)

    return P


def pbc_eigenvalues(omega, D_r, kx, ky):
    """
    Eigenvalues of P(k) at a single k-point.
    Returns sorted array of 4 eigenvalues (complex in general).
    """
    P = build_Pk(omega, D_r, kx, ky)
    evals = np.linalg.eigvals(P)
    # Sort by real part descending
    idx = np.argsort(-evals.real)
    return evals[idx]


def pbc_band_structure(omega, D_r, Nk=200):
    """
    Compute 4-band structure along high-symmetry path:
      Γ(0,0) → X(π,0) → M(π,π) → Γ(0,0)

    Returns:
      k_path: 1D array of cumulative k-distance
      bands:  (4, len(k_path)) array of eigenvalues
      ticks:  list of (position, label) for high-symmetry points
    """
    # High-symmetry points
    Gamma = np.array([0.0, 0.0])
    X = np.array([np.pi, 0.0])
    M = np.array([np.pi, np.pi])

    segments = [(Gamma, X, '$\\Gamma$', '$X$'),
                (X, M, '$X$', '$M$'),
                (M, Gamma, '$M$', '$\\Gamma$')]

    k_path = []
    bands = []
    tick_positions = [0.0]
    tick_labels = ['$\\Gamma$']
    cumulative = 0.0

    for start, end, label_s, label_e in segments:
        dk = end - start
        seg_length = np.linalg.norm(dk)
        t_vals = np.linspace(0, 1, Nk, endpoint=False)

        for t in t_vals:
            kpt = start + t * dk
            evals = pbc_eigenvalues(omega, D_r, kpt[0], kpt[1])
            bands.append(evals)
            k_path.append(cumulative + t * seg_length)

        cumulative += seg_length
        tick_positions.append(cumulative)
        tick_labels.append(label_e)

    k_path = np.array(k_path)
    bands = np.array(bands).T  # shape (4, N_total)

    # Sort bands by continuity (tracking by nearest neighbor)
    bands_sorted = np.empty_like(bands)
    bands_sorted[:, 0] = bands[:, 0]
    for i in range(1, bands.shape[1]):
        # Hungarian-style: assign each new eigenvalue to nearest previous
        prev = bands_sorted[:, i - 1]
        curr = bands[:, i]
        used = [False] * 4
        assignment = [0] * 4
        for b in range(4):
            dists = np.abs(curr - prev[b])
            order = np.argsort(dists)
            for j in order:
                if not used[j]:
                    assignment[b] = j
                    used[j] = True
                    break
        for b in range(4):
            bands_sorted[b, i] = curr[assignment[b]]

    ticks = list(zip(tick_positions, tick_labels))
    return k_path, bands_sorted, ticks


def pbc_full_bz(omega, D_r, Nk=50):
    """
    All eigenvalues on a Nk×Nk grid in the full BZ [-π,π)².
    Returns (4, Nk, Nk) array of eigenvalues.
    """
    kx_arr = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    ky_arr = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    evals = np.zeros((4, Nk, Nk), dtype=complex)

    for i, kx in enumerate(kx_arr):
        for j, ky in enumerate(ky_arr):
            evals[:, i, j] = pbc_eigenvalues(omega, D_r, kx, ky)

    return kx_arr, ky_arr, evals


# ============================================================
# 4B: OBC spectrum
# ============================================================

def obc_spectrum(omega, D_r, L):
    """
    Full eigenvalue spectrum of the (4L²)×(4L²) OBC transition matrix.

    Returns:
      eigenvalues: array of all eigenvalues
      edge_weights: for each eigenvector, fraction of |ψ|² on edge sites
    """
    P = build_transition_matrix(omega, D_r, L)
    P_dense = P.toarray()

    eigenvalues, eigenvectors = np.linalg.eig(P_dense)

    # Edge weight: fraction of eigenvector weight on edge sites
    edge_mask = np.zeros(4 * L * L, dtype=bool)
    for d in range(4):
        offset = d * L * L
        for x in range(L):
            for y in range(L):
                if x == 0 or x == L - 1 or y == 0 or y == L - 1:
                    edge_mask[offset + y * L + x] = True

    edge_weights = np.zeros(len(eigenvalues))
    for n in range(len(eigenvalues)):
        v = eigenvectors[:, n]
        weight = np.abs(v) ** 2
        total = weight.sum()
        if total > 0:
            edge_weights[n] = weight[edge_mask].sum() / total

    return eigenvalues, edge_weights


# ============================================================
# Figures
# ============================================================

def classify_eigenvalues(evals):
    """
    Classify 4 eigenvalues into real pair (±a) and imaginary pair (±ib).

    Bipartite structure guarantees eigenvalues come in ±λ pairs.
    Empirically, one pair is real and the other purely imaginary.

    Returns:
      a: positive real eigenvalue magnitude (a >= 0)
      b: positive imaginary eigenvalue magnitude (b >= 0)
    """
    # Sort by |Re| descending — real pair has largest |Re|
    idx = np.argsort(-np.abs(evals.real))
    real_pair = evals[idx[:2]]
    imag_pair = evals[idx[2:]]
    a = np.max(np.abs(real_pair.real))
    b = np.max(np.abs(imag_pair.imag))
    return a, b


def pbc_band_structure_classified(omega, D_r, Nk=200):
    """
    Compute band structure along Γ→X→M→Γ, classifying into
    real pair ±a(k) and imaginary pair ±ib(k).

    Returns:
      k_path, a_vals, b_vals, ticks
    """
    Gamma = np.array([0.0, 0.0])
    X = np.array([np.pi, 0.0])
    M = np.array([np.pi, np.pi])

    segments = [(Gamma, X, '$\\Gamma$', '$X$'),
                (X, M, '$X$', '$M$'),
                (M, Gamma, '$M$', '$\\Gamma$')]

    k_path = []
    a_vals = []
    b_vals = []
    tick_positions = [0.0]
    tick_labels = ['$\\Gamma$']
    cumulative = 0.0

    for start, end, label_s, label_e in segments:
        dk = end - start
        seg_length = np.linalg.norm(dk)
        t_vals = np.linspace(0, 1, Nk, endpoint=False)

        for t in t_vals:
            kpt = start + t * dk
            evals = pbc_eigenvalues(omega, D_r, kpt[0], kpt[1])
            a, b = classify_eigenvalues(evals)
            a_vals.append(a)
            b_vals.append(b)
            k_path.append(cumulative + t * seg_length)

        cumulative += seg_length
        tick_positions.append(cumulative)
        tick_labels.append(label_e)

    return (np.array(k_path), np.array(a_vals), np.array(b_vals),
            list(zip(tick_positions, tick_labels)))


def fig4b_band_structure():
    """
    Fig 4(b): Band structure for ω = 0.35, 0.5, 0.65.

    Eigenvalues come in two pairs: ±a(k) (real) and ±ib(k) (imaginary).
    Top row: Re(λ) — shows ±a(k) dispersive, 0 flat line for imaginary pair.
    Bottom row: full complex plane loci as k sweeps the BZ.
    Gap = min_k a(k); closes at ω = 0.5.
    """
    D_r = 0.1
    omegas = [0.35, 0.5, 0.65]
    titles = ['$\\omega = 0.35$ (topological)', '$\\omega = 0.5$ (critical)',
              '$\\omega = 0.65$ (trivial)']

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    for col, (omega, title) in enumerate(zip(omegas, titles)):
        # --- Top row: Re(λ) band structure ---
        k_path, a_vals, b_vals, ticks = pbc_band_structure_classified(
            omega, D_r, Nk=300)

        ax = axes[0, col]
        ax.fill_between(k_path, a_vals, -a_vals, alpha=0.15, color='#1f77b4')
        ax.plot(k_path, a_vals, color='#1f77b4', lw=1.8, label='$+a(k)$')
        ax.plot(k_path, -a_vals, color='#1f77b4', lw=1.8, label='$-a(k)$')
        ax.axhline(0, color='gray', lw=0.8, ls='-', alpha=0.4)

        # High-symmetry ticks
        for pos, label in ticks:
            ax.axvline(x=pos, color='gray', ls='--', lw=0.5, alpha=0.5)
        ax.set_xticks([pos for pos, _ in ticks])
        ax.set_xticklabels([label for _, label in ticks], fontsize=11)
        ax.set_ylabel('Re($\\lambda$)', fontsize=12)
        ax.set_title(title, fontsize=12)
        ax.set_ylim(-1.1, 1.1)

        # Annotate gap
        gap = np.min(a_vals)
        ax.axhline(y=gap, color='red', ls=':', lw=1, alpha=0.6)
        ax.axhline(y=-gap, color='red', ls=':', lw=1, alpha=0.6)
        ax.text(0.05, 0.05, f'gap = {gap:.4f}', transform=ax.transAxes,
                fontsize=10, fontweight='bold', color='red')

        # --- Bottom row: complex plane loci ---
        ax = axes[1, col]
        _, _, pbc_evals = pbc_full_bz(omega, D_r, Nk=50)
        all_evals = pbc_evals.ravel()
        ax.scatter(all_evals.real, all_evals.imag, s=1, alpha=0.3,
                   color='#1f77b4', edgecolors='none')
        ax.set_xlabel('Re($\\lambda$)', fontsize=11)
        ax.set_ylabel('Im($\\lambda$)', fontsize=11)
        ax.set_title(f'Complex plane ($\\omega={omega}$)', fontsize=11)
        ax.set_aspect('equal')
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.15, 1.15)
        ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
        ax.axvline(0, color='gray', lw=0.5, alpha=0.3)

        # Draw unit circle
        theta = np.linspace(0, 2 * np.pi, 200)
        ax.plot(np.cos(theta), np.sin(theta), 'k--', lw=0.5, alpha=0.2)

    axes[0, 0].legend(fontsize=9, loc='upper right')
    plt.suptitle(f'Fig 4(b): PBC spectrum ($D_r = {D_r}$)',
                 fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig('tcrw_fig4b_pbc_bands.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig4b_pbc_bands.png")


def fig4_obc_complex_plane():
    """
    Fig 4(f)/(g): OBC spectrum in the complex plane.
    Color by edge-localization weight of each eigenvector.

    Panel layout (matches paper Fig 4(f) and 4(g)):
      Row 1 = Fig 4(g): ω ∈ {0.5, 0.7, 0.9} at fixed D_r = 0.1
      Row 2 = Fig 4(f): ω = 1.0 at D_r ∈ {0.65, 0.5, 0.35}
    """
    L = 10

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # Row 1 = Fig 4(g): varying omega at fixed D_r
    D_r = 0.1
    omegas = [0.5, 0.7, 0.9]
    for col, omega in enumerate(omegas):
        print(f"  OBC spectrum: omega={omega}, D_r={D_r}, L={L}...")
        evals, ew = obc_spectrum(omega, D_r, L)
        ax = axes[0, col]
        sc = ax.scatter(evals.real, evals.imag, c=ew, cmap='hot_r',
                        s=8, alpha=0.7, vmin=0, vmax=1,
                        edgecolors='none')
        ax.set_xlabel('Re($\\lambda$)', fontsize=11)
        ax.set_ylabel('Im($\\lambda$)', fontsize=11)
        ax.set_title(f'$\\omega={omega}$, $D_r={D_r}$', fontsize=11)
        ax.set_aspect('equal')
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-0.6, 0.6)
        ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
        ax.axvline(0, color='gray', lw=0.5, alpha=0.3)

    plt.colorbar(sc, ax=axes[0, :].tolist(), shrink=0.6, label='Edge weight',
                 pad=0.02)

    # Row 2 = Fig 4(f): varying D_r at fixed omega = 1 (paper's values)
    omega = 1.0
    D_r_values = [0.65, 0.5, 0.35]
    for col, D_r in enumerate(D_r_values):
        print(f"  OBC spectrum: omega={omega}, D_r={D_r}, L={L}...")
        evals, ew = obc_spectrum(omega, D_r, L)
        ax = axes[1, col]
        sc2 = ax.scatter(evals.real, evals.imag, c=ew, cmap='hot_r',
                         s=8, alpha=0.7, vmin=0, vmax=1,
                         edgecolors='none')
        ax.set_xlabel('Re($\\lambda$)', fontsize=11)
        ax.set_ylabel('Im($\\lambda$)', fontsize=11)
        ax.set_title(f'$\\omega={omega}$, $D_r={D_r}$', fontsize=11)
        ax.set_aspect('equal')
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-0.6, 0.6)
        ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
        ax.axvline(0, color='gray', lw=0.5, alpha=0.3)

    plt.colorbar(sc2, ax=axes[1, :].tolist(), shrink=0.6, label='Edge weight',
                 pad=0.02)

    plt.suptitle(f'Fig 4: OBC eigenvalue spectrum ($L={L}$)',
                 fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig('tcrw_fig4_obc_spectrum.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig4_obc_spectrum.png")


def fig4_gap_closing():
    """
    Spectral gap vs ω at several D_r values.

    Gap = min_k a(k), where ±a(k) is the real eigenvalue pair.
    Gap closes (a→0) at ω = 0.5 for D_r < 0.5.
    This is the topological phase transition.
    """
    omega_scan = np.linspace(0.0, 1.0, 101)
    D_r_values = [0.05, 0.1, 0.2, 0.3, 0.5]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(D_r_values)))

    fig, ax = plt.subplots(figsize=(8, 5))

    for D_r, color in zip(D_r_values, colors):
        gaps = []
        for omega in omega_scan:
            min_a = np.inf
            for kx in np.linspace(-np.pi, np.pi, 50, endpoint=False):
                for ky in np.linspace(-np.pi, np.pi, 50, endpoint=False):
                    evals = pbc_eigenvalues(omega, D_r, kx, ky)
                    a, b = classify_eigenvalues(evals)
                    if a < min_a:
                        min_a = a
            gaps.append(min_a)

        ax.plot(omega_scan, gaps, lw=1.8, color=color, label=f'$D_r={D_r}$')

    ax.axvline(x=0.5, color='red', ls=':', lw=1.5, alpha=0.7, label='$\\omega=0.5$')
    ax.set_xlabel('$\\omega$', fontsize=13)
    ax.set_ylabel('Spectral gap $\\min_k a(k)$', fontsize=13)
    ax.set_title('Gap closing at the topological phase transition', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_ylim(bottom=-0.01)
    plt.tight_layout()
    plt.savefig('tcrw_fig4_gap_closing.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig4_gap_closing.png")


def fig4_pbc_vs_obc():
    """
    PBC bulk eigenvalues overlaid with OBC spectrum.
    Shows that OBC has extra complex eigenvalues (edge modes)
    inside the PBC bulk gap.
    """
    L = 10
    D_r = 0.1
    omegas = [0.35, 0.65]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    for ax, omega in zip(axes, omegas):
        # PBC: sample full BZ
        _, _, pbc_evals = pbc_full_bz(omega, D_r, Nk=40)
        pbc_flat = pbc_evals.ravel()

        # OBC
        obc_evals, ew = obc_spectrum(omega, D_r, L)

        # Plot PBC as gray background
        ax.scatter(pbc_flat.real, pbc_flat.imag, s=2, color='gray',
                   alpha=0.3, label='PBC bulk', zorder=1)

        # OBC colored by edge weight
        sc = ax.scatter(obc_evals.real, obc_evals.imag, c=ew, cmap='hot_r',
                        s=12, alpha=0.8, vmin=0, vmax=1, edgecolors='none',
                        zorder=2)

        ax.set_xlabel('Re($\\lambda$)', fontsize=12)
        ax.set_ylabel('Im($\\lambda$)', fontsize=12)
        phase = 'topological' if omega < 0.5 else 'trivial'
        ax.set_title(f'$\\omega={omega}$ ({phase}), $D_r={D_r}$', fontsize=12)
        ax.set_aspect('equal')
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-0.6, 0.6)
        ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
        ax.legend(fontsize=9, loc='upper left')

    plt.colorbar(sc, ax=axes.tolist(), shrink=0.7, label='Edge weight', pad=0.02)
    plt.suptitle(f'PBC bulk vs OBC spectrum ($L={L}$, $D_r={D_r}$)',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('tcrw_fig4_pbc_vs_obc.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig4_pbc_vs_obc.png")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("=" * 60)
    print("TCRW Phase 4: Transition Matrix Spectrum")
    print("=" * 60)

    print("\n--- Fig 4(b): PBC band structure ---")
    fig4b_band_structure()

    print("\n--- Fig 4(c)/(f)/(g): OBC spectrum in complex plane ---")
    fig4_obc_complex_plane()

    print("\n--- Gap closing diagnostic ---")
    fig4_gap_closing()

    print("\n--- PBC vs OBC overlay ---")
    fig4_pbc_vs_obc()

    print("\n" + "=" * 60)
    print("Phase 4 complete.")
    print("=" * 60)
