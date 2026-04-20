"""
TCRW Fig 4(h): Spectrum with hybrid periodic-open boundary conditions
=====================================================================

Hybrid boundary conditions:
  x-direction: OPEN boundary conditions (hard walls at x=0 and x=L-1)
  y-direction: PERIODIC boundary conditions (wraps around)

Because y is periodic, we can Fourier transform in y → momentum k_y.
The result is a (4L)×(4L) matrix for each k_y value (L sites in x, 4 director states).

State indexing: state_idx = d*L + x, where x ∈ {0,...,L-1}, d ∈ {0,1,2,3}

Matrix construction for each k_y:
  Noise step (prob D_r): stay at x, rotate d.
    - CCW: (x, (d-1)%4) with prob omega*D_r
    - CW:  (x, (d+1)%4) with prob (1-omega)*D_r

  Chiral step (prob 1-D_r): translate in direction d, then rotate opposite.
    - d=0 (↑): y translation → multiply by e^{iky}
    - d=2 (↓): y translation → multiply by e^{-iky}
    - d=1 (→): x translation to x+1 (blocked if x+1 >= L)
    - d=3 (←): x translation to x-1 (blocked if x-1 < 0)

Direction encoding:
  d=0: ↑ (north, y+1)  DY[0] = 1
  d=1: → (east,  x+1)  DX[1] = 1
  d=2: ↓ (south, y-1)  DY[2] = -1
  d=3: ← (west,  x-1)  DX[3] = -1

Results in Re(λ) vs k_y spectrum, colored by edge localization weight.

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def build_hpbc_matrix(omega, D_r, L, ky):
    """
    Build the (4L)×(4L) transition matrix M(k_y) for hybrid PBC/OBC.

    Parameters
    ----------
    omega : float
        Chirality parameter (0 to 1).
    D_r : float
        Noise probability (0 to 1).
    L : int
        System size in x-direction (open BCs).
    ky : float
        Momentum in y-direction (periodic BC).

    Returns
    -------
    M : complex (4L)×(4L) numpy array
        M[j, i] = amplitude for state i → state j at wavevector k_y.

    States indexed by: state_idx = d*L + x
      where x ∈ {0,...,L-1}, d ∈ {0,1,2,3}
    """
    N = 4 * L
    M = np.zeros((N, N), dtype=complex)

    # Direction vectors: d=0(↑), 1(→), 2(↓), 3(←)
    DX = np.array([0, 1, 0, -1])
    DY = np.array([1, 0, -1, 0])

    for x in range(L):
        for d in range(4):
            i = d * L + x  # Current state index

            # ==========================================
            # Noise step: prob D_r, stay at x, rotate d
            # ==========================================

            # CCW: d → (d-1) mod 4, prob omega*D_r
            d_ccw = (d - 1) % 4
            j_ccw = d_ccw * L + x
            M[j_ccw, i] += omega * D_r

            # CW: d → (d+1) mod 4, prob (1-omega)*D_r
            d_cw = (d + 1) % 4
            j_cw = d_cw * L + x
            M[j_cw, i] += (1 - omega) * D_r

            # ==========================================
            # Chiral step: prob (1-D_r)
            # Translate in direction d, then rotate opposite
            # ==========================================
            nx = x + DX[d]
            ny_disp = DY[d]  # y-displacement (±1 for vertical, 0 for horizontal)

            if 0 <= nx < L:
                # Valid x-move: translate in (x, d), rotate opposite
                # Phase from y-translation: e^{i*k_y*ny_disp}
                phase = np.exp(1j * ky * ny_disp)

                # After translation, rotate with OPPOSITE chirality:
                # CW: d → (d+1) mod 4, prob omega*(1-D_r)
                d_cw_c = (d + 1) % 4
                j1 = d_cw_c * L + nx
                M[j1, i] += omega * (1 - D_r) * phase

                # CCW: d → (d-1) mod 4, prob (1-omega)*(1-D_r)
                d_ccw_c = (d - 1) % 4
                j2 = d_ccw_c * L + nx
                M[j2, i] += (1 - omega) * (1 - D_r) * phase
            else:
                # Blocked at boundary: stay at (x, d), no rotation
                j_stay = d * L + x
                M[j_stay, i] += (1 - D_r)

    return M


def compute_hpbc_spectrum(omega, D_r, L, ky):
    """
    Diagonalize M(k_y) and compute edge localization weights.

    Parameters
    ----------
    omega, D_r, L, ky : as in build_hpbc_matrix

    Returns
    -------
    eigenvalues : array of shape (4L,), complex
    edge_weights : array of shape (4L,), float
        Fraction of |eigenvector|^2 on boundary sites (x=0 or x=L-1).
    """
    M = build_hpbc_matrix(omega, D_r, L, ky)
    eigenvalues, eigenvectors = np.linalg.eig(M)

    # Edge localization: x=0 or x=L-1 (all directors)
    edge_indices = set()
    for d in range(4):
        edge_indices.add(d * L + 0)      # x=0
        edge_indices.add(d * L + (L - 1))  # x=L-1

    edge_weights = np.zeros(len(eigenvalues))
    for n in range(len(eigenvalues)):
        v = eigenvectors[:, n]
        weight = np.abs(v) ** 2
        total = weight.sum()
        if total > 1e-14:
            edge_weights[n] = sum(weight[idx] for idx in edge_indices) / total
        else:
            edge_weights[n] = 0.0

    return eigenvalues, edge_weights


def fig4h_hpbc_spectrum():
    """
    Fig 4(h): Spectrum with hybrid PBC/OBC.

    Sweep k_y from -π to π, plot Re(λ) vs k_y.
    Color by edge localization weight.
    """
    # Parameters
    omega = 1.0
    D_r = 0.1
    L = 10
    Nky = 200

    print(f"Computing HPBC spectrum: omega={omega}, D_r={D_r}, L={L}")
    print(f"Sweeping k_y from -π to π with {Nky} points...")

    # Sweep in k_y
    ky_vals = np.linspace(-np.pi, np.pi, Nky, endpoint=False)
    re_lambda_vals = []
    edge_weight_vals = []

    for ky in ky_vals:
        evals, ew = compute_hpbc_spectrum(omega, D_r, L, ky)
        re_lambda_vals.append(evals.real)
        edge_weight_vals.append(ew)

    # Reshape for plotting
    re_lambda_array = np.array(re_lambda_vals)  # shape (Nky, 4L)
    edge_weight_array = np.array(edge_weight_vals)  # shape (Nky, 4L)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))

    # Flatten k_y for scatter plot
    ky_plot = np.repeat(ky_vals, 4 * L)
    re_lambda_plot = re_lambda_array.ravel()
    edge_weight_plot = edge_weight_array.ravel()

    # Scatter plot: Re(λ) vs k_y, colored by edge weight
    norm = Normalize(vmin=0, vmax=1)
    cmap = plt.colormaps['hot_r']
    sc = ax.scatter(ky_plot, re_lambda_plot, c=edge_weight_plot, cmap=cmap,
                    norm=norm, s=12, alpha=0.6, edgecolors='none')

    ax.set_xlabel(r'$k_y$', fontsize=14, fontweight='bold')
    ax.set_ylabel(r'Re($\lambda$)', fontsize=14, fontweight='bold')
    ax.set_title(f'Fig 4(h): HPBC Spectrum ($\\omega={omega}$, $D_r={D_r}$, $L={L}$)',
                 fontsize=13, fontweight='bold')

    # k_y ticks
    ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
    ax.set_xticklabels([r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'],
                       fontsize=12)

    ax.axhline(0, color='gray', lw=0.8, ls='--', alpha=0.5)
    ax.axvline(0, color='gray', lw=0.8, ls='--', alpha=0.5)
    ax.grid(True, alpha=0.2)

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, label='Edge localization weight')
    cbar.set_label('Edge localization weight', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig('tcrw_fig4h_hpbc_spectrum.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig4h_hpbc_spectrum.png")


if __name__ == '__main__':
    print("=" * 70)
    print("TCRW Fig 4(h): Spectrum with Hybrid Periodic-Open Boundary Conditions")
    print("=" * 70)

    fig4h_hpbc_spectrum()

    print("\n" + "=" * 70)
    print("Fig 4(h) complete.")
    print("=" * 70)
