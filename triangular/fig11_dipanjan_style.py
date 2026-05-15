"""
Reproduce the visual style of Dipanjan's draft Fig. 11.

The draft shows a near-continuous heatmap where hexagonal Voronoi cells
tile seamlessly — no visible gaps, no black background peeking through.

We achieve this by using matplotlib RegularPolygon patches sized in data
coordinates (radius = 1/√3 of the NN distance) so they tile perfectly.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import RegularPolyCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize

from triangular_jmvr_corrected import (
    L_DEFAULT,
    T_DEFAULT,
    build_Mk_corrected,
    load_kmc_counts,
    theory_probability_on_lattice,
)

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)

L = L_DEFAULT
t_final = T_DEFAULT
gamma = 0.01
epsilon = 0.15


# ---------------------------------------------------------------------------
# Draft-matching colormap: dark purple → magenta → red → orange → yellow
# ---------------------------------------------------------------------------
# The draft's dark regions are purple (not black), matching 'plasma'
draft_cmap = plt.cm.plasma


def periodic_tiled_lattice(
    L: int,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return all periodic-image sites (x, y, n1, n2) inside a rectangular window."""
    sqrt3 = np.sqrt(3.0)
    # Supercell translation vectors (isotropic display: a=1, b=√3)
    T1x, T1y = 2.0 * L, 0.0
    T2x, T2y = float(L), sqrt3 * L
    n_images = 3

    xs, ys, n1s, n2s = [], [], [], []
    for n1 in range(L):
        for n2 in range(L):
            x_base = 2.0 * n1 + n2
            y_base = sqrt3 * n2
            for i1 in range(-n_images, n_images + 1):
                for i2 in range(-n_images, n_images + 1):
                    x = x_base + i1 * T1x + i2 * T2x
                    y = y_base + i1 * T1y + i2 * T2y
                    if x_range[0] <= x <= x_range[1] and y_range[0] <= y <= y_range[1]:
                        xs.append(x)
                        ys.append(y)
                        n1s.append(n1)
                        n2s.append(n2)
    return (
        np.array(xs), np.array(ys),
        np.array(n1s, dtype=int), np.array(n2s, dtype=int),
    )


def min_image_coordinates(L: int) -> tuple[np.ndarray, np.ndarray]:
    """Return X[n2,n1], Y[n2,n1] for the closest periodic image to origin."""
    X = np.zeros((L, L))
    Y = np.zeros((L, L))
    sqrt3 = np.sqrt(3.0)
    for n2 in range(L):
        for n1 in range(L):
            best = (0.0, 0.0, np.inf)
            for dn1 in (-1, 0, 1):
                for dn2 in (-1, 0, 1):
                    nn1 = n1 + dn1 * L
                    nn2 = n2 + dn2 * L
                    x = 2.0 * nn1 + nn2
                    y = sqrt3 * nn2
                    r2 = x * x + y * y
                    if r2 < best[2]:
                        best = (x, y, r2)
            X[n2, n1] = best[0]
            Y[n2, n1] = best[1]
    return X, Y


def plot_hex_heatmap(ax, x, y, values, cmap, norm, title):
    """
    Plot a heatmap using large overlapping circular scatter markers.

    This matches the draft's "pointillist" style: individual dots are
    visible but large enough to overlap, creating a soft mosaic look.
    """
    sc = ax.scatter(
        x, y, c=values,
        cmap=cmap, norm=norm,
        marker='o',
        s=155,          # circles overlap smoothly at this size
        edgecolors='none',
        linewidths=0,
        zorder=2,
    )

    ax.set_xlim(0, 58)
    ax.set_ylim(0, 52)
    ax.set_aspect('equal')
    ax.set_facecolor(cmap(0.0))
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('x', fontsize=9)
    ax.set_ylabel('y', fontsize=9)
    ax.tick_params(labelsize=8)

    return sc


def main():
    # Load data
    counts, n_walkers = load_kmc_counts(HERE / "kmc_triangular_counts.txt", L=L)
    P_kmc = counts / n_walkers

    print(f"KMC: N = {n_walkers:,}")
    print("Computing corrected theory...")
    P_theory = theory_probability_on_lattice(
        build_Mk_corrected, gamma, epsilon, L=L, t_final=t_final,
    )

    # Periodic tiling in isotropic display coordinates
    # Tile slightly beyond the view window to avoid edge artifacts
    x_win = (-3.0, 62.0)
    y_win = (-3.0, 56.0)
    X_tile, Y_tile, n1_tile, n2_tile = periodic_tiled_lattice(L, x_win, y_win)

    P_theory_tile = P_theory[n2_tile, n1_tile]
    P_kmc_tile = P_kmc[n2_tile, n1_tile]

    # Shared colorbar normalization
    vmin = 0.0
    vmax = max(P_theory_tile.max(), P_kmc_tile.max())
    norm = Normalize(vmin=vmin, vmax=vmax)

    # ------- Figure -------
    fig = plt.figure(figsize=(12, 8.5), facecolor='white')

    # Top row: (a) Theory, (b) KMC — tight layout with dedicated colorbar axes
    ax1 = fig.add_axes([0.04, 0.44, 0.37, 0.50])
    cax1 = fig.add_axes([0.415, 0.44, 0.012, 0.50])
    ax2 = fig.add_axes([0.50, 0.44, 0.37, 0.50])
    cax2 = fig.add_axes([0.875, 0.44, 0.012, 0.50])

    print(f"Plotting {len(X_tile)} hex patches per panel...")

    pc1 = plot_hex_heatmap(ax1, X_tile, Y_tile, P_theory_tile,
                           draft_cmap, norm,
                           f"(a)\n$\\gamma$={gamma}, $\\varepsilon$={epsilon}, t={int(t_final)}")
    pc2 = plot_hex_heatmap(ax2, X_tile, Y_tile, P_kmc_tile,
                           draft_cmap, norm,
                           f"(b)\n$\\gamma$={gamma}, $\\varepsilon$={epsilon}, N={n_walkers:,}")

    # Dedicated colorbar axes
    cb1 = fig.colorbar(pc1, cax=cax1)
    cb2 = fig.colorbar(pc2, cax=cax2)
    cb1.ax.tick_params(labelsize=7)
    cb2.ax.tick_params(labelsize=7)

    # Bottom: cross-section at n2=0 (y=0)
    ax3 = fig.add_axes([0.15, 0.06, 0.72, 0.30])
    X_min, Y_min = min_image_coordinates(L)

    # Cross-section: sites with Y ≈ 0
    mask = np.isclose(Y_min, 0.0)
    x_cross = X_min[mask]
    p_kmc_cross = P_kmc[mask]
    p_theory_cross = P_theory[mask]
    order = np.argsort(x_cross)
    x_cross = x_cross[order]
    p_kmc_cross = p_kmc_cross[order]
    p_theory_cross = p_theory_cross[order]

    # Shift to match the tiled window (center at ~29)
    x_cross_shifted = x_cross + 29.0

    ax3.plot(x_cross_shifted, p_theory_cross, '-', color='navy', lw=1.5,
             label='Theory', zorder=2)
    ax3.plot(x_cross_shifted, p_kmc_cross, 'o', color='darkred', ms=3.5,
             mfc='none', label='KMC', zorder=3)
    ax3.set_xlabel('x', fontsize=10)
    ax3.set_ylabel('Probability', fontsize=10)
    ax3.set_title(f'$\\gamma$={gamma}, $\\varepsilon$={epsilon}, $t$={int(t_final)}',
                  fontsize=10)
    ax3.legend(fontsize=9, framealpha=0.9)
    ax3.set_xlim(0, 60)

    outpath = OUT / "fig11_dipanjan_style_v2.png"
    fig.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"Saved: {outpath}")
    # Also overwrite the old name
    outpath2 = OUT / "fig11_dipanjan_style.png"
    fig.savefig(outpath2, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"Saved: {outpath2}")
    plt.close(fig)


if __name__ == "__main__":
    main()
