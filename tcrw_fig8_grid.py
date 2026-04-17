"""
Fig 8(d)-(l): nine OBC spectrum vignettes at specific (omega, D_r) points.

Covers the (omega, D_r) plane spanning topological, trivial, and transition
regimes. Topological phase: D_r < 0.5 AND omega != 0.5 with Zak phases pi.
Each panel: full OBC spectrum in complex plane, colored by edge weight.

Following the paper's sampling: rows span D_r values, columns span omega.
We choose a 3x3 grid that clearly shows: topological region (edges complex),
gap-closing line (omega=0.5), and trivial region (D_r large).

Author: Prashant Bisht, TIFR Hyderabad
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from tcrw_spectrum import obc_spectrum


def fig8_spectrum_grid(L=8, outfile='tcrw_fig8dl_spectrum_grid.png'):
    # 3x3 grid of parameter points — matches paper Fig 8(c) sampling.
    # ω ∈ {0.25, 0.5, 0.75} crosses ω=0.5 gap-closing line.
    # D_r ∈ {0.25, 0.5, 0.75} crosses D_r=0.5 topological-to-trivial transition.
    omega_vals = [0.25, 0.5, 0.75]
    D_r_vals   = [0.25, 0.5, 0.75]

    fig, axes = plt.subplots(3, 3, figsize=(11, 11))

    t0 = time.time()
    idx = 0
    panel_labels = ['(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)']

    for r, D_r in enumerate(D_r_vals):
        for c, om in enumerate(omega_vals):
            ax = axes[r, c]
            print(f"  Panel {panel_labels[idx]}: omega={om}, D_r={D_r}, L={L}...")
            evals, edge_w = obc_spectrum(om, D_r, L)

            # Normalize edge weight for colouring
            norm_w = edge_w / edge_w.max() if edge_w.max() > 0 else edge_w
            sc = ax.scatter(evals.real, evals.imag, c=norm_w, cmap='viridis',
                             s=10, alpha=0.85, edgecolors='none',
                             norm=Normalize(0, 1))

            # Unit circle for reference
            theta = np.linspace(0, 2*np.pi, 200)
            ax.plot(np.cos(theta), np.sin(theta), 'k--', lw=0.5, alpha=0.4)

            # Mark lambda=1
            ax.plot([1], [0], 'r*', markersize=10, mec='k', mew=0.5, zorder=10)

            # Axis cosmetics
            ax.axhline(0, color='k', lw=0.3, alpha=0.4)
            ax.axvline(0, color='k', lw=0.3, alpha=0.4)
            ax.set_xlim(-1.1, 1.1)
            ax.set_ylim(-1.1, 1.1)
            ax.set_aspect('equal')

            # Regime label
            if D_r < 0.5 and abs(om - 0.5) > 0.05:
                regime = 'topological'
                color = 'tab:green'
            elif abs(om - 0.5) < 0.05:
                regime = 'gap-closing'
                color = 'tab:red'
            else:
                regime = 'trivial'
                color = 'tab:gray'

            ax.set_title(f"{panel_labels[idx]} $\\omega={om}$, $D_r={D_r}$"
                          f"\n[{regime}]", fontsize=10, color=color)

            if r == 2:
                ax.set_xlabel('$\\mathrm{Re}(\\lambda)$', fontsize=11)
            if c == 0:
                ax.set_ylabel('$\\mathrm{Im}(\\lambda)$', fontsize=11)

            idx += 1

    # Single colorbar for the whole figure
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(sc, cax=cbar_ax, label='edge weight')

    fig.suptitle(f'Fig 8(d)-(l): OBC spectrum at 9 points in $(\\omega, D_r)$ plane, '
                  f'$L={L}$', fontsize=13, y=0.995)
    plt.subplots_adjust(right=0.90, hspace=0.3, wspace=0.25, top=0.94)
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outfile}  (total {time.time()-t0:.1f}s)")


if __name__ == '__main__':
    fig8_spectrum_grid()
