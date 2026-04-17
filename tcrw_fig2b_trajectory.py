"""
Fig 2(b): single long OBC trajectory of a CRW, 10^6 steps.

Reproduces Fig 2(b) of Osat et al. (arXiv:2602.12020).
Walker with omega=1 (fully chiral CW), D_r=1e-3, L=10, T=10^6 steps.
Trajectory rendered as a LineCollection colored by time.

Author: Prashant Bisht, TIFR Hyderabad
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, Normalize
from tcrw_core import simulate_tcrw_obc


def fig2b_trajectory(omega=1.0, D_r=1e-3, L=10, T_steps=10**6, seed=12345,
                      outfile='tcrw_fig2b_trajectory.png'):
    t0 = time.time()
    print(f"  Running OBC walker: omega={omega}, D_r={D_r}, L={L}, T={T_steps:.0e}")
    res = simulate_tcrw_obc(omega, D_r, L, T_steps, N_traj=1,
                             seed=seed, record_traj=True, record_interval=1)
    print(f"  Walk finished in {time.time()-t0:.1f}s")

    x = res['traj_x'].astype(float)
    y = res['traj_y'].astype(float)
    # Jitter points slightly off the integer lattice so overlapping visits
    # are visible (paper uses a similar visualisation)
    rng = np.random.default_rng(seed + 1)
    xj = x + rng.uniform(-0.3, 0.3, size=len(x))
    yj = y + rng.uniform(-0.3, 0.3, size=len(y))

    # LineCollection coloured by time (log scale like paper)
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    pts = np.column_stack([xj, yj]).reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    tvec = np.arange(1, len(xj))
    lc = LineCollection(segs, cmap='inferno',
                         norm=LogNorm(vmin=1, vmax=tvec.max()),
                         linewidth=0.35, alpha=0.7)
    lc.set_array(tvec)
    ax.add_collection(lc)

    # Start / End markers
    ax.plot(xj[0], yj[0], 'o', color='lime', markersize=8, mec='k',
            mew=0.8, label='start', zorder=10)
    ax.plot(xj[-1], yj[-1], 'o', color='red', markersize=8, mec='k',
            mew=0.8, label=f'end (t={T_steps:.0e})', zorder=10)

    # Show the hard-wall box
    ax.add_patch(plt.Rectangle((-0.5, -0.5), L, L, fill=False,
                                edgecolor='k', linewidth=1.5))

    ax.set_xlim(-0.8, L - 0.2)
    ax.set_ylim(-0.8, L - 0.2)
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_aspect('equal')
    ax.set_title(f'Fig 2(b): OBC trajectory, $\\omega={omega}$, $D_r={D_r}$, $L={L}$, $T={T_steps:.0e}$',
                 fontsize=11)
    cbar = fig.colorbar(lc, ax=ax, label='time step $t$', shrink=0.85, pad=0.02)
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outfile}")


if __name__ == '__main__':
    fig2b_trajectory()
