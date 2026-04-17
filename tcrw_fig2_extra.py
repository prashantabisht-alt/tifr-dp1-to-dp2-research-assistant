"""
tcrw_fig2_extra.py — Reproduce Fig 2(b) and Fig 2(i)-(o) from arXiv:2602.12020

Figure 1: Fig 2(b) — Single long OBC trajectory
  - Run TCRW OBC simulator for ω=1, D_r=1e-3, L=10, T=10^6 steps
  - Plot trajectory colored by time with grid overlay
  - Save as tcrw_fig2b_trajectory.png

Figure 2: Fig 2(i)-(o) — Defect boundary panels (6 panels)
  - 2 defect geometries: edge notch + internal hole
  - For each: compute P(x,y), J total current, J_ω chiral current
  - Layout: 2 rows × 3 columns (P, J, J_ω for each geometry)
  - Save as tcrw_fig2_defects.png

Author: Claude Code Agent
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm, Normalize

# Add path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tcrw_core import simulate_tcrw_obc
from tcrw_geometry import RectangleWithDefects, DX, DY
from tcrw_obc import build_transition_matrix_generic, exact_steady_state_generic


def figure1_single_trajectory():
    """
    Fig 2(b): Single long OBC trajectory, ω=1, D_r=1e-3, L=10, T=10^6
    """
    print("Generating Figure 1: Single OBC trajectory...")

    # Run simulation
    res = simulate_tcrw_obc(
        omega=1.0, D_r=1e-3, L=10, T_steps=1000000, N_traj=1,
        seed=42, record_traj=True, record_interval=1
    )

    traj_x = res['traj_x']
    traj_y = res['traj_y']
    times = res['times']

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot grid
    for i in range(11):
        ax.axhline(y=i - 0.5, color='gray', linewidth=0.5, alpha=0.3)
        ax.axvline(x=i - 0.5, color='gray', linewidth=0.5, alpha=0.3)

    # Build 2D visit histogram for background heatmap
    L = 10
    visit_count = np.zeros((L, L))
    for xi, yi in zip(traj_x, traj_y):
        visit_count[int(xi), int(yi)] += 1
    # Log-scale visit density
    visit_log = np.log10(1 + visit_count)
    ax.imshow(visit_log.T, origin='lower', cmap='Blues', alpha=0.5,
              extent=(-0.5, L-0.5, -0.5, L-0.5), interpolation='nearest')

    # Plot trajectory as line segments (subsample for visibility)
    # Use LineCollection for efficiency
    from matplotlib.collections import LineCollection
    # Subsample: take every 100th point to show path structure
    step = max(1, len(traj_x) // 10000)
    tx = traj_x[::step]
    ty = traj_y[::step]
    tt = times[::step]

    points = np.array([tx, ty]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(tt.min(), tt.max())
    lc = LineCollection(segments, cmap='viridis', norm=norm, alpha=0.3, linewidths=0.5)
    lc.set_array(tt[:-1])
    ax.add_collection(lc)

    # Overlay scatter at subsampled points
    scatter = ax.scatter(tx, ty, c=tt, cmap='viridis', s=2, alpha=0.4, norm=norm)

    # Labels and formatting
    ax.set_xlim(-0.5, 9.5)
    ax.set_ylim(-0.5, 9.5)
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(f'Fig 2(b): TCRW OBC Trajectory (ω=1, D_r=1e-3, L=10, T=10^6)', fontsize=14)

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Time (steps)', fontsize=11)

    plt.tight_layout()
    plt.savefig('tcrw_fig2b_trajectory.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved tcrw_fig2b_trajectory.png")
    print(f"    Trajectory length: {len(traj_x)} points")
    return fig


def build_split_matrices_generic(omega, D_r, mask):
    """
    Build noise-only and chiral-only transition sub-matrices for arbitrary geometry.

    P_full = P_noise + P_chiral (both are sub-stochastic; they sum to stochastic).

    Returns
    -------
    P_noise : sparse matrix (noise step contributions only)
    P_chiral : sparse matrix (chiral step contributions only)
    """
    import scipy.sparse as sp

    N = mask.n_states
    rows_n, cols_n, vals_n = [], [], []
    rows_c, cols_c, vals_c = [], [], []

    for x, y in mask.valid_sites:
        for d in range(4):
            i = mask.state_index(x, y, d)

            # --- Noise step (prob D_r): stay put, rotate ---
            d_ccw = (d - 1) % 4
            j_ccw = mask.state_index(x, y, d_ccw)
            rows_n.append(j_ccw); cols_n.append(i); vals_n.append(omega * D_r)

            d_cw = (d + 1) % 4
            j_cw = mask.state_index(x, y, d_cw)
            rows_n.append(j_cw); cols_n.append(i); vals_n.append((1 - omega) * D_r)

            # --- Chiral step (prob 1-D_r): translate + rotate ---
            nb = mask.neighbor(x, y, d)
            if nb is not None:
                nx, ny = nb
                d_cw_c = (d + 1) % 4
                j1 = mask.state_index(nx, ny, d_cw_c)
                rows_c.append(j1); cols_c.append(i); vals_c.append(omega * (1 - D_r))

                d_ccw_c = (d - 1) % 4
                j2 = mask.state_index(nx, ny, d_ccw_c)
                rows_c.append(j2); cols_c.append(i); vals_c.append((1 - omega) * (1 - D_r))
            else:
                # Blocked: stay at (x,y,d)
                j_stay = mask.state_index(x, y, d)
                rows_c.append(j_stay); cols_c.append(i); vals_c.append(1 - D_r)

    P_noise = sp.coo_matrix((vals_n, (rows_n, cols_n)), shape=(N, N)).tocsc()
    P_chiral = sp.coo_matrix((vals_c, (rows_c, cols_c)), shape=(N, N)).tocsc()
    return P_noise, P_chiral


def compute_currents_from_steady_state(omega, D_r, mask):
    """
    Compute exact currents J, J_Dr, J_omega for arbitrary lattice geometry.

    Decomposition follows the paper (Osat et al.):
      J_total(x,y) = J_Dr(x,y) + J_omega(x,y)

    where J_Dr weights by probability of ARRIVING via noise step (pi_N),
    and J_omega weights by probability of arriving via chiral step (pi_C).

    Returns
    -------
    Jx, Jy : dict mapping (x,y) -> total current component
    Jx_omega, Jy_omega : dict mapping (x,y) -> chiral current component
    Jx_Dr, Jy_Dr : dict mapping (x,y) -> noise-induced current component
    Pxy : dict mapping (x,y) -> occupation probability
    """
    import numpy as np

    # Get steady state
    Pxy, pi_vec = exact_steady_state_generic(omega, D_r, mask)

    # Split transition matrix
    P_noise, P_chiral = build_split_matrices_generic(omega, D_r, mask)

    # Arrival probabilities: who arrived via noise vs chiral
    pi_N = P_noise.toarray() @ pi_vec   # arrived via noise step
    pi_C = P_chiral.toarray() @ pi_vec  # arrived via chiral step

    n_sites = mask.n_sites

    Jx, Jy = {}, {}
    Jx_omega, Jy_omega = {}, {}
    Jx_Dr, Jy_Dr = {}, {}

    for si, (x, y) in enumerate(mask.valid_sites):
        jx_t, jy_t = 0.0, 0.0
        jx_om, jy_om = 0.0, 0.0
        jx_dr, jy_dr = 0.0, 0.0

        for d in range(4):
            state_i = d * n_sites + si
            nb = mask.neighbor(x, y, d)
            if nb is not None:
                dx_move = float(DX[d])
                dy_move = float(DY[d])

                # Total spatial current from (x,y,d)
                flux = pi_vec[state_i] * (1 - D_r)
                jx_t += flux * dx_move
                jy_t += flux * dy_move

                # Decomposed by arrival type
                flux_dr = pi_N[state_i] * (1 - D_r)
                jx_dr += flux_dr * dx_move
                jy_dr += flux_dr * dy_move

                flux_om = pi_C[state_i] * (1 - D_r)
                jx_om += flux_om * dx_move
                jy_om += flux_om * dy_move

        Jx[(x, y)] = jx_t
        Jy[(x, y)] = jy_t
        Jx_omega[(x, y)] = jx_om
        Jy_omega[(x, y)] = jy_om
        Jx_Dr[(x, y)] = jx_dr
        Jy_Dr[(x, y)] = jy_dr

    return Jx, Jy, Jx_omega, Jy_omega, Jx_Dr, Jy_Dr, Pxy


def dict_to_array(data_dict, L):
    """Convert dict mapping (x,y) -> value to L×L array."""
    arr = np.zeros((L, L))
    for (x, y), val in data_dict.items():
        if 0 <= x < L and 0 <= y < L:
            arr[x, y] = val
    return arr


def figure2_defect_boundaries():
    """
    Fig 2(i)-(o): Defect boundary panels (6 panels)
    Two geometries × three observables (P, J, J_ω)
    """
    print("Generating Figure 2: Defect boundary panels...")

    L = 10
    omega = 1.0
    D_r = 1e-3

    # Define geometries
    # Geometry 1: Edge notch (deformed boundary)
    # Remove a 3×2 rectangular notch from right edge
    blocked_notch = [(8, 4), (8, 5), (9, 4), (9, 5), (8, 3), (9, 3)]
    mask_notch = RectangleWithDefects(L, blocked_notch)

    # Geometry 2: Internal hole (3×3 block at center)
    blocked_hole = []
    cx, cy = L // 2, L // 2
    for dx in range(-1, 2):
        for dy in range(-1, 2):
            blocked_hole.append((cx + dx, cy + dy))
    mask_hole = RectangleWithDefects(L, blocked_hole)

    # Compute currents for both geometries (properly decomposed)
    print("  Computing currents for edge notch...")
    Jx_n, Jy_n, Jx_om_n, Jy_om_n, Jx_dr_n, Jy_dr_n, Pxy_n = \
        compute_currents_from_steady_state(omega, D_r, mask_notch)

    print("  Computing currents for internal hole...")
    Jx_h, Jy_h, Jx_om_h, Jy_om_h, Jx_dr_h, Jy_dr_h, Pxy_h = \
        compute_currents_from_steady_state(omega, D_r, mask_hole)

    # Convert to arrays
    P_notch = dict_to_array(Pxy_n, L)
    P_hole = dict_to_array(Pxy_h, L)

    # Set blocked sites to NaN for visualization
    for (x, y) in blocked_notch:
        P_notch[x, y] = np.nan
    for (x, y) in blocked_hole:
        P_hole[x, y] = np.nan

    # Helper: auto-scaled quiver plot for currents
    def plot_currents(ax, Jx_dict, Jy_dict, L, title, color, blocked_sites=None):
        """Plot current vector field with auto-scaling."""
        # Collect into arrays
        xs, ys, us, vs = [], [], [], []
        for (x, y), jx in Jx_dict.items():
            jy = Jy_dict[(x, y)]
            xs.append(x); ys.append(y)
            us.append(jx); vs.append(jy)
        xs, ys = np.array(xs), np.array(ys)
        us, vs = np.array(us), np.array(vs)
        mags = np.sqrt(us**2 + vs**2)
        max_mag = mags.max() if len(mags) > 0 else 1e-10

        # Background: show lattice + blocked sites
        bg = np.ones((L, L, 3)) * 0.95
        if blocked_sites:
            for (bx, by) in blocked_sites:
                if 0 <= bx < L and 0 <= by < L:
                    bg[bx, by] = [0.3, 0.3, 0.3]
        ax.imshow(bg.transpose(1, 0, 2), origin='lower',
                  extent=(-0.5, L-0.5, -0.5, L-0.5), interpolation='nearest')

        if max_mag > 1e-20:
            # Normalize arrows and color by magnitude
            norm_mags = mags / max_mag
            ax.quiver(xs, ys, us / max_mag, vs / max_mag,
                      norm_mags, cmap='coolwarm',
                      scale=15, scale_units='width', width=0.004,
                      alpha=0.85, clim=[0, 1])
        ax.set_xlim(-0.5, L - 0.5)
        ax.set_ylim(-0.5, L - 0.5)
        ax.set_aspect('equal')
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.grid(alpha=0.2)

    # Create figure: 2 rows × 3 columns
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # --- Row 1: Edge notch geometry ---

    # Panel (i): P(x,y) for notch
    ax = axes[0, 0]
    im = ax.imshow(P_notch.T, origin='lower', cmap='YlOrRd', interpolation='nearest')
    ax.set_title('(i) Edge Notch: P(x,y)', fontsize=12, fontweight='bold')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(im, ax=ax)

    # Panel (j): J(x,y) for notch
    plot_currents(axes[0, 1], Jx_n, Jy_n, L,
                  '(j) Edge Notch: J(x,y) total', 'blue', blocked_notch)

    # Panel (k): J_Dr(x,y) for notch
    plot_currents(axes[0, 2], Jx_dr_n, Jy_dr_n, L,
                  '(k) Edge Notch: J_Dr(x,y) noise', 'green', blocked_notch)

    # --- Row 2: Internal hole geometry ---

    # Panel (l): P(x,y) for hole
    ax = axes[1, 0]
    im = ax.imshow(P_hole.T, origin='lower', cmap='YlOrRd', interpolation='nearest')
    ax.set_title('(l) Internal Hole: P(x,y)', fontsize=12, fontweight='bold')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(im, ax=ax)

    # Panel (m): J(x,y) for hole
    plot_currents(axes[1, 1], Jx_h, Jy_h, L,
                  '(m) Internal Hole: J(x,y) total', 'blue', blocked_hole)

    # Panel (n): J_Dr(x,y) for hole
    plot_currents(axes[1, 2], Jx_dr_h, Jy_dr_h, L,
                  '(n) Internal Hole: J_Dr(x,y) noise', 'green', blocked_hole)

    plt.suptitle(f'Fig 2(i)-(o): Defect Boundaries (ω={omega}, D_r={D_r})',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig('tcrw_fig2_defects.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved tcrw_fig2_defects.png")
    return fig


if __name__ == '__main__':
    print("\n" + "="*70)
    print("TCRW Figure Reproduction: Fig 2(b) and Fig 2(i)-(o)")
    print("="*70)

    # Generate both figures
    figure1_single_trajectory()
    figure2_defect_boundaries()

    print("\n" + "="*70)
    print("All figures generated successfully!")
    print("="*70)
    print("\nOutput files:")
    print("  - tcrw_fig2b_trajectory.png")
    print("  - tcrw_fig2_defects.png")
