"""
TCRW Phase 1: PBC simulations reproducing Fig 1 of Speck et al.
================================================================

Fig 1(b): Sample trajectories for different omega, D_r = 10^-3
Fig 1(c): MSD vs t for different omega, D_r = 10^-3
Fig 1(d): Diffusion coefficient D vs omega for different D_r

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from tcrw_core import simulate_tcrw_pbc, measure_diffusion_coeff

# ============================================================
# Fig 1(b): Sample trajectories — PBC, D_r = 10^-3
# ============================================================
def fig1b():
    """
    Single long trajectory colored by time for different omega values.
    Paper uses omega = 0.0, 0.5, 0.75, 1.0 with D_r = 10^-3.
    """
    D_r = 1e-3
    T_steps = 1_000_000
    L = 1000  # large L so we never wrap (effectively infinite plane)
    omegas = [0.0, 0.5, 0.75, 1.0]

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    fig.suptitle(f'Fig 1(b): Sample trajectories (PBC, $D_r = {D_r}$)',
                 fontsize=13, y=1.02)

    for ax, omega in zip(axes, omegas):
        res = simulate_tcrw_pbc(omega, D_r, L, T_steps, N_traj=1,
                                 seed=123, record_traj=True, record_interval=1)
        tx = res['traj_x'].astype(float)
        ty = res['traj_y'].astype(float)

        # Subsample for plotting (every 100 steps)
        skip = 100
        tx = tx[::skip]
        ty = ty[::skip]
        t_color = np.arange(len(tx))

        # Color by time
        points = np.column_stack([tx, ty]).reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap='plasma', linewidth=0.5)
        lc.set_array(t_color[:-1].astype(float))
        ax.add_collection(lc)

        # Mark start and end
        ax.plot(tx[0], ty[0], 'go', ms=6, zorder=5, label='Start')
        ax.plot(tx[-1], ty[-1], 'rs', ms=6, zorder=5, label='End')

        ax.set_aspect('equal')
        ax.autoscale()
        ax.set_title(f'$\\omega = {omega}$', fontsize=12)
        if omega == 0.0:
            ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig('tcrw_fig1b_trajectories.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Fig 1(b) saved: tcrw_fig1b_trajectories.png")


# ============================================================
# Fig 1(c): MSD vs t — PBC, D_r = 10^-3
# ============================================================
def fig1c():
    """
    MSD vs t on log-log axes for omega = 0.5, 0.7, 1.0.
    Paper also shows MSD ~ t reference line.
    D_r = 10^-3.
    """
    D_r = 1e-3
    T_steps = 2_000_000
    N_traj = 1000
    L = 2000
    omegas = [0.5, 0.7, 1.0]
    colors = ['tab:blue', 'tab:orange', 'tab:red']
    record_interval = 100

    fig, ax = plt.subplots(figsize=(6, 5))

    for omega, color in zip(omegas, colors):
        print(f"  MSD: omega={omega}, N_traj={N_traj}, T={T_steps}...")
        res = simulate_tcrw_pbc(omega, D_r, L, T_steps, N_traj,
                                 seed=42, record_interval=record_interval)
        t = res['times'].astype(float)
        msd = res['msd']

        # Skip t=0
        mask = t > 0
        ax.loglog(t[mask], msd[mask], color=color, lw=1.5,
                  label=f'$\\omega = {omega}$')

    # Reference line: MSD ~ t
    t_ref = np.logspace(0, np.log10(T_steps), 100)
    ax.loglog(t_ref, t_ref, 'k--', lw=1, alpha=0.5, label='$\\mathrm{MSD} \\propto t$')

    ax.set_xlabel('$t$', fontsize=13)
    ax.set_ylabel('MSD', fontsize=13)
    ax.set_title(f'Fig 1(c): MSD vs $t$ ($D_r = {D_r}$)', fontsize=13)
    ax.legend(fontsize=11)
    ax.set_xlim(1, T_steps)

    plt.tight_layout()
    plt.savefig('tcrw_fig1c_MSD.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Fig 1(c) saved: tcrw_fig1c_MSD.png")


# ============================================================
# Fig 1(d): Diffusion coefficient D vs omega
# ============================================================
def fig1d():
    """
    D vs omega for different D_r values.
    Paper shows D_r = 0.0001, 0.001, 0.01.
    D should decrease linearly with chirality (away from omega=0.5),
    and be independent of D_r.
    """
    D_r_values = [0.0001, 0.001, 0.01]
    markers = ['o', 's', '^']
    colors = ['tab:blue', 'tab:orange', 'tab:green']

    # omega grid: 0 to 1 in steps of 0.05
    omega_grid = np.arange(0.0, 1.01, 0.05)

    fig, ax = plt.subplots(figsize=(6, 5))

    for D_r, marker, color in zip(D_r_values, markers, colors):
        D_vals = []
        for omega in omega_grid:
            print(f"  D(omega={omega:.2f}, D_r={D_r})...")
            D = measure_diffusion_coeff(omega, D_r, L=500, T_steps=500000,
                                         N_traj=300, seed=42)
            D_vals.append(D)
        ax.plot(omega_grid, D_vals, marker=marker, color=color, ms=5, lw=1.5,
                label=f'$D_r = {D_r}$')

    ax.set_xlabel('$\\omega$', fontsize=13)
    ax.set_ylabel('$D$', fontsize=13)
    ax.set_title('Fig 1(d): Diffusion coefficient vs chirality', fontsize=13)
    ax.legend(fontsize=11)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig('tcrw_fig1d_D_vs_omega.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Fig 1(d) saved: tcrw_fig1d_D_vs_omega.png")


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    print("=" * 60)
    print("TCRW Phase 1: Reproducing Fig 1 (PBC)")
    print("=" * 60)

    print("\n--- Fig 1(b): Trajectories ---")
    fig1b()

    print("\n--- Fig 1(c): MSD ---")
    fig1c()

    print("\n--- Fig 1(d): D vs omega ---")
    fig1d()

    print("\n" + "=" * 60)
    print("Phase 1 complete.")
    print("=" * 60)
