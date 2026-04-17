#!/usr/bin/env python3
"""
jcABP_simulation.py
====================
Euler-Maruyama simulation of the jerky chiral ABP (jcABP).

7-variable state per trajectory: (x, vx, ax, y, vy, ay, theta)

EOM (from Eq. 3, first-order form):
  dx/dt   = vx
  dvx/dt  = ax
  dax/dt  = -ax/tau_J - vx/tau_F^2 + v0*cos(theta)/tau_F^2 + sqrt(2D)/tau_F^2 * eta_x(t)
  dy/dt   = vy
  dvy/dt  = ay
  day/dt  = -ay/tau_J - vy/tau_F^2 + v0*sin(theta)/tau_F^2 + sqrt(2D)/tau_F^2 * eta_y(t)
  dtheta/dt = omega_0 + sqrt(2*D_r) * eta_r(t)

where omega_0 = 1/tau_C, D_r = 1/tau_P (or 0 for D_r=0 case).

All trajectories start at rest: x=vx=ax=y=vy=ay=theta=0.

Vectorised over N_traj trajectories using NumPy.
Accumulates ensemble averages on-the-fly: <x>, <y>, <x^2+y^2>.
"""

import numpy as np
import time as pytime


def run_jcABP_sim(tau_J, tau_F, tau_P, tau_C, D=0.5, v0=1.0, m=1.0,
                  T_max=20.0, dt=0.01, N_traj=10000, save_every=10,
                  seed=42, verbose=True):
    """
    Run N_traj trajectories of the jcABP model.

    Parameters
    ----------
    tau_J, tau_F : jerk and friction timescales
    tau_P : persistence time (= 1/D_r). Use tau_P > 1e10 for D_r = 0.
    tau_C : chirality time (= 1/omega_0). Use tau_C > 1e10 for omega_0 = 0.
    D     : translational diffusion coefficient
    v0    : self-propulsion speed
    m     : mass (default 1)
    T_max : total simulation time
    dt    : time step
    N_traj: number of ensemble trajectories
    save_every : save statistics every this many steps
    seed  : RNG seed

    Returns
    -------
    dict with keys:
      't'      : (N_save,) time array
      'mean_x' : (N_save,) ensemble-averaged <x(t)>
      'mean_y' : (N_save,) ensemble-averaged <y(t)>
      'msd'    : (N_save,) ensemble-averaged <x^2 + y^2>
    """
    rng = np.random.default_rng(seed)

    N_steps = int(T_max / dt)
    N_save  = N_steps // save_every + 1
    sqrt_dt = np.sqrt(dt)

    # Derived constants
    omega0 = 1.0 / tau_C if tau_C < 1e10 else 0.0
    D_r    = 1.0 / tau_P if tau_P < 1e10 else 0.0
    inv_tJ = 1.0 / tau_J
    inv_tF2 = 1.0 / tau_F**2
    noise_a = np.sqrt(2.0 * D) * inv_tF2   # prefactor for translational noise in a-equation
    noise_r = np.sqrt(2.0 * D_r)            # prefactor for rotational noise

    # State: all start at zero
    x  = np.zeros(N_traj)
    vx = np.zeros(N_traj)
    ax = np.zeros(N_traj)
    y  = np.zeros(N_traj)
    vy = np.zeros(N_traj)
    ay = np.zeros(N_traj)
    th = np.zeros(N_traj)

    # Output arrays
    t_out    = np.zeros(N_save)
    mean_x   = np.zeros(N_save)
    mean_y   = np.zeros(N_save)
    mean_r2  = np.zeros(N_save)

    # Save initial condition
    t_out[0] = 0.0
    save_idx = 1

    t0 = pytime.time()

    for step in range(1, N_steps + 1):
        # Random kicks
        Wx = rng.standard_normal(N_traj) * sqrt_dt
        Wy = rng.standard_normal(N_traj) * sqrt_dt
        Wr = rng.standard_normal(N_traj) * sqrt_dt

        cos_th = np.cos(th)
        sin_th = np.sin(th)

        # Euler-Maruyama update
        # Acceleration (jerk) equation
        dax = (-inv_tJ * ax - inv_tF2 * vx + v0 * inv_tF2 * cos_th) * dt + noise_a * Wx
        day = (-inv_tJ * ay - inv_tF2 * vy + v0 * inv_tF2 * sin_th) * dt + noise_a * Wy

        # Update in order: positions, velocities, accelerations, angle
        x  += vx * dt
        y  += vy * dt
        vx += ax * dt
        vy += ay * dt
        ax += dax
        ay += day
        th += omega0 * dt + noise_r * Wr

        # Save statistics
        if step % save_every == 0 and save_idx < N_save:
            t_out[save_idx]   = step * dt
            mean_x[save_idx]  = np.mean(x)
            mean_y[save_idx]  = np.mean(y)
            mean_r2[save_idx] = np.mean(x**2 + y**2)
            save_idx += 1

    elapsed = pytime.time() - t0
    if verbose:
        print(f"  Simulation done: {N_steps} steps x {N_traj} traj in {elapsed:.1f}s")

    return dict(t=t_out[:save_idx], mean_x=mean_x[:save_idx],
                mean_y=mean_y[:save_idx], msd=mean_r2[:save_idx])


# =============================================================================
# VERIFICATION: Section 3 (no chirality) — Figs 1 & 2
# =============================================================================

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from scipy.signal import fftconvolve

    plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12,
                         'legend.fontsize': 10, 'lines.linewidth': 2})

    # --- Re-use analytic functions from Section 3 ---
    def get_poles(tau_J, tau_F):
        disc = 1.0/tau_F**2 - 1.0/(4.0*tau_J**2)
        sq = np.sqrt(complex(disc))
        return (-1j/(2*tau_J) + sq), (-1j/(2*tau_J) - sq)

    def G_vec(t_arr, tau_J, tau_F, m=1.0):
        lam = tau_J * m
        o1, o2 = get_poles(tau_J, tau_F)
        denom = lam * (o2 - o1)
        t = np.asarray(t_arr, dtype=float)
        tc = t.astype(complex)
        term1 = (1 - np.exp(-1j*o1*tc)) / o1
        term2 = (1 - np.exp(-1j*o2*tc)) / o2
        G = ((-term1 + term2) / denom).real
        G[t <= 0] = 0.0
        return G

    def compute_section3_analytic(tau_J, tau_F, tau_P=1.0, D=0.5, v0=1.0, m=1.0,
                                   T_max=30.0, N=15000):
        dt = T_max / N
        t = np.linspace(0.0, T_max, N+1)
        gamma = tau_J * m / tau_F**2
        G = G_vec(t, tau_J, tau_F, m)
        E = np.exp(-t / tau_P)
        conv_GE = fftconvolve(G, E)[:N+1] * dt
        mean_x = gamma * v0 * conv_GE
        msd_act = 2.0 * gamma**2 * v0**2 * np.cumsum(G * conv_GE) * dt
        msd_thm = 4.0 * D * gamma**2 * np.cumsum(G**2) * dt
        return t, mean_x, msd_act + msd_thm

    # ---------------------------------------------------------------
    # Section 3 panel (a): tau_J=0.2, tau_F=0.28  (smooth relaxation)
    # ---------------------------------------------------------------
    print("=== Section 3 panel (a): tau_J=0.2, tau_F=0.28 ===")
    tJ_a, tF_a = 0.2, 0.28
    tau_P, v0, D = 1.0, 1.0, 0.5
    T_max_a = 20.0

    # Analytic
    t_an, mx_an, msd_an = compute_section3_analytic(tJ_a, tF_a, tau_P, D, v0,
                                                      T_max=T_max_a, N=20000)

    # Simulation — dt must be << tau_J = 0.2
    sim_a = run_jcABP_sim(tJ_a, tF_a, tau_P, tau_C=1e15, D=D, v0=v0,
                           T_max=T_max_a, dt=0.002, N_traj=20000,
                           save_every=25, seed=123)

    # ---------------------------------------------------------------
    # Section 3 panel (b): tau_J=20, tau_F=1.0  (oscillatory)
    # ---------------------------------------------------------------
    print("\n=== Section 3 panel (b): tau_J=20, tau_F=1.0 ===")
    tJ_b, tF_b = 20.0, 1.0
    T_max_b = 200.0

    t_bn, mx_bn, msd_bn = compute_section3_analytic(tJ_b, tF_b, tau_P, D, v0,
                                                      T_max=T_max_b, N=60000)

    sim_b = run_jcABP_sim(tJ_b, tF_b, tau_P, tau_C=1e15, D=D, v0=v0,
                           T_max=T_max_b, dt=0.02, N_traj=20000,
                           save_every=50, seed=456)

    # ---------------------------------------------------------------
    # OVERLAY PLOTS
    # ---------------------------------------------------------------
    lP = v0 * tau_P

    # Fig 1 overlay: mean displacement
    fig1, axes1 = plt.subplots(1, 2, figsize=(11, 4.5))

    for ax, label, t_a, mx_a, sim in [
        (axes1[0], '(a)', t_an, mx_an, sim_a),
        (axes1[1], '(b)', t_bn, mx_bn, sim_b),
    ]:
        # Analytic
        mask_a = t_a > 0.01
        ax.loglog(t_a[mask_a]/tau_P, mx_a[mask_a]/lP, 'r-', lw=2.5, label='Analytic')
        # Simulation
        mask_s = sim['t'] > 0.01
        ax.loglog(sim['t'][mask_s]/tau_P, sim['mean_x'][mask_s]/lP,
                  'ko', ms=2, alpha=0.5, label=f"MC ({sim['t'].size} pts)")
        ax.set_xlabel(r'$t/\tau_P$')
        ax.set_ylabel(r'$\langle x(t)\rangle / l_P$')
        ax.set_title(f'{label}')
        ax.legend()
        ax.grid(True, which='both', alpha=0.25)

    fig1.suptitle('Fig 1 — Mean displacement: Analytic vs MC', fontsize=12)
    plt.tight_layout()
    out1 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig1_overlay_MC.png'
    plt.savefig(out1, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {out1}")

    # Fig 2 overlay: MSD
    fig2, axes2 = plt.subplots(1, 2, figsize=(11, 4.5))

    for ax, label, t_a, msd_a, sim in [
        (axes2[0], '(a)', t_an, msd_an, sim_a),
        (axes2[1], '(b)', t_bn, msd_bn, sim_b),
    ]:
        mask_a = (t_a > 0.01) & (msd_a > 1e-10)
        ax.loglog(t_a[mask_a]/tau_P, msd_a[mask_a]/lP**2, 'r-', lw=2.5, label='Analytic')
        mask_s = (sim['t'] > 0.01) & (sim['msd'] > 1e-10)
        ax.loglog(sim['t'][mask_s]/tau_P, sim['msd'][mask_s]/lP**2,
                  'ko', ms=2, alpha=0.5, label='MC')
        ax.set_xlabel(r'$t/\tau_P$')
        ax.set_ylabel(r'$\mathrm{MSD}(t) / l_P^2$')
        ax.set_title(f'{label}')
        ax.legend()
        ax.grid(True, which='both', alpha=0.25)

    fig2.suptitle('Fig 2 — MSD: Analytic vs MC', fontsize=12)
    plt.tight_layout()
    out2 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig2_overlay_MC.png'
    plt.savefig(out2, dpi=150, bbox_inches='tight')
    print(f"Saved: {out2}")

    # ---------------------------------------------------------------
    # Section 4.1: D_r = 0, finite chirality — Figs 3 & 4
    # ---------------------------------------------------------------
    print("\n=== Section 4.1 panel (e): tau_J=0.2, tau_F=0.2 (D_r=0) ===")
    tJ_41, tF_41, tC_41 = 0.2, 0.2, 1.0
    T_max_41 = 25.0

    sim_41 = run_jcABP_sim(tJ_41, tF_41, tau_P=1e15, tau_C=tC_41,
                            D=D, v0=v0, T_max=T_max_41, dt=0.002,
                            N_traj=10000, save_every=25, seed=789)

    # Analytic (from section 4.1 code)
    def compute_section4p1_analytic(tau_J, tau_F, tau_C, D, v0, T_max, N):
        dt = T_max / N
        t = np.linspace(0.0, T_max, N+1)
        gamma = tau_J / tau_F**2
        G = G_vec(t, tau_J, tau_F)
        C = np.cos(t / tau_C)
        S = np.sin(t / tau_C)
        conv_GC = fftconvolve(G, C)[:N+1] * dt
        conv_GS = fftconvolve(G, S)[:N+1] * dt
        mean_x = gamma * v0 * conv_GC
        mean_y = gamma * v0 * conv_GS
        mean_r2 = mean_x**2 + mean_y**2
        msd_thm = 4.0 * D * gamma**2 * np.cumsum(G**2) * dt
        return t, mean_x, mean_y, mean_r2 + msd_thm

    t_41a, mx_41a, my_41a, msd_41a = compute_section4p1_analytic(
        tJ_41, tF_41, tC_41, D, v0, T_max_41, 30000)

    # Trajectory overlay (Fig 3 panel e style)
    fig3o, ax3 = plt.subplots(1, 1, figsize=(6, 5))
    ax3.plot(my_41a, mx_41a, 'g-', lw=1.5, label='Analytic')
    ax3.plot(sim_41['mean_y'], sim_41['mean_x'], 'ko', ms=2, alpha=0.4, label='MC')
    ax3.set_xlabel(r'$\langle y\rangle$')
    ax3.set_ylabel(r'$\langle x\rangle$')
    ax3.set_title(r'Fig 3(e) overlay: $\tau_J=0.2, \tau_F=0.2, D_r=0$')
    ax3.legend()
    ax3.grid(True, alpha=0.25)
    plt.tight_layout()
    out3 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig3e_overlay_MC.png'
    plt.savefig(out3, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {out3}")

    # MSD overlay (Fig 4 panel a style)
    fig4o, ax4 = plt.subplots(1, 1, figsize=(6, 5))
    mask_a = (t_41a > 0.05) & (msd_41a > 1e-10)
    ax4.loglog(t_41a[mask_a], msd_41a[mask_a], 'r-', lw=2.5, label='Analytic')
    mask_s = (sim_41['t'] > 0.05) & (sim_41['msd'] > 1e-10)
    ax4.loglog(sim_41['t'][mask_s], sim_41['msd'][mask_s], 'ko', ms=2, alpha=0.4, label='MC')
    ax4.set_xlabel(r'$t/\tau_C$')
    ax4.set_ylabel(r'MSD')
    ax4.set_title(r'Fig 4(a) overlay: $\tau_J=0.2, \tau_F=0.2, D_r=0$')
    ax4.legend()
    ax4.grid(True, which='both', alpha=0.25)
    plt.tight_layout()
    out4 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig4a_overlay_MC.png'
    plt.savefig(out4, dpi=150, bbox_inches='tight')
    print(f"Saved: {out4}")

    print("\nAll simulations complete.")
    plt.show()
