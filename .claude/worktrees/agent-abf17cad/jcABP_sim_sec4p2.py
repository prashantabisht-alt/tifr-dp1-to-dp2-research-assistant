#!/usr/bin/env python3
"""
MC verification for Section 4.2: Full model (D_r != 0, omega_0 != 0).
Overlays Euler-Maruyama simulation on analytic MSD (Fig 8).
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
import time as pytime

# ── simulation engine (from jcABP_simulation.py) ──
def run_jcABP_sim(tau_J, tau_F, tau_P, tau_C, D=0.5, v0=1.0, m=1.0,
                  T_max=20.0, dt=0.01, N_traj=10000, save_every=10,
                  seed=42, verbose=True):
    rng = np.random.default_rng(seed)
    N_steps = int(T_max / dt)
    N_save  = N_steps // save_every + 1
    sqrt_dt = np.sqrt(dt)
    omega0 = 1.0 / tau_C if tau_C < 1e10 else 0.0
    D_r    = 1.0 / tau_P if tau_P < 1e10 else 0.0
    inv_tJ = 1.0 / tau_J
    inv_tF2 = 1.0 / tau_F**2
    noise_a = np.sqrt(2.0 * D) * inv_tF2
    noise_r = np.sqrt(2.0 * D_r)

    x = np.zeros(N_traj); vx = np.zeros(N_traj); ax = np.zeros(N_traj)
    y = np.zeros(N_traj); vy = np.zeros(N_traj); ay = np.zeros(N_traj)
    th = np.zeros(N_traj)

    t_out = np.zeros(N_save); mean_x = np.zeros(N_save)
    mean_y = np.zeros(N_save); mean_r2 = np.zeros(N_save)
    save_idx = 1
    t0 = pytime.time()

    for step in range(1, N_steps + 1):
        Wx = rng.standard_normal(N_traj) * sqrt_dt
        Wy = rng.standard_normal(N_traj) * sqrt_dt
        Wr = rng.standard_normal(N_traj) * sqrt_dt
        cos_th = np.cos(th); sin_th = np.sin(th)
        dax = (-inv_tJ*ax - inv_tF2*vx + v0*inv_tF2*cos_th)*dt + noise_a*Wx
        day = (-inv_tJ*ay - inv_tF2*vy + v0*inv_tF2*sin_th)*dt + noise_a*Wy
        x += vx*dt; y += vy*dt
        vx += ax*dt; vy += ay*dt
        ax += dax; ay += day
        th += omega0*dt + noise_r*Wr
        if step % save_every == 0 and save_idx < N_save:
            t_out[save_idx] = step * dt
            mean_x[save_idx] = np.mean(x)
            mean_y[save_idx] = np.mean(y)
            mean_r2[save_idx] = np.mean(x**2 + y**2)
            save_idx += 1

    elapsed = pytime.time() - t0
    if verbose:
        print(f"  Done: {N_steps} steps x {N_traj} traj in {elapsed:.1f}s")
    return dict(t=t_out[:save_idx], mean_x=mean_x[:save_idx],
                mean_y=mean_y[:save_idx], msd=mean_r2[:save_idx])

# ── analytic functions ──
def get_poles(tau_J, tau_F):
    disc = 1.0/tau_F**2 - 1.0/(4.0*tau_J**2)
    return (-1j/(2*tau_J) + np.sqrt(complex(disc)),
            -1j/(2*tau_J) - np.sqrt(complex(disc)))

def G_vec(t_arr, tau_J, tau_F, m=1.0):
    lam = tau_J * m
    o1, o2 = get_poles(tau_J, tau_F)
    denom = lam * (o2 - o1)
    t = np.asarray(t_arr, dtype=float).astype(complex)
    G = ((-(1-np.exp(-1j*o1*t))/o1 + (1-np.exp(-1j*o2*t))/o2) / denom).real
    G[np.asarray(t_arr) <= 0] = 0.0
    return G

def compute_msd_analytic(tau_J, tau_F, tau_P, tau_C, D, v0, T_max, N):
    dt = T_max / N
    t = np.linspace(0.0, T_max, N+1)
    gamma = tau_J / tau_F**2
    G = G_vec(t, tau_J, tau_F)
    K = np.cos(t/tau_C) * np.exp(-t/tau_P)
    h = fftconvolve(G, K)[:N+1] * dt
    msd_act = 2.0 * gamma**2 * v0**2 * np.cumsum(G * h) * dt
    msd_thm = 4.0 * D * gamma**2 * np.cumsum(G**2) * dt
    return t, msd_act + msd_thm

# ── Fig 8 parameters ──
tau_P = 10.0; v0 = 0.05; tau_C = 0.5; D = 0.5
lP = v0 * tau_P  # = 0.5
D_c = D + v0**2 / (2.0*(1.0/tau_P + tau_P/tau_C**2))
print(f"D_c = {D_c:.6f}, 4*D_c = {4*D_c:.6f}\n")

panels = [
    dict(tau_J=0.1, tau_F=0.1, T_max=20.0,  dt=0.001, N_traj=20000,
         save_every=50, N_an=60000, label='(a)', x0=1e-3, x1=10),
    dict(tau_J=5.0, tau_F=1.0, T_max=200.0, dt=0.02,  N_traj=20000,
         save_every=25, N_an=100000, label='(b)', x0=0.01, x1=200),
]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
plt.rcParams.update({'font.size': 11})

for i, (ax, p) in enumerate(zip(axes, panels)):
    tJ, tF = p['tau_J'], p['tau_F']
    print(f"--- Panel {p['label']}: tau_J={tJ}, tau_F={tF} ---")

    # Analytic
    t_an, msd_an = compute_msd_analytic(tJ, tF, tau_P, tau_C, D, v0,
                                         p['T_max'], p['N_an'])
    # Simulation
    sim = run_jcABP_sim(tJ, tF, tau_P, tau_C, D=D, v0=v0,
                         T_max=p['T_max'], dt=p['dt'], N_traj=p['N_traj'],
                         save_every=p['save_every'], seed=42+i*1000)

    # Plot
    mask_a = (t_an > p['x0']*0.5) & (msd_an > 1e-20)
    ax.loglog(t_an[mask_a]/tau_P, msd_an[mask_a]/lP**2, 'r-', lw=2.5,
              label='Analytic')
    mask_s = (sim['t'] > p['x0']*0.5) & (sim['msd'] > 1e-20)
    ax.loglog(sim['t'][mask_s]/tau_P, sim['msd'][mask_s]/lP**2,
              'ko', ms=2, alpha=0.5, label=f"MC (N={p['N_traj']})")

    # Guides
    c5 = D / (5.0 * tF**4)
    t_g5 = np.logspace(np.log10(p['x0']), np.log10(min(1.0, p['x1']*0.05)), 25)
    ax.loglog(t_g5/tau_P, c5*t_g5**5/lP**2, 'b--', lw=1, alpha=0.6, label=r'$t^5$')
    t_g1 = np.logspace(np.log10(max(tau_P*0.3, 1)), np.log10(p['T_max']*0.9), 25)
    ax.loglog(t_g1/tau_P, 4*D_c*t_g1/lP**2, 'b:', lw=1, alpha=0.6,
              label=rf'$4D_c t$')

    ax.set_xlabel(r'$t/\tau_P$')
    ax.set_ylabel(r'$\mathrm{MSD}/l_P^2$')
    ax.set_title(rf"{p['label']} $\tau_J={tJ},\ \tau_F={tF}$")
    ax.set_xlim([p['x0']/tau_P, p['x1']/tau_P])
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.2)

    # Slope check
    i70, i90 = int(0.7*len(t_an)), int(0.9*len(t_an))
    slope_an = (msd_an[i90]-msd_an[i70])/(t_an[i90]-t_an[i70])
    t_s = sim['t']; msd_s = sim['msd']
    i70s, i90s = int(0.7*len(t_s)), int(0.9*len(t_s))
    slope_mc = (msd_s[i90s]-msd_s[i70s])/(t_s[i90s]-t_s[i70s])
    print(f"  Analytic slope = {slope_an:.6f}")
    print(f"  MC slope       = {slope_mc:.6f}")
    print(f"  Expected 4*D_c = {4*D_c:.6f}")
    print(f"  MC/analytic     = {slope_mc/slope_an:.4f}\n")

fig.suptitle(r'Fig 8 — MSD full model: Analytic vs MC', fontsize=13)
plt.tight_layout()
out = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig8_overlay_MC.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved: {out}")
plt.show()
