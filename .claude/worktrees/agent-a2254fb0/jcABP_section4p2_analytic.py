#!/usr/bin/env python3
"""
jcABP_section4p2_analytic.py
=============================
Section 4.2: Full model — finite chirality AND rotational diffusion (D_r != 0).
Reproduces Figs 5, 6, 7 (mean 2D trajectory) and Fig 8 (MSD).

KEY PHYSICS:
  theta(t) = omega_0*t + sqrt(2*D_r)*W(t)  => stochastic orientation.
  <n_hat(t)> = e^{-t/tau_P} (cos(t/tau_C), sin(t/tau_C))    [Eq. 48]
  <n_hat(t').n_hat(t'')> = cos[(t'-t'')/tau_C] e^{-|t'-t''|/tau_P}  [Eq. 5]

  Mean trajectory (Eq. 16):
    <x(t)> = gamma*v0 * conv(G, e^{-t/tau_P} cos(t/tau_C))
    <y(t)> = gamma*v0 * conv(G, e^{-t/tau_P} sin(t/tau_C))

  Two decaying modes (Eq. 39-40):
    A-mode:  freq 1/tau_C,  decay 1/tau_P  ->  circular (logarithmic) spiral
    B-mode:  freq alpha,    decay 1/(2*tau_J)  ->  distorted elliptical spiral
    tau_P >> 2*tau_J => A dominates at large times => spira mirabilis
    tau_P << 2*tau_J => B dominates at large times => Lissajous pattern

  Spiral center (Eq. 38): (x_c, y_c) = (v0/(tau_P*sigma), v0/(tau_C*sigma))
    where sigma = 1/tau_P^2 + 1/tau_C^2

  MSD (Eq. 17):
    K(u) = cos(u/tau_C) * e^{-u/tau_P}    [orientation correlation kernel]
    h(s) = conv(G, K)[:N+1] * dt           [convolution of G with K]
    MSD_act = 2*gamma^2*v0^2 * cumsum(G * h) * dt
    MSD_thm = 4*D*gamma^2   * cumsum(G^2) * dt
    MSD = MSD_act + MSD_thm

  Long-time MSD (Eq. 42):
    MSD -> 4*D_c*t  where  D_c = D + v0^2 / (2*(1/tau_P + tau_P/tau_C^2))

Dimensionless: tau_P as time unit, l_P = v0*tau_P as length unit.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

plt.rcParams.update({'font.size': 10, 'axes.labelsize': 11,
                     'legend.fontsize': 9, 'lines.linewidth': 0.8})

# =============================================================================
# CORE FUNCTIONS
# =============================================================================

def get_poles(tau_J, tau_F):
    disc = 1.0/tau_F**2 - 1.0/(4.0*tau_J**2)
    sq = np.sqrt(complex(disc))
    return (-1j/(2*tau_J) + sq), (-1j/(2*tau_J) - sq)


def G_vec(t_arr, tau_J, tau_F, m=1.0):
    """Green's function G(t), Eq.(14). Handles tau_J < 0 (exploding case)."""
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


def compute_mean_traj(tau_J, tau_F, tau_P, tau_C, v0=1.0, m=1.0,
                      T_max=100.0, N=50000):
    """
    Mean 2D trajectory for full jcABP model (D_r != 0).
    Returns (t, mean_x, mean_y).
    """
    dt    = T_max / N
    t     = np.linspace(0.0, T_max, N+1)
    gamma = tau_J * m / tau_F**2

    G  = G_vec(t, tau_J, tau_F, m)
    Ex = np.exp(-t / tau_P) * np.cos(t / tau_C)
    Ey = np.exp(-t / tau_P) * np.sin(t / tau_C)

    mean_x = gamma * v0 * fftconvolve(G, Ex)[:N+1] * dt
    mean_y = gamma * v0 * fftconvolve(G, Ey)[:N+1] * dt

    return t, mean_x, mean_y


def compute_msd_full(tau_J, tau_F, tau_P, tau_C, D=0.5, v0=1.0, m=1.0,
                     T_max=100.0, N=50000):
    """
    MSD for full jcABP model (Eq. 17).
    K(u) = cos(u/tau_C)*exp(-u/tau_P), h = conv(G,K),
    MSD = 2*gamma^2*v0^2*cumsum(G*h)*dt + 4*D*gamma^2*cumsum(G^2)*dt.
    Returns (t, msd).
    """
    dt    = T_max / N
    t     = np.linspace(0.0, T_max, N+1)
    gamma = tau_J * m / tau_F**2

    G = G_vec(t, tau_J, tau_F, m)
    K = np.cos(t / tau_C) * np.exp(-t / tau_P)

    h       = fftconvolve(G, K)[:N+1] * dt
    msd_act = 2.0 * gamma**2 * v0**2 * np.cumsum(G * h) * dt
    msd_thm = 4.0 * D * gamma**2     * np.cumsum(G**2) * dt

    return t, msd_act + msd_thm


def subsample_logspace(t, y, n_pts=500, t_min=None, t_max=None):
    t_min = t_min or max(t[1], 1e-6)
    t_max = t_max or t[-1]
    t_log = np.logspace(np.log10(t_min), np.log10(t_max), n_pts)
    return t_log, np.interp(t_log, t, y)


# =============================================================================
# FIGURE 5: Damped Lissajous patterns
# =============================================================================
# tau_P = 10, D = 0.5, v0 = 0.05, tau_C = 0.05*tau_P = 0.5
# Row 1 (tau_P ~ 2*tau_J): tau_J = 0.5*tau_P = 5
# Row 2 (tau_P << 2*tau_J): tau_J = 2*tau_P = 20
# tau_F values: 0.001*tau_P, 0.005*tau_P, 0.01*tau_P, 0.2*tau_P
# =============================================================================

print("=== Figure 5: Damped Lissajous patterns ===\n")

tau_P5 = 10.0; D5 = 0.5; v05 = 0.05; tau_C5 = 0.05 * tau_P5  # = 0.5
lP5 = v05 * tau_P5  # = 0.5

# Spiral center (Eq. 38)
sigma5 = 1.0/tau_P5**2 + 1.0/tau_C5**2
xc5 = v05 / (tau_P5 * sigma5)
yc5 = v05 / (tau_C5 * sigma5)
print(f"  Spiral center: ({xc5/lP5:.6f}, {yc5/lP5:.6f}) in l_P units")

tF_vals = [0.001*tau_P5, 0.005*tau_P5, 0.01*tau_P5, 0.2*tau_P5]
rows5 = [
    dict(tau_J=0.5*tau_P5, label_prefix='Row1', row_title=r'$\tau_P \approx 2\tau_J$'),
    dict(tau_J=2.0*tau_P5, label_prefix='Row2', row_title=r'$\tau_P \ll 2\tau_J$'),
]
labels5 = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)']

fig5, axes5 = plt.subplots(2, 4, figsize=(16, 8))
idx = 0

for row_idx, rp in enumerate(rows5):
    tJ = rp['tau_J']
    T_max = max(10*tau_P5, 5*2*abs(tJ))
    for col_idx, tF in enumerate(tF_vals):
        ax = axes5[row_idx, col_idx]
        lab = labels5[idx]; idx += 1
        gamma = tJ / tF**2
        # Adaptive grid: dt must resolve min(tF, tau_C)
        dt_target = min(tF, tau_C5) / 5.0
        N = max(int(T_max / dt_target), 30000)

        print(f"  {lab}: tau_J={tJ:.1f}, tau_F={tF:.4f}, gamma={gamma:.1f}, "
              f"T_max={T_max:.0f}, N={N}")

        t, mx, my = compute_mean_traj(tJ, tF, tau_P5, tau_C5, v0=v05,
                                       T_max=T_max, N=N)
        ax.plot(my/lP5, mx/lP5, 'g-', lw=0.3)
        ax.set_xlabel(r'$\langle y\rangle/l_P$')
        ax.set_ylabel(r'$\langle x\rangle/l_P$')
        ax.set_title(rf"{lab} $\tau_J={tJ:.0f}, \tau_F={tF:.3f}$", fontsize=9)
        ax.grid(True, alpha=0.2)
        ax.plot(yc5/lP5, xc5/lP5, 'r+', ms=6, mew=1.5)

fig5.suptitle(r'Fig 5 — Damped Lissajous ($\tau_P=10, \tau_C=0.5, v_0=0.05$)',
              fontsize=12, y=1.01)
plt.tight_layout()
out5 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig5_Lissajous.png'
plt.savefig(out5, dpi=150, bbox_inches='tight')
print(f"\n  Saved: {out5}\n")


# =============================================================================
# FIGURE 6: Spira mirabilis vs Lissajous at large times
# =============================================================================
# Same global params as Fig 5.
# Row 1 (tau_P >> 2*tau_J): tau_J = 0.1*tau_P = 1
# Row 2 (tau_P << 2*tau_J): tau_J = 2*tau_P = 20
# =============================================================================

print("=== Figure 6: Spira mirabilis vs Lissajous ===\n")

rows6 = [
    dict(tau_J=0.1*tau_P5, label_prefix='Row1', row_title=r'$\tau_P \gg 2\tau_J$'),
    dict(tau_J=2.0*tau_P5, label_prefix='Row2', row_title=r'$\tau_P \ll 2\tau_J$'),
]
labels6 = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)']

fig6, axes6 = plt.subplots(2, 4, figsize=(16, 8))
idx = 0

for row_idx, rp in enumerate(rows6):
    tJ = rp['tau_J']
    T_max = max(10*tau_P5, 5*2*abs(tJ))
    for col_idx, tF in enumerate(tF_vals):
        ax = axes6[row_idx, col_idx]
        lab = labels6[idx]; idx += 1
        gamma = tJ / tF**2
        dt_target = min(tF, tau_C5) / 5.0
        N = max(int(T_max / dt_target), 30000)

        print(f"  {lab}: tau_J={tJ:.1f}, tau_F={tF:.4f}, gamma={gamma:.1f}, "
              f"T_max={T_max:.0f}, N={N}")

        t, mx, my = compute_mean_traj(tJ, tF, tau_P5, tau_C5, v0=v05,
                                       T_max=T_max, N=N)
        ax.plot(my/lP5, mx/lP5, 'g-', lw=0.3)
        ax.set_xlabel(r'$\langle y\rangle/l_P$')
        ax.set_ylabel(r'$\langle x\rangle/l_P$')
        ax.set_title(rf"{lab} $\tau_J={tJ:.1f}, \tau_F={tF:.3f}$", fontsize=9)
        ax.grid(True, alpha=0.2)
        ax.plot(yc5/lP5, xc5/lP5, 'r+', ms=6, mew=1.5)

fig6.suptitle(r'Fig 6 — Spira mirabilis vs Lissajous ($\tau_P=10, \tau_C=0.5, v_0=0.05$)',
              fontsize=12, y=1.01)
plt.tight_layout()
out6 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig6_spira_vs_lissajous.png'
plt.savefig(out6, dpi=150, bbox_inches='tight')
print(f"\n  Saved: {out6}\n")


# =============================================================================
# FIGURE 7: Damped vs exploding spirals (positive/negative tau_J)
# =============================================================================
# tau_P = 0.1, D = 0.5, v0 = 0.05, tau_F = 0.1*tau_P = 0.01
# Row 1: |tau_J| = tau_P = 0.1
# Row 2: |tau_J| = 0.5*tau_P = 0.05
# Columns: (pos, tau_C=0.1*tP), (neg, tau_C=0.1*tP),
#           (pos, tau_C=0.05*tP), (neg, tau_C=0.05*tP)
# =============================================================================

print("=== Figure 7: Damped vs exploding spirals ===\n")

tau_P7 = 0.1; D7 = 0.5; v07 = 0.05; tau_F7 = 0.01
lP7 = v07 * tau_P7  # = 0.005

params_fig7 = [
    # Row 1: |tau_J| = tau_P
    dict(tau_J=0.1,  tau_C=0.01,  label='(a)', kind='Damped spira mirabilis'),
    dict(tau_J=-0.1, tau_C=0.01,  label='(b)', kind='Exploding spira mirabilis'),
    dict(tau_J=0.1,  tau_C=0.005, label='(c)', kind='Damped Lissajous'),
    dict(tau_J=-0.1, tau_C=0.005, label='(d)', kind='Exploding Lissajous'),
    # Row 2: |tau_J| = 0.5*tau_P
    dict(tau_J=0.05,  tau_C=0.01,  label='(e)', kind='Damped spira mirabilis'),
    dict(tau_J=-0.05, tau_C=0.01,  label='(f)', kind='Exploding spira mirabilis'),
    dict(tau_J=0.05,  tau_C=0.005, label='(g)', kind='Damped Lissajous'),
    dict(tau_J=-0.05, tau_C=0.005, label='(h)', kind='Exploding Lissajous'),
]

fig7, axes7 = plt.subplots(2, 4, figsize=(16, 8))
ax_flat7 = axes7.flatten()

for i, p in enumerate(params_fig7):
    ax = ax_flat7[i]
    tJ, tC = p['tau_J'], p['tau_C']
    gamma = tJ / tau_F7**2

    # For negative tau_J: trajectory explodes.  Limit T_max.
    if tJ > 0:
        T_max = max(5*tau_P7, 5*2*tJ)
        N_lim = 60000
    else:
        # Growth rate ~ 1/(2|tau_J|). Want exp(T_max/(2|tJ|)) ~ 1e4 max.
        T_max = min(2*abs(tJ) * np.log(1e4), 10*tau_P7)
        N_lim = 80000

    dt_target = min(abs(tJ), tau_F7, tC) / 10.0
    N = min(max(int(T_max / dt_target), 40000), N_lim)

    print(f"  {p['label']}: tau_J={tJ:+.3f}, tau_C={tC:.3f}, gamma={gamma:.1f}, "
          f"T_max={T_max:.4f}, N={N}, ({p['kind']})")

    t, mx, my = compute_mean_traj(tJ, tau_F7, tau_P7, tC, v0=v07,
                                   T_max=T_max, N=N)

    # Clip extreme values for plotting stability
    clip = 100 * lP7
    mx_c = np.clip(mx, -clip, clip)
    my_c = np.clip(my, -clip, clip)

    color = 'g' if tJ > 0 else 'r'
    ax.plot(my_c/lP7, mx_c/lP7, color + '-', lw=0.3)
    ax.set_xlabel(r'$\langle y\rangle/l_P$')
    ax.set_ylabel(r'$\langle x\rangle/l_P$')
    ax.set_title(rf"{p['label']} $\tau_J={tJ:+.2f}, \tau_C={tC}$", fontsize=8)
    ax.grid(True, alpha=0.2)

    # Mark spiral center
    sigma7 = 1.0/tau_P7**2 + 1.0/tC**2
    yc7 = v07 / (tC * sigma7)
    xc7 = v07 / (tau_P7 * sigma7)
    ax.plot(yc7/lP7, xc7/lP7, 'k+', ms=6, mew=1.5)

fig7.suptitle(r'Fig 7 — Damped vs exploding spirals ($\tau_P=0.1, \tau_F=0.01, v_0=0.05$)',
              fontsize=12, y=1.01)
plt.tight_layout()
out7 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig7_damped_exploding.png'
plt.savefig(out7, dpi=150, bbox_inches='tight')
print(f"\n  Saved: {out7}\n")


# =============================================================================
# FIGURE 8: MSD with chirality and rotational diffusion
# =============================================================================
# tau_P = 10, v0 = 0.05, tau_C = 0.5, D = 0.5
# (a) tau_J = 0.01*tau_P = 0.1,  tau_F = 0.01*tau_P = 0.1
# (b) tau_J = 0.5*tau_P = 5.0,   tau_F = 0.1*tau_P = 1.0
# =============================================================================

print("=== Figure 8: MSD (full model) ===\n")

tau_P8 = 10.0; v08 = 0.05; tau_C8 = 0.5; D8 = 0.5
lP8 = v08 * tau_P8  # = 0.5

# Effective diffusion constant (Eq. 42)
D_c = D8 + v08**2 / (2.0 * (1.0/tau_P8 + tau_P8/tau_C8**2))
print(f"  D_eff = {D_c:.6f},  4*D_eff = {4*D_c:.6f}")
print(f"  4*D = {4*D8:.2f}")

params_fig8 = [
    dict(tau_J=0.1,  tau_F=0.1, T_max=20.0,  N=60000,  label='(a)', x0=0.001, x1=10),
    dict(tau_J=5.0,  tau_F=1.0, T_max=200.0, N=100000, label='(b)', x0=0.01,  x1=200),
]

fig8, axes8 = plt.subplots(1, 2, figsize=(11, 5))

for ax, p in zip(axes8, params_fig8):
    tJ, tF = p['tau_J'], p['tau_F']
    gamma = tJ / tF**2
    print(f"  {p['label']}: tau_J={tJ}, tau_F={tF}, gamma={gamma:.4f}")

    t, msd = compute_msd_full(tJ, tF, tau_P8, tau_C8, D=D8, v0=v08,
                               T_max=p['T_max'], N=p['N'])

    t_log, msd_log = subsample_logspace(t, msd, t_min=p['x0']*0.5, t_max=p['T_max'])
    mask = msd_log > 1e-20 * lP8**2

    ax.loglog(t_log[mask]/tau_P8, msd_log[mask]/lP8**2, 'r-', lw=2.5,
              label='MSD (Eq. 17)')

    # Short-time guide: t^5
    c5 = D8 / (5.0 * tF**4)
    t_g5 = np.logspace(np.log10(p['x0']), np.log10(min(1.0, p['x1']*0.1)), 30)
    ax.loglog(t_g5/tau_P8, c5*t_g5**5/lP8**2, 'k--', lw=1.2, label=r'$\propto t^5$')

    # Long-time guide: 4*D_c*t
    t_g1 = np.logspace(np.log10(max(tau_P8*0.5, 1.0)), np.log10(p['T_max']*0.9), 25)
    ax.loglog(t_g1/tau_P8, 4*D_c*t_g1/lP8**2, 'k:', lw=1.2,
              label=rf'$4D_c t\ (D_c={D_c:.4f})$')

    ax.set_xlabel(r'$t/\tau_P$')
    ax.set_ylabel(r'$\mathrm{MSD}(t)/l_P^2$')
    ax.set_title(rf"{p['label']} $\tau_J={tJ}, \tau_F={tF}$", fontsize=10)
    ax.set_xlim([p['x0']/tau_P8, p['x1']/tau_P8])
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.25)

fig8.suptitle(r'Fig 8 — MSD, full model ($\tau_P=10, \tau_C=0.5, v_0=0.05, D=0.5$)',
              fontsize=12)
plt.tight_layout()
out8 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig8_MSD_full.png'
plt.savefig(out8, dpi=150, bbox_inches='tight')
print(f"\n  Saved: {out8}\n")


# =============================================================================
# SANITY CHECKS
# =============================================================================
print("=== Sanity checks ===\n")

# Check 1: long-time MSD slope vs 4*D_c
print(f"[MSD] Long-time slope (expect 4*D_c = {4*D_c:.6f}):")
for p in params_fig8:
    tJ, tF = p['tau_J'], p['tau_F']
    t, msd = compute_msd_full(tJ, tF, tau_P8, tau_C8, D=D8, v0=v08,
                               T_max=p['T_max'], N=p['N'])
    i70 = int(0.70 * len(t))
    i90 = int(0.90 * len(t))
    slope = (msd[i90] - msd[i70]) / (t[i90] - t[i70])
    print(f"  {p['label']}: slope = {slope:.6f}, ratio = {slope/(4*D_c):.4f}")

# Check 2: mean trajectory converges to spiral center
print(f"\n[Traj] Spiral center convergence (expect ({xc5/lP5:.4f}, {yc5/lP5:.4f}) in l_P units):")
for tJ_val, tF_val in [(5.0, 0.01), (1.0, 0.1), (20.0, 2.0)]:
    T_max = max(10*tau_P5, 5*2*abs(tJ_val)) + 50
    dt_t = min(tF_val, tau_C5) / 5.0
    N_t = max(int(T_max / dt_t), 40000)
    t, mx, my = compute_mean_traj(tJ_val, tF_val, tau_P5, tau_C5, v0=v05,
                                   T_max=T_max, N=N_t)
    # Average last bit where both modes have decayed
    n_end = min(len(t)//20, 5000)
    mx_end = np.mean(mx[-n_end:]) / lP5
    my_end = np.mean(my[-n_end:]) / lP5
    print(f"  tau_J={tJ_val:5.1f}, tau_F={tF_val:.3f}: "
          f"x_c/lP={mx_end:.6f}, y_c/lP={my_end:.6f}")

# Check 3: short-time MSD scaling t^5
print(f"\n[MSD] Short-time: MSD ~ D*t^5/(5*tau_F^4)")
for p in params_fig8:
    tJ, tF = p['tau_J'], p['tau_F']
    t, msd = compute_msd_full(tJ, tF, tau_P8, tau_C8, D=D8, v0=v08,
                               T_max=p['T_max'], N=p['N'])
    c5 = D8 / (5.0 * tF**4)
    for t_check in [0.001, 0.005, 0.01]:
        idx = int(t_check / (p['T_max']/p['N']))
        if 0 < idx < len(t):
            msd_num  = msd[idx]
            msd_anal = c5 * t[idx]**5
            ratio = msd_num / msd_anal if msd_anal > 0 else float('nan')
            print(f"  {p['label']}: t={t[idx]:.4f}: ratio={ratio:.4f}")

print("\nDone.")
plt.show()
