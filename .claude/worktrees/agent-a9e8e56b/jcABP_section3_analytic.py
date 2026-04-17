#!/usr/bin/env python3
"""
jcABP_section3_analytic.py
===========================
Analytic evaluation of mean displacement and MSD for the jerky chiral ABP
in the VANISHING CHIRALITY limit (Section 3 of Jose & Löwen, arXiv 2025).

Reproduces:
  Fig 1 — Mean displacement <x(t)>/l_P  vs  t/tau_P  (log-log)
  Fig 2 — MSD(t)/l_P^2               vs  t/tau_P  (log-log)

KEY PHYSICS (Section 3, omega_0 = 0):
  Full EOM:  lambda * x''' + m * x'' + gamma * x' = gamma*v0*n_x + sqrt(2D)*gamma*eta
  Orientation: theta_dot = sqrt(2*D_r)*eta_r  (no chirality => mean n_y = 0 always)
  Orientation correlation: <n_x(t')> = exp(-t'/tau_P)

  G(t) satisfies:  lambda*G''' + m*G'' + gamma*G' = delta(t)
  G(t) -> 1/gamma as t -> inf  [crucial: DC response, NOT zero!]

  Mean displacement:
    <x(t)> = gamma*v0 * int_0^t G(t-t') exp(-t'/tau_P) dt'

  MSD (Eq. 23):
    MSD(t) = 2*gamma^2*v0^2 * int_0^t G(s)*h(s) ds     [activity]
           + 4*D*gamma^2    * int_0^t G^2(s) ds          [thermal]
    where h(s) = int_0^s G(s-u) exp(-u/tau_P) du  [same convolution as <x>/(gamma*v0)]

  Long-time checks:
    G(t) -> 1/gamma  =>  h(t) -> tau_P/gamma
    MSD_act -> 2*gamma^2*v0^2*(tau_P/gamma^2)*t = 2*v0^2*tau_P*t
    MSD_therm -> 4*D*gamma^2*(1/gamma^2)*t = 4*D*t
    Total: (4D + 2*v0^2*tau_P)*t  ✓  [Eq. 25]

NUMERICAL APPROACH:
  Use linear time grid + scipy.signal.fftconvolve.
  Avoids catastrophic cancellation that plagues the split
  exp(-t/tau_P) * cumtrapz[G(s)*exp(s/tau_P)] approach for large t.
  Then subsample to log-spaced points for log-log plots.

Dimensionless units: tau_P=1, v0=1, m=1  =>  lambda=tau_J, gamma=tau_J/tau_F^2.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

plt.rcParams.update({'font.size': 12, 'axes.labelsize': 13,
                     'legend.fontsize': 11, 'lines.linewidth': 2})

# =============================================================================
# 1.  GREEN'S FUNCTION  G(t)  — Eq. (14)  [VECTORIZED]
# =============================================================================

def get_poles(tau_J, tau_F):
    """
    Eq. (12): omega_{1,2} = -i/(2*tau_J) ± sqrt(1/tau_F^2 - 1/(4*tau_J^2))
    discriminant > 0: real sqrt  => oscillations in G (case tau_F < 2*tau_J)
    discriminant < 0: imag sqrt  => pure exponential decay (case tau_F > 2*tau_J)
    """
    disc = 1.0/tau_F**2 - 1.0/(4.0*tau_J**2)
    sq   = np.sqrt(complex(disc))
    return (-1j/(2*tau_J) + sq), (-1j/(2*tau_J) - sq)


def G_vec(t_arr, tau_J, tau_F, m=1.0):
    """
    G(t) from Eq. (14), vectorised over t_arr.

    G(t) = Theta(t)/[lambda*(omega2-omega1)] * sum_{j=1,2} (-1)^j*(1-exp(-i*omega_j*t))/omega_j

    lambda = m*tau_J  (m=1 in dimensionless units => lambda = tau_J)

    Short-time:  G(t) ~ t^2/(2*lambda)
    Long-time:   G(t) -> 1/gamma = tau_F^2/tau_J
    """
    lam = tau_J * m
    o1, o2 = get_poles(tau_J, tau_F)
    denom = lam * (o2 - o1)

    t = np.asarray(t_arr, dtype=float)
    tc = t.astype(complex)
    term1 = (1 - np.exp(-1j*o1*tc)) / o1   # j=1, sign = -1
    term2 = (1 - np.exp(-1j*o2*tc)) / o2   # j=2, sign = +1
    G = ((-term1 + term2) / denom).real
    G[t <= 0] = 0.0
    return G


# =============================================================================
# 2.  MEAN DISPLACEMENT AND MSD via FFT CONVOLUTION on a linear grid
# =============================================================================
#
#  For a linear grid t_k = k*dt (k=0,...,N):
#
#  <x(t_k)>  = gamma*v0 * (G ⊛ E)[k] * dt       where E[j] = exp(-t_j/tau_P)
#
#  h(t_k)    = (G ⊛ E)[k] * dt                  [same convolution!]
#
#  MSD_act(t_k) = 2*gamma^2*v0^2 * cumsum(G*h)[k] * dt
#  MSD_thm(t_k) = 4*D*gamma^2    * cumsum(G^2)[k] * dt
#
#  fftconvolve(G, E)[:N+1] gives the one-sided causal convolution.
# =============================================================================

def compute_section3(tau_J, tau_F, tau_P=1.0, D=0.5, v0=1.0, m=1.0,
                     T_max=30.0, N=15000):
    """
    Compute <x(t)> and MSD(t) on a linear grid of N+1 points in [0, T_max].
    Returns (t, mean_x, msd).
    """
    dt    = T_max / N
    t     = np.linspace(0.0, T_max, N+1)
    gamma = tau_J * m / tau_F**2

    G = G_vec(t, tau_J, tau_F, m)
    E = np.exp(-t / tau_P)                           # exp(-t/tau_P)

    # One-sided convolution: (G ⊛ E)[k] = sum_{j=0}^{k} G[k-j] * E[j]
    conv_GE = fftconvolve(G, E)[:N+1] * dt           # h(t_k)

    mean_x   = gamma * v0 * conv_GE                  # <x(t)>
    msd_act  = 2.0 * gamma**2 * v0**2 * np.cumsum(G * conv_GE) * dt
    msd_thm  = 4.0 * D * gamma**2    * np.cumsum(G**2) * dt
    msd      = msd_act + msd_thm

    return t, mean_x, msd


def subsample_logspace(t, y, n_pts=400, t_min=None, t_max=None):
    """Return (t, y) subsampled at log-spaced t values for log-log plots."""
    t_min = t_min or max(t[1], 1e-4)
    t_max = t_max or t[-1]
    t_log = np.logspace(np.log10(t_min), np.log10(t_max), n_pts)
    y_log = np.interp(t_log, t, y)
    return t_log, y_log


# =============================================================================
# 3.  FIGURE 1 — Mean displacement  <x(t)>/l_P  vs  t/tau_P
# =============================================================================
#  Panel (a): tau_J=0.2*tau_P, tau_F=0.28*tau_P  (smooth, alpha~2.55/tau_P)
#  Panel (b): tau_J=20*tau_P,  tau_F=1.0*tau_P   (oscillatory, alpha~1/tau_P)
# =============================================================================

tau_P, v0, m = 1.0, 1.0, 1.0
lP = v0 * tau_P

print("=== Figure 1: Mean displacement ===")

params_fig1 = [
    dict(tau_J=0.2,  tau_F=0.28, T_max=20.0,  N=20000, label='(a)', x0=1e-2, x1=15),
    dict(tau_J=20.0, tau_F=1.0,  T_max=200.0, N=60000, label='(b)', x0=1e-1, x1=150),
]

fig1, axes = plt.subplots(1, 2, figsize=(10, 4.5))

for ax, p in zip(axes, params_fig1):
    tJ, tF = p['tau_J'], p['tau_F']
    disc = 1/tF**2 - 1/(4*tJ**2)
    alpha = np.sqrt(abs(disc))
    print(f"  {p['label']}: tau_J={tJ}, tau_F={tF}")
    print(f"    disc={disc:.4f}, |alpha|={alpha:.4f}/tau_P, "
          f"{'oscillatory' if disc>0 else 'smooth'}")

    t_lin, mx_lin, _ = compute_section3(tJ, tF, tau_P, v0=v0, D=0.5,
                                         T_max=p['T_max'], N=p['N'])
    t_plot, mx_plot = subsample_logspace(t_lin, mx_lin,
                                          t_min=p['x0']*0.5, t_max=p['T_max'])

    mask = mx_plot > 1e-8 * lP
    ax.loglog(t_plot[mask]/tau_P, mx_plot[mask]/lP, 'r-', lw=2.5,
              label=r'$\langle x(t)\rangle$ (Eq. 20)')

    # Short-time t^3 guide: <x> ~ v0*t^3/(6*tau_F^2)
    c3 = v0 / (6.0 * tF**2 * lP)
    t_g = np.logspace(np.log10(p['x0']), -0.5, 30)
    ax.loglog(t_g, c3 * t_g**3, 'k--', lw=1.2, label=r'$\propto t^3$')

    ax.axhline(1.0, color='gray', ls=':', lw=1.2, label=r'$l_P$')
    ax.set_xlabel(r'$t/\tau_P$')
    ax.set_ylabel(r'$\langle x(t)\rangle / l_P$')
    ax.set_title(rf"{p['label']} $\tau_J={tJ}\tau_P$, $\tau_F={tF}\tau_P$")
    ax.set_xlim([p['x0'], p['x1']])
    ax.set_ylim([1e-5, 2])
    ax.legend(fontsize=10)
    ax.grid(True, which='both', alpha=0.25)

plt.suptitle('Fig 1 — Mean displacement, no chirality (Sec. 3)', fontsize=12)
plt.tight_layout()
out1 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig1_mean_disp.png'
plt.savefig(out1, dpi=150, bbox_inches='tight')
print(f"  Saved: {out1}\n")

# =============================================================================
# 4.  FIGURE 2 — MSD(t)/l_P^2  vs  t/tau_P
# =============================================================================
#  D = 0.5 fixed.
#  Panel (a): tau_J=0.2,  tau_F=0.28  (same as Fig1a)
#  Panel (b): tau_J=20,   tau_F=0.2   [NOTE: tau_F=0.2, not 1.0]
# =============================================================================

D = 0.5

params_fig2 = [
    dict(tau_J=0.2,  tau_F=0.28, T_max=40.0,  N=20000, label='(a)', x0=0.05, x1=25),
    dict(tau_J=20.0, tau_F=0.2,  T_max=250.0, N=80000, label='(b)', x0=0.05, x1=200),
]

print("=== Figure 2: MSD ===")
print(f"  D={D}, expected long-time slope = 4D + 2*v0^2*tau_P = {4*D + 2*v0**2*tau_P:.2f}")

fig2, axes = plt.subplots(1, 2, figsize=(10, 4.5))

for ax, p in zip(axes, params_fig2):
    tJ, tF = p['tau_J'], p['tau_F']
    gamma = tJ / tF**2
    G_inf = 1.0/gamma         # G(t->inf) = 1/gamma = tau_F^2/tau_J
    print(f"  {p['label']}: tau_J={tJ}, tau_F={tF}, gamma={gamma:.4f}, G(inf)={G_inf:.4f}")

    t_lin, _, msd_lin = compute_section3(tJ, tF, tau_P, D=D, v0=v0,
                                          T_max=p['T_max'], N=p['N'])
    t_plot, msd_plot = subsample_logspace(t_lin, msd_lin,
                                           t_min=p['x0']*0.5, t_max=p['T_max'])

    mask = msd_plot > 1e-14 * lP**2
    ax.loglog(t_plot[mask]/tau_P, msd_plot[mask]/lP**2, 'r-', lw=2.5,
              label='MSD (Eq. 24)')

    # Guide: t^5 (short time; thermal dominates: MSD ~ D*t^5/(5*tau_F^4))
    c5 = D / (5.0 * tF**4 * lP**2)
    t_g5 = np.logspace(np.log10(p['x0']), 0.0, 30)
    ax.loglog(t_g5, c5 * t_g5**5, 'k--', lw=1.2, label=r'$\propto t^5$')

    # Guide: t^1 (long time: slope = 4D + 2*v0^2*tau_P)
    slope_lt = 4.0*D + 2.0*v0**2*tau_P
    t_g1 = np.logspace(0.8, np.log10(p['x1']*0.8), 20)
    ax.loglog(t_g1, slope_lt * t_g1 / lP**2 * 0.5, 'k:', lw=1.2, label=r'$\propto t$')

    ax.set_xlabel(r'$t/\tau_P$')
    ax.set_ylabel(r'$\mathrm{MSD}(t)/l_P^2$')
    ax.set_title(rf"{p['label']} $\tau_J={tJ}\tau_P$, $\tau_F={tF}\tau_P$, $D={D}$")
    ax.set_xlim([p['x0'], p['x1']])
    ax.legend(fontsize=10)
    ax.grid(True, which='both', alpha=0.25)

plt.suptitle('Fig 2 — MSD, no chirality (Sec. 3)', fontsize=12)
plt.tight_layout()
out2 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig2_MSD.png'
plt.savefig(out2, dpi=150, bbox_inches='tight')
print(f"  Saved: {out2}\n")

# =============================================================================
# 5.  SANITY CHECKS
# =============================================================================
print("=== Sanity checks ===")

# G(t) short-time: G(t) ~ t^2/(2*tau_J)
print("\n[G(t)] Short-time: G(t) ~ t^2/(2*tau_J)")
for tJ, tF in [(0.2, 0.28), (20.0, 1.0)]:
    t_test = np.array([0.001, 0.005, 0.01])
    G_num  = G_vec(t_test, tJ, tF)
    G_anal = t_test**2 / (2.0 * tJ)
    for tt, gn, ga in zip(t_test, G_num, G_anal):
        print(f"  tJ={tJ:5.1f}, t={tt:.3f}: ratio={gn/ga:.6f}")

# G(t) long-time: G -> 1/gamma = tau_F^2/tau_J
print("\n[G(t)] Long-time: G(t) -> tau_F^2/tau_J")
for tJ, tF in [(0.2, 0.28), (20.0, 1.0)]:
    G_inf = tF**2/tJ
    t_check = np.array([5*tJ, 10*tJ, 20*tJ])
    G_num = G_vec(t_check, tJ, tF)
    for tc, gn in zip(t_check, G_num):
        print(f"  tJ={tJ:5.1f}, t={tc:7.2f}: G_num={gn:.6f}, G_inf={G_inf:.6f}, ratio={gn/G_inf:.6f}")

# Long-time MSD slope
print("\n[MSD] Long-time slope vs expected")
expected_slope = 4.0*D + 2.0*v0**2*tau_P
for p in params_fig2:
    tJ, tF = p['tau_J'], p['tau_F']
    t_lin, _, msd_lin = compute_section3(tJ, tF, tau_P, D=D, v0=v0,
                                          T_max=p['T_max'], N=p['N'])
    # Estimate slope from last 10% of time
    i80 = int(0.80 * len(t_lin))
    i90 = int(0.90 * len(t_lin))
    slope = (msd_lin[i90] - msd_lin[i80]) / (t_lin[i90] - t_lin[i80])
    print(f"  {p['label']}: tJ={tJ:5.1f}, tF={tF:.3f}: "
          f"numerical slope = {slope:.4f}, expected = {expected_slope:.4f}, "
          f"ratio = {slope/expected_slope:.4f}")

# Mean displacement long-time: <x> -> l_P = v0*tau_P
print("\n[Mean disp] Long-time: <x>/l_P -> 1")
for p in params_fig1:
    tJ, tF = p['tau_J'], p['tau_F']
    t_lin, mx_lin, _ = compute_section3(tJ, tF, tau_P, v0=v0, D=0.5,
                                         T_max=p['T_max'], N=p['N'])
    mx_end = mx_lin[-1] / lP
    print(f"  {p['label']}: tJ={tJ:5.1f}: <x(T_max)>/l_P = {mx_end:.6f}  (expect 1.0)")

print("\nDone. Figures saved.")
plt.show()
