#!/usr/bin/env python3
"""
jcABP_section4p1_analytic.py
=============================
Section 4.1: Finite chirality, D_r = 0 (zero rotational noise).
Reproduces Fig 3 (mean 2D trajectory) and Fig 4 (MSD) of Jose & Löwen arXiv 2025.

KEY PHYSICS (Section 4.1, D_r = 0  =>  tau_P -> inf):
  Orientation is deterministic:  theta(t) = omega_0 * t  =>  n_hat = (cos(t/tau_C), sin(t/tau_C))
  No orientational noise => activity force is fully deterministic.

  Mean trajectory (Eq. 26):
    <r(t)> = gamma*v0 * int_0^t G(t-t') n_hat(t') dt'
    <x(t)> = gamma*v0 * [G ⊛ cos(t/tau_C)]
    <y(t)> = gamma*v0 * [G ⊛ sin(t/tau_C)]

  Long-time: circular orbit, center (0, v0*tau_C), radius r (Eq. 32), freq 1/tau_C.
    r = v0*tau_C / sqrt(tau_F^4/(tau_J^2 * tau_C^2) + (tau_F^2/tau_C^2 - 1)^2)
    When tau_F^2/tau_C^2 = 1  =>  denominator -> tau_F^2/(tau_J*tau_C)  =>  r -> v0*tau_J
    This is the RESONANCE condition: alpha ≈ 1/tau_C = omega_0.

  MSD (Eq. 34, with D_r = 0  =>  no orientational averaging):
    MSD(t) = |<r(t)>|^2  +  4*D*gamma^2 * int_0^t G^2(s) ds
    Activity part: deterministic, bounded (circular orbit at long times)
    Thermal part: -> 4Dt (same as Section 3 thermal piece)

  Limiting behaviors (Eq. 35):
    MSD(t) ->  D*t^5/(5*tau_F^4) + v0^2*t^6/(36*tau_F^4) + ...   [short time]
    MSD(t) ->  4Dt                                                  [long time]

NUMERICAL APPROACH:
  FFT convolution on a uniform time grid:
    <x(t)> = gamma*v0 * fftconvolve(G, cos(t/tau_C))[:N+1] * dt
    <y(t)> = gamma*v0 * fftconvolve(G, sin(t/tau_C))[:N+1] * dt
  Then MSD = <x>^2 + <y>^2  +  4*D*gamma^2 * cumsum(G^2)*dt.

Dimensionless units: tau_C = 1, v0 = 1, m = 1  =>  lambda = tau_J, gamma = tau_J/tau_F^2.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12,
                     'legend.fontsize': 10, 'lines.linewidth': 1.5})

# =============================================================================
# 1.  GREEN'S FUNCTION  G(t)  — same as Section 3 code
# =============================================================================

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


# =============================================================================
# 2.  MEAN TRAJECTORY AND MSD  (Section 4.1: D_r = 0)
# =============================================================================

def compute_section4p1(tau_J, tau_F, tau_C=1.0, D=0.5, v0=1.0, m=1.0,
                       T_max=100.0, N=50000):
    """
    Mean 2D trajectory and MSD for D_r = 0 (deterministic orientation).
    Returns (t, mean_x, mean_y, msd).
    """
    dt    = T_max / N
    t     = np.linspace(0.0, T_max, N+1)
    gamma = tau_J * m / tau_F**2

    G = G_vec(t, tau_J, tau_F, m)
    C = np.cos(t / tau_C)
    S = np.sin(t / tau_C)

    # Mean trajectory:  <x> = gamma*v0 * (G ⊛ cos), <y> = gamma*v0 * (G ⊛ sin)
    conv_GC = fftconvolve(G, C)[:N+1] * dt
    conv_GS = fftconvolve(G, S)[:N+1] * dt

    mean_x = gamma * v0 * conv_GC
    mean_y = gamma * v0 * conv_GS

    # MSD = |<r>|^2 + thermal
    mean_r_sq = mean_x**2 + mean_y**2
    msd_thm   = 4.0 * D * gamma**2 * np.cumsum(G**2) * dt
    msd        = mean_r_sq + msd_thm

    return t, mean_x, mean_y, msd


def orbit_radius(tau_J, tau_F, tau_C=1.0, v0=1.0):
    """Eq. (32): long-time circular orbit radius."""
    num = v0 * tau_C
    den = np.sqrt(tau_F**4 / (tau_J**2 * tau_C**2) + (tau_F**2/tau_C**2 - 1)**2)
    return num / den


def subsample_logspace(t, y, n_pts=500, t_min=None, t_max=None):
    t_min = t_min or max(t[1], 1e-4)
    t_max = t_max or t[-1]
    t_log = np.logspace(np.log10(t_min), np.log10(t_max), n_pts)
    return t_log, np.interp(t_log, t, y)


# =============================================================================
# 3.  FIGURE 3 — Mean 2D trajectory  <x(t)> vs <y(t)>
# =============================================================================
#  Axes:  x-axis = <y(t)>/(v0*tau_C),  y-axis = <x(t)>/(v0*tau_C)
#  5 panels with different timescale ratios.
# =============================================================================

tau_C, v0, m = 1.0, 1.0, 1.0

print("=== Figure 3: Mean 2D trajectory (D_r = 0) ===\n")

params_fig3 = [
    dict(tau_J=1.0,  tau_F=2*np.sqrt(2), T_max=80,   N=50000,  label='(a)'),
    dict(tau_J=6.0,  tau_F=2*np.sqrt(2), T_max=200,  N=80000,  label='(b)'),
    dict(tau_J=10.0, tau_F=2.0,          T_max=250,  N=100000, label='(c)'),
    dict(tau_J=10.0, tau_F=1.0,          T_max=300,  N=120000, label='(d)'),
    dict(tau_J=0.2,  tau_F=0.2,          T_max=25,   N=30000,  label='(e)'),
]

# Layout: 3 rows x 2 cols, panel (d) larger
fig3 = plt.figure(figsize=(12, 14))
gs = fig3.add_gridspec(3, 2, hspace=0.35, wspace=0.3)

ax_map = {
    '(a)': fig3.add_subplot(gs[0, 0]),
    '(b)': fig3.add_subplot(gs[1, 0]),
    '(c)': fig3.add_subplot(gs[2, 0]),
    '(d)': fig3.add_subplot(gs[0:2, 1]),   # tall panel for resonance
    '(e)': fig3.add_subplot(gs[2, 1]),
}

for p in params_fig3:
    tJ, tF = p['tau_J'], p['tau_F']
    r_orb  = orbit_radius(tJ, tF, tau_C, v0)
    disc   = 1.0/tF**2 - 1.0/(4.0*tJ**2)
    alpha  = np.sqrt(abs(disc))
    kind   = 'oscillatory' if disc > 0 else 'overdamped'

    print(f"  {p['label']}: tau_J={tJ}, tau_F={tF:.4f}")
    print(f"    disc={disc:.4f}, alpha={alpha:.4f}/tau_C ({kind}), r_orbit={r_orb:.4f}")

    t, mx, my, _ = compute_section4p1(tJ, tF, tau_C, D=0.5, v0=v0,
                                       T_max=p['T_max'], N=p['N'])

    ax = ax_map[p['label']]
    ax.plot(my/(v0*tau_C), mx/(v0*tau_C), 'g-', lw=0.4)
    ax.set_xlabel(r'$\langle y(t)\rangle / (v_0\tau_C)$')
    ax.set_ylabel(r'$\langle x(t)\rangle / (v_0\tau_C)$')
    ax.set_title(rf"{p['label']} $\tau_J={tJ}\tau_C,\ \tau_F={tF:.2f}\tau_C$",
                 fontsize=10)
    ax.grid(True, alpha=0.25)

    # Mark long-time orbit center
    ax.plot(v0*tau_C/(v0*tau_C), 0, 'r+', ms=8, mew=1.5, label='orbit center')

fig3.suptitle('Fig 3 — Mean trajectory, $D_r = 0$ (Sec. 4.1)', fontsize=13, y=0.99)
out3 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig3_mean_traj_Dr0.png'
plt.savefig(out3, dpi=150, bbox_inches='tight')
print(f"\n  Saved: {out3}\n")


# =============================================================================
# 4.  FIGURE 4 — MSD(t)/(v0*tau_C)^2  vs  t/tau_C   (log-log)
# =============================================================================
#  D = 0.5 fixed.
#  (a) tau_J = 0.2*tau_C, tau_F = 0.2*tau_C
#  (b) tau_J = 10*tau_C,  tau_F = 0.4*tau_C
# =============================================================================

D = 0.5
print("=== Figure 4: MSD (D_r = 0) ===")
print(f"  D={D}, expected long-time slope = 4D = {4*D:.2f}\n")

params_fig4 = [
    dict(tau_J=0.2,  tau_F=0.2, T_max=200,  N=80000,  label='(a)',
         x0=0.1, x1=200),
    dict(tau_J=10.0, tau_F=0.4, T_max=200,  N=100000, label='(b)',
         x0=0.1, x1=200),
]

fig4, axes4 = plt.subplots(1, 2, figsize=(11, 5))

for ax, p in zip(axes4, params_fig4):
    tJ, tF = p['tau_J'], p['tau_F']
    gamma  = tJ / tF**2
    r_orb  = orbit_radius(tJ, tF, tau_C, v0)

    print(f"  {p['label']}: tau_J={tJ}, tau_F={tF}, gamma={gamma:.4f}, r_orbit={r_orb:.4f}")

    t, mx, my, msd = compute_section4p1(tJ, tF, tau_C, D=D, v0=v0,
                                          T_max=p['T_max'], N=p['N'])

    lC = v0 * tau_C   # length unit
    t_log, msd_log = subsample_logspace(t, msd, t_min=p['x0']*0.5, t_max=p['T_max'])

    mask = msd_log > 1e-14 * lC**2
    ax.loglog(t_log[mask]/tau_C, msd_log[mask]/lC**2, 'r-', lw=2.5,
              label='MSD (Eq. 34)')

    # Guide: t^5 (short time: thermal dominates: MSD ~ D*t^5/(5*tau_F^4))
    c5 = D / (5.0 * tF**4)
    t_g5 = np.logspace(np.log10(p['x0']), 0.0, 30)
    ax.loglog(t_g5, c5 * t_g5**5 / lC**2, 'k--', lw=1.2, label=r'$\propto t^5$')

    # Guide: t (long time: MSD -> 4Dt)
    t_g1 = np.logspace(0.5, np.log10(p['x1']*0.8), 25)
    # Offset to visually match the actual curve
    ax.loglog(t_g1, 4*D * t_g1 / lC**2, 'k:', lw=1.2, label=r'$\propto t$')

    ax.set_xlabel(r'$t/\tau_C$')
    ax.set_ylabel(r'$\mathrm{MSD}(t) / (v_0 \tau_C)^2$')
    ax.set_title(rf"{p['label']} $\tau_J={tJ}\tau_C$, $\tau_F={tF}\tau_C$, $D={D}$",
                 fontsize=10)
    ax.set_xlim([p['x0'], p['x1']])
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.25)

fig4.suptitle('Fig 4 — MSD, $D_r = 0$ (Sec. 4.1)', fontsize=12)
plt.tight_layout()
out4 = '/sessions/youthful-inspiring-planck/mnt/TIFR DP1 to DP2 Research Assistant/fig4_MSD_Dr0.png'
plt.savefig(out4, dpi=150, bbox_inches='tight')
print(f"\n  Saved: {out4}\n")


# =============================================================================
# 5.  SANITY CHECKS
# =============================================================================
print("=== Sanity checks ===\n")

# Check 1: long-time orbit center = (0, v0*tau_C)
print("[Trajectory] Long-time orbit center (expect (0, 1)):")
for p in params_fig3:
    tJ, tF = p['tau_J'], p['tau_F']
    t, mx, my, _ = compute_section4p1(tJ, tF, tau_C, D=0.5, v0=v0,
                                       T_max=p['T_max'], N=p['N'])
    # Average over last full orbit period (2*pi*tau_C)
    n_period = int(2*np.pi*tau_C / (p['T_max']/p['N']))
    if n_period < 10:
        n_period = p['N'] // 10
    mx_avg = np.mean(mx[-n_period:]) / (v0*tau_C)
    my_avg = np.mean(my[-n_period:]) / (v0*tau_C)
    print(f"  {p['label']}: <x>_avg/(v0*tC) = {mx_avg:.4f}  (expect 0),  "
          f"<y>_avg/(v0*tC) = {my_avg:.4f}  (expect 1)")

# Check 2: long-time orbit radius
print(f"\n[Trajectory] Long-time orbit radius:")
for p in params_fig3:
    tJ, tF = p['tau_J'], p['tau_F']
    r_theory = orbit_radius(tJ, tF, tau_C, v0)
    t, mx, my, _ = compute_section4p1(tJ, tF, tau_C, D=0.5, v0=v0,
                                       T_max=p['T_max'], N=p['N'])
    n_period = int(2*np.pi*tau_C / (p['T_max']/p['N']))
    if n_period < 10:
        n_period = p['N'] // 10
    # Amplitude of oscillation in the last period
    mx_tail = mx[-n_period:]
    my_tail = my[-n_period:]
    mx_amp = (np.max(mx_tail) - np.min(mx_tail)) / 2 / (v0*tau_C)
    my_amp = (np.max(my_tail) - np.min(my_tail)) / 2 / (v0*tau_C)
    r_num  = (mx_amp + my_amp) / 2   # rough estimate
    print(f"  {p['label']}: r_theory={r_theory:.4f}, r_numerical≈{r_num:.4f}, "
          f"ratio={r_num/r_theory:.4f}")

# Check 3: MSD long-time slope
print(f"\n[MSD] Long-time slope (expect 4D = {4*D:.2f}):")
for p in params_fig4:
    tJ, tF = p['tau_J'], p['tau_F']
    t, mx, my, msd = compute_section4p1(tJ, tF, tau_C, D=D, v0=v0,
                                          T_max=p['T_max'], N=p['N'])
    # Use the thermal-only part to estimate slope (activity is bounded)
    gamma = tJ / tF**2
    msd_thm = 4.0 * D * gamma**2 * np.cumsum(G_vec(t, tJ, tF)**2) * (t[1]-t[0])
    # Slope from last 20% of time
    i70, i90 = int(0.70*len(t)), int(0.90*len(t))
    slope_thm = (msd_thm[i90] - msd_thm[i70]) / (t[i90] - t[i70])
    slope_tot = (msd[i90] - msd[i70]) / (t[i90] - t[i70])
    print(f"  {p['label']}: thermal slope = {slope_thm:.4f}, total slope ≈ {slope_tot:.4f}  "
          f"(thermal ratio = {slope_thm/(4*D):.4f})")

# Check 4: short-time MSD scaling
print(f"\n[MSD] Short-time: MSD ~ D*t^5/(5*tau_F^4)")
for p in params_fig4:
    tJ, tF = p['tau_J'], p['tau_F']
    t, _, _, msd = compute_section4p1(tJ, tF, tau_C, D=D, v0=v0,
                                       T_max=p['T_max'], N=p['N'])
    c5 = D / (5.0 * tF**4)
    # Check at a few early times
    for t_check in [0.01, 0.05, 0.1]:
        idx = int(t_check / (p['T_max']/p['N']))
        if idx > 0 and idx < len(t):
            msd_num  = msd[idx]
            msd_anal = c5 * t[idx]**5
            ratio = msd_num / msd_anal if msd_anal > 0 else float('nan')
            print(f"  {p['label']}: t={t[idx]:.3f}: MSD_num={msd_num:.6e}, "
                  f"D*t^5/(5*tF^4)={msd_anal:.6e}, ratio={ratio:.4f}")

print("\nDone. Figures saved.")
plt.show()
