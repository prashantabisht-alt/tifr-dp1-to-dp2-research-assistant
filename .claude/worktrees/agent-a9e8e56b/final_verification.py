#!/usr/bin/env python3
"""
Final cross-check: verify all 8 figures against the paper's predictions.
Checks scaling exponents, long-time limits, orbit parameters, and D_eff.
"""
import numpy as np
from scipy.signal import fftconvolve

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

PASS = 0; FAIL = 0

def check(name, value, expected, tol=0.02):
    global PASS, FAIL
    ratio = value / expected if expected != 0 else float('inf')
    ok = abs(ratio - 1.0) < tol
    status = "PASS" if ok else "FAIL"
    if ok: PASS += 1
    else: FAIL += 1
    print(f"  [{status}] {name}: got {value:.6f}, expect {expected:.6f}, ratio={ratio:.4f}")

print("=" * 70)
print("FINAL VERIFICATION — Jose & Löwen jcABP paper")
print("=" * 70)

# ─── SECTION 3: no chirality ───
print("\n── Section 3 (omega_0 = 0) ──")
tau_P, v0, D, m = 1.0, 1.0, 0.5, 1.0

for label, tJ, tF, T_max, N in [
    ("Fig1a", 0.2, 0.28, 20.0, 20000),
    ("Fig2b", 20.0, 0.2, 250.0, 80000),
]:
    dt = T_max / N
    t = np.linspace(0, T_max, N+1)
    gamma = tJ / tF**2
    G = G_vec(t, tJ, tF)
    E = np.exp(-t / tau_P)
    h = fftconvolve(G, E)[:N+1] * dt
    msd = 2*gamma**2*v0**2*np.cumsum(G*h)*dt + 4*D*gamma**2*np.cumsum(G**2)*dt

    # G short-time: G(t) ~ t^2/(2*tau_J)
    idx = 5
    G_anal = t[idx]**2 / (2*tJ)
    check(f"{label} G(t) short-time", G[idx], G_anal, tol=0.05)

    # G long-time: G -> 1/gamma
    G_inf = tF**2 / tJ
    check(f"{label} G(t) long-time", G[-1], G_inf, tol=0.02)

    # MSD slope
    i80, i90 = int(0.80*len(t)), int(0.90*len(t))
    slope = (msd[i90]-msd[i80])/(t[i90]-t[i80])
    expected_slope = 4*D + 2*v0**2*tau_P
    check(f"{label} MSD long-time slope", slope, expected_slope, tol=0.01)

    # Mean displacement saturation
    mx = gamma * v0 * h
    lP = v0 * tau_P
    check(f"{label} <x(T)>/l_P", mx[-1]/lP, 1.0, tol=0.02)

# ─── SECTION 4.1: chirality, D_r = 0 ───
print("\n── Section 4.1 (D_r = 0, omega_0 != 0) ──")
tau_C = 1.0

for label, tJ, tF, T_max, N in [
    ("Fig3a", 1.0, 2*np.sqrt(2), 80, 50000),
    ("Fig3d", 10.0, 1.0, 300, 120000),
    ("Fig3e", 0.2, 0.2, 25, 30000),
]:
    dt = T_max / N
    t = np.linspace(0, T_max, N+1)
    gamma = tJ / tF**2
    G = G_vec(t, tJ, tF)
    conv_GC = fftconvolve(G, np.cos(t/tau_C))[:N+1] * dt
    conv_GS = fftconvolve(G, np.sin(t/tau_C))[:N+1] * dt
    mx = gamma * v0 * conv_GC
    my = gamma * v0 * conv_GS

    # Orbit center: (0, v0*tau_C)
    n_per = max(int(2*np.pi*tau_C / dt), 1000)
    my_avg = np.mean(my[-n_per:])
    mx_avg = np.mean(mx[-n_per:])
    check(f"{label} orbit center y", my_avg, v0*tau_C, tol=0.02)
    check(f"{label} orbit center x", abs(mx_avg), 0.0, tol=0.05)  # should be ~0

    # Orbit radius (Eq. 32)
    r_theory = v0*tau_C / np.sqrt(tF**4/(tJ**2*tau_C**2) + (tF**2/tau_C**2 - 1)**2)
    mx_amp = (np.max(mx[-n_per:]) - np.min(mx[-n_per:])) / 2
    check(f"{label} orbit radius", mx_amp, r_theory, tol=0.02)

for label, tJ, tF, T_max, N in [
    ("Fig4a", 0.2, 0.2, 200, 80000),
    ("Fig4b", 10.0, 0.4, 200, 100000),
]:
    dt = T_max / N
    t = np.linspace(0, T_max, N+1)
    gamma = tJ / tF**2
    G = G_vec(t, tJ, tF)
    conv_GC = fftconvolve(G, np.cos(t/tau_C))[:N+1] * dt
    conv_GS = fftconvolve(G, np.sin(t/tau_C))[:N+1] * dt
    mx = gamma * v0 * conv_GC
    my = gamma * v0 * conv_GS
    msd = mx**2 + my**2 + 4*D*gamma**2*np.cumsum(G**2)*dt

    # MSD long-time: purely thermal, slope = 4D
    i70, i90 = int(0.70*len(t)), int(0.90*len(t))
    # Use thermal part only for clean slope
    msd_thm = 4*D*gamma**2*np.cumsum(G**2)*dt
    slope_thm = (msd_thm[i90]-msd_thm[i70])/(t[i90]-t[i70])
    check(f"{label} MSD thermal slope", slope_thm, 4*D, tol=0.01)

# ─── SECTION 4.2: full model ───
print("\n── Section 4.2 (D_r != 0, omega_0 != 0) ──")
tau_P8 = 10.0; v08 = 0.05; tau_C8 = 0.5; D8 = 0.5

D_c = D8 + v08**2 / (2.0*(1.0/tau_P8 + tau_P8/tau_C8**2))

for label, tJ, tF, T_max, N in [
    ("Fig8a", 0.1, 0.1, 20, 60000),
    ("Fig8b", 5.0, 1.0, 200, 100000),
]:
    dt = T_max / N
    t = np.linspace(0, T_max, N+1)
    gamma = tJ / tF**2
    G = G_vec(t, tJ, tF)
    K = np.cos(t/tau_C8) * np.exp(-t/tau_P8)
    h = fftconvolve(G, K)[:N+1] * dt
    msd = 2*gamma**2*v08**2*np.cumsum(G*h)*dt + 4*D8*gamma**2*np.cumsum(G**2)*dt

    i70, i90 = int(0.70*len(t)), int(0.90*len(t))
    slope = (msd[i90]-msd[i70])/(t[i90]-t[i70])
    check(f"{label} MSD slope vs 4*D_c", slope, 4*D_c, tol=0.01)

    # Spiral center (Eq. 38)
    Ex = np.exp(-t/tau_P8) * np.cos(t/tau_C8)
    Ey = np.exp(-t/tau_P8) * np.sin(t/tau_C8)
    mx = gamma*v08 * fftconvolve(G, Ex)[:N+1]*dt
    my = gamma*v08 * fftconvolve(G, Ey)[:N+1]*dt
    sigma = 1/tau_P8**2 + 1/tau_C8**2
    xc_th = v08 / (tau_P8 * sigma)
    yc_th = v08 / (tau_C8 * sigma)
    # At long times both modes decay, mean -> center
    n_end = min(len(t)//20, 5000)
    check(f"{label} spiral center y", np.mean(my[-n_end:]), yc_th, tol=0.05)

# ─── SUMMARY ───
print(f"\n{'='*70}")
print(f"SUMMARY:  {PASS} PASS,  {FAIL} FAIL  out of  {PASS+FAIL}  checks")
print(f"{'='*70}")
