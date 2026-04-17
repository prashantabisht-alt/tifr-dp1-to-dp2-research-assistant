"""
Run high-statistics assembly (Fig 6d + 6e) using Python backend.
Optimized: fewer extreme points, lower max_steps, suppressed verbose output.
"""
import sys, os, time, contextlib, io
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import tcrw_accuracy_boost as tab
tab.HAS_NUMBA = False
from tcrw_accuracy_boost import run_assembly_fast

N_TRIALS = 15
L_ARENA = 40
MAX_STEPS = 200000  # reduced from 500000

def run_silent(omega, D_r, n_trials):
    """Run assembly suppressing per-tile progress prints."""
    with contextlib.redirect_stdout(io.StringIO()):
        tau, npl = run_assembly_fast(omega, D_r, n_trials=n_trials,
                                      L_arena=L_ARENA, max_steps_per_tile=MAX_STEPS)
    return tau, npl

# ============================================================
# Fig 6(d): τ_SA vs D_r for several ω
# ============================================================
print("=" * 70)
print(f"FIG 6(d): τ_SA vs D_r ({N_TRIALS} trials/point)")
print("=" * 70)

D_r_vals = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5])
omega_vals = [0.0, 0.5, 0.8, 1.0]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
markers = ['o', 's', '^', 'D']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

results_d = {}
t_total = time.time()
for oi, omega in enumerate(omega_vals):
    means = []
    stds = []
    completion = []

    for D_r in D_r_vals:
        t0 = time.time()
        tau_arr, npl_arr = run_silent(omega, D_r, N_TRIALS)
        dt = time.time() - t0

        completed = npl_arr == 25
        frac_completed = completed.mean()
        if frac_completed > 0:
            tau_completed = tau_arr[completed]
            m = tau_completed.mean()
            s = tau_completed.std() / np.sqrt(len(tau_completed))
        else:
            m = MAX_STEPS * 24
            s = 0

        means.append(m)
        stds.append(s)
        completion.append(frac_completed)
        print(f"  ω={omega:.1f} D_r={D_r:.3f}: τ={m:.0f}±{s:.0f}, "
              f"complete={frac_completed:.0%} ({dt:.1f}s)")

    means = np.array(means)
    stds = np.array(stds)
    completion = np.array(completion)
    results_d[omega] = (means, stds, completion)

    ax1.errorbar(D_r_vals, means, yerr=stds,
                 marker=markers[oi], ms=7, lw=2, capsize=4,
                 color=colors[oi], label=f'ω={omega}')

    ax2.plot(D_r_vals, completion, marker=markers[oi], ms=7, lw=2,
             color=colors[oi], label=f'ω={omega}')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('$D_r$', fontsize=13)
ax1.set_ylabel(r'Assembly time $\tau_{SA}$ (steps)', fontsize=13)
ax1.set_title(f'Fig 6(d): τ_SA vs D_r ({N_TRIALS} trials/point)', fontsize=13)
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3, which='both')

ax2.set_xscale('log')
ax2.set_xlabel('$D_r$', fontsize=13)
ax2.set_ylabel('Fraction completed', fontsize=13)
ax2.set_title(f'Completion fraction ({N_TRIALS} trials/point)', fontsize=13)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(-0.05, 1.05)

plt.tight_layout()
plt.savefig('tcrw_fig6d_highstats.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"\nSaved: tcrw_fig6d_highstats.png ({time.time()-t_total:.0f}s total)")

# ============================================================
# Fig 6(e): τ_SA vs ω for several D_r
# ============================================================
print("\n" + "=" * 70)
print(f"FIG 6(e): τ_SA vs ω ({N_TRIALS} trials/point)")
print("=" * 70)

omega_sweep = np.linspace(0, 1, 9)  # 9 points instead of 11
D_r_vals_e = [0.1, 0.2, 0.3]
colors_e = ['#2166ac', '#d6604d', '#4dac26']

fig, ax = plt.subplots(figsize=(10, 7))

t_total2 = time.time()
for di, D_r in enumerate(D_r_vals_e):
    means = []
    stds = []

    for omega in omega_sweep:
        t0 = time.time()
        tau_arr, npl_arr = run_silent(omega, D_r, N_TRIALS)
        dt = time.time() - t0

        completed = npl_arr == 25
        if completed.any():
            tau_c = tau_arr[completed]
            m = tau_c.mean()
            s = tau_c.std() / np.sqrt(len(tau_c))
        else:
            m = MAX_STEPS * 24
            s = 0

        means.append(m)
        stds.append(s)
        print(f"  D_r={D_r:.2f} ω={omega:.2f}: τ={m:.0f}±{s:.0f} ({dt:.1f}s)")

    ax.errorbar(omega_sweep, means, yerr=stds,
                marker='o', ms=7, lw=2, capsize=4,
                color=colors_e[di], label=f'$D_r$={D_r}')

ax.set_xlabel('Chirality $\\omega$', fontsize=13)
ax.set_ylabel(r'Assembly time $\tau_{SA}$ (steps)', fontsize=13)
ax.set_title(f'Fig 6(e): τ_SA vs ω ({N_TRIALS} trials/point)', fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_xlim(-0.05, 1.05)

plt.tight_layout()
plt.savefig('tcrw_fig6e_highstats.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"\nSaved: tcrw_fig6e_highstats.png ({time.time()-t_total2:.0f}s total)")

# ============================================================
# Chiral advantage ratio
# ============================================================
print("\n" + "=" * 70)
print("CHIRAL ADVANTAGE SUMMARY")
print("=" * 70)
if 0.0 in results_d and 1.0 in results_d:
    tau_achiral = results_d[0.0][0]
    tau_chiral = results_d[1.0][0]
    ratio = tau_achiral / tau_chiral
    for i, D_r in enumerate(D_r_vals):
        r = ratio[i]
        label = "CHIRAL ADVANTAGE" if r > 1 else "no advantage"
        print(f"  D_r={D_r:.3f}: τ_achiral/τ_chiral = {r:.3f}  [{label}]")

print("\nDone!")
