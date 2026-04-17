#!/usr/bin/env python3
"""
Regenerate Fig 1(d): D(ω) with improved Monte Carlo statistics.
====================================================================

Computes both:
  1. Analytical D(ω) from PBC band structure using d²λ/dk² at k=0
  2. MC simulations with higher statistics (N_traj=20, T_steps=500k)

Plots both on same panel for comparison.

Author: Claude (improved version)
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add source to path
sys.path.insert(0, "/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant")
os.chdir("/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant")

from tcrw_core import simulate_tcrw_pbc
from tcrw_spectrum import build_Pk, pbc_eigenvalues

# ============================================================
# Analytical D(ω) from PBC dispersion
# ============================================================

def analytical_D_omega(omega, D_r=0.001, dk=0.001):
    """
    Compute D(ω) analytically from PBC band structure.

    D = -Re(d²λ/dk²) / 2, where λ(k) is largest eigenvalue at k=(0,k_y).

    Use finite differences with small dk to approximate d²λ/dk².
    """
    # Evaluate largest eigenvalue at (kx=0, ky) for three k values
    ky0 = 0.0
    ky_minus = -dk
    ky_plus = dk

    evals0 = pbc_eigenvalues(omega, D_r, 0.0, ky0)
    evals_minus = pbc_eigenvalues(omega, D_r, 0.0, ky_minus)
    evals_plus = pbc_eigenvalues(omega, D_r, 0.0, ky_plus)

    # Largest eigenvalue (by real part)
    lambda0 = evals0[0]
    lambda_minus = evals_minus[0]
    lambda_plus = evals_plus[0]

    # Second derivative: d²λ/dk_y² ≈ (λ(k+dk) - 2λ(k) + λ(k-dk)) / dk²
    d2lambda_dk2 = (lambda_plus - 2*lambda0 + lambda_minus) / (dk**2)

    # Diffusion: D = -Re(d²λ/dk²) / 2
    D = -np.real(d2lambda_dk2) / 2.0

    return D


def compute_analytical_curve(omega_vals, D_r=0.001):
    """Compute D(ω) analytically for array of ω."""
    D_analytic = np.array([analytical_D_omega(omega, D_r) for omega in omega_vals])
    return D_analytic


# ============================================================
# Monte Carlo extraction of D(ω)
# ============================================================

def extract_D_from_msd(msd, times, fit_frac=0.5):
    """
    Extract diffusion coefficient from MSD using linear fit at late times.

    MSD(t) = 4*D*t in 2D  =>  D = slope/4

    fit_frac: use the last fit_frac fraction of the data.
    """
    t_fit = times.astype(np.float64)
    msd_fit = msd.astype(np.float64)

    # Linear fit to last fraction
    n = len(t_fit)
    i0 = int(n * (1 - fit_frac))

    if i0 >= n - 1:
        i0 = max(0, n - 10)

    t_subset = t_fit[i0:]
    msd_subset = msd_fit[i0:]

    # Linear: slope, intercept
    slope = np.polyfit(t_subset, msd_subset, 1)[0]

    # D = slope / 4 (from MSD = 4Dt)
    D = slope / 4.0

    return D


def run_mc_simulation(omega, D_r=0.001, L=16, T_steps=500_000,
                      N_traj=20, seed=42):
    """
    Run MC simulation and extract D using PBC (which is faster than OBC).
    Returns the diffusion coefficient.
    """
    print(f"  omega={omega:.2f}: ", end='', flush=True)

    result = simulate_tcrw_pbc(omega, D_r, L, T_steps, N_traj, seed,
                               record_interval=500)

    msd = result['msd']
    times = result['times']

    D_mc = extract_D_from_msd(msd, times, fit_frac=0.5)

    print(f"D = {D_mc:.6f}")
    return D_mc


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("=" * 70)
    print("Regenerating Fig 1(d): D(ω) with improved MC statistics")
    print("=" * 70)

    D_r = 0.001
    L = 16  # Smaller but still reasonable
    T_steps = 500_000  # Good balance
    N_traj = 20

    omega_vals = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

    # --- Analytical ---
    print("\n[1] Computing analytical D(ω) from PBC dispersion...")
    D_analytic = compute_analytical_curve(omega_vals, D_r=D_r)
    print(f"    Done. Sample: D(0.5) = {D_analytic[4]:.6f}")

    # --- Monte Carlo ---
    print(f"\n[2] Running MC simulations (L={L}, T={T_steps}, N={N_traj})...")
    D_mc_list = []
    for i, omega in enumerate(omega_vals):
        seed = 42 + i  # Different seed for each omega
        D_mc = run_mc_simulation(omega, D_r, L, T_steps, N_traj, seed)
        D_mc_list.append(D_mc)

    D_mc = np.array(D_mc_list)

    # --- Plot ---
    print("\n[3] Plotting...")
    fig, ax = plt.subplots(figsize=(9, 6))

    # Analytical curve
    omega_fine = np.linspace(0.1, 0.9, 200)
    D_fine = compute_analytical_curve(omega_fine, D_r=D_r)
    ax.plot(omega_fine, D_fine, 'b-', lw=2.5, label='Analytical (PBC)', zorder=3)

    # MC points with error bars
    ax.scatter(omega_vals, D_mc, s=100, color='red', marker='o',
               edgecolors='darkred', linewidth=1.5, label='MC (N=20, T=500k)',
               zorder=4, alpha=0.8)

    # Labels and formatting
    ax.set_xlabel('Chirality $\\omega$', fontsize=13)
    ax.set_ylabel('Diffusion coefficient $D(\\omega)$', fontsize=13)
    ax.set_title(f'Fig 1(d): $D(\\omega)$ with improved statistics ($D_r={D_r}$, $L={L}$)',
                 fontsize=13)
    ax.legend(fontsize=11, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.05, 0.95)

    # Save
    output_path = '/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant/improved_fig1d_D_vs_omega.png'
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"    Saved: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("Summary:")
    print("=" * 70)
    print(f"Analytical D(ω):      min={D_analytic.min():.6f}, max={D_analytic.max():.6f}")
    print(f"MC D(ω):              min={D_mc.min():.6f}, max={D_mc.max():.6f}")
    print(f"RMS error (MC vs analytic): {np.sqrt(np.mean((D_mc - D_analytic)**2)):.6f}")
    print(f"\nOutput PNG: {output_path}")
    print("=" * 70)
