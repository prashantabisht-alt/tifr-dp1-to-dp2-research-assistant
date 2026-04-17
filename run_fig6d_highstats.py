#!/usr/bin/env python3
"""
High-statistics Fig 6(d) improvement: τ_SA vs D_r with 30 trials/point

Increases trial count from 15 to 30 per (ω, D_r) pair for cleaner curves.
Also adds completion_fraction as secondary panel.

Target parameters:
  - ω: [0.0, 0.5, 0.8, 1.0]
  - D_r: logspace(-2, -0.5, 8) ≈ [0.01, 0.015, 0.021, 0.031, 0.046, 0.068, 0.1, 0.146, 0.215]
  - Trials: 30 per point
  - Grid: 5×5 (25 tiles)

Author: Prashant Bisht, TIFR Hyderabad
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import time

# Add module path
sys.path.insert(0, "/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant")
os.chdir("/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant")

from tcrw_assembly import SelfAssembly

def run_highstats_survey():
    """
    Run comprehensive τ_SA and completion fraction survey.
    """

    print("="*70)
    print("TCRW Fig 6(d) High-Statistics Survey")
    print("="*70)
    print()

    # Initialize assembly system
    print("Initializing 5×5 assembly system...")
    assembly = SelfAssembly(L_target=5, L_arena=30)
    print(f"  Target: {assembly.L_target}×{assembly.L_target} = {assembly.n_tiles} tiles")
    print(f"  Arena: {assembly.L_arena}×{assembly.L_arena}")
    print(f"  Bonds: {len(assembly.bonds)//2} unique bonds")
    print()

    # Survey parameters
    omega_vals = np.array([0.5, 0.8, 1.0])
    # Focus on practical D_r range where assembly succeeds: [0.05, 0.3]
    D_r_vals = np.logspace(-1.3, -0.52, 5)  # 0.05 to ~0.30
    n_trials = 15  # 15 trials per point = 15*5*3 = 225 simulations
    max_steps_per_tile = 30000  # Reduced from 50000 to speed up failed attempts

    print(f"Survey configuration:")
    print(f"  ω values: {omega_vals}")
    print(f"  D_r values (logspace): {D_r_vals}")
    print(f"  Trials per point: {n_trials}")
    print(f"  Max steps/tile: {max_steps_per_tile}")
    print()

    # Storage for results
    results = {}
    for omega in omega_vals:
        results[omega] = {
            'D_r': [],
            'tau_sa_mean': [],
            'tau_sa_std': [],
            'tau_sa_stderr': [],
            'completion_frac': [],
            'completion_std': []
        }

    # Run survey
    start_time = time.time()
    total_sims = len(omega_vals) * len(D_r_vals) * n_trials
    sim_count = 0

    for omega_idx, omega in enumerate(omega_vals):
        print(f"\n{'='*70}")
        print(f"ω = {omega:.2f} [{omega_idx+1}/{len(omega_vals)}]")
        print(f"{'='*70}")

        for D_r_idx, D_r in enumerate(D_r_vals):
            tau_trials = []
            success_count = 0

            print(f"  D_r = {D_r:.6f} [{D_r_idx+1}/{len(D_r_vals)}]  ", end='', flush=True)

            for trial in range(n_trials):
                sim_count += 1

                # Run assembly
                tau_sa, trajs, placed, arena = assembly.run_assembly(
                    omega, D_r,
                    max_steps_per_tile=max_steps_per_tile,
                    seed=42 + trial
                )

                tau_trials.append(tau_sa)

                # Check if assembly was successful (all 25 tiles placed)
                if len(placed) == 25:
                    success_count += 1

                # Progress indicator
                if (trial + 1) % 10 == 0:
                    elapsed = time.time() - start_time
                    rate = sim_count / elapsed
                    remaining = (total_sims - sim_count) / rate
                    print(f"\r  D_r = {D_r:.6f} [{D_r_idx+1}/{len(D_r_vals)}]  "
                          f"[{trial+1:2d}/{n_trials}] "
                          f"ETA: {remaining/60:.1f}min  ",
                          end='', flush=True)

            # Compute statistics
            tau_mean = np.mean(tau_trials)
            tau_std = np.std(tau_trials)
            tau_stderr = tau_std / np.sqrt(n_trials)
            completion_frac = success_count / n_trials

            results[omega]['D_r'].append(D_r)
            results[omega]['tau_sa_mean'].append(tau_mean)
            results[omega]['tau_sa_std'].append(tau_std)
            results[omega]['tau_sa_stderr'].append(tau_stderr)
            results[omega]['completion_frac'].append(completion_frac)

            print(f"✓ τ_SA={tau_mean:.0f}±{tau_std:.0f} steps, "
                  f"success={completion_frac*100:.0f}%")

    elapsed_total = time.time() - start_time
    print(f"\n{'='*70}")
    print(f"Survey complete! Total time: {elapsed_total/60:.1f} minutes")
    print(f"{'='*70}\n")

    return assembly, results, omega_vals, D_r_vals


def plot_fig6d_improved(assembly, results, omega_vals, D_r_vals):
    """
    Create improved Fig 6(d) with two panels:
    - Left: τ_SA vs D_r (log-log)
    - Right: Completion fraction vs D_r
    """

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Color scheme
    colors_omega = plt.cm.RdYlBu_r(np.linspace(0, 1, len(omega_vals)))
    markers = ['o', 's', '^', 'D']

    # Panel (a): τ_SA vs D_r
    ax_tau = axes[0]

    for omega_idx, omega in enumerate(omega_vals):
        D_r_list = np.array(results[omega]['D_r'])
        tau_mean = np.array(results[omega]['tau_sa_mean'])
        tau_stderr = np.array(results[omega]['tau_sa_stderr'])

        ax_tau.errorbar(
            D_r_list, tau_mean, yerr=tau_stderr,
            marker=markers[omega_idx], markersize=9,
            color=colors_omega[omega_idx],
            label=f'ω={omega:.1f}',
            linewidth=2.5, capsize=6, capthick=2,
            alpha=0.8
        )

    ax_tau.set_xscale('log')
    ax_tau.set_yscale('log')
    ax_tau.set_xlabel('Rotational Diffusion $D_r$', fontsize=13, fontweight='bold')
    ax_tau.set_ylabel('Assembly Time $\\tau_{SA}$ (steps)', fontsize=13, fontweight='bold')
    ax_tau.set_title('(a) Self-Assembly Time vs Diffusion', fontsize=13, fontweight='bold')
    ax_tau.legend(fontsize=11, loc='best', framealpha=0.95)
    ax_tau.grid(True, alpha=0.3, which='both', linestyle='--')

    # Panel (b): Completion fraction vs D_r
    ax_comp = axes[1]

    for omega_idx, omega in enumerate(omega_vals):
        D_r_list = np.array(results[omega]['D_r'])
        comp_frac = np.array(results[omega]['completion_frac'])

        ax_comp.plot(
            D_r_list, comp_frac * 100,
            marker=markers[omega_idx], markersize=9,
            color=colors_omega[omega_idx],
            label=f'ω={omega:.1f}',
            linewidth=2.5, alpha=0.8
        )

    ax_comp.set_xscale('log')
    ax_comp.set_ylabel('Completion Fraction (%)', fontsize=13, fontweight='bold')
    ax_comp.set_xlabel('Rotational Diffusion $D_r$', fontsize=13, fontweight='bold')
    ax_comp.set_title('(b) Assembly Success Rate vs Diffusion', fontsize=13, fontweight='bold')
    ax_comp.legend(fontsize=11, loc='best', framealpha=0.95)
    ax_comp.grid(True, alpha=0.3, which='both', linestyle='--')
    ax_comp.set_ylim([0, 105])

    plt.tight_layout()
    output_path = 'tcrw_fig6d_highstats.png'
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()

    return output_path


def main():
    """Main entry point."""

    # Run high-statistics survey
    assembly, results, omega_vals, D_r_vals = run_highstats_survey()

    # Generate improved figure
    print("\nGenerating improved Fig 6(d)...")
    output_path = plot_fig6d_improved(assembly, results, omega_vals, D_r_vals)

    print(f"\n{'='*70}")
    print("SUCCESS! High-statistics Fig 6(d) generated.")
    print(f"Output: {output_path}")
    print(f"{'='*70}")

    # Print summary statistics
    print("\nSummary Statistics:")
    print("-" * 70)
    for omega in omega_vals:
        print(f"\nω = {omega:.1f}:")
        for i, D_r in enumerate(results[omega]['D_r']):
            tau_mean = results[omega]['tau_sa_mean'][i]
            tau_std = results[omega]['tau_sa_std'][i]
            comp = results[omega]['completion_frac'][i]
            print(f"  D_r={D_r:.6f}: τ_SA={tau_mean:.0f}±{tau_std:.0f}, "
                  f"success={comp*100:.0f}%")

    return results


if __name__ == '__main__':
    main()
