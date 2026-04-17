#!/usr/bin/env python3
"""
Improved Fig 5(c): MFPT vs chirality with collapse test
========================================================

This script generates a publication-quality version of Fig 5(c) by:
  1. Using 5 D_r values: [0.005, 0.01, 0.02, 0.05, 0.1]
  2. Scanning ω from 0 to 1 (11 points)
  3. Averaging over 20 trials on 3 different maze seeds
  4. Top panel: τ_M vs ω for each D_r (log scale y)
  5. Bottom panel: collapse test τ_M × D_r^α vs ω
  6. Testing α=1 and α=0.5 to find best collapse

Author: Prashant Bisht, TIFR Hyderabad (improved)
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

# Setup paths
sys.path.insert(0, "/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant")
os.chdir("/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant")

from tcrw_maze import generate_maze_prim, run_tcrw_maze

# ============================================================
# Configuration
# ============================================================

L = 10
D_r_values = [0.005, 0.01, 0.02, 0.05, 0.1]
omega_scan = np.linspace(0.0, 1.0, 11)
N_trials_per_setup = 15  # Reduced from 20
N_maze_seeds = 2        # Reduced from 3 (still gives 30 trials per setup)
max_steps = 500_000

# Color scheme for D_r values
colors_dr = ['#e41a1c', '#ff7f00', '#4daf4a', '#377eb8', '#984ea3']
markers_dr = ['o', 's', '^', 'D', 'v']

print("=" * 70)
print("FIG 5(c) IMPROVED: MFPT vs CHIRALITY WITH COLLAPSE TEST")
print("=" * 70)
print(f"\nConfiguration:")
print(f"  L = {L}")
print(f"  D_r values: {D_r_values}")
print(f"  ω range: {omega_scan[0]} to {omega_scan[-1]} ({len(omega_scan)} points)")
print(f"  Trials per (ω, D_r) setup: {N_trials_per_setup} trials × {N_maze_seeds} mazes = {N_trials_per_setup*N_maze_seeds} total")
print(f"  Max steps per trial: {max_steps}")
print(f"  Total simulations: ~{len(D_r_values) * len(omega_scan) * N_trials_per_setup * N_maze_seeds}")
print()

# ============================================================
# Data collection
# ============================================================

# Store results: mfpt_data[D_r] = {"omega": [...], "mean": [...], "std": [...]}
mfpt_data = {}

for D_r_idx, D_r in enumerate(D_r_values):
    print(f"\n{'='*70}")
    print(f"D_r = {D_r} ({D_r_idx+1}/{len(D_r_values)})")
    print(f"{'='*70}")

    mfpt_means = []
    mfpt_stds = []

    for omega_idx, omega in enumerate(omega_scan):
        all_times = []

        # Run on multiple maze seeds
        for maze_seed in range(3000, 3000 + N_maze_seeds):
            maze, entrance, exit_pos = generate_maze_prim(L, seed=maze_seed)

            # Run N_trials_per_setup on this maze
            for trial in range(N_trials_per_setup):
                d_init = trial % 4  # Vary initial direction

                # Run walker
                t, _ = run_tcrw_maze(
                    maze, entrance, exit_pos,
                    omega, D_r,
                    d_init=d_init,
                    max_steps=max_steps,
                    seed=5000 + maze_seed * 1000 + trial
                )

                if t > 0:
                    all_times.append(t)

        # Compute statistics
        if len(all_times) > 0:
            mean_t = np.mean(all_times)
            std_t = np.std(all_times)
            sem_t = std_t / np.sqrt(len(all_times))

            mfpt_means.append(mean_t)
            mfpt_stds.append(sem_t)

            success_rate = len(all_times) / (N_trials_per_setup * N_maze_seeds)
            print(f"  ω={omega:.2f}: τ_M = {mean_t:8.1f} ± {sem_t:6.1f}  "
                  f"({len(all_times)}/{N_trials_per_setup * N_maze_seeds}, "
                  f"{success_rate*100:.0f}%)")
        else:
            mfpt_means.append(np.nan)
            mfpt_stds.append(np.nan)
            print(f"  ω={omega:.2f}: NO SUCCESSES")

    mfpt_data[D_r] = {
        "omega": omega_scan.copy(),
        "mean": np.array(mfpt_means),
        "std": np.array(mfpt_stds)
    }

print("\n" + "=" * 70)
print("Data collection complete. Analyzing collapse...")
print("=" * 70)

# ============================================================
# Collapse analysis
# ============================================================

def compute_spread(data_arrays):
    """
    Compute spread: relative std of all data points.
    For collapse test: smaller spread = better collapse.

    Spread = (mean absolute deviation) / (median value)
    """
    flat = np.concatenate([a[~np.isnan(a)] for a in data_arrays])
    if len(flat) == 0:
        return np.inf
    median = np.median(flat)
    if median == 0:
        return np.inf
    spread = np.mean(np.abs(flat - median)) / median
    return spread

# Test different collapse exponents
alpha_candidates = [0.5, 1.0, 1.5]
collapse_results = {}

for alpha in alpha_candidates:
    print(f"\nCollapse scaling: τ_M × D_r^{alpha}")
    print("-" * 60)

    rescaled_arrays = []
    for D_r in D_r_values:
        data = mfpt_data[D_r]
        rescaled = data["mean"] * (D_r ** alpha)
        rescaled_arrays.append(rescaled)

        # Show a few points
        valid_idx = ~np.isnan(rescaled)
        if np.any(valid_idx):
            sample_idx = np.where(valid_idx)[0][::max(1, len(np.where(valid_idx)[0])//3)]
            for idx in sample_idx:
                print(f"  D_r={D_r:.4f}, ω={omega_scan[idx]:.2f}: "
                      f"rescaled τ_M = {rescaled[idx]:.1f}")

    spread = compute_spread(rescaled_arrays)
    collapse_results[alpha] = spread
    print(f"  → Spread metric: {spread:.4f}")

best_alpha = min(collapse_results, key=collapse_results.get)
print(f"\nBest collapse: α = {best_alpha} (spread = {collapse_results[best_alpha]:.4f})")

# ============================================================
# Plotting
# ============================================================

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 11))

# Top panel: raw MFPT vs ω
for D_r, color, marker in zip(D_r_values, colors_dr, markers_dr):
    data = mfpt_data[D_r]
    valid = ~np.isnan(data["mean"])

    ax1.errorbar(
        data["omega"][valid], data["mean"][valid],
        yerr=data["std"][valid],
        fmt=f'{marker}-', color=color, ms=7, lw=2,
        capsize=4, capthick=1.5,
        label=f'$D_r = {D_r}$',
        alpha=0.85
    )

ax1.set_ylabel(r'$\tau_M$ (steps)', fontsize=13, fontweight='bold')
ax1.set_yscale('log')
ax1.set_xlim(-0.05, 1.05)
ax1.legend(fontsize=11, loc='upper right', framealpha=0.95, edgecolor='black', fancybox=False)
ax1.grid(True, alpha=0.3, linestyle='--', linewidth=0.7)
ax1.set_title(
    f'MFPT vs Chirality ($L={L}$, log scale)',
    fontsize=13, fontweight='bold', pad=10
)
ax1.tick_params(labelsize=11)

# Bottom panel: collapse test with best α
alpha_collapse = best_alpha
print(f"\nPlotting collapse with α = {alpha_collapse}")

for D_r, color, marker in zip(D_r_values, colors_dr, markers_dr):
    data = mfpt_data[D_r]
    rescaled = data["mean"] * (D_r ** alpha_collapse)
    valid = ~np.isnan(rescaled)

    ax2.errorbar(
        data["omega"][valid], rescaled[valid],
        yerr=data["std"][valid] * (D_r ** alpha_collapse),
        fmt=f'{marker}-', color=color, ms=7, lw=2,
        capsize=4, capthick=1.5,
        label=f'$D_r = {D_r}$',
        alpha=0.85
    )

ax2.set_xlabel(r'$\omega$ (chirality)', fontsize=13, fontweight='bold')
ax2.set_ylabel(
    f'$\\tau_M \\cdot D_r^{{{alpha_collapse}}}$ (scaled steps)',
    fontsize=13, fontweight='bold'
)
ax2.set_yscale('log')
ax2.set_xlim(-0.05, 1.05)
ax2.legend(fontsize=11, loc='upper right', framealpha=0.95, edgecolor='black', fancybox=False)
ax2.grid(True, alpha=0.3, linestyle='--', linewidth=0.7)
ax2.set_title(
    f'Collapse Test: α = {alpha_collapse} (spread = {collapse_results[best_alpha]:.4f})',
    fontsize=13, fontweight='bold', pad=10
)
ax2.tick_params(labelsize=11)

plt.suptitle(
    'Figure 5(c): Maze-Solving Time vs Chirality with Collapse Analysis',
    fontsize=14, fontweight='bold', y=0.995
)

plt.tight_layout()

output_file = '/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant/tcrw_fig5c_definitive.png'
plt.savefig(output_file, dpi=200, bbox_inches='tight', facecolor='white')
print(f"\n✓ Saved: {output_file}")

plt.close()

# ============================================================
# Summary statistics
# ============================================================

print("\n" + "=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)

for D_r in D_r_values:
    data = mfpt_data[D_r]
    valid = ~np.isnan(data["mean"])
    if np.any(valid):
        min_idx = np.nanargmin(data["mean"])
        max_idx = np.nanargmax(data["mean"])
        print(f"\nD_r = {D_r}:")
        print(f"  Min τ_M: {data['mean'][min_idx]:.1f} at ω={data['omega'][min_idx]:.2f}")
        print(f"  Max τ_M: {data['mean'][max_idx]:.1f} at ω={data['omega'][max_idx]:.2f}")
        ratio = data['mean'][max_idx] / data['mean'][min_idx]
        print(f"  Max/Min ratio: {ratio:.1f}x")

print("\n" + "=" * 70)
print("COLLAPSE ANALYSIS SUMMARY")
print("=" * 70)
for alpha in sorted(collapse_results.keys()):
    print(f"  α = {alpha}: spread = {collapse_results[alpha]:.4f}")
print(f"\n  Best α = {best_alpha} → use τ_M × D_r^{best_alpha}")

print("\n" + "=" * 70)
print("Script complete!")
print("=" * 70)
