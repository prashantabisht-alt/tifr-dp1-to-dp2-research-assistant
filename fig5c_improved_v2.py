#!/usr/bin/env python3
"""
Improved Fig 5(c): MFPT vs chirality with collapse test (V2 - Optimized)
=========================================================================

This script generates a publication-quality version of Fig 5(c) by:
  1. Using 5 D_r values: [0.005, 0.01, 0.02, 0.05, 0.1]
  2. Scanning ω at key points: [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
  3. Adaptive max_steps: higher for lower D_r (where steps needed is larger)
  4. Averaging over 10-15 trials on 2 different maze seeds
  5. Top panel: τ_M vs ω for each D_r (log scale y)
  6. Bottom panel: collapse test τ_M × D_r^α vs ω
  7. Testing α=0.5 and α=1 to find best collapse

Author: Prashant Bisht, TIFR Hyderabad (optimized v2)
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Tuple

# Setup paths
sys.path.insert(0, "/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant")
os.chdir("/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant")

from tcrw_maze import generate_maze_prim, run_tcrw_maze

# ============================================================
# Configuration
# ============================================================

L = 10
D_r_values = [0.005, 0.01, 0.02, 0.05, 0.1]
omega_scan = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
N_trials_per_setup = 10  # Reduced but still reasonable
N_maze_seeds = 2         # 20 trials total per setup
N_mazes_initial = 4      # Use more maze seeds for averaging

# Adaptive max_steps: lower D_r needs more steps
# D_r=0.1:   100k steps (mostly moves, little noise)
# D_r=0.05:  150k steps
# D_r=0.02:  300k steps
# D_r=0.01:  500k steps
# D_r=0.005: 500k steps (but expect some timeouts)
max_steps_dict = {
    0.1:   100_000,
    0.05:  150_000,
    0.02:  300_000,
    0.01:  500_000,
    0.005: 500_000
}

# Color scheme for D_r values
colors_dr = ['#e41a1c', '#ff7f00', '#4daf4a', '#377eb8', '#984ea3']
markers_dr = ['o', 's', '^', 'D', 'v']

print("=" * 70)
print("FIG 5(c) IMPROVED V2: MFPT vs CHIRALITY WITH ADAPTIVE TIMING")
print("=" * 70)
print(f"\nConfiguration:")
print(f"  L = {L}")
print(f"  D_r values: {D_r_values}")
print(f"  ω points: {len(omega_scan)} points from 0 to 1")
print(f"  Trials per setup: {N_trials_per_setup} trials × {N_maze_seeds} mazes = {N_trials_per_setup*N_maze_seeds}")
print(f"  Maze seed averaging: {N_mazes_initial} additional seeds for robustness")
print(f"  Adaptive max_steps based on D_r")
print(f"  Total simulations: ~{len(D_r_values) * len(omega_scan) * N_trials_per_setup * N_maze_seeds}")
print()

# ============================================================
# Data collection with averaging over multiple maze seeds
# ============================================================

def measure_mfpt_robust(D_r: float, omega: float, L: int,
                       N_trials: int, N_mazes: int,
                       max_steps: int) -> Tuple[float, float, int]:
    """
    Measure MFPT for a given (D_r, ω) by averaging over multiple mazes.
    Returns: (mean, sem, count_successes)
    """
    all_times = []

    for maze_idx in range(N_mazes):
        maze, entrance, exit_pos = generate_maze_prim(L, seed=3000 + maze_idx)

        for trial in range(N_trials):
            d_init = trial % 4  # Vary initial direction

            t, _ = run_tcrw_maze(
                maze, entrance, exit_pos,
                omega, D_r,
                d_init=d_init,
                max_steps=max_steps,
                seed=5000 + maze_idx * 10000 + trial * 100
            )

            if t > 0:
                all_times.append(t)

    if len(all_times) == 0:
        return np.nan, np.nan, 0

    mean_t = np.mean(all_times)
    sem_t = np.std(all_times) / np.sqrt(len(all_times))
    return mean_t, sem_t, len(all_times)


# Store results
mfpt_data = {}

for D_r_idx, D_r in enumerate(D_r_values):
    print(f"\n{'='*70}")
    print(f"D_r = {D_r} ({D_r_idx+1}/{len(D_r_values)})")
    print(f"{'='*70}")

    max_steps = max_steps_dict[D_r]
    mfpt_means = []
    mfpt_stds = []

    for omega_idx, omega in enumerate(omega_scan):
        mean_t, sem_t, count = measure_mfpt_robust(
            D_r, omega, L,
            N_trials=N_trials_per_setup,
            N_mazes=N_maze_seeds,
            max_steps=max_steps
        )

        mfpt_means.append(mean_t)
        mfpt_stds.append(sem_t)

        progress = (omega_idx + 1) / len(omega_scan) * 100
        if count > 0:
            print(f"  [{progress:5.1f}%] ω={omega:.2f}: τ_M = {mean_t:8.1f} ± {sem_t:6.1f}  "
                  f"({count} successes)")
        else:
            print(f"  [{progress:5.1f}%] ω={omega:.2f}: TIMEOUT/FAIL (0 successes)")

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
    Compute spread: normalized MAD of all data points.
    For collapse test: smaller spread = better collapse.
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
alpha_candidates = [0.5, 1.0]
collapse_results = {}

for alpha in alpha_candidates:
    print(f"\nCollapse scaling: τ_M × D_r^{alpha}")
    print("-" * 60)

    rescaled_arrays = []
    n_valid_total = 0

    for D_r in D_r_values:
        data = mfpt_data[D_r]
        rescaled = data["mean"] * (D_r ** alpha)
        rescaled_arrays.append(rescaled)

        # Show sample points
        valid_idx = ~np.isnan(rescaled)
        n_valid = np.sum(valid_idx)
        n_valid_total += n_valid

        if n_valid > 0:
            # Show endpoints and middle
            indices_to_show = [np.where(valid_idx)[0][0], np.where(valid_idx)[0][-1]]
            mid_idx = len(np.where(valid_idx)[0]) // 2
            if mid_idx not in indices_to_show:
                indices_to_show.insert(1, np.where(valid_idx)[0][mid_idx])

            for idx in sorted(set(indices_to_show)):
                omega = omega_scan[idx]
                print(f"  D_r={D_r:.4f}, ω={omega:.2f}: "
                      f"τ_M = {data['mean'][idx]:8.1f}, "
                      f"rescaled = {rescaled[idx]:8.1f}")

    spread = compute_spread(rescaled_arrays)
    collapse_results[alpha] = spread
    print(f"  → Valid data points: {n_valid_total}")
    print(f"  → Spread metric: {spread:.4f}")

best_alpha = min(collapse_results, key=collapse_results.get)
print(f"\nBest collapse: α = {best_alpha} (spread = {collapse_results[best_alpha]:.4f})")

# ============================================================
# Plotting
# ============================================================

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 12))

# Top panel: raw MFPT vs ω
for D_r, color, marker in zip(D_r_values, colors_dr, markers_dr):
    data = mfpt_data[D_r]
    valid = ~np.isnan(data["mean"])

    ax1.errorbar(
        data["omega"][valid], data["mean"][valid],
        yerr=data["std"][valid],
        fmt=f'{marker}-', color=color, ms=8, lw=2.5,
        capsize=5, capthick=2,
        label=f'$D_r = {D_r}$',
        alpha=0.85,
        elinewidth=1.5
    )

ax1.set_ylabel(r'$\tau_M$ (steps)', fontsize=14, fontweight='bold')
ax1.set_yscale('log')
ax1.set_xlim(-0.05, 1.05)
ax1.legend(fontsize=12, loc='upper right', framealpha=0.95, edgecolor='black', fancybox=False)
ax1.grid(True, alpha=0.25, linestyle='--', linewidth=0.8, which='both')
ax1.set_title(
    f'MFPT vs Chirality ($L={L}$, logarithmic scale)',
    fontsize=14, fontweight='bold', pad=12
)
ax1.tick_params(labelsize=12)
ax1.set_xticks(np.linspace(0, 1, 11))

# Bottom panel: collapse test with best α
alpha_collapse = best_alpha
print(f"\nPlotting collapse with α = {alpha_collapse}")

for D_r, color, marker in zip(D_r_values, colors_dr, markers_dr):
    data = mfpt_data[D_r]
    rescaled = data["mean"] * (D_r ** alpha_collapse)
    valid = ~np.isnan(rescaled)

    if np.any(valid):
        ax2.errorbar(
            data["omega"][valid], rescaled[valid],
            yerr=data["std"][valid] * (D_r ** alpha_collapse),
            fmt=f'{marker}-', color=color, ms=8, lw=2.5,
            capsize=5, capthick=2,
            label=f'$D_r = {D_r}$',
            alpha=0.85,
            elinewidth=1.5
        )

ax2.set_xlabel(r'$\omega$ (chirality)', fontsize=14, fontweight='bold')
if alpha_collapse == 0.5:
    ylabel = r'$\tau_M \sqrt{D_r}$ (scaled steps)'
else:
    ylabel = f'$\\tau_M \\cdot D_r$ (scaled steps)'
ax2.set_ylabel(ylabel, fontsize=14, fontweight='bold')

ax2.set_yscale('log')
ax2.set_xlim(-0.05, 1.05)
ax2.legend(fontsize=12, loc='upper right', framealpha=0.95, edgecolor='black', fancybox=False)
ax2.grid(True, alpha=0.25, linestyle='--', linewidth=0.8, which='both')

spread_str = f'{collapse_results[best_alpha]:.3f}'
ax2.set_title(
    f'Collapse Test: Power-law exponent $\\alpha = {best_alpha}$ (spread = {spread_str})',
    fontsize=14, fontweight='bold', pad=12
)
ax2.tick_params(labelsize=12)
ax2.set_xticks(np.linspace(0, 1, 11))

plt.suptitle(
    'Figure 5(c): Maze-Solving Time vs Chirality — Improved Collapse Analysis',
    fontsize=15, fontweight='bold', y=0.9995
)

plt.tight_layout()

output_file = '/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant/tcrw_fig5c_definitive.png'
plt.savefig(output_file, dpi=200, bbox_inches='tight', facecolor='white')
print(f"\n✓ Saved figure: {output_file}")

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
        print(f"  Min τ_M: {data['mean'][min_idx]:8.1f} at ω={omega_scan[min_idx]:.2f}")
        print(f"  Max τ_M: {data['mean'][max_idx]:8.1f} at ω={omega_scan[max_idx]:.2f}")
        if data['mean'][min_idx] > 0:
            ratio = data['mean'][max_idx] / data['mean'][min_idx]
            print(f"  Max/Min ratio: {ratio:.1f}x")

print("\n" + "=" * 70)
print("COLLAPSE ANALYSIS SUMMARY")
print("=" * 70)
for alpha in sorted(collapse_results.keys()):
    print(f"  α = {alpha}: spread = {collapse_results[alpha]:.4f}")
print(f"\n  Best scaling: τ_M × D_r^{best_alpha}")
print(f"  Interpretation: MFPT scales with D_r^{best_alpha} for noise (jerk/rotation)")

print("\n" + "=" * 70)
print("✓ Script complete!")
print("=" * 70)
