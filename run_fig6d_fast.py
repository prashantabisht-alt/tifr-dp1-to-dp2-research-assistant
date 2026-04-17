#!/usr/bin/env python3
"""
Optimized Fig 6(d) reproduction: τ_SA vs D_r with 20 trials/point
Removes verbose print statements for speed.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from collections import deque
import time

sys.path.insert(0, "/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant")
os.chdir("/sessions/optimistic-epic-mayer/mnt/TIFR DP1 to DP2 Research Assistant")

# Import just the structure and bond system
from tcrw_assembly import SelfAssembly

# Direction vectors: d = 0(↑), 1(→), 2(↓), 3(←)
DX = np.array([0, 1, 0, -1], dtype=np.int64)
DY = np.array([1, 0, -1, 0], dtype=np.int64)


def fast_run_assembly(assembly, omega, D_r, max_steps_per_tile=30000, seed=42):
    """
    Run assembly without verbose output (faster).
    """
    rng = np.random.default_rng(seed)
    L_arena = assembly.L_arena
    arena = np.full((L_arena, L_arena), -1, dtype=int)

    # Place seed
    cx, cy = L_arena // 2, L_arena // 2
    seed_tid = 12  # center tile for 5x5
    arena[cx, cy] = seed_tid
    placed = {seed_tid: (cx, cy)}

    # Build release order (BFS from seed)
    L_target = assembly.L_target
    target = assembly.target
    bonds = assembly.bonds

    tid_to_target = {}
    for r in range(L_target):
        for c in range(L_target):
            tid_to_target[target[r, c]] = (r, c)

    visited_tiles = {seed_tid}
    queue = deque([seed_tid])
    release_order = []

    while queue:
        tid = queue.popleft()
        r, c = tid_to_target[tid]
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < L_target and 0 <= nc < L_target:
                nb_tid = target[nr, nc]
                if nb_tid not in visited_tiles:
                    visited_tiles.add(nb_tid)
                    release_order.append(nb_tid)
                    queue.append(nb_tid)

    total_steps = 0

    for tile_id in release_order:
        # Release tile
        placed_tiles = list(placed.values())
        anchor = placed_tiles[rng.integers(0, len(placed_tiles))]

        x, y = None, None
        for _ in range(100):
            angle = rng.random() * 2 * np.pi
            r_release = rng.integers(4, 12)
            sx = int(anchor[0] + r_release * np.cos(angle))
            sy = int(anchor[1] + r_release * np.sin(angle))

            if (0 <= sx < L_arena and 0 <= sy < L_arena and arena[sx, sy] == -1):
                x, y = sx, sy
                break

        if x is None:
            for _ in range(100):
                sx = rng.integers(0, L_arena)
                sy = rng.integers(0, L_arena)
                if arena[sx, sy] == -1:
                    x, y = sx, sy
                    break

        if x is None:
            total_steps += max_steps_per_tile
            continue

        d = rng.integers(0, 4)
        bound = False

        for step in range(max_steps_per_tile):
            r_step = rng.random()
            r_rot = rng.random()

            if r_step < D_r:
                # NOISE STEP
                if r_rot < omega:
                    d = (d - 1) % 4
                else:
                    d = (d + 1) % 4
            else:
                # CHIRAL STEP
                nx = x + DX[d]
                ny = y + DY[d]

                if (0 <= nx < L_arena and 0 <= ny < L_arena and arena[nx, ny] == -1):
                    x, y = nx, ny
                    if r_rot < omega:
                        d = (d + 1) % 4
                    else:
                        d = (d - 1) % 4

            # Check binding
            for dd in range(4):
                ax = x + DX[dd]
                ay = y + DY[dd]

                if 0 <= ax < L_arena and 0 <= ay < L_arena:
                    neighbor_tid = arena[ax, ay]
                    if neighbor_tid >= 0:
                        facing_edge = dd
                        if (tile_id, facing_edge) in bonds:
                            if bonds[(tile_id, facing_edge)] == neighbor_tid:
                                arena[x, y] = tile_id
                                placed[tile_id] = (x, y)
                                bound = True
                                break

            if bound:
                total_steps += step + 1
                break

        if not bound:
            total_steps += max_steps_per_tile

    return total_steps, len(placed)


def main():
    print("="*70)
    print("TCRW Fig 6(d) - Fast High-Statistics Survey")
    print("="*70)
    print()

    assembly = SelfAssembly(L_target=5, L_arena=30)
    print(f"System initialized: {assembly.L_target}×{assembly.L_target}, "
          f"arena {assembly.L_arena}×{assembly.L_arena}")
    print()

    # Parameters
    omega_vals = np.array([0.5, 0.8, 1.0])
    D_r_vals = np.logspace(-1.3, -0.52, 5)  # 0.05 to 0.30
    n_trials = 20
    max_steps_per_tile = 30000

    print(f"Survey config:")
    print(f"  ω values: {omega_vals}")
    print(f"  D_r values: {D_r_vals}")
    print(f"  Trials: {n_trials} per point")
    print(f"  Max steps/tile: {max_steps_per_tile}")
    print()

    # Results storage
    results = {omega: {'D_r': [], 'tau_sa_mean': [], 'tau_sa_std': [],
                       'tau_sa_stderr': [], 'completion_frac': []}
               for omega in omega_vals}

    total_sims = len(omega_vals) * len(D_r_vals) * n_trials
    sim_count = 0
    start_time = time.time()

    for omega_idx, omega in enumerate(omega_vals):
        print(f"\nω = {omega:.1f}:")

        for D_r_idx, D_r in enumerate(D_r_vals):
            tau_trials = []
            success_count = 0

            for trial in range(n_trials):
                sim_count += 1

                tau_sa, n_placed = fast_run_assembly(
                    assembly, omega, D_r,
                    max_steps_per_tile=max_steps_per_tile,
                    seed=42 + trial
                )

                tau_trials.append(tau_sa)
                if n_placed == 25:
                    success_count += 1

                if (trial + 1) % 5 == 0:
                    elapsed = time.time() - start_time
                    rate = sim_count / max(elapsed, 0.1)
                    remaining = (total_sims - sim_count) / max(rate, 0.1)
                    pct = 100 * sim_count / total_sims
                    print(f"  D_r={D_r:.4f}: [{trial+1:2d}/{n_trials}] ({pct:.0f}%) "
                          f"ETA {remaining/60:.1f}min", end='\r', flush=True)

            # Statistics
            tau_mean = np.mean(tau_trials)
            tau_std = np.std(tau_trials)
            tau_stderr = tau_std / np.sqrt(n_trials)
            comp_frac = success_count / n_trials

            results[omega]['D_r'].append(D_r)
            results[omega]['tau_sa_mean'].append(tau_mean)
            results[omega]['tau_sa_std'].append(tau_std)
            results[omega]['tau_sa_stderr'].append(tau_stderr)
            results[omega]['completion_frac'].append(comp_frac)

            print(f"  D_r={D_r:.4f}: τ_SA={tau_mean:.0f}±{tau_std:.0f} steps, "
                  f"success={comp_frac*100:.0f}%")

    elapsed_total = time.time() - start_time
    print(f"\nComputation complete: {elapsed_total/60:.1f} minutes\n")

    # Generate figure
    print("Generating Fig 6(d)...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    colors_omega = plt.cm.RdYlBu_r(np.linspace(0, 1, len(omega_vals)))
    markers = ['o', 's', '^']

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
            linewidth=2.5, capsize=6, capthick=2, alpha=0.8
        )

    ax_tau.set_xscale('log')
    ax_tau.set_yscale('log')
    ax_tau.set_xlabel('Rotational Diffusion $D_r$', fontsize=13, fontweight='bold')
    ax_tau.set_ylabel('Assembly Time $\\tau_{SA}$ (steps)', fontsize=13, fontweight='bold')
    ax_tau.set_title('(a) Self-Assembly Time vs Diffusion', fontsize=13, fontweight='bold')
    ax_tau.legend(fontsize=11, loc='best')
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
    ax_comp.set_title('(b) Assembly Success Rate', fontsize=13, fontweight='bold')
    ax_comp.legend(fontsize=11, loc='best')
    ax_comp.grid(True, alpha=0.3, which='both', linestyle='--')
    ax_comp.set_ylim([0, 105])

    plt.tight_layout()
    output_path = 'tcrw_fig6d_highstats.png'
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()

    # Summary
    print("\n" + "="*70)
    print("RESULTS SUMMARY")
    print("="*70)
    for omega in omega_vals:
        print(f"\nω = {omega:.1f}:")
        for i, D_r in enumerate(results[omega]['D_r']):
            tau = results[omega]['tau_sa_mean'][i]
            tau_err = results[omega]['tau_sa_stderr'][i]
            comp = results[omega]['completion_frac'][i]
            print(f"  D_r={D_r:.4f}: τ_SA={tau:.0f}±{tau_err:.0f}, "
                  f"success={comp*100:.0f}%")


if __name__ == '__main__':
    main()
