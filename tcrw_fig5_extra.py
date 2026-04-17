"""
TCRW Phase 6C: Extended Figure 5 Analysis
==========================================

This module generates three additional figures extending the TCRW maze results:

Fig 1: Fig 5(b) — MFPT heatmap by starting position
  A single maze (L=10, seed=42) with MFPT measured from each passage cell.
  Walls are gray; passage cells colored by MFPT (blue=fast, red=slow).

Fig 2: Fig 5(e) — Disconnected (non-simply-connected) mazes
  Standard mazes are trees (Prim's). We create loops by removing extra walls.
  Test effect of loops on MFPT across different D_r values.

Fig 3: Fig 5(f) — Wide-passage mazes
  Standard mazes have passage width=1. We scale them up (passage_width=3).
  Show visit-frequency heatmaps for ω=0, 0.5, 1 (like Fig 5a but wide).

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from collections import deque
import sys


# ============================================================
# Import functions from tcrw_maze.py
# ============================================================

def generate_maze_prim(L, seed=42):
    """
    Generate a random maze using randomized Prim's algorithm.
    """
    rng = np.random.default_rng(seed)
    N = 2 * L + 1
    maze = np.ones((N, N), dtype=np.int32)
    visited = np.zeros((L, L), dtype=bool)

    start_i = rng.integers(0, L)
    start_j = rng.integers(0, L)
    visited[start_i, start_j] = True
    maze[2 * start_i + 1, 2 * start_j + 1] = 0

    frontier = []

    def add_frontier(i, j):
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ni, nj = i + di, j + dj
            if 0 <= ni < L and 0 <= nj < L and not visited[ni, nj]:
                wi = 2 * i + 1 + di
                wj = 2 * j + 1 + dj
                frontier.append((ni, nj, wi, wj))

    add_frontier(start_i, start_j)

    while frontier:
        idx = rng.integers(0, len(frontier))
        ni, nj, wi, wj = frontier[idx]
        frontier[idx] = frontier[-1]
        frontier.pop()

        if visited[ni, nj]:
            continue

        visited[ni, nj] = True
        maze[2 * ni + 1, 2 * nj + 1] = 0
        maze[wi, wj] = 0

        add_frontier(ni, nj)

    entrance = (1, 0)
    maze[entrance[0], entrance[1]] = 0

    exit_pos = (2 * L - 1, 2 * L)
    maze[exit_pos[0], exit_pos[1]] = 0

    return maze, entrance, exit_pos


def run_tcrw_maze(maze, entrance, exit_pos, omega, D_r, d_init=1,
                  max_steps=10_000_000, seed=42, record_traj=False):
    """
    Run a single TCRW walker through a maze until it reaches the exit.
    """
    rng = np.random.default_rng(seed)
    DX = [0, 1, 0, -1]
    DY = [1, 0, -1, 0]

    x, y = entrance
    d = d_init

    traj = [(x, y)] if record_traj else None

    for t in range(1, max_steps + 1):
        r_step = rng.random()
        r_rot = rng.random()

        if r_step < D_r:
            # Noise step: stay, rotate
            if r_rot < omega:
                d = (d - 1) % 4   # CCW
            else:
                d = (d + 1) % 4   # CW
        else:
            # Chiral step: try translate, then rotate
            nx = x + DX[d]
            ny = y + DY[d]

            if (0 <= nx < maze.shape[0] and 0 <= ny < maze.shape[1]
                    and maze[nx, ny] == 0):
                x, y = nx, ny
                if r_rot < omega:
                    d = (d + 1) % 4   # CW
                else:
                    d = (d - 1) % 4   # CCW

        if record_traj:
            traj.append((x, y))

        if (x, y) == exit_pos:
            return t, traj

    return -1, traj


# ============================================================
# Figure 5(b): MFPT heatmap by starting position
# ============================================================

def fig5b_mfpt_heatmap():
    """
    For a single maze (L=10, seed=42), measure MFPT from each passage cell.

    Create heatmap: passage cells colored by MFPT, walls gray.
    """
    print("\n--- Fig 5(b): MFPT heatmap by starting position ---")

    L = 10
    omega = 1.0
    D_r = 0.01
    max_steps = 5_000_000

    # Generate maze once
    maze, entrance, exit_pos = generate_maze_prim(L, seed=42)
    N = maze.shape[0]

    # Find all passage cells
    passage_cells = []
    for i in range(N):
        for j in range(N):
            if maze[i, j] == 0:
                passage_cells.append((i, j))

    print(f"  Maze size: {N}×{N}, {len(passage_cells)} passage cells")
    print(f"  Parameters: ω={omega}, D_r={D_r}")
    print(f"  Measuring MFPT from each passage cell...")

    # For each passage cell, measure MFPT by averaging over director states and seeds
    mfpt_map = np.full((N, N), np.nan)
    count_completed = 0

    for idx, (start_x, start_y) in enumerate(passage_cells):
        if (idx + 1) % 10 == 0:
            print(f"    [{idx + 1}/{len(passage_cells)}] Processing...")

        times = []

        # Run multiple trials: average over 4 director states × 5 random seeds
        for d_init in range(4):
            for seed_offset in range(5):
                t, _ = run_tcrw_maze(
                    maze, (start_x, start_y), exit_pos,
                    omega, D_r, d_init=d_init,
                    max_steps=max_steps,
                    seed=42 + seed_offset * 1000 + idx * 10
                )
                if t > 0:
                    times.append(t)

        if len(times) > 0:
            mfpt_map[start_x, start_y] = np.mean(times)
            count_completed += 1

    print(f"  Completed: {count_completed}/{len(passage_cells)} cells")

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))

    # Build composite RGB image: walls gray, passages colored by MFPT
    rgb = np.zeros((N, N, 3))

    # Walls: dark gray
    wall_mask = maze == 1
    rgb[wall_mask] = [0.35, 0.35, 0.35]

    # Passages: color by MFPT
    passage_mask = maze == 0
    valid_mfpt = mfpt_map[passage_mask & ~np.isnan(mfpt_map)]

    if len(valid_mfpt) > 0:
        min_mfpt = np.min(valid_mfpt)
        max_mfpt = np.max(valid_mfpt)
        norm_mfpt = (valid_mfpt - min_mfpt) / (max_mfpt - min_mfpt) if max_mfpt > min_mfpt else np.zeros_like(valid_mfpt)

        # Color map: blue (low/fast) → white → red (high/slow)
        # Using a simple colormap: blue=fast, red=slow
        cmap = mcolors.LinearSegmentedColormap.from_list(
            'mfpt', ['#0000FF', '#FFFFFF', '#FF0000']
        )

        # Apply colormap to valid passage cells
        valid_mask = passage_mask & ~np.isnan(mfpt_map)
        normalized = (mfpt_map[valid_mask] - min_mfpt) / (max_mfpt - min_mfpt) if max_mfpt > min_mfpt else np.zeros(np.sum(valid_mask))
        colors = cmap(normalized)
        rgb[valid_mask] = colors[:, :3]

        # Unvisited passages (due to max_steps cutoff): light gray
        unvisited_mask = passage_mask & np.isnan(mfpt_map)
        rgb[unvisited_mask] = [0.95, 0.95, 0.95]

    # Display image
    ax.imshow(rgb.transpose(1, 0, 2), origin='lower',
              extent=(-0.5, N-0.5, -0.5, N-0.5), interpolation='nearest')

    # Mark entrance and exit
    ax.plot(*entrance, 'o', color='#00AA00', ms=10, zorder=5,
            markeredgecolor='white', markeredgewidth=1, label='Entrance')
    ax.plot(*exit_pos, '*', color='#FFFF00', ms=15, zorder=5,
            markeredgecolor='black', markeredgewidth=1, label='Exit')

    ax.set_aspect('equal')
    ax.set_xlim(-1, N)
    ax.set_ylim(-1, N)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(fontsize=10, loc='upper right')

    ax.set_title(
        f'Fig 5(b): MFPT heatmap ($L={L}$, $\\omega={omega}$, $D_r={D_r}$)',
        fontsize=13, pad=15
    )

    # Colorbar for MFPT
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min_mfpt, vmax=max_mfpt))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r'MFPT $\tau_M$ (steps)', fontsize=11)

    plt.tight_layout()
    plt.savefig('tcrw_fig5b_mfpt_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: tcrw_fig5b_mfpt_heatmap.png\n")


# ============================================================
# Figure 5(e): Disconnected (non-simply-connected) mazes
# ============================================================

def create_loopy_maze(maze_base, n_loops=5, seed=42):
    """
    Start with a Prim's tree maze and create loops by removing extra walls.

    Parameters
    ----------
    maze_base : 2D array
        Base maze from Prim's algorithm (0=passage, 1=wall)
    n_loops : int
        Number of wall segments to remove (creating loops)
    seed : int
        Random seed

    Returns
    -------
    maze_loopy : 2D array
        Modified maze with loops
    """
    maze_loopy = maze_base.copy()
    rng = np.random.default_rng(seed)
    N = maze_base.shape[0]

    # Find all wall cells (candidates for removal)
    wall_cells = []
    for i in range(N):
        for j in range(N):
            if maze_base[i, j] == 1:
                # Check if removing this wall would create a loop
                # (i.e., it connects two passage cells)
                neighbors = 0
                for di, dj in [(-1,0), (1,0), (0,-1), (0,1)]:
                    ni, nj = i + di, j + dj
                    if 0 <= ni < N and 0 <= nj < N and maze_base[ni, nj] == 0:
                        neighbors += 1
                if neighbors >= 2:  # Wall between 2+ passage cells
                    wall_cells.append((i, j))

    # Randomly remove n_loops walls
    if len(wall_cells) > 0:
        n_to_remove = min(n_loops, len(wall_cells))
        indices_to_remove = rng.choice(len(wall_cells), size=n_to_remove, replace=False)
        for idx in indices_to_remove:
            i, j = wall_cells[idx]
            maze_loopy[i, j] = 0

    return maze_loopy


def fig5e_disconnected_mazes():
    """
    Compare MFPT in tree mazes vs loopy mazes across D_r values.

    Key result: loops trap random walkers → large MFPT at small D_r.
    """
    print("\n--- Fig 5(e): Disconnected (loopy) mazes ---")

    L = 8
    omega = 1.0
    D_r_values = [0.001, 0.01, 0.1]
    max_steps = 5_000_000
    N_trials = 10
    N_mazes = 5

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # --- Test 1: Tree maze (standard Prim's) ---
    print("  Tree maze (standard Prim's):")
    mfpt_tree = []
    mfpt_tree_err = []

    for D_r in D_r_values:
        all_times = []
        for maze_idx in range(N_mazes):
            maze, entrance, exit_pos = generate_maze_prim(L, seed=1000 + maze_idx)
            for trial in range(N_trials):
                d_init = trial % 4
                t, _ = run_tcrw_maze(
                    maze, entrance, exit_pos, omega, D_r,
                    d_init=d_init, max_steps=max_steps,
                    seed=2000 + maze_idx * 100 + trial
                )
                if t > 0:
                    all_times.append(t)

        if len(all_times) > 0:
            m = np.mean(all_times)
            e = np.std(all_times) / np.sqrt(len(all_times))
            mfpt_tree.append(m)
            mfpt_tree_err.append(e)
            print(f"    D_r={D_r}: τ_M = {m:.0f} ± {e:.0f}")
        else:
            mfpt_tree.append(np.nan)
            mfpt_tree_err.append(np.nan)
            print(f"    D_r={D_r}: DNF")

    # --- Test 2: Loopy maze (Prim's + removed walls) ---
    print("  Loopy maze (Prim's + removed walls):")
    mfpt_loopy = []
    mfpt_loopy_err = []

    for D_r in D_r_values:
        all_times = []
        for maze_idx in range(N_mazes):
            maze_tree, entrance, exit_pos = generate_maze_prim(L, seed=1000 + maze_idx)
            maze = create_loopy_maze(maze_tree, n_loops=5, seed=3000 + maze_idx)
            for trial in range(N_trials):
                d_init = trial % 4
                t, _ = run_tcrw_maze(
                    maze, entrance, exit_pos, omega, D_r,
                    d_init=d_init, max_steps=max_steps,
                    seed=4000 + maze_idx * 100 + trial
                )
                if t > 0:
                    all_times.append(t)

        if len(all_times) > 0:
            m = np.mean(all_times)
            e = np.std(all_times) / np.sqrt(len(all_times))
            mfpt_loopy.append(m)
            mfpt_loopy_err.append(e)
            print(f"    D_r={D_r}: τ_M = {m:.0f} ± {e:.0f}")
        else:
            mfpt_loopy.append(np.nan)
            mfpt_loopy_err.append(np.nan)
            print(f"    D_r={D_r}: DNF")

    # Plot comparison
    x_pos = np.arange(len(D_r_values))
    width = 0.35

    ax1.bar(x_pos - width/2, mfpt_tree, width, label='Tree maze', color='#2166ac', alpha=0.8)
    ax1.bar(x_pos + width/2, mfpt_loopy, width, label='Loopy maze', color='#b2182b', alpha=0.8)

    ax1.set_xlabel(r'$D_r$', fontsize=12)
    ax1.set_ylabel(r'MFPT $\tau_M$ (steps)', fontsize=12)
    ax1.set_yscale('log')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels([f'{D_r}' for D_r in D_r_values])
    ax1.legend(fontsize=11)
    ax1.set_title(f'MFPT comparison ($L={L}$, $\\omega={omega}$)', fontsize=12)
    ax1.grid(axis='y', alpha=0.3)

    # Plot ratio (loopy / tree)
    ratio = []
    ratio_err = []
    for m_tree, e_tree, m_loopy, e_loopy in zip(mfpt_tree, mfpt_tree_err, mfpt_loopy, mfpt_loopy_err):
        if m_tree > 0 and m_loopy > 0:
            r = m_loopy / m_tree
            # Error propagation
            e_r = r * np.sqrt((e_tree/m_tree)**2 + (e_loopy/m_loopy)**2)
            ratio.append(r)
            ratio_err.append(e_r)
        else:
            ratio.append(np.nan)
            ratio_err.append(np.nan)

    valid = ~np.isnan(ratio)
    ax2.errorbar(np.array(D_r_values)[valid], np.array(ratio)[valid],
                 yerr=np.array(ratio_err)[valid],
                 fmt='o-', color='#d8b365', ms=8, lw=2, capsize=5)
    ax2.axhline(y=1, color='k', linestyle='--', alpha=0.5, label='Tree = Loopy')
    ax2.set_xlabel(r'$D_r$', fontsize=12)
    ax2.set_ylabel(r'MFPT ratio $\tau_M^{\text{loopy}} / \tau_M^{\text{tree}}$', fontsize=12)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.legend(fontsize=10)
    ax2.set_title('Loop effect on MFPT', fontsize=12)
    ax2.grid(alpha=0.3)

    plt.suptitle('Fig 5(e): Loops in mazes trap walkers', fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig('tcrw_fig5e_disconnected.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: tcrw_fig5e_disconnected.png\n")


# ============================================================
# Figure 5(f): Wide-passage mazes
# ============================================================

def scale_maze_wide(L_logical, passage_width=3, seed=42):
    """
    Generate a logical maze (standard Prim's) and scale it up by passage_width.

    Each logical passage cell becomes a passage_width×passage_width block.
    Walls remain 1 cell thick.

    Parameters
    ----------
    L_logical : int
        Size parameter for logical maze
    passage_width : int
        Scaling factor for passage cells
    seed : int
        Random seed

    Returns
    -------
    maze_wide : 2D array
        Scaled maze
    entrance_wide : (int, int)
        Entrance in scaled coordinates
    exit_wide : (int, int)
        Exit in scaled coordinates
    """
    # Generate logical maze
    maze_logical, entrance_logical, exit_logical = generate_maze_prim(L_logical, seed=seed)

    N_logical = maze_logical.shape[0]

    # Build scaled maze by expanding each cell
    # Logical coordinate (i, j) in original maze...
    # In logical maze: odd i,j = passage, even i,j = wall
    # Strategy: for each logical cell, determine if it's wall or passage
    # If passage: expand to passage_width×passage_width block
    # If wall: expand to 1×(passage_width or 1) blocks

    # Actually simpler approach:
    # New maze size: for each row/col, compute its type in logical maze
    # Type is: odd=passage, even=wall
    # For passage row/col: expand to passage_width cells
    # For wall row/col: keep as 1 cell

    new_rows = []
    new_cols = []

    for i in range(N_logical):
        if i % 2 == 1:  # Passage row
            new_rows.extend([i] * passage_width)
        else:  # Wall row
            new_rows.append(i)

    for j in range(N_logical):
        if j % 2 == 1:  # Passage column
            new_cols.extend([j] * passage_width)
        else:  # Wall column
            new_cols.append(j)

    N_wide = len(new_rows)
    maze_wide = np.zeros((N_wide, N_wide), dtype=np.int32)

    for i_new, i_logical in enumerate(new_rows):
        for j_new, j_logical in enumerate(new_cols):
            maze_wide[i_new, j_new] = maze_logical[i_logical, j_logical]

    # Map entrance and exit to new coordinates
    # entrance_logical in maze_logical → find index in new_rows/new_cols
    entrance_i = entrance_logical[0]
    entrance_j = entrance_logical[1]
    entrance_wide = (new_rows.index(entrance_i), new_cols.index(entrance_j))

    exit_i = exit_logical[0]
    exit_j = exit_logical[1]
    exit_wide = (new_rows.index(exit_i), new_cols.index(exit_j))

    return maze_wide, entrance_wide, exit_wide


def fig5f_wide_maze():
    """
    Generate wide-passage mazes and show visit frequency for ω=0, 0.5, 1.
    """
    print("\n--- Fig 5(f): Wide-passage mazes ---")

    L_logical = 5
    passage_width = 3
    D_r = 1e-3
    omegas = [0.0, 0.5, 1.0]
    labels = [r'$\omega=0$', r'$\omega=0.5$', r'$\omega=1$']

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

    for ax, omega, label in zip(axes, omegas, labels):
        print(f"  Processing ω={omega}...")

        # Generate wide maze
        maze, entrance, exit_pos = scale_maze_wide(L_logical, passage_width=passage_width, seed=42)
        N = maze.shape[0]

        # Run walker (allow longer max_steps for wide mazes)
        max_s = 20_000_000 if omega == 0.5 else 10_000_000
        steps, traj = run_tcrw_maze(maze, entrance, exit_pos, omega, D_r,
                                    d_init=1, max_steps=max_s, seed=123,
                                    record_traj=True)

        # Build visit-frequency heatmap
        visit_count = np.zeros((N, N), dtype=float)
        if traj is not None:
            for (tx, ty) in traj:
                visit_count[tx, ty] += 1

        # Normalize for display
        passage_mask = maze == 0
        visit_display = np.zeros((N, N))
        if np.max(visit_count[passage_mask]) > 0:
            visit_display[passage_mask] = (
                np.log10(1 + visit_count[passage_mask]) /
                np.log10(1 + np.max(visit_count[passage_mask]))
            )

        # Create composite RGB: walls dark gray, passages warm colored
        rgb = np.zeros((N, N, 3))
        wall_mask = maze == 1
        rgb[wall_mask] = [0.35, 0.35, 0.35]

        unvisited = passage_mask & (visit_count == 0)
        rgb[unvisited] = [1.0, 1.0, 1.0]

        visited = passage_mask & (visit_count > 0)
        if np.any(visited):
            v = visit_display[visited]
            rgb[visited, 0] = 1.0 - 0.6 * v
            rgb[visited, 1] = 0.9 - 0.7 * v
            rgb[visited, 2] = 0.8 - 0.75 * v

        ax.imshow(rgb.transpose(1, 0, 2), origin='lower',
                  extent=(-0.5, N-0.5, -0.5, N-0.5), interpolation='nearest')

        # Mark entrance and exit
        ax.plot(*entrance, 'o', color='#2266CC', ms=8, zorder=5,
                markeredgecolor='white', markeredgewidth=0.8, label='Start')
        ax.plot(*exit_pos, '*', color='#22CC22', ms=12, zorder=5,
                markeredgecolor='white', markeredgewidth=0.5, label='Exit')

        result_str = f'{steps:,} steps' if steps > 0 else 'DNF'
        ax.set_title(f'{label}\n{result_str}', fontsize=13)
        ax.set_aspect('equal')
        ax.set_xlim(-1, N)
        ax.set_ylim(-1, N)
        ax.set_xticks([])
        ax.set_yticks([])

    axes[0].legend(fontsize=9, loc='upper right', framealpha=0.9)

    plt.suptitle(f'Fig 5(f): Wide-passage mazes ($L_{{\\rm logical}}={L_logical}$, '
                 f'$w={passage_width}$, $D_r={D_r}$)',
                 fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig('tcrw_fig5f_wide_maze.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  Saved: tcrw_fig5f_wide_maze.png\n")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("=" * 70)
    print("TCRW Phase 6C: Extended Figure 5 Analysis")
    print("=" * 70)

    # Quick sanity check
    print("\n--- Sanity check: maze generation ---")
    L = 5
    maze, entrance, exit_pos = generate_maze_prim(L, seed=42)
    print(f"  Standard maze (L={L}): {maze.shape}")
    print(f"  Entrance: {entrance}, Exit: {exit_pos}")

    # Loopy maze
    maze_loopy = create_loopy_maze(maze, n_loops=2, seed=42)
    print(f"  Loopy maze: {np.sum(maze)} → {np.sum(maze_loopy)} wall cells")

    # Wide maze
    maze_wide, ent_wide, exit_wide = scale_maze_wide(3, passage_width=3, seed=42)
    print(f"  Wide maze (logical L=3, width=3): {maze_wide.shape}")

    # Generate figures
    print("\n" + "=" * 70)
    fig5b_mfpt_heatmap()
    fig5e_disconnected_mazes()
    fig5f_wide_maze()

    print("=" * 70)
    print("Phase 6C complete. Generated 3 figures:")
    print("  - tcrw_fig5b_mfpt_heatmap.png")
    print("  - tcrw_fig5e_disconnected.png")
    print("  - tcrw_fig5f_wide_maze.png")
    print("=" * 70)
