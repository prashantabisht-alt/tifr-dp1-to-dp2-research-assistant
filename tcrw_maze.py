"""
TCRW Phase 6B: Maze Solving
============================

A topological chiral random walker can efficiently solve mazes using the
"hand-on-the-wall" strategy: it finds an edge and follows it, thereby
systematically traversing any simply connected maze.

Algorithm:
  1. Generate random maze using Prim's algorithm on an L×L grid.
     The maze is encoded on a (2L+1) × (2L+1) array where:
       - Cells at odd coordinates (2i+1, 2j+1) are rooms
       - Cells between rooms can be walls (1) or passages (0)
       - Even-coordinate cells are always walls (pillars/borders)
  2. Place TCRW walker at entrance, run until it reaches exit.
  3. Measure mean first-passage time (MFPT) τ_M.

Key predictions (Fig 5):
  (a) Chiral walkers (ω=0,1) follow edges systematically → fast
  (b) Achiral walker (ω=0.5) does random walk → slow
  (c) τ_M vs ω shows minimum at ω=0,1 (parabolic shape)
  (d) Scaling: τ_M ~ L² (chiral) vs τ_M ~ L³ (achiral)

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from collections import deque


# ============================================================
# Maze generation: Prim's algorithm
# ============================================================

def generate_maze_prim(L, seed=42):
    """
    Generate a random maze using randomized Prim's algorithm.

    The maze lives on a (2L+1) × (2L+1) grid.
    Rooms are at positions (2i+1, 2j+1) for i,j = 0..L-1.
    Walls are everywhere else initially; Prim's algorithm carves passages.

    Parameters
    ----------
    L : int
        Number of rooms per side (L × L rooms).
    seed : int
        Random seed.

    Returns
    -------
    maze : (2L+1, 2L+1) int array
        0 = open (passage/room), 1 = wall
    entrance : (int, int)
        Coordinates of entrance cell in maze array.
    exit_pos : (int, int)
        Coordinates of exit cell in maze array.
    """
    rng = np.random.default_rng(seed)
    N = 2 * L + 1
    maze = np.ones((N, N), dtype=np.int32)

    # Mark rooms as open
    visited = np.zeros((L, L), dtype=bool)

    # Start from random room
    start_i = rng.integers(0, L)
    start_j = rng.integers(0, L)
    visited[start_i, start_j] = True
    maze[2 * start_i + 1, 2 * start_j + 1] = 0

    # Frontier: walls between visited and unvisited rooms
    frontier = []

    def add_frontier(i, j):
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ni, nj = i + di, j + dj
            if 0 <= ni < L and 0 <= nj < L and not visited[ni, nj]:
                # Wall between (i,j) and (ni,nj) in maze coords
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
        maze[2 * ni + 1, 2 * nj + 1] = 0  # open the room
        maze[wi, wj] = 0                    # open the wall between

        add_frontier(ni, nj)

    # Entrance: top-left corner room (0,0) → maze coord (1,1)
    # Open an entrance on the left border
    entrance = (1, 0)
    maze[entrance[0], entrance[1]] = 0

    # Exit: bottom-right corner room (L-1, L-1) → maze coord (2L-1, 2L-1)
    # Open an exit on the right border
    exit_pos = (2 * L - 1, 2 * L)
    maze[exit_pos[0], exit_pos[1]] = 0

    return maze, entrance, exit_pos


# ============================================================
# TCRW walker in a maze
# ============================================================

def run_tcrw_maze(maze, entrance, exit_pos, omega, D_r, d_init=1,
                  max_steps=10_000_000, seed=42, record_traj=False):
    """
    Run a single TCRW walker through a maze until it reaches the exit.

    The maze grid encodes walls (1) and passages (0). The walker moves
    on passage cells. Walls block translation (same as OBC blocking rule).

    Parameters
    ----------
    maze : 2D int array
        0=passage, 1=wall.
    entrance : (int, int)
        Starting position.
    exit_pos : (int, int)
        Target position.
    omega : float
        Chirality.
    D_r : float
        Rotational noise probability.
    d_init : int
        Initial director (0=↑, 1=→, 2=↓, 3=←). Default: → (pointing into maze).
    max_steps : int
        Safety cutoff.
    seed : int
    record_traj : bool
        If True, record the full trajectory.

    Returns
    -------
    steps : int
        Number of steps to reach exit (-1 if max_steps exceeded).
    traj : list of (x,y) or None
        Trajectory if recorded.
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

            # Check bounds and wall
            if (0 <= nx < maze.shape[0] and 0 <= ny < maze.shape[1]
                    and maze[nx, ny] == 0):
                x, y = nx, ny
                # Rotate after successful move
                if r_rot < omega:
                    d = (d + 1) % 4   # CW
                else:
                    d = (d - 1) % 4   # CCW
            # else: blocked (wall or out of bounds), no move, no rotation

        if record_traj:
            traj.append((x, y))

        # Check if exit reached
        if (x, y) == exit_pos:
            return t, traj

    return -1, traj


def measure_mfpt(maze, entrance, exit_pos, omega, D_r,
                 N_trials=100, max_steps=5_000_000, seed=42):
    """
    Measure mean first-passage time over N_trials.
    Average over random initial directors.
    """
    rng = np.random.default_rng(seed)
    times = []

    for trial in range(N_trials):
        d_init = rng.integers(0, 4)
        t, _ = run_tcrw_maze(maze, entrance, exit_pos, omega, D_r,
                             d_init=d_init,
                             max_steps=max_steps,
                             seed=seed + trial + 1)
        if t > 0:
            times.append(t)

    if len(times) == 0:
        return np.nan, np.nan

    return np.mean(times), np.std(times) / np.sqrt(len(times))


# ============================================================
# Figures
# ============================================================

def fig5a_trajectories():
    """
    Fig 5(a): Sample trajectories through a maze for ω=0, 0.5, 1.

    Matches paper style:
      - Warm gradient coloring (dark orange → light yellow for time)
      - Maze walls in gray/black
      - Blue circle = entrance, green star = exit
      - L=20 maze for visual clarity
    """
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LinearSegmentedColormap

    L = 12
    D_r = 1e-3
    omegas = [0.0, 0.5, 1.0]
    labels = [r'$\omega=0$', r'$\omega=0.5$', r'$\omega=1$']

    maze, entrance, exit_pos = generate_maze_prim(L, seed=42)

    # Paper-matching warm colormap: dark (early) → light (late)
    cmap_traj = LinearSegmentedColormap.from_list(
        'traj_warm', ['#552200', '#AA4400', '#DD6600', '#FF9933', '#FFCC66', '#FFEEAA'])

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

    for ax, omega, label in zip(axes, omegas, labels):
        # Run walker
        max_s = 10_000_000 if omega == 0.5 else 5_000_000
        steps, traj = run_tcrw_maze(maze, entrance, exit_pos, omega, D_r,
                                    d_init=1, max_steps=max_s, seed=123,
                                    record_traj=True)

        N = maze.shape[0]

        # Build visit-frequency heatmap (paper style)
        visit_count = np.zeros((N, N), dtype=float)
        if traj is not None:
            for (tx, ty) in traj:
                visit_count[tx, ty] += 1

        # Normalize per passage cell: log scale for better visibility
        # Wall cells stay at 0
        passage_mask = maze == 0
        visit_display = np.zeros((N, N))
        if np.max(visit_count[passage_mask]) > 0:
            visit_display[passage_mask] = (
                np.log10(1 + visit_count[passage_mask]) /
                np.log10(1 + np.max(visit_count[passage_mask]))
            )

        # Create composite image: walls=dark gray, passages colored by visits
        # RGB array: walls=[0.35,0.35,0.35], passages colored warm
        rgb = np.zeros((N, N, 3))
        # Walls
        wall_mask = maze == 1
        rgb[wall_mask] = [0.35, 0.35, 0.35]
        # Unvisited passages = white
        unvisited = passage_mask & (visit_count == 0)
        rgb[unvisited] = [1.0, 1.0, 1.0]
        # Visited passages: warm colormap (white → orange → dark brown)
        visited = passage_mask & (visit_count > 0)
        if np.any(visited):
            v = visit_display[visited]
            # Interpolate: 0=light peach, 1=dark brown
            rgb[visited, 0] = 1.0 - 0.6 * v      # R: 1.0 → 0.4
            rgb[visited, 1] = 0.9 - 0.7 * v      # G: 0.9 → 0.2
            rgb[visited, 2] = 0.8 - 0.75 * v     # B: 0.8 → 0.05

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

    plt.suptitle(f'Fig 5(a): Maze trajectories ($L={L}$, $D_r={D_r}$)',
                 fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig('tcrw_fig5a_maze_traj.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig5a_maze_traj.png")


def fig5c_mfpt_vs_omega():
    """
    Fig 5(c): MFPT τ_M vs ω for different D_r.

    Top: raw τ_M. Bottom: τ_M rescaled by 1/D_r.
    Average over multiple random mazes.
    """
    L = 10
    N_mazes = 10
    N_trials_per_maze = 4
    max_steps = 2_000_000

    D_r_values = [0.1, 0.01, 0.001]
    omega_scan = np.linspace(0.0, 1.0, 9)
    colors = ['#e41a1c', '#377eb8', '#4daf4a']
    markers = ['o', 's', '^']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

    for D_r, color, marker in zip(D_r_values, colors, markers):
        print(f"\n  MFPT vs ω for D_r = {D_r}:")
        mfpt_mean = []
        mfpt_err = []

        for omega in omega_scan:
            all_times = []
            for maze_idx in range(N_mazes):
                maze, entrance, exit_pos = generate_maze_prim(L, seed=1000 + maze_idx)
                for trial in range(N_trials_per_maze):
                    d_init = trial % 4
                    t, _ = run_tcrw_maze(maze, entrance, exit_pos, omega, D_r,
                                         d_init=d_init, max_steps=max_steps,
                                         seed=2000 + maze_idx * 100 + trial)
                    if t > 0:
                        all_times.append(t)

            if len(all_times) > 0:
                m = np.mean(all_times)
                e = np.std(all_times) / np.sqrt(len(all_times))
                mfpt_mean.append(m)
                mfpt_err.append(e)
                print(f"    ω={omega:.2f}: τ_M = {m:.0f} ± {e:.0f} "
                      f"({len(all_times)} successes)")
            else:
                mfpt_mean.append(np.nan)
                mfpt_err.append(np.nan)
                print(f"    ω={omega:.2f}: DNF")

        mfpt_mean = np.array(mfpt_mean)
        mfpt_err = np.array(mfpt_err)

        # Top: raw
        ax1.errorbar(omega_scan, mfpt_mean, yerr=mfpt_err,
                     fmt=f'{marker}-', color=color, ms=6, lw=1.5,
                     capsize=3, label=f'$D_r = {D_r}$')

        # Bottom: rescaled by 1/D_r
        valid = ~np.isnan(mfpt_mean)
        ax2.errorbar(omega_scan[valid], mfpt_mean[valid] * D_r,
                     yerr=mfpt_err[valid] * D_r,
                     fmt=f'{marker}-', color=color, ms=6, lw=1.5,
                     capsize=3, label=f'$D_r = {D_r}$')

    ax1.set_ylabel(r'$\tau_M$', fontsize=13)
    ax1.set_yscale('log')
    ax1.legend(fontsize=10)
    ax1.set_title(f'MFPT vs chirality ($L={L}$, {N_mazes} mazes)', fontsize=12)

    ax2.set_xlabel(r'$\omega$', fontsize=13)
    ax2.set_ylabel(r'$\tau_M \cdot D_r$', fontsize=13)
    ax2.set_yscale('log')
    ax2.legend(fontsize=10)
    ax2.set_title('Rescaled by $D_r$', fontsize=12)

    plt.suptitle(r'Fig 5(c): Maze-solving time $\tau_M$ vs $\omega$',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('tcrw_fig5c_mfpt_vs_omega.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("\nSaved: tcrw_fig5c_mfpt_vs_omega.png")


def fig5d_mfpt_vs_L():
    """
    Fig 5(d): MFPT τ_M vs maze size L.

    Key result: τ_M ~ L² for chiral (ω→0,1), τ_M ~ L³ for achiral (ω=0.5).
    """
    D_r = 0.01
    N_mazes = 8
    N_trials_per_maze = 3
    max_steps = 2_000_000

    L_values = [5, 7, 10, 15]
    omega_values = [0.5, 0.8, 1.0]
    colors = ['#377eb8', '#ff7f00', '#e41a1c']

    fig, ax = plt.subplots(figsize=(8, 6))

    for omega, color in zip(omega_values, colors):
        print(f"\n  τ_M vs L for ω = {omega}:")
        means = []
        errs = []
        L_good = []

        for L in L_values:
            all_times = []
            for maze_idx in range(N_mazes):
                maze, entrance, exit_pos = generate_maze_prim(L, seed=3000 + maze_idx)
                for trial in range(N_trials_per_maze):
                    d_init = trial % 4
                    t, _ = run_tcrw_maze(maze, entrance, exit_pos, omega, D_r,
                                         d_init=d_init, max_steps=max_steps,
                                         seed=4000 + maze_idx * 100 + trial)
                    if t > 0:
                        all_times.append(t)

            if len(all_times) > 5:
                m = np.mean(all_times)
                e = np.std(all_times) / np.sqrt(len(all_times))
                means.append(m)
                errs.append(e)
                L_good.append(L)
                print(f"    L={L}: τ_M = {m:.0f} ± {e:.0f} "
                      f"({len(all_times)}/{N_mazes * N_trials_per_maze})")
            else:
                print(f"    L={L}: too few successes ({len(all_times)})")

        if len(L_good) > 1:
            ax.errorbar(L_good, means, yerr=errs, fmt='o-', color=color,
                        ms=6, lw=1.5, capsize=3, label=f'$\\omega={omega}$')

    # Reference lines — anchor to actual data if available
    L_ref = np.logspace(np.log10(4), np.log10(25), 50)
    ax.plot(L_ref, 3e3 * (L_ref / 10)**2, 'k--', lw=1.5, alpha=0.5,
            label=r'$\propto L^2$')
    ax.plot(L_ref, 8e4 * (L_ref / 10)**3, 'k:', lw=1.5, alpha=0.5,
            label=r'$\propto L^3$')

    ax.set_xlabel('$L$ (maze size)', fontsize=13)
    ax.set_ylabel(r'$\tau_M$', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=9, ncol=2)
    ax.set_title(f'Fig 5(d): MFPT scaling with maze size ($D_r={D_r}$)',
                 fontsize=13)

    plt.tight_layout()
    plt.savefig('tcrw_fig5d_mfpt_vs_L.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("\nSaved: tcrw_fig5d_mfpt_vs_L.png")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("=" * 60)
    print("TCRW Phase 6B: Maze Solving")
    print("=" * 60)

    # Quick maze sanity check
    print("\n--- Sanity check: maze generation ---")
    L = 5
    maze, entrance, exit_pos = generate_maze_prim(L, seed=42)
    print(f"  Maze size: {maze.shape}")
    print(f"  Entrance: {entrance}, Exit: {exit_pos}")
    print(f"  Open cells: {np.sum(maze == 0)}, Wall cells: {np.sum(maze == 1)}")

    # Check maze is solvable (BFS)
    from collections import deque
    visited = set()
    queue = deque([entrance])
    visited.add(entrance)
    while queue:
        cx, cy = queue.popleft()
        if (cx, cy) == exit_pos:
            print("  Maze is solvable: YES")
            break
        for dx, dy in [(-1,0),(1,0),(0,-1),(0,1)]:
            nx, ny = cx+dx, cy+dy
            if (0 <= nx < maze.shape[0] and 0 <= ny < maze.shape[1]
                    and maze[nx, ny] == 0 and (nx, ny) not in visited):
                visited.add((nx, ny))
                queue.append((nx, ny))
    else:
        print("  Maze is solvable: NO — BUG!")

    # Quick MFPT test
    print("\n--- Quick MFPT test ---")
    for omega in [0.0, 0.5, 1.0]:
        t, _ = run_tcrw_maze(maze, entrance, exit_pos, omega, D_r=0.01,
                             d_init=1, max_steps=1_000_000, seed=42)
        print(f"  ω={omega}: {t} steps")

    # Figures
    print("\n--- Fig 5(a): Maze trajectories ---")
    fig5a_trajectories()

    print("\n--- Fig 5(c): MFPT vs ω ---")
    fig5c_mfpt_vs_omega()

    print("\n--- Fig 5(d): MFPT vs L ---")
    fig5d_mfpt_vs_L()

    print("\n" + "=" * 60)
    print("Phase 6B complete.")
    print("=" * 60)
