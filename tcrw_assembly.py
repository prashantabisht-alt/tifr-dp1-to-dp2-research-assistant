"""
Self-Assembly of Patchy Tiles via TCRW
=======================================

Simulates the assembly of 25 unique patchy tiles into a 5×5 target structure
using Topological Chiral Random Walkers (TCRW).

Model:
  - 25 unique tiles arranged in target 5×5 grid with identities 0-24
  - Tiles interact via Hebbian (identity-specific) bonds
  - Each free tile performs TCRW until it reaches an adjacent binding site
  - Binding is irreversible
  - Total assembly time τ_SA is measured in TCRW steps

Key physics: Chiral walkers (ω≠0.5) follow cluster edge → faster assembly.

Reference: Osat, Meyberg, Metson & Speck, arXiv:2602.12020, Figure 6

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.collections import PatchCollection
import sys
import os

# Direction vectors: d = 0(↑), 1(→), 2(↓), 3(←)
DX = np.array([0, 1, 0, -1], dtype=np.int64)
DY = np.array([1, 0, -1, 0], dtype=np.int64)


class SelfAssembly:
    """
    Self-assembly of 25 unique tiles into 5×5 target structure.

    Target layout:
      00 01 02 03 04
      05 06 07 08 09
      10 11 12 13 14
      15 16 17 18 19
      20 21 22 23 24

    Tile 12 is at the center (2,2) and serves as the seed.
    """

    def __init__(self, L_target=5, L_arena=30):
        """
        Initialize assembly system.

        Parameters
        ----------
        L_target : int
            Size of target structure (L_target × L_target tiles)
        L_arena : int
            Size of simulation arena
        """
        self.L_target = L_target
        self.L_arena = L_arena
        self.n_tiles = L_target * L_target

        # Target structure: tile_id at each position
        # tile_id = row * L_target + col
        self.target = np.arange(self.n_tiles).reshape(L_target, L_target)

        # Build bond compatibility matrix
        # bonds[(tile_a, edge_d)] = tile_b means tile_a's edge d bonds to tile_b
        # Edge d matches arena direction d: 0=↑(y+1), 1=→(x+1), 2=↓(y-1), 3=←(x-1)
        #
        # Target grid rows map to arena y via: y = cy + (seed_r - r)
        # So row-1 → y+1 (north in arena), row+1 → y-1 (south in arena)
        # Target cols map to arena x via: x = cx + (c - seed_c)
        # So col+1 → x+1 (east in arena)
        self.bonds = {}

        for r in range(L_target):
            for c in range(L_target):
                tid = self.target[r, c]

                # North (d=0, ↑, y+1 in arena): neighbor at row-1 in target grid
                if r - 1 >= 0:
                    nb_tid = self.target[r - 1, c]
                    self.bonds[(tid, 0)] = nb_tid      # tile's north edge
                    self.bonds[(nb_tid, 2)] = tid      # neighbor's south edge

                # East (d=1, →, x+1 in arena): neighbor at col+1 in target grid
                if c + 1 < L_target:
                    nb_tid = self.target[r, c + 1]
                    self.bonds[(tid, 1)] = nb_tid      # tile's east edge
                    self.bonds[(nb_tid, 3)] = tid      # neighbor's west edge

    def get_interaction_matrix(self):
        """
        Return 25×25 binary interaction matrix I[i,j] = 1 if tiles i,j are neighbors.
        """
        I = np.zeros((self.n_tiles, self.n_tiles), dtype=int)
        for (i, _), j in self.bonds.items():
            I[i, j] = 1
        # Make symmetric
        I = I + I.T
        return I

    def run_assembly(self, omega, D_r, max_steps_per_tile=50000, seed=42):
        """
        Run complete assembly simulation.

        Parameters
        ----------
        omega : float
            Chirality parameter (0 <= omega <= 1)
        D_r : float
            Rotational diffusion (noise step probability)
        max_steps_per_tile : int
            Maximum TCRW steps allowed per tile before timeout
        seed : int
            Random seed

        Returns
        -------
        tau_SA : int
            Total TCRW steps for full assembly
        tile_trajectories : dict
            tile_id -> list of (x,y) positions during its walk
        placed : dict
            tile_id -> (arena_x, arena_y) placement position
        arena : array
            Final arena state
        """
        rng = np.random.default_rng(seed)

        # Arena: -1 = empty, >=0 = tile_id
        arena = np.full((self.L_arena, self.L_arena), -1, dtype=int)

        # Place seed tile (tile 12) at arena center
        cx, cy = self.L_arena // 2, self.L_arena // 2
        seed_tid = self.target[self.L_target//2, self.L_target//2]  # tile 12
        arena[cx, cy] = seed_tid

        placed = {seed_tid: (cx, cy)}
        # Target position of each tile in terms of offset from seed
        # seed_tid is at target position (L_target//2, L_target//2)
        seed_r, seed_c = self.L_target // 2, self.L_target // 2
        # Map tile_id -> (target_row, target_col)
        tid_to_target = {}
        for r in range(self.L_target):
            for c in range(self.L_target):
                tid_to_target[self.target[r, c]] = (r, c)

        # Build a queue: only release tiles that have at least one neighbor
        # already placed (BFS order from seed). This ensures binding is possible.
        from collections import deque
        release_order = []
        visited_tiles = {seed_tid}
        queue = deque([seed_tid])

        while queue:
            tid = queue.popleft()
            r, c = tid_to_target[tid]
            # Check all 4 target-grid neighbors
            for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < self.L_target and 0 <= nc < self.L_target:
                    nb_tid = self.target[nr, nc]
                    if nb_tid not in visited_tiles:
                        visited_tiles.add(nb_tid)
                        release_order.append(nb_tid)
                        queue.append(nb_tid)

        remaining = release_order  # BFS order guarantees bindable neighbors exist

        total_steps = 0
        tile_trajectories = {}

        for tile_idx, tile_id in enumerate(remaining):
            # Release tile near the cluster (proximity-based)
            placed_tiles = list(placed.values())
            anchor = placed_tiles[rng.integers(0, len(placed_tiles))]

            # Release within a ring around the cluster
            max_release_attempts = 100
            released = False
            for _ in range(max_release_attempts):
                angle = rng.random() * 2 * np.pi
                r_release = rng.integers(4, 12)
                sx = int(anchor[0] + r_release * np.cos(angle))
                sy = int(anchor[1] + r_release * np.sin(angle))

                if (0 <= sx < self.L_arena and 0 <= sy < self.L_arena and
                    arena[sx, sy] == -1):
                    x, y = sx, sy
                    released = True
                    break

            if not released:
                # Fallback: release at random empty position
                for _ in range(100):
                    sx = rng.integers(0, self.L_arena)
                    sy = rng.integers(0, self.L_arena)
                    if arena[sx, sy] == -1:
                        x, y = sx, sy
                        break

            d = rng.integers(0, 4)
            traj = [(x, y)]
            bound = False

            for step in range(max_steps_per_tile):
                # TCRW step
                r_step = rng.random()
                r_rot = rng.random()

                if r_step < D_r:
                    # NOISE STEP: rotate only
                    if r_rot < omega:
                        d = (d - 1) % 4  # CCW: ↑→←→↓→→→↑
                    else:
                        d = (d + 1) % 4  # CW:  ↑→→→↓→←→↑
                else:
                    # CHIRAL STEP: translate + rotate (opposite chirality)
                    nx = x + DX[d]
                    ny = y + DY[d]

                    # Check bounds and collision
                    if (0 <= nx < self.L_arena and 0 <= ny < self.L_arena and
                        arena[nx, ny] == -1):
                        x, y = nx, ny
                        # Rotate OPPOSITE to noise step
                        if r_rot < omega:
                            d = (d + 1) % 4  # CW (opposite of noise CCW)
                        else:
                            d = (d - 1) % 4  # CCW (opposite of noise CW)
                    # else: blocked (wall or occupied), no move no rotate

                traj.append((x, y))

                # Check all 4 neighbors for valid binding site
                for dd in range(4):
                    ax = x + DX[dd]
                    ay = y + DY[dd]

                    if 0 <= ax < self.L_arena and 0 <= ay < self.L_arena:
                        neighbor_tid = arena[ax, ay]
                        if neighbor_tid >= 0:
                            # There's a tile at (ax, ay)
                            # Our edge facing the neighbor is edge dd
                            # (dd points from us toward the neighbor)
                            facing_edge = dd

                            # Check if this is a valid bond
                            if (tile_id, facing_edge) in self.bonds:
                                if self.bonds[(tile_id, facing_edge)] == neighbor_tid:
                                    # Valid bond! Bind this tile
                                    arena[x, y] = tile_id
                                    placed[tile_id] = (x, y)
                                    bound = True
                                    break

                if bound:
                    total_steps += step + 1
                    break

            if not bound:
                total_steps += max_steps_per_tile

            tile_trajectories[tile_id] = traj

            if (tile_idx + 1) % 5 == 0:
                n_placed = len(placed)
                print(f"  Placed {n_placed}/25 tiles (τ_SA = {total_steps})")

        return total_steps, tile_trajectories, placed, arena


def plot_fig6a(assembly):
    """
    Panel (a): Tile catalog + target structure + interaction matrix
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Left: Tile catalog (5×5 grid of unique colors)
    ax = axes[0]
    L = assembly.L_target
    colors = plt.cm.tab20c(np.linspace(0, 1, assembly.n_tiles))

    # Create colored patches for each tile
    for r in range(L):
        for c in range(L):
            tid = assembly.target[r, c]
            rect = mpatches.Rectangle((c, L-1-r), 1, 1,
                                      facecolor=colors[tid],
                                      edgecolor='black', linewidth=1)
            ax.add_patch(rect)
            ax.text(c+0.5, L-1-r+0.5, str(tid), ha='center', va='center',
                   fontsize=8, fontweight='bold')

    ax.set_xlim(-0.5, L+0.5)
    ax.set_ylim(-0.5, L+0.5)
    ax.set_aspect('equal')
    ax.set_title('Tile Catalog', fontsize=12, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])

    # Center: Target structure with bonds
    ax = axes[1]

    # Draw tiles
    for r in range(L):
        for c in range(L):
            tid = assembly.target[r, c]
            rect = mpatches.Rectangle((c, L-1-r), 1, 1,
                                      facecolor=colors[tid],
                                      edgecolor='black', linewidth=1, alpha=0.7)
            ax.add_patch(rect)

    # Draw bonds
    for r in range(L):
        for c in range(L):
            tid = assembly.target[r, c]

            # Up neighbor
            if r + 1 < L:
                cx, cy = c + 0.5, L - r - 0.5
                ax.plot([cx, cx], [cy, cy - 1], 'r-', linewidth=2, alpha=0.6)

            # Right neighbor
            if c + 1 < L:
                cx, cy = c + 0.5, L - 1 - r + 0.5
                ax.plot([cx, cx + 1], [cy, cy], 'r-', linewidth=2, alpha=0.6)

    ax.set_xlim(-0.5, L+0.5)
    ax.set_ylim(-0.5, L+0.5)
    ax.set_aspect('equal')
    ax.set_title('Target Structure (5×5)', fontsize=12, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])

    # Right: Interaction matrix
    ax = axes[2]
    I = assembly.get_interaction_matrix()
    im = ax.imshow(I, cmap='Blues', aspect='auto', origin='upper')
    ax.set_title('Interaction Matrix', fontsize=12, fontweight='bold')
    ax.set_xlabel('Tile j')
    ax.set_ylabel('Tile i')
    plt.colorbar(im, ax=ax, label='Bond exists')

    plt.tight_layout()
    plt.savefig('tcrw_fig6a_tiles.png', dpi=150, bbox_inches='tight')
    print("Saved: tcrw_fig6a_tiles.png")
    plt.close()


def plot_fig6bc(assembly, omega_vals=[0.5, 1.0], D_r=0.01):
    """
    Panels (b)-(c): Assembly trajectories for achiral vs chiral
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    labels = ['Achiral (ω=0.5)', 'Chiral (ω=1.0)']

    for idx, omega in enumerate(omega_vals):
        ax = axes[idx]

        print(f"Simulating trajectory for ω={omega}, D_r={D_r}...")
        tau_sa, trajs, placed, arena = assembly.run_assembly(omega, D_r, seed=42)

        # Pick one interesting tile to visualize (not the seed)
        # Use tile 5 (near seed)
        tile_to_plot = 5
        if tile_to_plot in trajs:
            traj = trajs[tile_to_plot]
            xs, ys = zip(*traj)

            # Color by time
            times = np.arange(len(xs))
            colors_cmap = plt.cm.viridis(times / max(1, len(xs) - 1))

            ax.scatter(xs, ys, c=times, cmap='viridis', s=20, alpha=0.6)

            # Draw growing cluster (placed tiles)
            cluster_x = [pos[0] for tid, pos in placed.items()]
            cluster_y = [pos[1] for tid, pos in placed.items()]
            ax.scatter(cluster_x, cluster_y, c='red', s=100, marker='s',
                      alpha=0.3, label='Placed tiles')

        ax.set_xlim(-1, assembly.L_arena)
        ax.set_ylim(-1, assembly.L_arena)
        ax.set_aspect('equal')
        ax.set_title(labels[idx] + f'\nτ_SA={tau_sa}', fontsize=11, fontweight='bold')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        if idx == 0:
            ax.legend(loc='upper right', fontsize=9)

        cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='viridis',
                           norm=plt.Normalize(0, max(1, len(xs)-1))), ax=ax)
        cbar.set_label('Time step')

    plt.tight_layout()
    plt.savefig('tcrw_fig6bc_trajectories.png', dpi=150, bbox_inches='tight')
    print("Saved: tcrw_fig6bc_trajectories.png")
    plt.close()


def plot_fig6d(assembly, D_r_vals=None, omega_vals=None, n_trials=2):
    """
    Panel (d): τ_SA vs D_r for different ω (log-log)
    """
    if D_r_vals is None:
        D_r_vals = [0.01, 0.05, 0.1, 0.2]
    if omega_vals is None:
        omega_vals = [0.5, 0.8, 1.0]

    fig, ax = plt.subplots(figsize=(10, 7))

    colors_omega = plt.cm.RdYlBu_r(np.linspace(0, 1, len(omega_vals)))
    markers = ['o', 's', '^', 'D']

    for omega_idx, omega in enumerate(omega_vals):
        tau_sa_vals = []
        tau_sa_errs = []

        for D_r in D_r_vals:
            print(f"  Running τ_SA vs D_r: ω={omega}, D_r={D_r}")

            tau_trials = []
            for trial in range(n_trials):
                tau_sa, _, _, _ = assembly.run_assembly(omega, D_r, seed=42+trial)
                tau_trials.append(tau_sa)

            tau_avg = np.mean(tau_trials)
            tau_std = np.std(tau_trials)

            tau_sa_vals.append(tau_avg)
            tau_sa_errs.append(tau_std)

        tau_sa_vals = np.array(tau_sa_vals)
        tau_sa_errs = np.array(tau_sa_errs)

        ax.errorbar(D_r_vals, tau_sa_vals, yerr=tau_sa_errs,
                   marker=markers[omega_idx], markersize=8,
                   color=colors_omega[omega_idx],
                   label=f'ω={omega}', linewidth=2, capsize=5)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Rotational Diffusion $D_r$', fontsize=12, fontweight='bold')
    ax.set_ylabel('Assembly Time $\\tau_{SA}$ (steps)', fontsize=12, fontweight='bold')
    ax.set_title('Self-Assembly Time vs Diffusion (Fig 6d)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=11, loc='best')
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig('tcrw_fig6d_tau_vs_Dr.png', dpi=150, bbox_inches='tight')
    print("Saved: tcrw_fig6d_tau_vs_Dr.png")
    plt.close()


def plot_fig6e(assembly, omega_vals=None, D_r=0.01, n_trials=2):
    """
    Panel (e): τ_SA vs ω
    """
    if omega_vals is None:
        omega_vals = np.linspace(0, 1, 9)

    fig, ax = plt.subplots(figsize=(10, 7))

    print(f"Running τ_SA vs ω: D_r={D_r}")

    tau_sa_vals = []
    tau_sa_errs = []

    for omega in omega_vals:
        print(f"  ω={omega:.2f}")

        tau_trials = []
        for trial in range(n_trials):
            tau_sa, _, _, _ = assembly.run_assembly(omega, D_r, seed=42+trial)
            tau_trials.append(tau_sa)

        tau_avg = np.mean(tau_trials)
        tau_std = np.std(tau_trials)

        tau_sa_vals.append(tau_avg)
        tau_sa_errs.append(tau_std)

    tau_sa_vals = np.array(tau_sa_vals)
    tau_sa_errs = np.array(tau_sa_errs)

    ax.errorbar(omega_vals, tau_sa_vals, yerr=tau_sa_errs,
               marker='o', markersize=8, color='steelblue',
               linewidth=2, capsize=5)

    ax.set_xlabel('Chirality Parameter $\\omega$', fontsize=12, fontweight='bold')
    ax.set_ylabel('Assembly Time $\\tau_{SA}$ (steps)', fontsize=12, fontweight='bold')
    ax.set_title('Self-Assembly Time vs Chirality (Fig 6e)', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.05, 1.05)

    plt.tight_layout()
    plt.savefig('tcrw_fig6e_tau_vs_omega.png', dpi=150, bbox_inches='tight')
    print("Saved: tcrw_fig6e_tau_vs_omega.png")
    plt.close()


if __name__ == '__main__':
    print("="*70)
    print("TCRW Self-Assembly Simulation (Figure 6 from arXiv:2602.12020)")
    print("="*70)
    print()

    # Initialize assembly system
    print("Initializing 5×5 assembly system...")
    assembly = SelfAssembly(L_target=5, L_arena=30)
    print(f"  Target: {assembly.L_target}×{assembly.L_target} = {assembly.n_tiles} tiles")
    print(f"  Arena: {assembly.L_arena}×{assembly.L_arena}")
    print(f"  Bonds: {len(assembly.bonds)//2} unique bonds")
    print()

    # Generate Fig 6(a): Tiles + target + interaction matrix
    print("Generating Fig 6(a): Tile catalog, target structure, interaction matrix...")
    plot_fig6a(assembly)
    print()

    # Generate Fig 6(b)-(c): Trajectories
    print("Generating Fig 6(b)-(c): Assembly trajectories...")
    plot_fig6bc(assembly, omega_vals=[0.5, 1.0], D_r=0.01)
    print()

    # Generate Fig 6(d): τ_SA vs D_r
    print("Generating Fig 6(d): τ_SA vs D_r...")
    print("(This will take a few minutes...)")
    plot_fig6d(assembly,
               D_r_vals=[0.01, 0.05, 0.1, 0.2],
               omega_vals=[0.5, 0.8, 1.0],
               n_trials=2)
    print()

    # Generate Fig 6(e): τ_SA vs ω
    print("Generating Fig 6(e): τ_SA vs ω...")
    print("(This will take a few minutes...)")
    plot_fig6e(assembly,
               omega_vals=np.linspace(0, 1, 9),
               D_r=0.01,
               n_trials=2)
    print()

    print("="*70)
    print("All figures generated successfully!")
    print("  - tcrw_fig6a_tiles.png")
    print("  - tcrw_fig6bc_trajectories.png")
    print("  - tcrw_fig6d_tau_vs_Dr.png")
    print("  - tcrw_fig6e_tau_vs_omega.png")
    print("="*70)
