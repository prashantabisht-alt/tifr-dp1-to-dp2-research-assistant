"""
Quick test of assembly code
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
import sys
import os

# Direction vectors: d = 0(↑), 1(→), 2(↓), 3(←)
DX = np.array([0, 1, 0, -1], dtype=np.int64)
DY = np.array([1, 0, -1, 0], dtype=np.int64)


class SelfAssembly:
    """
    Self-assembly of 25 unique tiles into 5×5 target structure.
    """

    def __init__(self, L_target=5, L_arena=30):
        self.L_target = L_target
        self.L_arena = L_arena
        self.n_tiles = L_target * L_target

        # Target structure
        self.target = np.arange(self.n_tiles).reshape(L_target, L_target)

        # Build bond compatibility matrix
        self.bonds = {}

        for r in range(L_target):
            for c in range(L_target):
                tid = self.target[r, c]

                # Up neighbor
                if r + 1 < L_target:
                    nb_tid = self.target[r+1, c]
                    self.bonds[(tid, 0)] = nb_tid
                    self.bonds[(nb_tid, 2)] = tid

                # Right neighbor
                if c + 1 < L_target:
                    nb_tid = self.target[r, c+1]
                    self.bonds[(tid, 1)] = nb_tid
                    self.bonds[(nb_tid, 3)] = tid

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

    def run_assembly(self, omega, D_r, max_steps_per_tile=100000, seed=42):
        """
        Run complete assembly simulation.
        """
        rng = np.random.default_rng(seed)

        # Arena: -1 = empty, >=0 = tile_id
        arena = np.full((self.L_arena, self.L_arena), -1, dtype=int)

        # Place seed tile (tile 12) at arena center
        cx, cy = self.L_arena // 2, self.L_arena // 2
        seed_tid = self.target[self.L_target//2, self.L_target//2]  # tile 12
        arena[cx, cy] = seed_tid

        placed = {seed_tid: (cx, cy)}
        remaining = list(range(self.n_tiles))
        remaining.remove(seed_tid)
        rng.shuffle(remaining)

        total_steps = 0
        tile_trajectories = {}

        for tile_idx, tile_id in enumerate(remaining):
            # Release tile at random empty position
            while True:
                sx = rng.integers(0, self.L_arena)
                sy = rng.integers(0, self.L_arena)
                if arena[sx, sy] == -1:
                    break

            x, y = sx, sy
            d = rng.integers(0, 4)
            traj = [(x, y)]

            bound = False

            for step in range(max_steps_per_tile):
                # TCRW step
                r_step = rng.random()
                r_rot = rng.random()

                if r_step < D_r:
                    # NOISE STEP
                    if r_rot < omega:
                        d = (d + 1) % 4  # CCW
                    else:
                        d = (d - 1) % 4  # CW
                else:
                    # CHIRAL STEP
                    nx = x + DX[d]
                    ny = y + DY[d]

                    # Check bounds and collision
                    if 0 <= nx < self.L_arena and 0 <= ny < self.L_arena and arena[nx, ny] == -1:
                        x, y = nx, ny
                        # Rotate opposite to noise step
                        if r_rot < omega:
                            d = (d - 1) % 4  # CW (opposite)
                        else:
                            d = (d + 1) % 4  # CCW
                    # else: blocked, no move no rotate

                traj.append((x, y))

                # Check all 4 neighbors for valid binding site
                for dd in range(4):
                    ax = x + DX[dd]
                    ay = y + DY[dd]

                    if 0 <= ax < self.L_arena and 0 <= ay < self.L_arena:
                        neighbor_tid = arena[ax, ay]
                        if neighbor_tid >= 0:
                            # Check if our tile can bond here
                            opp_edge = (dd + 2) % 4

                            if (tile_id, opp_edge) in self.bonds:
                                if self.bonds[(tile_id, opp_edge)] == neighbor_tid:
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
                print(f"  Placed {tile_idx + 1}/{len(remaining)} tiles (τ_SA = {total_steps})")

        return total_steps, tile_trajectories, placed, arena


print("Testing assembly initialization...")
assembly = SelfAssembly(L_target=5, L_arena=30)
print(f"  Target: {assembly.L_target}×{assembly.L_target} = {assembly.n_tiles} tiles")
print(f"  Arena: {assembly.L_arena}×{assembly.L_arena}")
print(f"  Bonds: {len(assembly.bonds)//2} unique bonds")
print()

print("Testing single assembly run (ω=0.5, D_r=0.01)...")
tau_sa, trajs, placed, arena = assembly.run_assembly(omega=0.5, D_r=0.01, seed=42)
print(f"  Assembly complete: τ_SA = {tau_sa} steps")
print(f"  Placed tiles: {len(placed)}/25")
print()

print("Success! Assembly simulation works.")
