import random
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
from matplotlib.colors import LogNorm
import argparse
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import os


# use these codes for visual arrows!
#directions = ['\u2191', '\u2192', '\u2193', '\u2190']

# Directions and their mapping
directions = ['↑', '→', '↓', '←']
dir_to_vec = {
    '↑': (0, 1),
    '→': (1, 0),
    '↓': (0, -1),
    '←': (-1, 0)
}


def rotate_direction(current, clockwise=True):
    idx = directions.index(current)
    if clockwise:
        return directions[(idx + 1) % 4]
    else:
        return directions[(idx - 1) % 4]

class ChiralWalker:
    def __init__(self, i=0, j=0, direction='↑', D_r=0.1, omega=1.0, L=6, defects=None, record_path=False):
        self.i = i
        self.j = j
        self.d = direction
        self.D_r = D_r
        self.omega = omega
        self.L = L
        self.record_path = record_path
        self.path = [(self.i, self.j)]
        self.P = np.zeros((L + 1, L + 1), dtype=np.float64)
        self.J1 = np.zeros((L + 1, L + 1, 2), dtype=np.float64)
        self.J2 = np.zeros((L + 1, L + 1, 2), dtype=np.float64)
        self.noise_step=False
        self.steps = 0

        # Initialize all sites as accessible
        self.accessible = np.ones((L + 1, L + 1), dtype=bool)

        # mark defects as inaccessible
        if defects is not None:
            for (x, y) in defects:
                if 0 <= x <= L and 0 <= y <= L:
                    self.accessible[x, y] = False
        
    def in_bounds(self, i, j):
        return 0 <= i <= self.L and 0 <= j <= self.L and self.accessible[i, j]

    def calculate_currents (self, new_i, new_j):
        dx, dy = (new_i - self.i), (new_j - self.j)
        if self.noise_step:
            self.J1[self.i][self.j][0] += dx
            self.J1[self.i][self.j][1] += dy
        else:
            self.J2[self.i][self.j][0] += dx
            self.J2[self.i][self.j][1] += dy
        self.steps+=1
        
    def step(self):
        if random.random() < self.D_r:
            # Rotational noise step: only rotate (in ooposite chirality)
            clockwise = random.random() > self.omega
            self.d = rotate_direction(self.d, clockwise=clockwise)
            self.P[self.i,self.j]+=1
            self.noise_step=True
        else:
            # Chiral move: move AND rotation (if move is blocked then rotation will not happen)
            di, dj = dir_to_vec[self.d]
            new_i = self.i + di
            new_j = self.j + dj

            if self.in_bounds(new_i, new_j):

                self.calculate_currents(new_i, new_j)
                
                # Successful move
                self.i = new_i
                self.j = new_j
                
                if self.record_path:
                    self.path.append((self.i, self.j))

                # Rotate after move
                clockwise = random.random() < self.omega  # Chiral bias
                self.d = rotate_direction(self.d, clockwise=clockwise)
                self.noise_step=False
            else:
                # Move blocked by boundary: NO move AND NO rotation
                if self.record_path:
                    self.path.append((self.i, self.j))  # Still record time step

            self.P[self.i,self.j]+=1

######################################################################
# Steady States Analysis
######################################################################

def index(i, j, d, L):
    #Map (i,j,direction) to a unique state index.
    return (i * (L + 1) + j) * len(directions) + directions.index(d)

def build_sparse_transition_matrix(L, omega, D_r, defects=None):
    """
    Constructs the Discrete Transition Matrix W with Defects.
    W[i, j] represents the Probability of State j -> State i.
    (Columns sum to 1).
    """
    num_states = (L + 1) * (L + 1) * len(directions)
    defect_set = set(defects) if defects else set()
    
    # We build the matrix in (Destination, Source) format
    row = [] # Destinations
    col = [] # Sources
    data = []
    
    for i in range(L + 1):
        for j in range(L + 1):
            
            # If the SOURCE site is a defect, it's a dummy state.
            # We give it a self-loop of 1.0 to keep the matrix stochastic.
            # It will be disconnected from the rest.
            if (i, j) in defect_set:
                for d in directions:
                    idx = index(i, j, d, L)
                    row.append(idx)
                    col.append(idx)
                    data.append(1.0)
                continue

            for d in directions:
                src_idx = index(i, j, d, L)

                ### 1. Rotational Noise Step (Prob D_r)
                # Simulation: clockwise = random.random() > omega
                # P(CW) = 1 - omega, P(CCW) = omega
                
                cw_dir = rotate_direction(d, clockwise=True)
                ccw_dir = rotate_direction(d, clockwise=False)

                # Source -> CW (Dest)
                row.append(index(i, j, cw_dir, L))
                col.append(src_idx)
                data.append(D_r * (1 - omega))

                # Source -> CCW (Dest)
                row.append(index(i, j, ccw_dir, L))
                col.append(src_idx)
                data.append(D_r * omega)

                ### 2. Chiral Move Step (Prob 1 - D_r)
                # Simulation: clockwise = random.random() < omega
                # P(CW) = omega, P(CCW) = 1 - omega
                
                di, dj = dir_to_vec[d]
                new_i, new_j = i + di, j + dj
                
                # Check if move is VALID (In bounds AND not a defect)
                is_blocked = False
                if not (0 <= new_i <= L and 0 <= new_j <= L):
                    is_blocked = True
                elif (new_i, new_j) in defect_set:
                    is_blocked = True

                if not is_blocked:
                    # Valid Move: Move THEN Rotate
                    dest_cw = index(new_i, new_j, cw_dir, L)
                    dest_ccw = index(new_i, new_j, ccw_dir, L)
                    
                    # Source -> Move -> CW
                    row.append(dest_cw)
                    col.append(src_idx)
                    data.append((1 - D_r) * omega)
                    
                    # Source -> Move -> CCW
                    row.append(dest_ccw)
                    col.append(src_idx)
                    data.append((1 - D_r) * (1 - omega))
                else:
                    # Blocked (Boundary or Defect): NO Move AND NO Rotation
                    # The walker stays in state 'src_idx'
                    row.append(src_idx)
                    col.append(src_idx)
                    data.append(1 - D_r)

    # Create COO sparse matrix (Dest, Source)
    W = sp.coo_matrix((data, (row, col)), shape=(num_states, num_states))
    
    # Convert to CSR for efficiency
    W = W.tocsr()
    
    return W

def solve_steady_state_sparse(W, L, sigma=1, defects=None):
    """
    Compute steady-state distribution P such that W P = P.
    """
    # Solve for Eigenvalue 1
    vals, vecs = spla.eigs(W, k=1, sigma=sigma+1e-8, which='LM')
    
    steady_state = np.real(vecs[:, 0])
    
    # --- Defect Cleanup ---
    # Ensure defect sites have exactly 0 probability 
    if defects:
        defect_set = set(defects)
        for i, j in defect_set:
            if 0 <= i <= L and 0 <= j <= L:
                for d in directions:
                    idx = index(i, j, d, L)
                    steady_state[idx] = 0.0

    # Normalize
    steady_state = steady_state / np.sum(steady_state)
    
    # Ensure no negative small values
    steady_state[steady_state < 0] = 0
    steady_state = steady_state / np.sum(steady_state)
    
    return steady_state, np.abs(vals[0])

def calculate_J1_J2_with_boundaries(P_ss, W_sparse, L, D_r, omega, defects=None):
    
    #Decomposes currents into J_Dr (Rotation origin) and J_omega (Move origin).
    
    num_states = (L + 1) * (L + 1) * 4
    defect_set = set(defects) if defects else set()
    
    P_rot = np.zeros(num_states, dtype=np.float64)
    P_chiral = np.zeros(num_states, dtype=np.float64)

    W_coo = W_sparse.tocoo()
    
    # Step 1: Infer Flux flavor (Rotation or Chiral)
    for idx_dest, idx_src, prob in zip(W_coo.row, W_coo.col, W_coo.data):
        if idx_src == idx_dest: continue # Ignore self-loops (Blocked moves)
        
        i_dest = (idx_dest // 4) // (L + 1)
        j_dest = (idx_dest // 4) % (L + 1)
        i_src  = (idx_src  // 4) // (L + 1)
        j_src  = (idx_src  // 4) % (L + 1)

        # Skip if source or dest is a defect (Should represent 0 flux anyway)
        if (i_src, j_src) in defect_set or (i_dest, j_dest) in defect_set:
            continue

        flux = P_ss[idx_src] * prob

        if i_src == i_dest and j_src == j_dest:
            # Spatial position didn't change -> Pure Rotation
            P_rot[idx_dest] += flux
        else:
            # Spatial position changed -> Chiral Move
            P_chiral[idx_dest] += flux

    # Normalize flavors
    total_flux_in = P_rot + P_chiral
    mask = total_flux_in > 1e-12
    
    P_rot[mask] = P_ss[mask] * (P_rot[mask] / total_flux_in[mask])
    P_chiral[mask] = P_ss[mask] * (P_chiral[mask] / total_flux_in[mask])
    
    P_rot[~mask] = 0
    P_chiral[~mask] = 0

    # Step 2: Compute Outgoing Currents
    J1 = np.zeros((L + 1, L + 1, 2), dtype=np.float64)
    J2 = np.zeros((L + 1, L + 1, 2), dtype=np.float64)

    for idx_dest, idx_src, prob in zip(W_coo.row, W_coo.col, W_coo.data):
        if idx_src == idx_dest: continue

        i_src  = (idx_src  // 4) // (L + 1)
        j_src  = (idx_src  // 4) % (L + 1)
        i_dest = (idx_dest // 4) // (L + 1)
        j_dest = (idx_dest // 4) % (L + 1)
        
        # Skip non-spatial transitions and defects
        if (i_src == i_dest and j_src == j_dest): continue
        if (i_src, j_src) in defect_set: continue

        dx, dy = i_dest - i_src, j_dest - j_src
        disp = np.array([dx, dy], dtype=np.float64)

        J1[i_src, j_src] += P_rot[idx_src] * prob * disp
        J2[i_src, j_src] += P_chiral[idx_src] * prob * disp

    return J1, J2

########################################################

def compute_P_bulk(P_ss, L):
    # Sum probabilities for all interior (i, j) positions, excluding the boundaries
    bulk_prob = 0
    num_bulk_sites = 0
    
    for i in range(1, L):
        for j in range(1, L):
            # Sum over all directions for each (i, j) position
            bulk_prob += sum(P_ss[index(i, j, d, L)] for d in directions)
            num_bulk_sites += 1
    
    # Average probability at bulk sites
    P_bulk = bulk_prob / num_bulk_sites if num_bulk_sites > 0 else 0
    return P_bulk


def compute_P_edge(P_ss, L):
    # Sum probabilities for all edge (i, j) positions
    edge_prob = 0
    num_edge_sites = 0
    
    # Check all boundary positions
    for i in [0, L]:
        for j in range(L + 1):
            # Sum over all directions for each (i, j) position
            edge_prob += sum(P_ss[index(i, j, d, L)] for d in directions)
            num_edge_sites += 1
    
    for j in [0, L]:
        for i in range(1, L):  # Avoid double-counting the corners
            edge_prob += sum(P_ss[index(i, j, d, L)] for d in directions)
            num_edge_sites += 1
    
    # Average probability at edge sites
    P_edge = edge_prob / num_edge_sites if num_edge_sites > 0 else 0
    return P_edge

################### Visualization ######################

def visualize_steady_state(P_ss, L, defects=None):
    """
    Visualize steady-state probability (summed over directions)
    with optional defect sites masked out.
    """
    P_grid = np.zeros((L + 1, L + 1))
    defect_set = set(defects or [])

    # Build 2D grid (sum over directions, skip defects)
    for i in range(L + 1):
        for j in range(L + 1):
            if (i, j) in defect_set:
                P_grid[i, j] = np.nan  # mark defects for visualization
            else:
                P_grid[i, j] = sum(P_ss[index(i, j, d, L)] for d in directions)

    # Avoid log(0) by masking zeros and defects
    valid = np.isfinite(P_grid) & (P_grid > 0)
    vmin = np.min(P_grid[valid])
    vmax = np.max(P_grid[valid])

    plt.figure(figsize=(6, 5))
    plt.imshow(P_grid.T, origin="lower", cmap="hot",
               interpolation="nearest", norm=LogNorm(vmin=np.min(P_grid[P_grid > 0]), vmax=np.max(P_grid[P_grid>0])))

    plt.colorbar(label="Log(Probability)")
    
    # Mark defect regions (optional visual)
    if defects:
        for (i, j) in defect_set:
            plt.scatter(i, j, color="r", s=1)

    
    plt.title(f"Steady State Distribution (L={L})")
    plt.xlabel("i")
    plt.ylabel("j")
    plt.show()