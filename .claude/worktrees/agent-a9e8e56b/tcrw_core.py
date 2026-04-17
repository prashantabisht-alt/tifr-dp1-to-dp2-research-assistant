"""
TCRW Core Simulation Engine
============================
Topological Chiral Random Walker on a 2D square lattice.

Model rules (Osat, Meyberg, Metson & Speck, arXiv:2602.12020):
  State: (x, y, d) where d ∈ {0,1,2,3} = {↑, →, ↓, ←}  (CCW ordering)

  At each discrete time step:
    With prob D_r:  NOISE STEP
      - Walker stays at (x,y)
      - Director rotates: CCW (d -> d+1 mod 4) with prob ω
                          CW  (d -> d-1 mod 4) with prob (1-ω)

    With prob (1-D_r):  CHIRAL STEP
      - Walker translates one step in direction d
      - Then director rotates: CW  (d -> d-1 mod 4) with prob ω
                               CCW (d -> d+1 mod 4) with prob (1-ω)
      NOTE: rotation chirality in chiral step is OPPOSITE to noise step!

      OBC special rule: if translation would move walker off-grid,
      the entire chiral step is blocked (no translation AND no rotation).

  Direction encoding:
    d=0: ↑  (y+1)
    d=1: →  (x+1)
    d=2: ↓  (y-1)
    d=3: ←  (x-1)

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np

# Direction vectors: dx[d], dy[d] for d = 0(↑), 1(→), 2(↓), 3(←)
DX = np.array([0, 1, 0, -1], dtype=np.int64)
DY = np.array([1, 0, -1, 0], dtype=np.int64)


def simulate_tcrw_pbc(omega, D_r, L, T_steps, N_traj=1, seed=42,
                       record_traj=False, record_interval=1):
    """
    Simulate TCRW on L×L lattice with periodic boundary conditions.

    Parameters
    ----------
    omega : float
        Chirality parameter, 0 <= omega <= 1.
        omega=1: chiral step rotates CW, noise rotates CCW
        omega=0: chiral step rotates CCW, noise rotates CW
        omega=0.5: achiral
    D_r : float
        Probability of noise step (rotational diffusion), 0 <= D_r <= 1.
    L : int
        Lattice size (L×L grid, sites 0..L-1).
    T_steps : int
        Number of discrete time steps.
    N_traj : int
        Number of independent trajectories (vectorized).
    seed : int
        Random seed.
    record_traj : bool
        If True, record full (x,y) trajectory of the FIRST walker.
    record_interval : int
        Record trajectory every this many steps.

    Returns
    -------
    result : dict with keys:
        'msd' : array (T_steps+1,) — mean square displacement averaged over trajectories
                 MSD(t) = <(x(t)-x(0))^2 + (y(t)-y(0))^2>
        'traj_x', 'traj_y' : arrays of recorded trajectory (if record_traj=True)
        'final_x', 'final_y', 'final_d' : final positions and directors
    """
    rng = np.random.default_rng(seed)

    # Initial state: all walkers start at (L//2, L//2) with random director
    x = np.full(N_traj, L // 2, dtype=np.int64)
    y = np.full(N_traj, L // 2, dtype=np.int64)
    d = rng.integers(0, 4, size=N_traj)

    # Track unwrapped displacement for MSD (don't mod these)
    ux = np.zeros(N_traj, dtype=np.int64)  # unwrapped x
    uy = np.zeros(N_traj, dtype=np.int64)  # unwrapped y

    # MSD storage
    n_records = T_steps // record_interval + 1
    msd = np.zeros(n_records, dtype=np.float64)
    times = np.zeros(n_records, dtype=np.int64)

    # Trajectory storage (first walker only)
    if record_traj:
        traj_x = np.zeros(n_records, dtype=np.int64)
        traj_y = np.zeros(n_records, dtype=np.int64)
        traj_x[0] = ux[0]
        traj_y[0] = uy[0]

    rec_idx = 1  # next record slot

    for t in range(1, T_steps + 1):
        # Draw random numbers for all walkers
        r_step = rng.random(N_traj)       # noise vs chiral
        r_rot  = rng.random(N_traj)       # rotation direction

        is_noise = r_step < D_r
        is_chiral = ~is_noise

        # --- Noise step: stay put, rotate ---
        # With encoding d=0(↑),1(→),2(↓),3(←):
        #   CW  = d+1 mod 4  (angle decreases: ↑→→→↓→←)
        #   CCW = d-1 mod 4  (angle increases: ↑→←→↓→→)
        # Paper: noise rotates CCW with prob omega, CW with prob (1-omega)
        noise_ccw = is_noise & (r_rot < omega)
        noise_cw  = is_noise & (r_rot >= omega)
        d[noise_ccw] = (d[noise_ccw] - 1) % 4   # CCW = d-1
        d[noise_cw]  = (d[noise_cw] + 1) % 4    # CW  = d+1

        # --- Chiral step: translate in direction d, then rotate ---
        # Translation
        step_dx = DX[d]  # displacement for each walker's current director
        step_dy = DY[d]

        # PBC: just translate and wrap
        x[is_chiral] = (x[is_chiral] + step_dx[is_chiral]) % L
        y[is_chiral] = (y[is_chiral] + step_dy[is_chiral]) % L
        ux[is_chiral] += step_dx[is_chiral]
        uy[is_chiral] += step_dy[is_chiral]

        # Rotation after chiral step: OPPOSITE chirality to noise
        # Paper: chiral rotates CW with prob omega, CCW with prob (1-omega)
        chiral_cw  = is_chiral & (r_rot < omega)
        chiral_ccw = is_chiral & (r_rot >= omega)
        d[chiral_cw]  = (d[chiral_cw] + 1) % 4   # CW  = d+1
        d[chiral_ccw] = (d[chiral_ccw] - 1) % 4   # CCW = d-1

        # Record
        if t % record_interval == 0:
            r2 = ux.astype(np.float64)**2 + uy.astype(np.float64)**2
            msd[rec_idx] = np.mean(r2)
            times[rec_idx] = t
            if record_traj:
                traj_x[rec_idx] = ux[0]
                traj_y[rec_idx] = uy[0]
            rec_idx += 1

    result = {
        'msd': msd[:rec_idx],
        'times': times[:rec_idx],
        'final_x': x, 'final_y': y, 'final_d': d
    }
    if record_traj:
        result['traj_x'] = traj_x[:rec_idx]
        result['traj_y'] = traj_y[:rec_idx]
    return result


def simulate_tcrw_obc(omega, D_r, L, T_steps, N_traj=1, seed=42,
                       record_traj=False, record_interval=1,
                       track_visits=False, track_currents=False,
                       track_step_type=False):
    """
    Simulate TCRW on L×L lattice with open (hard-wall) boundary conditions.

    OBC rule: lattice sites are 0..L-1. Boundary = sites that the walker
    cannot cross. If a chiral step would move the walker to x<0, x>=L,
    y<0, or y>=L, the ENTIRE chiral step is blocked (no move, no rotation).

    Parameters
    ----------
    (same as PBC, plus:)
    track_visits : bool
        If True, accumulate visit histogram P(x,y) over all trajectories.
    track_currents : bool
        If True, track net displacement current J(x,y) and its decomposition.
    track_step_type : bool
        If True, track whether previous step was noise (needed for J decomposition).

    Returns
    -------
    result : dict with various arrays depending on flags.
    """
    rng = np.random.default_rng(seed)

    # Initial state: random position inside the grid, random director
    x = rng.integers(0, L, size=N_traj)
    y = rng.integers(0, L, size=N_traj)
    d = rng.integers(0, 4, size=N_traj)

    # For MSD: track displacement from start
    x0 = x.copy()
    y0 = y.copy()

    # Visit histogram
    if track_visits:
        visits = np.zeros((L, L), dtype=np.int64)

    # Current tracking: J_x(x,y), J_y(x,y) and decomposition
    if track_currents:
        Jx = np.zeros((L, L), dtype=np.float64)
        Jy = np.zeros((L, L), dtype=np.float64)
        # Decomposition
        Jx_Dr = np.zeros((L, L), dtype=np.float64)
        Jy_Dr = np.zeros((L, L), dtype=np.float64)
        Jx_omega = np.zeros((L, L), dtype=np.float64)
        Jy_omega = np.zeros((L, L), dtype=np.float64)
        prev_was_noise = np.zeros(N_traj, dtype=bool)

    # MSD storage
    n_records = T_steps // record_interval + 1
    msd = np.zeros(n_records, dtype=np.float64)
    times = np.zeros(n_records, dtype=np.int64)

    if record_traj:
        traj_x = np.zeros(n_records, dtype=np.int64)
        traj_y = np.zeros(n_records, dtype=np.int64)
        traj_x[0] = x[0]
        traj_y[0] = y[0]

    rec_idx = 1

    for t in range(1, T_steps + 1):
        r_step = rng.random(N_traj)
        r_rot  = rng.random(N_traj)

        is_noise = r_step < D_r
        is_chiral = ~is_noise

        # Accumulate visits BEFORE the step
        if track_visits:
            for i in range(N_traj):
                visits[x[i], y[i]] += 1

        # --- Noise step: CCW (d-1) with prob omega, CW (d+1) with prob (1-omega) ---
        noise_ccw = is_noise & (r_rot < omega)
        noise_cw  = is_noise & (r_rot >= omega)
        d[noise_ccw] = (d[noise_ccw] - 1) % 4   # CCW = d-1
        d[noise_cw]  = (d[noise_cw] + 1) % 4    # CW  = d+1

        # --- Chiral step with OBC ---
        if np.any(is_chiral):
            new_x = x + DX[d]
            new_y = y + DY[d]

            # Check which chiral walkers can actually move
            can_move = is_chiral & (new_x >= 0) & (new_x < L) & (new_y >= 0) & (new_y < L)
            blocked = is_chiral & ~can_move

            # Record current BEFORE moving (outflow from current position)
            if track_currents:
                for i in np.where(can_move)[0]:
                    ddx = DX[d[i]]
                    ddy = DY[d[i]]
                    Jx[x[i], y[i]] += ddx
                    Jy[x[i], y[i]] += ddy
                    if prev_was_noise[i]:
                        Jx_Dr[x[i], y[i]] += ddx
                        Jy_Dr[x[i], y[i]] += ddy
                    else:
                        Jx_omega[x[i], y[i]] += ddx
                        Jy_omega[x[i], y[i]] += ddy

            # Move
            x[can_move] = new_x[can_move]
            y[can_move] = new_y[can_move]

            # Rotate after chiral step: CW (d+1) with prob omega, CCW (d-1) with prob (1-omega)
            # Only for walkers that actually moved (blocked → no rotation)
            chiral_cw  = can_move & (r_rot < omega)
            chiral_ccw = can_move & (r_rot >= omega)
            d[chiral_cw]  = (d[chiral_cw] + 1) % 4   # CW  = d+1
            d[chiral_ccw] = (d[chiral_ccw] - 1) % 4   # CCW = d-1

        # Track previous step type for current decomposition
        if track_currents:
            prev_was_noise[:] = False
            prev_was_noise[is_noise] = True

        # Record
        if t % record_interval == 0:
            r2 = (x.astype(np.float64) - x0.astype(np.float64))**2 + \
                 (y.astype(np.float64) - y0.astype(np.float64))**2
            msd[rec_idx] = np.mean(r2)
            times[rec_idx] = t
            if record_traj:
                traj_x[rec_idx] = x[0]
                traj_y[rec_idx] = y[0]
            rec_idx += 1

    # Final visit count
    if track_visits:
        for i in range(N_traj):
            visits[x[i], y[i]] += 1

    result = {
        'msd': msd[:rec_idx],
        'times': times[:rec_idx],
        'final_x': x, 'final_y': y, 'final_d': d,
    }
    if record_traj:
        result['traj_x'] = traj_x[:rec_idx]
        result['traj_y'] = traj_y[:rec_idx]
    if track_visits:
        result['visits'] = visits
    if track_currents:
        total_steps = T_steps * N_traj
        result['Jx'] = Jx / total_steps
        result['Jy'] = Jy / total_steps
        result['Jx_Dr'] = Jx_Dr / total_steps
        result['Jy_Dr'] = Jy_Dr / total_steps
        result['Jx_omega'] = Jx_omega / total_steps
        result['Jy_omega'] = Jy_omega / total_steps
    return result


def measure_diffusion_coeff(omega, D_r, L=200, T_steps=500000, N_traj=500,
                             seed=42, fit_frac=0.5):
    """
    Measure the diffusion coefficient D from long-time MSD slope.
    MSD = 4Dt in 2D.  Returns D.

    fit_frac: use the last fit_frac fraction of the data for the linear fit.
    """
    res = simulate_tcrw_pbc(omega, D_r, L, T_steps, N_traj, seed,
                             record_interval=100)
    t = res['times'].astype(np.float64)
    msd = res['msd']

    # Linear fit to last half
    n = len(t)
    i0 = int(n * (1 - fit_frac))
    t_fit = t[i0:]
    msd_fit = msd[i0:]

    # D = slope / 4 from MSD = 4Dt
    slope = np.polyfit(t_fit, msd_fit, 1)[0]
    D = slope / 4.0
    return D


if __name__ == '__main__':
    # Quick sanity check
    print("=== TCRW Core Sanity Checks ===\n")

    # Check 1: omega=0.5 should give ordinary diffusion
    D = measure_diffusion_coeff(omega=0.5, D_r=0.001, L=200, T_steps=200000,
                                 N_traj=200, seed=42)
    print(f"omega=0.5, D_r=0.001: D = {D:.4f}")

    # Check 2: D_r=0, omega=0 => deterministic CCW orbit, period 4 on PBC
    res = simulate_tcrw_pbc(omega=0.0, D_r=0.0, L=10, T_steps=8, N_traj=1,
                             seed=42, record_traj=True, record_interval=1)
    print(f"\nD_r=0, omega=0 (det. CCW chiral rotor):")
    print(f"  traj_x: {res['traj_x']}")
    print(f"  traj_y: {res['traj_y']}")
    print(f"  Should trace CCW orbit with period 4")

    # Check 3: D_r=1 => pure spinor, no translation
    res = simulate_tcrw_pbc(omega=1.0, D_r=1.0, L=10, T_steps=1000, N_traj=100,
                             seed=42)
    print(f"\nD_r=1, omega=1 (pure spinor): final MSD = {res['msd'][-1]:.6f}")
    print(f"  Should be 0 (no translation)")

    print("\nDone.")
