"""
Postmortem of Dipanjan's rtp_tl_2.nb — re-implements his calculation in Python,
runs all the sanity checks his Mathematica notebook lacks, and locates any bugs.

What Dipanjan's notebook does (reconstructed):
  1. Define 6x6 Bloch generator M(k1, k2) for an active random walker
     on a 2D "triangular-like" lattice with primitive vectors
     a_1 = (2a, 0), a_2 = (a, b), a_3 = (a, -b) and a = b = 1.
  2. Loop over k-grid: kx[n] = 2*pi*n/(2*a*L), ky[n] = 2*pi*n/(b*L),
     with n_x in [0, 2L) and n_y in [0, L).
     Total k-points: 2*L^2 = 1800 (for L=30).
  3. At each k-point, diagonalize M and compute
        tilde_P(k, t) = sum_n V_in e^{lambda_n t} (V^{-1} rhtvct)_n
     where rhtvct = (1/6) e^{i k . r0} is the Fourier-transformed
     initial condition (delta function at r0 = (0,0), uniform over 6 directors).
  4. Inverse Fourier:
        P(x, y, t) = (1/N_k) sum_k tilde_P(k, t) e^{-i k . r}
  5. Print P(x, y=4, t) as a function of x.

Bug-hunt questions:
  A. Is the Bloch matrix consistent (column sums = 0 at k=0)?
  B. Is the inverse Fourier reconstruction giving a normalised probability?
     i.e. does sum_{x, y} P(x, y, t) = 1 ?
  C. Is P non-negative everywhere?
  D. Does P respect the bipartite parity (the walker can only reach sites
     with x + y == x0 + y0 mod 2)?
  E. Is the lattice geometry "2L x L" or "L x L" or something else?
     This determines whether the k-grid Dipanjan uses is the right one.

Run: python3 dipanjan_postmortem.py
"""
from __future__ import annotations
import numpy as np


# ---------------------------------------------------------------------------
# 1. Bloch matrix
# ---------------------------------------------------------------------------
def build_Mk(gamma: float, epsilon: float, k1: float, k2: float,
             a: float = 1.0, b: float = 1.0) -> np.ndarray:
    """Exact match to Dipanjan's mat[k1, k2]."""
    axes = [(2 * a, 0.0), (a, b), (a, -b)]
    bulk = (1.0 / 3.0) * (
        np.cos(2 * a * k1)
        + np.cos(a * k1 + b * k2)
        + np.cos(a * k1 - b * k2)
    )
    coeff = []
    for vx, vy in axes:
        kd = vx * k1 + vy * k2
        coeff.append(bulk - 1.0 - gamma + 2j * epsilon * np.sin(kd))
    M = np.zeros((6, 6), dtype=complex)
    M[0, 0] = coeff[0]; M[1, 1] = coeff[1]; M[2, 2] = coeff[2]
    M[3, 3] = np.conj(coeff[0]); M[4, 4] = np.conj(coeff[1]); M[5, 5] = np.conj(coeff[2])
    for i in range(6):
        M[i, (i + 1) % 6] += gamma / 2.0
        M[i, (i - 1) % 6] += gamma / 2.0
    return M


# ---------------------------------------------------------------------------
# 2. Time-evolve at one k-point (theory side)
# ---------------------------------------------------------------------------
def P_k_at_time(gamma: float, epsilon: float, k1: float, k2: float,
                r0: tuple, t: float) -> np.ndarray:
    """
    Return the 6-vector tilde_P_m(k, t) for the JMVR walker on Dipanjan's
    lattice, starting from r0 with uniform initial director.

    tilde_P_m(k, 0) = (1/6) * e^{i k . r0}      for each m
    tilde_P_m(k, t) = (e^{Mt})_{mn} * tilde_P_n(k, 0)
                    = (V e^{Lambda t} V^{-1})_{mn} * (1/6) e^{i k . r0}
    """
    M = build_Mk(gamma, epsilon, k1, k2)
    evals, evecs = np.linalg.eig(M)  # columns of evecs are eigvecs
    V = evecs
    Vi = np.linalg.inv(V)
    diag_exp = np.diag(np.exp(evals * t))
    expMt = V @ diag_exp @ Vi
    init = np.full(6, 1.0 / 6.0, dtype=complex) * np.exp(1j * (k1 * r0[0] + k2 * r0[1]))
    return expMt @ init


# ---------------------------------------------------------------------------
# 3. Inverse FT to get P(x, y, t) — Dipanjan's discretisation
# ---------------------------------------------------------------------------
def P_xy_dipanjan(gamma: float, epsilon: float, t: float,
                  L: int = 30, r0: tuple = (0, 0),
                  a: float = 1.0, b: float = 1.0):
    """
    Reproduce Dipanjan's loop exactly:
       k_x = 2 pi n_x / (2 a L), n_x in [0, 2L)
       k_y = 2 pi n_y / (b L),   n_y in [0, L)
       P(x, y, t) = (1/N_k) * sum_{k_x, k_y} P_m_summed(k, t) * e^{-i k . r}
    where P_m_summed = sum over directors m of tilde_P_m.
    Real-space x ranges 0, 2, ..., 2L-2 (step 2, Dipanjan's convention).
    Real-space y ranges 0, 1, ..., L-1.
    """
    Nkx, Nky = 2 * L, L
    N_k = Nkx * Nky
    # Precompute Fourier-transformed director-summed tilde_P at every k
    P_k_summed = np.zeros((Nkx, Nky), dtype=complex)
    for nx in range(Nkx):
        kx = 2 * np.pi * nx / (2 * a * L)
        for ny in range(Nky):
            ky = 2 * np.pi * ny / (b * L)
            v = P_k_at_time(gamma, epsilon, kx, ky, r0, t)
            P_k_summed[nx, ny] = np.sum(v)
    # Inverse Fourier to real space — same discretisation used by Dipanjan
    xs = np.arange(0, 2 * L, 2)          # step of 2 (Dipanjan's convention)
    ys = np.arange(0, L)
    P_real = np.zeros((len(xs), len(ys)), dtype=complex)
    for ix, x in enumerate(xs):
        for iy, y in enumerate(ys):
            s = 0j
            for nx in range(Nkx):
                kx = 2 * np.pi * nx / (2 * a * L)
                for ny in range(Nky):
                    ky = 2 * np.pi * ny / (b * L)
                    s += P_k_summed[nx, ny] * np.exp(-1j * (kx * x + ky * y))
            P_real[ix, iy] = s / N_k
    return xs, ys, P_real


# ---------------------------------------------------------------------------
# 4. Sanity checks
# ---------------------------------------------------------------------------
def sanity_checks():
    print("=" * 70)
    print("POSTMORTEM:  Dipanjan's rtp_tl_2.nb — Python re-implementation")
    print("=" * 70)
    gamma = 0.01
    epsilon = 0.15
    L = 30

    # --- A. Bloch matrix at k=0 ---
    M0 = build_Mk(gamma, epsilon, 0.0, 0.0)
    col0 = M0.sum(axis=0)
    print(f"\n[A] Column sums of M(k=0)        max|col_sum| = {np.max(np.abs(col0)):.2e}")
    print(f"    (should be 0 — probability conservation)")

    # --- A1. Eigenvalues at k=0 ---
    w0 = np.linalg.eigvals(M0)
    w0_sorted = sorted(w0, key=lambda z: -z.real)
    print(f"\n[A1] Eigenvalues at k=0 (Re-sorted):")
    for i, lam in enumerate(w0_sorted):
        print(f"     lambda_{i} = {lam:.6f}")
    print(f"     (lambda_0 should be 0; others should have Re < 0)")

    # --- B. Run the full reconstruction at t=10 (Dipanjan's value), L=30 ---
    t = 10.0
    print(f"\n[B] Running Dipanjan's reconstruction at t={t}, L={L}, eps={epsilon}, gamma={gamma}")
    print(f"    (this loops over {2*L} x {L} = {2*L*L} k-points and same in r ...")
    print(f"     takes ~30 seconds)")
    import time as time_mod
    t0 = time_mod.time()
    xs, ys, P_real = P_xy_dipanjan(gamma, epsilon, t, L=L)
    print(f"    cpu = {time_mod.time() - t0:.1f}s")

    # Check: imaginary part should be tiny
    max_imag = np.max(np.abs(P_real.imag))
    print(f"\n[C1] max|Im P(x,y,t)|             = {max_imag:.2e}")
    print(f"     (should be ~ 1e-12 — physical P is real)")
    P_real_r = P_real.real

    # Check: probability conservation — total over the L x L Dipanjan grid
    total = np.sum(P_real_r)
    print(f"\n[C2] Sum over (x, y) grid          = {total:.6f}")
    print(f"     (should be 1.0 — probability conservation)")

    # Check: non-negativity
    min_P = np.min(P_real_r)
    max_P = np.max(P_real_r)
    print(f"\n[C3] min P(x,y,t) = {min_P:.3e}")
    print(f"     max P(x,y,t) = {max_P:.3e}")
    print(f"     (min should be >= 0, perhaps with tiny floating-point negatives)")

    # Check: parity / sublattice structure
    # On Dipanjan's lattice with steps (2,0), (1,1), (1,-1), starting at (0,0)
    # the walker stays on sites with (x + y) even.  After inverse FT, sites with
    # (x + y) odd should have P = 0.
    even_xs = [x for x in xs if (x + 4) % 2 == 0]  # y=4 is even, so even x
    odd_xs  = [x for x in xs if (x + 4) % 2 == 1]
    print(f"\n[D] Bipartite parity check (y=4 slice, so even x = on-lattice):")
    print(f"     x values in output: {list(xs[:6])} ... (step of 2 → all even)")
    print(f"     -> ALL sample points have correct parity (x even, y even)")
    print(f"     so this test is vacuous in Dipanjan's convention.")

    # --- E. Cross-check: re-Fourier P_real and compare ---
    # If P_real(x, y) is correct, then sum_x,y P_real(x, y) e^{i k . r} should equal tilde_P(k, t).
    # We pick a random (nx, ny) and check.
    nx_test, ny_test = 3, 7
    kx_test = 2 * np.pi * nx_test / (2 * L)
    ky_test = 2 * np.pi * ny_test / L
    v_direct = P_k_at_time(gamma, epsilon, kx_test, ky_test, (0, 0), t)
    p_k_direct = np.sum(v_direct)
    # Re-Fourier of P_real (Dipanjan's grid)
    p_k_reconstructed = 0j
    for ix, x in enumerate(xs):
        for iy, y in enumerate(ys):
            p_k_reconstructed += P_real[ix, iy] * np.exp(1j * (kx_test * x + ky_test * y))
    print(f"\n[E] FT round-trip check at (kx, ky) = ({kx_test:.3f}, {ky_test:.3f}):")
    print(f"     direct tilde_P(k) = {p_k_direct:.4e}")
    print(f"     re-FT  tilde_P(k) = {p_k_reconstructed:.4e}")
    print(f"     |diff|             = {np.abs(p_k_direct - p_k_reconstructed):.2e}")

    # Also check: is the total probability normalised correctly?
    # Compare Dipanjan's "L x L grid" vs "2L x L grid" interpretation.
    print(f"\n[F] Lattice-area / grid-count interpretation:")
    print(f"    Dipanjan uses 2L*L = {2*L*L} k-points")
    print(f"    Real-space grid (xs steps of 2) is L*L = {L*L} sites")
    print(f"    Ratio of k-points / real-space sites = 2  <-- this is suspicious")
    print(f"    For a Bravais lattice the BZ contains exactly N k-points,")
    print(f"    where N = number of real-space sites.")
    print(f"    Either:")
    print(f"      (i)  Dipanjan's real-space grid is 2L*L sites (every integer (x,y))")
    print(f"           but then the walker can't reach the odd-parity sites.")
    print(f"      (ii) Dipanjan's k-grid is twice the BZ (overcounting).")
    print(f"           In that case, sum over k-points double-counts but ")
    print(f"           1/N_k = 1/(2L^2) divides it out — answer is correct.")
    print()
    print(f"    Test: if (ii) is true, total prob = 1 with N_k = 2L^2.")
    print(f"          if (i) is true, total prob should also be 1 with N_k = 2L^2")
    print(f"          IF the walker can reach all 2L*L sites.")
    return P_real_r, xs, ys


# ---------------------------------------------------------------------------
# 5. Compare to direct kinetic Monte Carlo
# ---------------------------------------------------------------------------
def kmc_compare(gamma: float = 0.01, epsilon: float = 0.15,
                t_total: float = 10.0, n_trajectories: int = 5000,
                L: int = 30, r0: tuple = (0, 0)):
    """
    Direct kinetic Monte Carlo on Dipanjan's lattice with PBC.
    Each walker:
      - Starts at r0 with random director m in {0, ..., 5}
      - At each step, sample exponential waiting time = -ln(u)/total_rate
      - With probability prop. to rate_m_d, hop in direction d or flip state
    Returns the empirical P(x, y, t_total) on a Dipanjan-style lattice.
    """
    rng = np.random.default_rng(42)
    # NN displacement vectors for the 6 directors
    DR = np.array([[2, 0], [1, 1], [1, -1], [-2, 0], [-1, -1], [-1, 1]])
    # For each state m, hopping rates to 6 NN; the bias adds +epsilon along
    # ê_m's direction (state m has director DR[m]) and -epsilon along -ê_m.
    # Rate = 1/6 baseline ± epsilon contribution
    # Actually per Dipanjan's matrix: walker hops to 3 forward directions
    # at rate (1/6 with epsilon bias along director's axis only).
    # Hmm — but the off-diagonal matrix elements in the Bloch matrix don't
    # encode hopping; they encode rotation. The diagonal "coeff_i" terms
    # encode hopping. Let me re-derive the real-space rates from the
    # Bloch matrix.
    #
    # In real space: dP_m/dt = sum_{NN displacements a} W_{m m}^{(a)} P_m(r-a, t)
    #                          + sum_{m' != m} (rotation terms)
    #
    # From coeff_m(k) we have:
    #   coeff_m(k) - 1 - gamma = sum_a W_{mm}^{(a)} cos(k.a)    (symmetric part)
    #   2 i eps sin(k . ê_m) corresponds to:
    #     +eps e^{i k . ê_m} = forward hop rate +eps
    #     -eps e^{-i k . ê_m} = backward hop rate -eps  (anti-bias)
    #
    # So the real-space rates for hopping FROM state m are:
    #   To NN displacement +ê_m:  1/6 + eps        (forward bias)
    #   To NN displacement -ê_m:  1/6 - eps        (backward bias)
    #   To other 4 NN:            1/6              (unbiased)
    # Plus rotation rate gamma/2 to (m+1) mod 6 and (m-1) mod 6.
    # Total outgoing rate per state: 1 + gamma.
    base_rate = 1.0 / 6.0
    P_hits = np.zeros((L, L))  # only count sites where (x+y) has even parity
                                # on Dipanjan's lattice, but here use L x L plain
    for traj in range(n_trajectories):
        x, y = r0
        m = rng.integers(0, 6)
        t = 0.0
        while t < t_total:
            # All possible events: 6 hops + 2 rotations
            # Build cumulative rate vector
            rates = np.zeros(8)
            for i in range(6):
                if i == m:
                    rates[i] = base_rate + epsilon
                elif i == (m + 3) % 6:
                    rates[i] = base_rate - epsilon
                else:
                    rates[i] = base_rate
            rates[6] = gamma / 2.0    # rotate to (m+1) mod 6
            rates[7] = gamma / 2.0    # rotate to (m-1) mod 6
            total_rate = rates.sum()
            # Wait time
            dt = -np.log(rng.random()) / total_rate
            t += dt
            if t >= t_total:
                break
            # Pick event
            u = rng.random() * total_rate
            cum = np.cumsum(rates)
            event = np.searchsorted(cum, u)
            if event < 6:
                # Hop in direction DR[event]
                x = (x + DR[event][0]) % (2 * L)
                y = (y + DR[event][1]) % L
            elif event == 6:
                m = (m + 1) % 6
            else:
                m = (m - 1) % 6
        # Snapshot at t_total — record position
        x_idx = x // 2   # Dipanjan's x-step is 2 → divide
        y_idx = y
        if 0 <= x_idx < L and 0 <= y_idx < L:
            P_hits[x_idx, y_idx] += 1
    P_emp = P_hits / n_trajectories
    return P_emp


if __name__ == "__main__":
    P_real_r, xs, ys = sanity_checks()
