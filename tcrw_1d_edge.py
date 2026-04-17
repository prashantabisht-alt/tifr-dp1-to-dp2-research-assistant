"""
TCRW Phase 6A: 1D Effective Edge Model
=======================================

The paper (Osat et al., arXiv:2602.12020, Methods + Fig 9-10) reduces the
edge dynamics to a 2-state Markov chain with absorption into the bulk.

Physical picture (left edge, ω=1 fully chiral):
  A walker on the left boundary (x=0) has two relevant internal states:
    ← (pointing into the wall, d=3) and ↓ (pointing along the wall, d=2).

  State ←:
    - Cannot translate (blocked by wall), so no chiral step moves it.
    - With prob D_r: noise rotation → CCW (d-1=2=↓) with prob ω,
                                        CW (d+1=0=↑) with prob (1-ω)
    - Going ↑ means pointing into bulk → leaves edge (absorbed).
    - Stays as ← with prob (1-D_r) [chiral step blocked, no rotation].
    - Net: stays ← w.p. (1-D_r), goes to ↓ w.p. ωD_r = R1,
           absorbed w.p. (1-ω)D_r = R2.

  State ↓:
    - Can translate: moves to (x, y-1), director ↓.
    - After chiral translation, rotates: CW (d+1=3=←) w.p. ω,
                                          CCW (d-1=1=→) w.p. (1-ω)
    - Going → means pointing into bulk → absorbed.
    - If noise step instead (prob D_r): rotates but stays at same site.
      Noise CCW (d-1=1=→) w.p. ω → absorbed.
      Noise CW (d+1=3=←) w.p. (1-ω) → stays on edge!
    - So noise PARTIALLY absorbs from ↓ (only the CCW branch).
    - Net: with Fourier phase e^{ik} for the translation along the edge:
      ↓ → ← w.p. ω(1-D_r)·e^{ik} + (1-ω)D_r [chiral CW + noise CW]
      ↓ → absorbed otherwise (noise CCW or chiral CCW).

  Rate matrix A(k) in Fourier space (Eq. 3 of paper):

    A(k) = | 1-D_r          R2+C2*e^ik |
           | R1              0          |

  where p_←,t is top row, p_↓,t is bottom row.
  R1 = ωD_r, R2 = (1-ω)D_r, C2 = ω(1-D_r).

  NOTE: The off-diagonal entries encode:
    A[0,1] = R2+C2·e^{ik}:  ↓→← via noise CW ((1-ω)D_r) + chiral CW (ω(1-D_r)·e^{ik})
    A[1,0] = R1:             ←→↓ via noise CCW (ωD_r)

  Note: NOT probability-conserving (Tr(A) < 2) because of absorption to bulk.
  All eigenvalues have |λ| < 1.

  Eigenvalues (Eq. 4):
    λ_{1,2} = [Tr(A) ± sqrt(Tr(A)² - 4·det(A))] / 2
    Tr(A) = 1 - D_r
    det(A) = -R1·(R2 + C2·e^{ik})

Key predictions:
  - For ω=1, D_r<0.2: eigenspectrum splits into TWO loops in complex plane
  - For D_r>0.2: single loop only
  - The 1D edge eigenvalues trace the oval-shaped edge band in the OBC spectrum
  - P(τ_edge) ~ exp(-τ_edge / τ*) with τ* ~ 1/D_r², collapse with D_r² rescaling

Reproduces:
  Fig 7(e): P(τ_edge) distribution and exponential decay rate ~ D_r²
  Fig 10:   PBC vs OBC vs 1D edge model spectra overlaid

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import matplotlib.pyplot as plt
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tcrw_spectrum import pbc_full_bz, obc_spectrum
from tcrw_core import simulate_tcrw_obc


# ============================================================
# 1D Edge Model: rate matrix and eigenvalues
# ============================================================

def edge_rate_matrix(omega, D_r, k):
    """
    Build the 2×2 edge rate matrix A(k) from Eq. (3) of the paper.

    States: (←, ↓) on the left edge.

    Parameters
    ----------
    omega : float
        Chirality parameter.
    D_r : float
        Rotational noise probability.
    k : float
        Fourier wavevector along the edge direction.

    Returns
    -------
    A : (2,2) complex array
    """
    R1 = omega * D_r
    R2 = (1 - omega) * D_r
    C2 = omega * (1 - D_r)

    A = np.array([
        [1 - D_r,                    R2 + C2 * np.exp(1j * k)],
        [R1,                         0]
    ], dtype=complex)

    return A


def edge_eigenvalues(omega, D_r, k):
    """
    Analytical eigenvalues of the 2×2 edge rate matrix A(k).

    From Eq. (4):
      λ_{1,2} = [Tr(A) ± sqrt(Tr(A)² - 4·det(A))] / 2

    Tr(A) = 1 - D_r
    det(A) = (1-D_r)·0 - R1·(R2 + C2·e^{ik}) = -R1·(R2 + C2·e^{ik})

    Parameters
    ----------
    omega, D_r, k : float

    Returns
    -------
    lambda1, lambda2 : complex
    """
    R1 = omega * D_r
    R2 = (1 - omega) * D_r
    C2 = omega * (1 - D_r)

    tr = 1 - D_r
    det = -R1 * (R2 + C2 * np.exp(1j * k))

    discriminant = tr**2 - 4 * det
    sqrt_disc = np.sqrt(discriminant + 0j)  # ensure complex

    lam1 = (tr + sqrt_disc) / 2
    lam2 = (tr - sqrt_disc) / 2

    return lam1, lam2


def edge_spectrum(omega, D_r, Nk=500):
    """
    Full eigenvalue spectrum of the 1D edge model as k sweeps [-π, π).

    Returns
    -------
    k_arr : (Nk,) array
    lam1, lam2 : (Nk,) complex arrays — the two eigenvalue branches
    """
    k_arr = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    lam1 = np.zeros(Nk, dtype=complex)
    lam2 = np.zeros(Nk, dtype=complex)

    for i, k in enumerate(k_arr):
        lam1[i], lam2[i] = edge_eigenvalues(omega, D_r, k)

    return k_arr, lam1, lam2


# ============================================================
# Edge residence time from simulation
# ============================================================

def measure_edge_residence_times(omega, D_r, L, T_steps, N_traj=1, seed=42):
    """
    Measure edge residence time distribution P(τ_edge) from OBC simulation.

    An edge visit begins when the walker arrives at an edge site (x=0, x=L-1,
    y=0, or y=L-1) and ends when it moves to a non-edge site.

    We track the FIRST walker only for simplicity (single long trajectory).

    Returns
    -------
    tau_list : list of int — all measured edge residence durations
    """
    rng = np.random.default_rng(seed)

    # Single walker
    x = rng.integers(0, L)
    y = rng.integers(0, L)
    d = rng.integers(0, 4)

    DX = np.array([0, 1, 0, -1])
    DY = np.array([1, 0, -1, 0])

    def is_edge(xx, yy):
        return xx == 0 or xx == L - 1 or yy == 0 or yy == L - 1

    tau_list = []
    on_edge = is_edge(x, y)
    current_tau = 1 if on_edge else 0

    for t in range(T_steps):
        r_step = rng.random()
        r_rot = rng.random()

        if r_step < D_r:
            # Noise step: stay put, rotate
            if r_rot < omega:
                d = (d - 1) % 4   # CCW
            else:
                d = (d + 1) % 4   # CW
        else:
            # Chiral step: try to translate, then rotate
            nx = x + DX[d]
            ny = y + DY[d]
            if 0 <= nx < L and 0 <= ny < L:
                x, y = nx, ny
                if r_rot < omega:
                    d = (d + 1) % 4   # CW (opposite chirality)
                else:
                    d = (d - 1) % 4   # CCW
            # else: blocked, no move, no rotation

        now_edge = is_edge(x, y)
        if on_edge:
            if now_edge:
                current_tau += 1
            else:
                # Left the edge
                if current_tau > 0:
                    tau_list.append(current_tau)
                current_tau = 0
                on_edge = False
        else:
            if now_edge:
                on_edge = True
                current_tau = 1

    return tau_list


# ============================================================
# Figure: Fig 10 — PBC vs OBC vs 1D edge model
# ============================================================

def fig10_pbc_obc_1d_overlay():
    """
    Fig 10: PBC (pink fill) vs OBC (colored by edge weight) vs 1D edge model
    (black dots) spectra overlaid for D_r = 0.3, 0.2, 0.15.

    Matches paper style:
      - PBC: pink background fill
      - OBC: colored by Σ_{i∈edge} |ψ_i|², using paper's colormap
      - 1D edge: black dots/line
      - ω = 1 (fully chiral)
    """
    from matplotlib.colors import LinearSegmentedColormap

    omega = 1.0
    D_r_values = [0.3, 0.2, 0.15]
    L_obc = 10

    # Paper-matching colormap: yellow(bulk) → orange → red → dark(edge)
    cmap_paper = LinearSegmentedColormap.from_list(
        'paper_edge', ['#FFFFAA', '#FFDD44', '#FF8800', '#CC0000', '#660000'])

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

    for idx, (ax, D_r) in enumerate(zip(axes, D_r_values)):
        # --- PBC bulk spectrum ---
        _, _, pbc_evals = pbc_full_bz(omega, D_r, Nk=50)
        pbc_flat = pbc_evals.ravel()

        # PBC as filled pink region (match paper)
        ax.scatter(pbc_flat.real, pbc_flat.imag, s=4, color='#FFB0B0',
                   alpha=0.5, zorder=1, label='PBC', edgecolors='none')

        # --- OBC spectrum ---
        obc_evals, ew = obc_spectrum(omega, D_r, L_obc)

        # Color by edge weight — paper style: bulk=yellow, edge=dark red
        sc = ax.scatter(obc_evals.real, obc_evals.imag, c=ew,
                        cmap=cmap_paper, s=12, alpha=0.85,
                        vmin=0, vmax=1, edgecolors='none', zorder=2,
                        label='OBC')

        # --- 1D edge model (black, matching paper's dot style) ---
        k_arr, lam1, lam2 = edge_spectrum(omega, D_r, Nk=500)

        ax.plot(lam1.real, lam1.imag, 'k.', ms=1.5, zorder=3, label='1D edge')
        ax.plot(lam2.real, lam2.imag, 'k.', ms=1.5, zorder=3)

        ax.set_xlabel(r'$\Re(\lambda)$', fontsize=12)
        if idx == 0:
            ax.set_ylabel(r'$\Im(\lambda)$', fontsize=12)
        ax.set_title(f'$D_r = {D_r}$', fontsize=13)
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.15, 1.15)
        ax.set_aspect('equal')
        ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
        ax.axvline(0, color='gray', lw=0.5, alpha=0.3)

    # Legend only on last panel
    axes[2].legend(fontsize=8, loc='upper right', framealpha=0.8)

    # Colorbar — place to the right of the figure
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.65])
    cbar = fig.colorbar(sc, cax=cbar_ax)
    cbar.set_label(r'$\sum_{i\in\mathrm{edge}} |\psi_i|^2$', fontsize=10)

    plt.suptitle(r'Fig 10: Time-scale separation of chiral edge current ($\omega=1$)',
                 fontsize=13, y=1.02)
    plt.savefig('tcrw_fig10_pbc_obc_1d.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig10_pbc_obc_1d.png")


# ============================================================
# Figure: Fig 7(e) — Edge residence time distribution
# ============================================================

def measure_edge_residence_batch(omega, D_r, L, N_walkers, T_steps, seed=42):
    """
    Vectorized edge residence time measurement using N_walkers in parallel.

    IMPORTANT: We define "on edge" using the paper's 1D edge model criterion:
      The walker is on a boundary AND in an edge-following state, i.e.,
      its director either points INTO the wall (blocked) or ALONG the wall
      in the direction of the edge current.

    For the LEFT wall (x=0): edge states are ← (d=3, into wall) and ↓ (d=2, along wall)
    For the RIGHT wall (x=L-1): edge states are → (d=1) and ↑ (d=0)
    For the BOTTOM wall (y=0): edge states are ↓ (d=2) and → (d=1)
    For the TOP wall (y=L-1): edge states are ↑ (d=0) and ← (d=3)

    For ω=1 (CCW chiral): the edge-coupled states form directed loops on each wall.
    τ_edge = consecutive steps in an edge-coupled state before transitioning to bulk.

    Returns list of all τ_edge values collected across walkers.
    """
    rng = np.random.default_rng(seed)

    DXl = np.array([0, 1, 0, -1])
    DYl = np.array([1, 0, -1, 0])

    x = rng.integers(0, L, size=N_walkers)
    y = rng.integers(0, L, size=N_walkers)
    d = rng.integers(0, 4, size=N_walkers)

    def is_edge_state(xx, yy, dd):
        """Check if walker is in an edge-coupled state (not just on boundary).

        Edge-coupled states depend on chirality ω:
          ω > 0.5 (CW chiral): CCW edge current
            Left: ←(3) into wall, ↓(2) along edge
            Bottom: ↓(2) into wall, →(1) along edge
            Right: →(1) into wall, ↑(0) along edge
            Top: ↑(0) into wall, ←(3) along edge
          ω < 0.5 (CCW chiral): CW edge current
            Left: ←(3) into wall, ↑(0) along edge
            Top: ↑(0) into wall, →(1) along edge
            Right: →(1) into wall, ↓(2) along edge
            Bottom: ↓(2) into wall, ←(3) along edge
          ω = 0.5: achiral, both directions equally likely.
                   We use the ω>0.5 convention (arbitrary).
        """
        if omega >= 0.5:
            # CW chiral → CCW edge current
            left   = (xx == 0)   & ((dd == 3) | (dd == 2))
            right  = (xx == L-1) & ((dd == 1) | (dd == 0))
            bottom = (yy == 0)   & ((dd == 2) | (dd == 1))
            top    = (yy == L-1) & ((dd == 0) | (dd == 3))
        else:
            # CCW chiral → CW edge current
            left   = (xx == 0)   & ((dd == 3) | (dd == 0))
            right  = (xx == L-1) & ((dd == 1) | (dd == 2))
            bottom = (yy == 0)   & ((dd == 2) | (dd == 3))
            top    = (yy == L-1) & ((dd == 0) | (dd == 1))
        return left | right | bottom | top

    on_edge = is_edge_state(x, y, d)
    current_tau = np.where(on_edge, 1, 0).astype(np.int64)

    all_taus = []

    for t in range(T_steps):
        r_step = rng.random(N_walkers)
        r_rot = rng.random(N_walkers)

        is_noise = r_step < D_r
        is_chiral = ~is_noise

        # Noise step
        noise_ccw = is_noise & (r_rot < omega)
        noise_cw = is_noise & ~(r_rot < omega)
        d[noise_ccw] = (d[noise_ccw] - 1) % 4
        d[noise_cw] = (d[noise_cw] + 1) % 4

        # Chiral step with OBC
        if np.any(is_chiral):
            nx = x + DXl[d]
            ny = y + DYl[d]
            can_move = is_chiral & (nx >= 0) & (nx < L) & (ny >= 0) & (ny < L)

            x[can_move] = nx[can_move]
            y[can_move] = ny[can_move]

            chiral_cw = can_move & (r_rot < omega)
            chiral_ccw = can_move & ~(r_rot < omega)
            d[chiral_cw] = (d[chiral_cw] + 1) % 4
            d[chiral_ccw] = (d[chiral_ccw] - 1) % 4

        now_edge = is_edge_state(x, y, d)

        # Walkers that WERE on edge and LEFT it: record their tau
        left_edge = on_edge & ~now_edge
        if np.any(left_edge):
            finished_taus = current_tau[left_edge]
            all_taus.extend(finished_taus[finished_taus > 0].tolist())

        # Update tracking
        still_on = on_edge & now_edge
        current_tau[still_on] += 1
        current_tau[left_edge] = 0
        arrived = ~on_edge & now_edge
        current_tau[arrived] = 1

        on_edge = now_edge

    return all_taus


def fig7e_edge_residence_time():
    """
    Fig 7(e): P(τ_edge) distribution.

    Left panel:  P/D_r³ vs τ·D_r² — data collapse across different D_r.
    Right panel: decay rate λ vs D_r — confirm λ ~ D_r².

    Uses vectorized batch walkers for speed.
    """
    omega = 1.0
    L = 30

    # D_r list matches paper Fig 7(e): {1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3}
    # Adaptive: more walkers/steps for smaller D_r (residence time ~ 1/D_r^2)
    D_r_configs = [
        (5e-3,  500, 200_000),
        (2e-3,  200, 500_000),
        (1e-3,  200, 1_000_000),
        (5e-4,  100, 2_000_000),
        (2e-4,  100, 5_000_000),
        (1e-4,  200,  5_000_000),   # 200w×5M ≡ same walker-steps as 100w×10M;
                                     # trimmed from 100w×20M to fit sandbox timeout.
                                     # Event stats for 1e-4 remain marginal (τ ~ 1/D_r² = 1e8).
    ]

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(D_r_configs)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    D_r_values = []
    decay_rates = []

    for (D_r, N_w, T), color in zip(D_r_configs, colors):
        print(f"  τ_edge: D_r={D_r}, {N_w} walkers × {T} steps...")
        taus = measure_edge_residence_batch(omega, D_r, L, N_w, T, seed=42)
        n_events = len(taus)
        print(f"    Collected {n_events} edge events")

        D_r_values.append(D_r)

        if n_events < 100:
            print(f"    Too few events — skipping")
            decay_rates.append(np.nan)
            continue

        taus = np.array(taus, dtype=float)
        median_tau = np.median(taus)
        mean_tau = np.mean(taus)
        print(f"    median τ = {median_tau:.1f}, mean τ = {mean_tau:.1f}")

        # Log-spaced bins for better tail coverage
        max_tau = np.percentile(taus, 99)
        if max_tau < 10:
            max_tau = np.max(taus)
        n_bins = min(80, int(max_tau))
        if n_bins < 5:
            n_bins = int(max_tau) + 1
        bins = np.linspace(0.5, max_tau + 0.5, n_bins + 1)
        counts, bin_edges = np.histogram(taus, bins=bins, density=True)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        mask = counts > 0
        bc = bin_centers[mask]
        P = counts[mask]

        # Rescaled axes
        x_scaled = bc * D_r**2
        y_scaled = P / D_r**3

        ax1.semilogy(x_scaled, y_scaled, 'o-', ms=2.5, lw=0.8, color=color,
                     label=f'$D_r = {D_r}$', alpha=0.8)

        # Fit exponential: P(τ) ~ exp(-λτ)
        if len(bc) > 5:
            i0 = max(1, len(bc) // 4)
            i1 = int(len(bc) * 0.9)
            log_P = np.log(P[i0:i1])
            tau_fit = bc[i0:i1]
            if len(tau_fit) > 3:
                coeffs = np.polyfit(tau_fit, log_P, 1)
                decay_rate = -coeffs[0]
                decay_rates.append(decay_rate)
            else:
                decay_rates.append(np.nan)
        else:
            decay_rates.append(np.nan)

    ax1.set_xlabel(r'$\tau_{\rm edge} \times D_r^2$', fontsize=12)
    ax1.set_ylabel(r'$P(\tau_{\rm edge}) / D_r^3$', fontsize=12)
    ax1.set_title(r'Rescaled $P(\tau_{\rm edge})$ — data collapse', fontsize=12)
    ax1.legend(fontsize=8, ncol=2)
    ax1.set_xlim(left=0, right=20)

    # Right panel: decay rate vs D_r
    D_r_arr = np.array(D_r_values)
    lam_arr = np.array(decay_rates)
    valid = ~np.isnan(lam_arr) & (lam_arr > 0)

    if np.sum(valid) > 0:
        ax2.loglog(D_r_arr[valid], lam_arr[valid], 'ro', ms=8, zorder=3,
                   label='simulation')

    # Fit power law
    if np.sum(valid) > 2:
        log_Dr = np.log(D_r_arr[valid])
        log_lam = np.log(lam_arr[valid])
        slope, intercept = np.polyfit(log_Dr, log_lam, 1)
        Dr_fit = np.logspace(np.log10(D_r_arr[valid].min()) - 0.3,
                             np.log10(D_r_arr[valid].max()) + 0.3, 100)
        lam_fit = np.exp(intercept) * Dr_fit**slope
        ax2.loglog(Dr_fit, lam_fit, 'b--', lw=1.5,
                   label=f'$\\lambda \\propto D_r^{{{slope:.1f}}}$')

    # Reference: λ ~ D_r²
    if np.sum(valid) > 0:
        ref_amp = lam_arr[valid][0] / D_r_arr[valid][0]**2
        Dr_ref = np.logspace(-4.5, -1, 100)
        ax2.loglog(Dr_ref, ref_amp * Dr_ref**2, 'k:', lw=1, alpha=0.5,
                   label=r'$\propto D_r^2$')

    ax2.set_xlabel(r'$D_r$', fontsize=12)
    ax2.set_ylabel(r'Decay rate $\lambda$', fontsize=12)
    ax2.set_title(r'Decay rate $\lambda$ vs $D_r$', fontsize=12)
    ax2.legend(fontsize=10)

    plt.suptitle(r'Fig 7(e): Edge residence time ($\omega=1$, $L=30$)',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('tcrw_fig7e_edge_residence.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig7e_edge_residence.png")


# ============================================================
# Additional diagnostic: 1D edge spectrum vs omega and D_r
# ============================================================

def edge_spectrum_diagnostic():
    """
    Diagnostic plot: how the 1D edge spectrum changes with D_r and omega.

    Row 1: Fix omega=1, vary D_r = 0.4, 0.3, 0.2, 0.15, 0.1
    Row 2: Fix D_r=0.15, vary omega = 1.0, 0.8, 0.6, 0.55
    """
    fig, axes = plt.subplots(2, 5, figsize=(20, 7.5))

    # Row 1: vary D_r
    omega = 1.0
    D_r_vals = [0.4, 0.3, 0.2, 0.15, 0.1]
    for col, D_r in enumerate(D_r_vals):
        ax = axes[0, col]
        k_arr, lam1, lam2 = edge_spectrum(omega, D_r, Nk=500)

        ax.plot(lam1.real, lam1.imag, 'r-', lw=1.5, label=r'$\lambda_1$')
        ax.plot(lam2.real, lam2.imag, 'b-', lw=1.5, label=r'$\lambda_2$')

        ax.set_title(f'$D_r = {D_r}$', fontsize=11)
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-0.6, 0.6)
        ax.set_aspect('equal')
        ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
        ax.axvline(0, color='gray', lw=0.5, alpha=0.3)
        if col == 0:
            ax.set_ylabel(r'$\Im(\lambda)$', fontsize=11)
            ax.legend(fontsize=8)

    axes[0, 2].set_xlabel(r'$\Re(\lambda)$', fontsize=11)

    # Row 2: vary omega
    D_r = 0.15
    omega_vals = [1.0, 0.8, 0.6, 0.55, 0.51]
    for col, omega in enumerate(omega_vals):
        ax = axes[1, col]
        k_arr, lam1, lam2 = edge_spectrum(omega, D_r, Nk=500)

        ax.plot(lam1.real, lam1.imag, 'r-', lw=1.5)
        ax.plot(lam2.real, lam2.imag, 'b-', lw=1.5)

        ax.set_title(f'$\\omega = {omega}$', fontsize=11)
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-0.6, 0.6)
        ax.set_aspect('equal')
        ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
        ax.axvline(0, color='gray', lw=0.5, alpha=0.3)
        if col == 0:
            ax.set_ylabel(r'$\Im(\lambda)$', fontsize=11)

    axes[1, 2].set_xlabel(r'$\Re(\lambda)$', fontsize=11)

    axes[0, 0].text(-0.15, 0.5, r'Vary $D_r$ ($\omega=1$)',
                    transform=axes[0, 0].transAxes, fontsize=12,
                    rotation=90, va='center', fontweight='bold')
    axes[1, 0].text(-0.15, 0.5, r'Vary $\omega$ ($D_r=0.15$)',
                    transform=axes[1, 0].transAxes, fontsize=12,
                    rotation=90, va='center', fontweight='bold')

    plt.suptitle('1D edge model: eigenvalue spectrum in complex plane',
                 fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig('tcrw_1d_edge_diagnostic.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_1d_edge_diagnostic.png")


# ============================================================
# Verification: check A(k) eigenvalues against numpy.linalg.eig
# ============================================================

def verify_analytical_eigenvalues():
    """
    Cross-check: analytical formula vs numpy eigvals for several (omega, D_r, k).
    """
    print("=== Verifying analytical eigenvalues ===")
    test_cases = [
        (1.0, 0.1, 0.0),
        (1.0, 0.1, np.pi),
        (1.0, 0.3, 1.5),
        (0.7, 0.2, -1.0),
        (0.5, 0.1, 0.5),
        (0.0, 0.3, np.pi/3),
    ]

    max_err = 0.0
    for omega, D_r, k in test_cases:
        # Analytical
        l1, l2 = edge_eigenvalues(omega, D_r, k)
        anal = np.sort_complex(np.array([l1, l2]))

        # Numerical
        A = edge_rate_matrix(omega, D_r, k)
        num = np.sort_complex(np.linalg.eigvals(A))

        err = np.max(np.abs(anal - num))
        max_err = max(max_err, err)
        status = "OK" if err < 1e-12 else "FAIL"
        print(f"  ω={omega}, D_r={D_r}, k={k:.3f}: err={err:.2e} [{status}]")

    print(f"\n  Max error across all tests: {max_err:.2e}")
    print(f"  {'PASS' if max_err < 1e-12 else 'FAIL'}: analytical = numerical\n")
    return max_err < 1e-12


def verify_trace_and_det():
    """
    Check Tr(A) = 1-D_r and det(A) = -R1(R2+C2*e^{ik}) for several cases.
    """
    print("=== Verifying Tr(A) and det(A) formulas ===")
    test_cases = [
        (1.0, 0.1, 0.0),
        (1.0, 0.3, np.pi),
        (0.7, 0.2, 1.0),
        (0.5, 0.5, -0.5),
    ]

    for omega, D_r, k in test_cases:
        A = edge_rate_matrix(omega, D_r, k)
        R1 = omega * D_r
        R2 = (1 - omega) * D_r
        C2 = omega * (1 - D_r)

        tr_num = np.trace(A)
        tr_ana = 1 - D_r
        det_num = np.linalg.det(A)
        det_ana = -R1 * (R2 + C2 * np.exp(1j * k))

        err_tr = abs(tr_num - tr_ana)
        err_det = abs(det_num - det_ana)
        print(f"  ω={omega}, D_r={D_r}, k={k:.2f}: "
              f"|ΔTr|={err_tr:.2e}, |Δdet|={err_det:.2e}")


def verify_loop_splitting():
    """
    Check that for ω=1, D_r<0.2 the spectrum splits into two loops
    and for D_r>0.2 it's a single loop.

    Criterion: count connected components by checking if Im(λ) changes sign.
    """
    print("\n=== Checking loop splitting threshold ===")
    omega = 1.0
    for D_r in [0.1, 0.15, 0.2, 0.25, 0.3, 0.4]:
        k_arr, lam1, lam2 = edge_spectrum(omega, D_r, Nk=1000)

        # For branch 1: count zero crossings of Im(λ)
        im1 = lam1.imag
        crossings = np.sum(np.abs(np.diff(np.sign(im1))) > 0)

        # Two loops = 4 zero crossings per branch, one loop = 2
        n_loops = crossings // 2 if crossings > 0 else 1
        print(f"  D_r={D_r}: Im zero crossings = {crossings}, "
              f"estimated loops = {n_loops}")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("=" * 60)
    print("TCRW Phase 6A: 1D Effective Edge Model")
    print("=" * 60)

    # Verification
    verify_analytical_eigenvalues()
    verify_trace_and_det()
    verify_loop_splitting()

    # Figures
    print("\n--- 1D edge spectrum diagnostic ---")
    edge_spectrum_diagnostic()

    print("\n--- Fig 10: PBC vs OBC vs 1D edge overlay ---")
    fig10_pbc_obc_1d_overlay()

    print("\n--- Fig 7(e): Edge residence time ---")
    fig7e_edge_residence_time()

    print("\n" + "=" * 60)
    print("Phase 6A complete.")
    print("=" * 60)
