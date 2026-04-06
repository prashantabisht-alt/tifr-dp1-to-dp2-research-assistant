"""
TCRW Phase 3: Edge currents and their decomposition
=====================================================

Exact computation of probability currents from the transition matrix.

Physics:
  The total current J(x,y) is the net probability flux out of each site.
  Only chiral steps produce spatial displacement; noise steps don't move.

  Decomposition J = J_Dr + J_omega:
    - J_Dr: current from chiral translations where the PREVIOUS step was noise
    - J_omega: current from chiral translations where the PREVIOUS step was chiral

  Key result: the edge current runs OPPOSITE to the walker's chirality
  (quantum Hall skipping orbit analogy). For omega=0.5 (achiral), the
  edge-circulating component of J_Dr vanishes (|J_Dr|/|J_om| ~ D_r on walls).

Method:
  1. Build transition matrix P and decompose: P = P_noise + P_chiral
  2. Find steady state pi: P @ pi = pi
  3. Compute arrival probabilities:
       pi_N(s) = [P_noise @ pi](s)   (arrived at s via noise)
       pi_C(s) = [P_chiral @ pi](s)  (arrived at s via chiral)
     Note: pi_N + pi_C = pi (exact)
  4. Current from state (x,y,d):
       displacement = (DX[d], DY[d]) if move is valid, else (0,0)
       J_Dr(x,y) += pi_N(x,y,d) * (1-D_r) * displacement * I[valid]
       J_omega(x,y) += pi_C(x,y,d) * (1-D_r) * displacement * I[valid]

Reproduces:
  Fig 2(c)-(e): J, J_omega, J_Dr vector fields for omega=1
  Fig 2(h)-(j): same for omega=0 (achiral)
  Fig 3(b): |J_Dr|/|J_omega| ratio vs D_r
  Fig 3(c)-(e): current vectors & angles along left edge vs D_r
  Fig 3(g)-(j): same vs omega

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.colors as mcolors
import sys
import os

# Import from Phase 2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tcrw_obc import (build_transition_matrix, exact_steady_state,
                       state_index, DX, DY)


# ============================================================
# Core: exact current computation
# ============================================================

def build_split_matrices(omega, D_r, L):
    """
    Build noise-only and chiral-only parts of the transition matrix.
    P_noise[j,i]  = prob of noise  step taking state i to state j
    P_chiral[j,i] = prob of chiral step taking state i to state j
    P = P_noise + P_chiral  (column-stochastic)
    """
    import scipy.sparse as sp

    N = 4 * L * L
    # Noise part
    rows_n, cols_n, vals_n = [], [], []
    # Chiral part
    rows_c, cols_c, vals_c = [], [], []

    for x in range(L):
        for y in range(L):
            for d in range(4):
                i = state_index(x, y, d, L)

                # --- Noise step (prob D_r) ---
                d_ccw = (d - 1) % 4
                j_ccw = state_index(x, y, d_ccw, L)
                rows_n.append(j_ccw); cols_n.append(i)
                vals_n.append(omega * D_r)

                d_cw = (d + 1) % 4
                j_cw = state_index(x, y, d_cw, L)
                rows_n.append(j_cw); cols_n.append(i)
                vals_n.append((1 - omega) * D_r)

                # --- Chiral step (prob 1-D_r) ---
                nx = x + DX[d]
                ny = y + DY[d]

                if 0 <= nx < L and 0 <= ny < L:
                    d_cw_c = (d + 1) % 4
                    j1 = state_index(nx, ny, d_cw_c, L)
                    rows_c.append(j1); cols_c.append(i)
                    vals_c.append(omega * (1 - D_r))

                    d_ccw_c = (d - 1) % 4
                    j2 = state_index(nx, ny, d_ccw_c, L)
                    rows_c.append(j2); cols_c.append(i)
                    vals_c.append((1 - omega) * (1 - D_r))
                else:
                    j_stay = state_index(x, y, d, L)
                    rows_c.append(j_stay); cols_c.append(i)
                    vals_c.append(1 - D_r)

    P_noise = sp.coo_matrix((vals_n, (rows_n, cols_n)), shape=(N, N)).tocsc()
    P_chiral = sp.coo_matrix((vals_c, (rows_c, cols_c)), shape=(N, N)).tocsc()
    return P_noise, P_chiral


def exact_currents(omega, D_r, L):
    """
    Compute exact steady-state probability currents J, J_Dr, J_omega.

    Returns dict with:
      'Pxy':   (L,L) spatial probability distribution
      'Jx', 'Jy':       (L,L) total current vector field
      'Jx_Dr', 'Jy_Dr': (L,L) noise-induced current (J_Dr)
      'Jx_om', 'Jy_om': (L,L) chiral current (J_omega)
      'pi':    full state vector
    """
    # 1. Steady state
    Pxy, pi = exact_steady_state(omega, D_r, L)

    # 2. Split transition matrix
    P_noise, P_chiral = build_split_matrices(omega, D_r, L)

    # 3. Arrival probabilities
    #    pi_N(s) = prob of being at s AND having arrived via noise
    #    pi_C(s) = prob of being at s AND having arrived via chiral
    if L <= 30:
        pi_N = (P_noise.toarray() @ pi)
        pi_C = (P_chiral.toarray() @ pi)
    else:
        pi_N = P_noise @ pi
        pi_C = P_chiral @ pi

    # Sanity: pi_N + pi_C should equal pi
    assert np.allclose(pi_N + pi_C, pi, atol=1e-12), \
        f"pi_N + pi_C != pi, max diff = {np.max(np.abs(pi_N + pi_C - pi))}"

    # 4. Compute currents
    #    For each state (x,y,d): if chiral step is valid, displacement is (DX[d], DY[d])
    #    Current contribution = pi_*(x,y,d) * (1-D_r) * displacement
    #    But wait — (1-D_r) is already included in pi_N and pi_C!
    #    No: pi_N(s) = sum_{s'} pi(s') * P_noise[s,s'] = the probability FLOWING
    #    into s via noise. The NEXT step from s is independent.
    #
    #    The current from s = (x,y,d) is the displacement from the CHIRAL part
    #    of the step starting from s. The probability of being at s is pi(s).
    #    The probability of taking a chiral step is (1-D_r).
    #    The displacement is DX[d], DY[d] if valid, else 0.
    #
    #    For the DECOMPOSITION: we weight by whether the walker ARRIVED at s
    #    via noise (pi_N) or via chiral (pi_C).

    Jx = np.zeros((L, L))
    Jy = np.zeros((L, L))
    Jx_Dr = np.zeros((L, L))
    Jy_Dr = np.zeros((L, L))
    Jx_om = np.zeros((L, L))
    Jy_om = np.zeros((L, L))

    for x in range(L):
        for y in range(L):
            for d in range(4):
                nx = x + DX[d]
                ny = y + DY[d]
                if 0 <= nx < L and 0 <= ny < L:
                    s = state_index(x, y, d, L)
                    dx = float(DX[d])
                    dy = float(DY[d])

                    # Prob of chiral step from this state = (1-D_r)
                    # Total current from (x,y,d)
                    flux = pi[s] * (1 - D_r)
                    Jx[x, y] += flux * dx
                    Jy[x, y] += flux * dy

                    # Decomposed: weight by arrival type
                    flux_Dr = pi_N[s] * (1 - D_r)
                    Jx_Dr[x, y] += flux_Dr * dx
                    Jy_Dr[x, y] += flux_Dr * dy

                    flux_om = pi_C[s] * (1 - D_r)
                    Jx_om[x, y] += flux_om * dx
                    Jy_om[x, y] += flux_om * dy

    return {
        'Pxy': Pxy, 'pi': pi,
        'Jx': Jx, 'Jy': Jy,
        'Jx_Dr': Jx_Dr, 'Jy_Dr': Jy_Dr,
        'Jx_om': Jx_om, 'Jy_om': Jy_om,
    }


# ============================================================
# Plotting helpers
# ============================================================

def current_magnitude(Jx, Jy):
    """Magnitude of current vector field."""
    return np.sqrt(Jx**2 + Jy**2)


def plot_vector_field(ax, Jx, Jy, L, title='', cmap='viridis',
                      log_color=True, scale_arrows=True):
    """
    Plot a current vector field as colored arrows on the lattice.
    Arrow direction = current direction; color = log10(magnitude).
    """
    mag = current_magnitude(Jx, Jy)
    X, Y = np.meshgrid(np.arange(L), np.arange(L), indexing='ij')

    # Normalize arrows to unit length for visibility
    mag_safe = np.where(mag > 0, mag, 1.0)
    Ux = Jx / mag_safe
    Uy = Jy / mag_safe

    if log_color and mag.max() > 0:
        mag_plot = mag.copy()
        mag_plot[mag_plot == 0] = mag_plot[mag_plot > 0].min() * 0.1
        colors = np.log10(mag_plot)
    else:
        colors = mag

    q = ax.quiver(X, Y, Ux, Uy, colors, cmap=cmap,
                   scale=L * 1.5, width=0.008, headwidth=4, headlength=5,
                   pivot='mid', clim=[colors[mag > 0].min() if (mag > 0).any() else -1,
                                      colors.max()])
    cbar = plt.colorbar(q, ax=ax, shrink=0.85, pad=0.02)
    if log_color:
        cbar.set_label('$\\log_{10} |\\mathbf{J}|$', fontsize=10)
    else:
        cbar.set_label('$|\\mathbf{J}|$', fontsize=10)

    ax.set_xlim(-0.5, L - 0.5)
    ax.set_ylim(-0.5, L - 0.5)
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=11)
    ax.set_xlabel('X', fontsize=10)
    ax.set_ylabel('Y', fontsize=10)
    return q


def left_edge_currents(result, L):
    """
    Extract current vectors along the left edge (x=0, y=0..L-1).
    Returns arrays of (Jx, Jy) for total, J_Dr, J_omega.
    """
    y_arr = np.arange(L)
    return {
        'y': y_arr,
        'Jx': result['Jx'][0, :],
        'Jy': result['Jy'][0, :],
        'Jx_Dr': result['Jx_Dr'][0, :],
        'Jy_Dr': result['Jy_Dr'][0, :],
        'Jx_om': result['Jx_om'][0, :],
        'Jy_om': result['Jy_om'][0, :],
    }


def current_angle(Jx, Jy):
    """Angle of current vector in radians, measured from +x axis."""
    return np.arctan2(Jy, Jx)


# ============================================================
# Fig 2(c)-(e): vector fields for omega=1 (chiral)
# Fig 2(h)-(j): vector fields for omega=0 (achiral)
# ============================================================

def fig2_vector_fields():
    """
    6-panel figure: top row omega=1, bottom row omega=0.5 (achiral).
    Columns: J_total, J_omega, J_Dr.
    L=10, D_r=10^-3.
    """
    L = 10
    D_r = 1e-3

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    for row, omega in enumerate([1.0, 0.5]):
        print(f"  Computing currents: omega={omega}, D_r={D_r}, L={L}...")
        res = exact_currents(omega, D_r, L)

        titles = [
            f'$\\vec{{J}}$ total',
            f'$\\vec{{J}}_\\omega$ (chiral)',
            f'$\\vec{{J}}_{{D_r}}$ (noise-induced)',
        ]
        Jx_list = [res['Jx'], res['Jx_om'], res['Jx_Dr']]
        Jy_list = [res['Jy'], res['Jy_om'], res['Jy_Dr']]

        for col in range(3):
            ax = axes[row, col]
            label = f"$\\omega={omega}$: {titles[col]}"
            plot_vector_field(ax, Jx_list[col], Jy_list[col], L, title=label)

        # Annotate rows
        row_label = 'CHIRAL ($\\omega=1$)' if omega == 1.0 else 'ACHIRAL ($\\omega=0.5$)'
        axes[row, 0].text(-0.25, 0.5, row_label, transform=axes[row, 0].transAxes,
                          fontsize=13, fontweight='bold', rotation=90,
                          va='center', ha='center')

    plt.suptitle(f'Fig 2: Current vector fields (OBC, L={L}, $D_r={D_r}$)',
                 fontsize=15, y=1.01)
    plt.tight_layout()
    plt.savefig('tcrw_fig2_currents.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig2_currents.png")


# ============================================================
# Fig 3(b): |J_Dr|/|J_omega| vs D_r on the left wall
# ============================================================

def fig3b_ratio_vs_Dr():
    """
    |J_Dr|/|J_omega| averaged over the left wall (x=0) vs D_r.
    Paper: ratio is ~constant for small D_r, diverges for large D_r.
    Uses multiple L values to show size-independence.
    """
    omega = 1.0
    D_r_values = np.logspace(-4, 0, 25)
    L_values = [4, 9, 19, 49]
    markers = ['o', 's', '^', 'D']
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(L_values)))

    fig, ax = plt.subplots(figsize=(7, 5))

    for L, marker, color in zip(L_values, markers, colors):
        ratios = []
        print(f"  Fig 3(b): L={L}...")
        for D_r in D_r_values:
            res = exact_currents(omega, D_r, L)
            # Average |J_Dr| and |J_omega| on left wall (x=0)
            mag_Dr = current_magnitude(res['Jx_Dr'][0, :], res['Jy_Dr'][0, :])
            mag_om = current_magnitude(res['Jx_om'][0, :], res['Jy_om'][0, :])
            mean_Dr = mag_Dr.mean()
            mean_om = mag_om.mean()
            if mean_om > 0:
                ratios.append(mean_Dr / mean_om)
            else:
                ratios.append(np.nan)
        ax.loglog(D_r_values, ratios, marker=marker, color=color, ms=5,
                  lw=1.2, label=f'$L={L}$', markerfacecolor='none',
                  markeredgewidth=1.2)

    ax.set_xlabel('$D_r$', fontsize=13)
    ax.set_ylabel('$|\\mathbf{J}_{D_r}| / |\\mathbf{J}_\\omega|$', fontsize=13)
    ax.set_title('Fig 3(b): Current ratio on left wall ($\\omega = 1$)', fontsize=13)
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig('tcrw_fig3b_Jratio_vs_Dr.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig3b_Jratio_vs_Dr.png")


# ============================================================
# Fig 3(c)-(e): currents along left edge vs D_r
# ============================================================

def fig3cde_left_edge_vs_Dr():
    """
    Left-edge current vectors and angles for omega=1, varying D_r.

    (c) |J_Dr| along left edge for different D_r (color-coded arrows)
    (d) |J_omega| along left edge
    (e) theta_{J_Dr} vs D_r (angle of noise-induced current)
    """
    omega = 1.0
    L = 10
    D_r_values = np.logspace(-4, -0.3, 8)
    cmap = plt.cm.plasma

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # --- (c) J_Dr vectors along left edge ---
    ax = axes[0]
    for idx, D_r in enumerate(D_r_values):
        res = exact_currents(omega, D_r, L)
        le = left_edge_currents(res, L)
        color = cmap(idx / (len(D_r_values) - 1))

        mag = current_magnitude(le['Jx_Dr'], le['Jy_Dr'])
        mag_safe = np.where(mag > 0, mag, 1.0)
        ux = le['Jx_Dr'] / mag_safe
        uy = le['Jy_Dr'] / mag_safe

        # Offset x position for visibility
        x_offset = idx * 0.12
        ax.quiver(np.full(L, x_offset), le['y'], ux, uy,
                  color=color, scale=15, width=0.006, headwidth=4,
                  headlength=5, pivot='mid')

    ax.set_xlim(-0.5, len(D_r_values) * 0.12 + 0.5)
    ax.set_ylim(-0.5, L - 0.5)
    ax.set_ylabel('Left Edge (Y)', fontsize=12)
    ax.set_title('$\\vec{J}_{D_r}$ along left edge', fontsize=12)
    ax.set_xticks([])
    # Add colorbar for D_r
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.LogNorm(
        vmin=D_r_values.min(), vmax=D_r_values.max()))
    cbar = plt.colorbar(sm, ax=ax, shrink=0.85)
    cbar.set_label('$D_r$', fontsize=11)

    # --- (d) J_omega vectors along left edge ---
    ax = axes[1]
    for idx, D_r in enumerate(D_r_values):
        res = exact_currents(omega, D_r, L)
        le = left_edge_currents(res, L)
        color = cmap(idx / (len(D_r_values) - 1))

        mag = current_magnitude(le['Jx_om'], le['Jy_om'])
        mag_safe = np.where(mag > 0, mag, 1.0)
        ux = le['Jx_om'] / mag_safe
        uy = le['Jy_om'] / mag_safe

        x_offset = idx * 0.12
        ax.quiver(np.full(L, x_offset), le['y'], ux, uy,
                  color=color, scale=15, width=0.006, headwidth=4,
                  headlength=5, pivot='mid')

    ax.set_xlim(-0.5, len(D_r_values) * 0.12 + 0.5)
    ax.set_ylim(-0.5, L - 0.5)
    ax.set_ylabel('Left Edge (Y)', fontsize=12)
    ax.set_title('$\\vec{J}_\\omega$ along left edge', fontsize=12)
    ax.set_xticks([])
    cbar2 = plt.colorbar(sm, ax=ax, shrink=0.85)
    cbar2.set_label('$D_r$', fontsize=11)

    # --- (e) theta_{J_Dr} vs D_r ---
    ax = axes[2]
    D_r_scan = np.logspace(-4, 0, 30)
    theta_Dr_avg = []
    theta_om_avg = []

    for D_r in D_r_scan:
        res = exact_currents(omega, D_r, L)
        le = left_edge_currents(res, L)

        # Average angle over interior edge sites (exclude corners)
        inner = slice(1, L - 1)
        mag_Dr = current_magnitude(le['Jx_Dr'][inner], le['Jy_Dr'][inner])
        mag_om = current_magnitude(le['Jx_om'][inner], le['Jy_om'][inner])

        if mag_Dr.sum() > 0:
            # Magnitude-weighted average angle
            angles_Dr = current_angle(le['Jx_Dr'][inner], le['Jy_Dr'][inner])
            theta_Dr_avg.append(np.average(angles_Dr, weights=mag_Dr))
        else:
            theta_Dr_avg.append(np.nan)

        if mag_om.sum() > 0:
            angles_om = current_angle(le['Jx_om'][inner], le['Jy_om'][inner])
            theta_om_avg.append(np.average(angles_om, weights=mag_om))
        else:
            theta_om_avg.append(np.nan)

    ax.semilogx(D_r_scan, np.array(theta_Dr_avg), 'o-', color='tab:red',
                ms=4, lw=1.5, label='$\\theta_{J_{D_r}}$', markerfacecolor='none')
    ax.semilogx(D_r_scan, np.array(theta_om_avg), 's-', color='tab:blue',
                ms=4, lw=1.5, label='$\\theta_{J_\\omega}$', markerfacecolor='none')

    ax.axhline(y=np.pi / 4, color='gray', ls='--', lw=0.8, alpha=0.5)
    ax.axhline(y=-np.pi / 2, color='gray', ls='--', lw=0.8, alpha=0.5)
    ax.axhline(y=0, color='gray', ls='-', lw=0.5, alpha=0.3)

    ax.set_yticks([-np.pi/2, -np.pi/4, 0, np.pi/4, np.pi/2])
    ax.set_yticklabels(['$-\\pi/2$', '$-\\pi/4$', '0', '$\\pi/4$', '$\\pi/2$'])
    ax.set_xlabel('$D_r$', fontsize=13)
    ax.set_ylabel('Current angle $\\theta$', fontsize=13)
    ax.set_title(f'Current angles vs $D_r$ ($\\omega={omega}$, L={L})', fontsize=12)
    ax.legend(fontsize=10)

    plt.suptitle('Fig 3(c)-(e): Left-edge currents vs $D_r$', fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig('tcrw_fig3cde_leftedge_vs_Dr.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig3cde_leftedge_vs_Dr.png")


# ============================================================
# Fig 3(g)-(j): currents along left edge vs omega
# ============================================================

def fig3ghij_left_edge_vs_omega():
    """
    Left-edge current vectors and angles for D_r=10^-3, varying omega.

    (g) |J_Dr|/|J_omega| vs omega
    (h) J_Dr vectors along left edge for different omega
    (i) J_omega vectors along left edge
    (j) theta of both currents vs omega
    """
    D_r = 1e-3
    L = 10
    omega_values = np.linspace(0.0, 1.0, 21)
    omega_arrows = np.linspace(0.0, 1.0, 9)
    cmap = plt.cm.coolwarm

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # --- (h) J_Dr vectors along left edge vs omega ---
    ax = axes[0]
    for idx, omega in enumerate(omega_arrows):
        res = exact_currents(omega, D_r, L)
        le = left_edge_currents(res, L)
        color = cmap(idx / (len(omega_arrows) - 1))

        mag = current_magnitude(le['Jx_Dr'], le['Jy_Dr'])
        mag_safe = np.where(mag > 0, mag, 1.0)
        ux = le['Jx_Dr'] / mag_safe
        uy = le['Jy_Dr'] / mag_safe

        x_offset = idx * 0.12
        ax.quiver(np.full(L, x_offset), le['y'], ux, uy,
                  color=color, scale=15, width=0.006, headwidth=4,
                  headlength=5, pivot='mid')

    ax.set_xlim(-0.5, len(omega_arrows) * 0.12 + 0.5)
    ax.set_ylim(-0.5, L - 0.5)
    ax.set_ylabel('Left Edge (Y)', fontsize=12)
    ax.set_title('$\\vec{J}_{D_r}$ along left edge', fontsize=12)
    ax.set_xticks([])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(
        vmin=0, vmax=1))
    cbar = plt.colorbar(sm, ax=ax, shrink=0.85)
    cbar.set_label('$\\omega$', fontsize=11)

    # --- (i) J_omega vectors along left edge vs omega ---
    ax = axes[1]
    for idx, omega in enumerate(omega_arrows):
        res = exact_currents(omega, D_r, L)
        le = left_edge_currents(res, L)
        color = cmap(idx / (len(omega_arrows) - 1))

        mag = current_magnitude(le['Jx_om'], le['Jy_om'])
        mag_safe = np.where(mag > 0, mag, 1.0)
        ux = le['Jx_om'] / mag_safe
        uy = le['Jy_om'] / mag_safe

        x_offset = idx * 0.12
        ax.quiver(np.full(L, x_offset), le['y'], ux, uy,
                  color=color, scale=15, width=0.006, headwidth=4,
                  headlength=5, pivot='mid')

    ax.set_xlim(-0.5, len(omega_arrows) * 0.12 + 0.5)
    ax.set_ylim(-0.5, L - 0.5)
    ax.set_ylabel('Left Edge (Y)', fontsize=12)
    ax.set_title('$\\vec{J}_\\omega$ along left edge', fontsize=12)
    ax.set_xticks([])
    cbar2 = plt.colorbar(sm, ax=ax, shrink=0.85)
    cbar2.set_label('$\\omega$', fontsize=11)

    # --- (j) theta of both currents vs omega ---
    ax = axes[2]
    theta_Dr_avg = []
    theta_om_avg = []

    for omega in omega_values:
        res = exact_currents(omega, D_r, L)
        le = left_edge_currents(res, L)
        inner = slice(1, L - 1)

        mag_Dr = current_magnitude(le['Jx_Dr'][inner], le['Jy_Dr'][inner])
        mag_om = current_magnitude(le['Jx_om'][inner], le['Jy_om'][inner])

        if mag_Dr.sum() > 1e-30:
            angles_Dr = current_angle(le['Jx_Dr'][inner], le['Jy_Dr'][inner])
            theta_Dr_avg.append(np.average(angles_Dr, weights=mag_Dr))
        else:
            theta_Dr_avg.append(np.nan)

        if mag_om.sum() > 1e-30:
            angles_om = current_angle(le['Jx_om'][inner], le['Jy_om'][inner])
            theta_om_avg.append(np.average(angles_om, weights=mag_om))
        else:
            theta_om_avg.append(np.nan)

    ax.plot(omega_values, theta_Dr_avg, 'o-', color='tab:red',
            ms=5, lw=1.5, label='$\\theta_{J_{D_r}}$', markerfacecolor='none')
    ax.plot(omega_values, theta_om_avg, 's-', color='tab:blue',
            ms=5, lw=1.5, label='$\\theta_{J_\\omega}$', markerfacecolor='none')

    ax.axhline(y=np.pi / 4, color='gray', ls='--', lw=0.8, alpha=0.5)
    ax.axhline(y=-np.pi / 4, color='gray', ls='--', lw=0.8, alpha=0.5)
    ax.axhline(y=0, color='gray', ls='-', lw=0.5, alpha=0.3)
    ax.axvline(x=0.5, color='gray', ls=':', lw=1, alpha=0.5)

    ax.set_yticks([-np.pi/2, -np.pi/4, 0, np.pi/4, np.pi/2])
    ax.set_yticklabels(['$-\\pi/2$', '$-\\pi/4$', '0', '$\\pi/4$', '$\\pi/2$'])
    ax.set_xlabel('$\\omega$', fontsize=13)
    ax.set_ylabel('Current angle $\\theta$', fontsize=13)
    ax.set_title(f'Current angles vs $\\omega$ ($D_r={D_r}$, L={L})', fontsize=12)
    ax.legend(fontsize=10)
    ax.set_xlim(-0.02, 1.02)

    plt.suptitle('Fig 3(g)-(j): Left-edge currents vs $\\omega$', fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig('tcrw_fig3ghij_leftedge_vs_omega.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig3ghij_leftedge_vs_omega.png")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("=" * 60)
    print("TCRW Phase 3: Edge Currents and Decomposition")
    print("=" * 60)

    print("\n--- Fig 2: Current vector fields ---")
    fig2_vector_fields()

    print("\n--- Fig 3(b): |J_Dr|/|J_omega| vs D_r ---")
    fig3b_ratio_vs_Dr()

    print("\n--- Fig 3(c)-(e): Left-edge currents vs D_r ---")
    fig3cde_left_edge_vs_Dr()

    print("\n--- Fig 3(g)-(j): Left-edge currents vs omega ---")
    fig3ghij_left_edge_vs_omega()

    print("\n" + "=" * 60)
    print("Phase 3 complete.")
    print("=" * 60)
