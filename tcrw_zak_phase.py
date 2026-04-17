"""
TCRW Phase 5: Topological invariant — the 2D Zak phase
========================================================

The Zak phase is the Berry phase accumulated by an eigenvector as
momentum winds around the Brillouin zone. For this non-Hermitian
stochastic system, it requires the BIORTHOGONAL framework:

  Left eigenvectors:   P(k)^† |φ_n> = λ_n^* |φ_n>
  Right eigenvectors:  P(k)   |ψ_n> = λ_n   |ψ_n>
  Biorthogonality:     <φ_m|ψ_n> = δ_{mn}

Berry connection (non-Hermitian):
  A_μ^n(k) = i <φ_n(k)| ∂_{k_μ} |ψ_n(k)>

Zak phase in direction μ at fixed perpendicular momentum:
  Φ_μ^n(k_perp) = ∮ dk_μ A_μ^n(k)

Discrete (gauge-invariant) Wilson loop:
  W_n = ∏_{j=0}^{N-1} <φ_n(k_j)|ψ_n(k_{j+1})>
  Φ_n = -arg(W_n)

Key result:
  Φ_x = Φ_y = π   for D_r < 0.5 and ω ≠ 0.5  (topological)
  Φ_x = Φ_y = 0   for D_r > 0.5                (trivial)
  Transition at ω = 0.5 where the spectral gap closes.

Reproduces:
  Fig 8(a): Φ_x in the (ω, D_r) plane
  Fig 8(b): Φ_y in the (ω, D_r) plane
  Fig 8(c): phase boundary overlaid with spectral gap closing

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
from scipy.linalg import eig as scipy_eig
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tcrw_spectrum import build_Pk, pbc_eigenvalues, classify_eigenvalues


# ============================================================
# Biorthogonal eigenvector machinery
# ============================================================

def biorthogonal_eig(P):
    """
    Compute biorthogonally normalised left and right eigenvectors of P.

    Returns:
      eigenvalues: (4,) array
      left:  (4,4) array — columns are left eigenvectors |φ_n>
      right: (4,4) array — columns are right eigenvectors |ψ_n>

    Convention:
      P @ right[:,n] = eigenvalues[n] * right[:,n]
      P^H @ left[:,n] = eigenvalues[n]^* * left[:,n]
      left[:,m].conj() @ right[:,n] = δ_{mn}   (biorthogonal)

    Uses scipy.linalg.eig with left=True.
    scipy returns vl such that vl[:,i].conj() @ P = w[i] * vl[:,i].conj()
    ⟹ P^H @ vl[:,i] = w[i]^* * vl[:,i]   ✓
    """
    eigenvalues, vl, vr = scipy_eig(P, left=True, right=True)

    # Biorthogonal normalisation: <φ_m|ψ_n> = δ_{mn}
    # For distinct eigenvalues, <φ_m|ψ_n> = 0 for m ≠ n automatically.
    # We just need to normalise the diagonal: <φ_n|ψ_n> = 1.
    for n in range(len(eigenvalues)):
        overlap = vl[:, n].conj() @ vr[:, n]
        if abs(overlap) > 1e-15:
            # Scale left eigenvector so overlap = 1
            vl[:, n] /= overlap.conj()
        else:
            # Degenerate or near-degenerate: leave as-is
            pass

    return eigenvalues, vl, vr


# ============================================================
# Wilson loop — single band
# ============================================================

def wilson_loop_single_band(omega, D_r, direction, k_perp,
                             band_selector, Nk=200):
    """
    Compute the Wilson loop (Berry phase) for a single band along
    one direction in the BZ.

    Parameters:
      direction:     'x' or 'y' — which k to wind
      k_perp:        fixed value of the perpendicular k
      band_selector: function(eigenvalues) → index of the band to track
                     at the FIRST k-point. Subsequent points use overlap tracking.
      Nk:            number of k-points in the loop

    Returns:
      zak_phase: the Zak phase Φ ∈ (-π, π]
      wilson_W:  the full Wilson loop product (complex number, |W| ≈ 1)
    """
    k_loop = np.linspace(-np.pi, np.pi, Nk, endpoint=False)

    # Step 1: Compute biorthogonal eigenvectors at all k-points
    all_evals = []
    all_left = []
    all_right = []
    for k_par in k_loop:
        if direction == 'x':
            P = build_Pk(omega, D_r, k_par, k_perp)
        else:
            P = build_Pk(omega, D_r, k_perp, k_par)

        evals, vl, vr = biorthogonal_eig(P)
        all_evals.append(evals)
        all_left.append(vl)
        all_right.append(vr)

    # Step 2: Select starting band
    band_idx_0 = band_selector(all_evals[0])

    # Step 3: Wilson loop with overlap-based tracking
    W = 1.0 + 0j
    prev_band = band_idx_0

    for j in range(Nk):
        j_next = (j + 1) % Nk

        # Left eigenvector at current k
        phi_j = all_left[j][:, prev_band]

        # Find best-matching band at k_{j+1}
        overlaps = np.array([abs(phi_j.conj() @ all_right[j_next][:, m])
                             for m in range(4)])
        best = np.argmax(overlaps)

        # Accumulate Wilson loop
        W *= phi_j.conj() @ all_right[j_next][:, best]

        prev_band = best

    zak_phase = -np.angle(W)
    return zak_phase, W


# ============================================================
# Band selectors
# ============================================================

def select_top_band(eigenvalues):
    """Select the band with eigenvalue closest to +1 (largest Re)."""
    return int(np.argmax(eigenvalues.real))


def select_bottom_band(eigenvalues):
    """Select the band with eigenvalue closest to -1 (smallest Re)."""
    return int(np.argmin(eigenvalues.real))


def select_imag_pos(eigenvalues):
    """Select the purely imaginary band with positive Im."""
    # The imaginary pair has Re ≈ 0
    imag_mask = np.abs(eigenvalues.real) < 0.5 * np.max(np.abs(eigenvalues.real))
    if imag_mask.sum() >= 2:
        candidates = np.where(imag_mask)[0]
        return int(candidates[np.argmax(eigenvalues[candidates].imag)])
    # Fallback: pick second largest by Re
    idx = np.argsort(-eigenvalues.real)
    return int(idx[1])


def select_imag_neg(eigenvalues):
    """Select the purely imaginary band with negative Im."""
    imag_mask = np.abs(eigenvalues.real) < 0.5 * np.max(np.abs(eigenvalues.real))
    if imag_mask.sum() >= 2:
        candidates = np.where(imag_mask)[0]
        return int(candidates[np.argmin(eigenvalues[candidates].imag)])
    idx = np.argsort(-eigenvalues.real)
    return int(idx[2])


# ============================================================
# Multi-band Wilson loop (non-Abelian)
# ============================================================

def wilson_loop_multiband(omega, D_r, direction, k_perp,
                           band_selectors, Nk=200):
    """
    Non-Abelian Wilson loop for a group of bands.

    The overlap matrix at each step is:
      O_{mn}(k_j, k_{j+1}) = <φ_m(k_j)|ψ_n(k_{j+1})>
    for m,n in the selected band group.

    The Wilson loop matrix is the ordered product of these overlap matrices.
    The Zak phase = -Im(ln det(W)).

    This is gauge-invariant and doesn't require band tracking.

    Parameters:
      band_selectors: list of functions, one per band in the group.
                      Each returns the band index at the first k-point.

    Returns:
      zak_phase: total Zak phase of the group
      W_matrix:  the Wilson loop matrix
    """
    k_loop = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    n_bands = len(band_selectors)

    # Compute eigenvectors at all k
    all_evals = []
    all_left = []
    all_right = []
    for k_par in k_loop:
        if direction == 'x':
            P = build_Pk(omega, D_r, k_par, k_perp)
        else:
            P = build_Pk(omega, D_r, k_perp, k_par)
        evals, vl, vr = biorthogonal_eig(P)
        all_evals.append(evals)
        all_left.append(vl)
        all_right.append(vr)

    # Select starting bands
    start_bands = [sel(all_evals[0]) for sel in band_selectors]

    # Track bands through the loop using overlap
    band_assignments = [start_bands.copy()]
    for j in range(1, Nk):
        prev_bands = band_assignments[-1]
        new_bands = []
        used = set()
        for m_idx, prev_b in enumerate(prev_bands):
            phi_prev = all_left[j - 1][:, prev_b]
            overlaps = np.array([abs(phi_prev.conj() @ all_right[j][:, n])
                                 for n in range(4)])
            # Find best unused match
            order = np.argsort(-overlaps)
            for candidate in order:
                if candidate not in used:
                    new_bands.append(int(candidate))
                    used.add(int(candidate))
                    break
        band_assignments.append(new_bands)

    # Wilson loop matrix product
    W = np.eye(n_bands, dtype=complex)
    for j in range(Nk):
        j_next = (j + 1) % Nk
        bands_j = band_assignments[j]
        bands_next = band_assignments[j_next] if j_next > 0 else start_bands

        # Overlap matrix
        O = np.zeros((n_bands, n_bands), dtype=complex)
        for m in range(n_bands):
            for n in range(n_bands):
                O[m, n] = (all_left[j][:, bands_j[m]].conj() @
                           all_right[j_next][:, bands_next[n]])
        W = W @ O

    zak_phase = -np.angle(np.linalg.det(W))
    return zak_phase, W


# ============================================================
# Convenience: compute Φ_x and Φ_y
# ============================================================

def compute_zak_phases(omega, D_r, Nk=300, k_perp=0.0, gap_threshold=0.005):
    """
    Compute Φ_x and Φ_y for the top band at fixed k_perp.

    If the spectral gap is below gap_threshold, returns (NaN, NaN)
    since the Wilson loop is ill-defined at the gap closing.

    Returns (Phi_x, Phi_y) in radians.
    """
    # Check spectral gap first
    min_a = np.inf
    for kx in np.linspace(-np.pi, np.pi, 30, endpoint=False):
        for ky in np.linspace(-np.pi, np.pi, 30, endpoint=False):
            evals = pbc_eigenvalues(omega, D_r, kx, ky)
            a, b = classify_eigenvalues(evals)
            if a < min_a:
                min_a = a
    if min_a < gap_threshold:
        return np.nan, np.nan

    Phi_x, _ = wilson_loop_single_band(omega, D_r, 'x', k_perp,
                                        select_top_band, Nk)
    Phi_y, _ = wilson_loop_single_band(omega, D_r, 'y', k_perp,
                                        select_top_band, Nk)
    return Phi_x, Phi_y


def compute_zak_phases_multiband(omega, D_r, Nk=300, k_perp=0.0):
    """
    Compute Φ_x and Φ_y using 2-band Wilson loop
    (top band + positive-imaginary band).
    """
    selectors = [select_top_band, select_imag_pos]
    Phi_x, _ = wilson_loop_multiband(omega, D_r, 'x', k_perp, selectors, Nk)
    Phi_y, _ = wilson_loop_multiband(omega, D_r, 'y', k_perp, selectors, Nk)
    return Phi_x, Phi_y


# ============================================================
# Figures
# ============================================================

def fig8_phase_diagram():
    """
    Fig 8: Zak phase diagram in the (ω, D_r) plane.

    (a) Φ_x — should be π for D_r < 0.5 and ω ≠ 0.5, else 0
    (b) Φ_y — same
    (c) spectral gap boundary overlaid

    Uses single-band Wilson loop for the top band.
    """
    N_omega = 61
    N_Dr = 61
    omega_arr = np.linspace(0.01, 0.99, N_omega)
    D_r_arr = np.linspace(0.01, 0.99, N_Dr)

    Phi_x = np.zeros((N_Dr, N_omega))
    Phi_y = np.zeros((N_Dr, N_omega))

    Nk = 200  # k-points for Wilson loop

    total = N_Dr * N_omega
    count = 0
    for i, D_r in enumerate(D_r_arr):
        for j, omega in enumerate(omega_arr):
            px, py = compute_zak_phases(omega, D_r, Nk=Nk, k_perp=0.0)
            Phi_x[i, j] = px
            Phi_y[i, j] = py
            count += 1
            if count % 200 == 0:
                print(f"  {count}/{total} ({100*count/total:.0f}%)")

    # Wrap to [0, 2π) then identify 0 and π
    # Zak phase should be quantized to 0 or ±π
    # Take abs to fold -π → π
    Phi_x_wrapped = np.abs(Phi_x)
    Phi_y_wrapped = np.abs(Phi_y)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # Custom colormap: blue for 0, red for π
    cmap = plt.cm.RdBu_r

    # (a) Phi_x
    ax = axes[0]
    im = ax.pcolormesh(omega_arr, D_r_arr, Phi_x_wrapped,
                        cmap=cmap, vmin=0, vmax=np.pi, shading='auto')
    ax.set_xlabel('$\\omega$', fontsize=13)
    ax.set_ylabel('$D_r$', fontsize=13)
    ax.set_title('$|\\Phi_x|$', fontsize=14)
    plt.colorbar(im, ax=ax, label='$|\\Phi_x|$ (rad)',
                 ticks=[0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
    ax.axhline(y=0.5, color='white', ls='--', lw=1.5, alpha=0.7)
    ax.axvline(x=0.5, color='white', ls='--', lw=1.5, alpha=0.7)

    # (b) Phi_y
    ax = axes[1]
    im2 = ax.pcolormesh(omega_arr, D_r_arr, Phi_y_wrapped,
                         cmap=cmap, vmin=0, vmax=np.pi, shading='auto')
    ax.set_xlabel('$\\omega$', fontsize=13)
    ax.set_ylabel('$D_r$', fontsize=13)
    ax.set_title('$|\\Phi_y|$', fontsize=14)
    plt.colorbar(im2, ax=ax, label='$|\\Phi_y|$ (rad)',
                 ticks=[0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
    ax.axhline(y=0.5, color='white', ls='--', lw=1.5, alpha=0.7)
    ax.axvline(x=0.5, color='white', ls='--', lw=1.5, alpha=0.7)

    # (c) Combined with gap boundary
    ax = axes[2]
    # Compute spectral gap on same grid
    from tcrw_spectrum import classify_eigenvalues
    gap = np.zeros((N_Dr, N_omega))
    for i, D_r in enumerate(D_r_arr):
        for j, omega in enumerate(omega_arr):
            min_a = np.inf
            for kx in np.linspace(-np.pi, np.pi, 30, endpoint=False):
                for ky in np.linspace(-np.pi, np.pi, 30, endpoint=False):
                    evals = pbc_eigenvalues(omega, D_r, kx, ky)
                    a, b = classify_eigenvalues(evals)
                    if a < min_a:
                        min_a = a
            gap[i, j] = min_a

    # Plot Phi_x + gap contour
    im3 = ax.pcolormesh(omega_arr, D_r_arr, Phi_x_wrapped,
                         cmap=cmap, vmin=0, vmax=np.pi, shading='auto')
    ax.contour(omega_arr, D_r_arr, gap, levels=[0.01],
               colors='lime', linewidths=2)
    ax.set_xlabel('$\\omega$', fontsize=13)
    ax.set_ylabel('$D_r$', fontsize=13)
    ax.set_title('$|\\Phi_x|$ + gap boundary (green)', fontsize=13)
    plt.colorbar(im3, ax=ax, label='$|\\Phi_x|$ (rad)',
                 ticks=[0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])

    plt.suptitle('Fig 8: Zak phase diagram (topological invariant)',
                 fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig('tcrw_fig8_zak_phase.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig8_zak_phase.png")

    # Also save raw data for debugging
    np.savez('tcrw_zak_phase_data.npz',
             omega_arr=omega_arr, D_r_arr=D_r_arr,
             Phi_x=Phi_x, Phi_y=Phi_y, gap=gap)
    print("Saved: tcrw_zak_phase_data.npz")


def fig8_linecuts():
    """
    Diagnostic 1D linecuts through the phase diagram.
    (a) Φ_x vs ω at fixed D_r = 0.1, 0.3
    (b) Φ_x vs D_r at fixed ω = 0.3, 0.8

    Uses compute_zak_phases (with gap check) to avoid divergence
    at the phase boundary.
    """
    Nk = 400

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (a) Φ_x vs ω
    ax = axes[0]
    D_r_vals = [0.1, 0.3]
    omega_scan = np.linspace(0.02, 0.98, 60)
    for D_r in D_r_vals:
        phi_x = []
        for omega in omega_scan:
            px, _ = compute_zak_phases(omega, D_r, Nk=Nk, k_perp=0.0,
                                        gap_threshold=0.01)
            phi_x.append(abs(px) if not np.isnan(px) else np.nan)
        ax.plot(omega_scan, phi_x, 'o-', ms=4, lw=1.5,
                label=f'$D_r={D_r}$')

    ax.axhline(y=np.pi, color='gray', ls='--', lw=0.8, alpha=0.5)
    ax.axhline(y=0, color='gray', ls='--', lw=0.8, alpha=0.5)
    ax.axvline(x=0.5, color='red', ls=':', lw=1, alpha=0.5,
               label='$\\omega=0.5$ (critical)')
    ax.set_xlabel('$\\omega$', fontsize=13)
    ax.set_ylabel('$|\\Phi_x|$', fontsize=13)
    ax.set_title('$|\\Phi_x|$ vs $\\omega$', fontsize=13)
    ax.set_yticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
    ax.set_yticklabels(['0', '$\\pi/4$', '$\\pi/2$', '$3\\pi/4$', '$\\pi$'])
    ax.legend(fontsize=10)
    ax.set_ylim(-0.2, np.pi + 0.3)

    # (b) Φ_x vs D_r
    ax = axes[1]
    omega_vals = [0.3, 0.8]
    D_r_scan = np.linspace(0.02, 0.98, 60)
    for omega in omega_vals:
        phi_x = []
        for D_r in D_r_scan:
            px, _ = compute_zak_phases(omega, D_r, Nk=Nk, k_perp=0.0,
                                        gap_threshold=0.01)
            phi_x.append(abs(px) if not np.isnan(px) else np.nan)
        ax.plot(D_r_scan, phi_x, 'o-', ms=4, lw=1.5,
                label=f'$\\omega={omega}$')

    ax.axhline(y=np.pi, color='gray', ls='--', lw=0.8, alpha=0.5)
    ax.axhline(y=0, color='gray', ls='--', lw=0.8, alpha=0.5)
    ax.axvline(x=0.5, color='red', ls=':', lw=1, alpha=0.5,
               label='$D_r=0.5$ (boundary)')
    ax.set_xlabel('$D_r$', fontsize=13)
    ax.set_ylabel('$|\\Phi_x|$', fontsize=13)
    ax.set_title('$|\\Phi_x|$ vs $D_r$', fontsize=13)
    ax.set_yticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
    ax.set_yticklabels(['0', '$\\pi/4$', '$\\pi/2$', '$3\\pi/4$', '$\\pi$'])
    ax.legend(fontsize=10)
    ax.set_ylim(-0.2, np.pi + 0.3)

    plt.suptitle('Zak phase linecuts (NaN at gap-closing boundary)',
                 fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig('tcrw_fig8_linecuts.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: tcrw_fig8_linecuts.png")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("=" * 60)
    print("TCRW Phase 5: Zak Phase (Topological Invariant)")
    print("=" * 60)

    print("\n--- Linecuts (diagnostic) ---")
    fig8_linecuts()

    print("\n--- Fig 8: Full phase diagram ---")
    fig8_phase_diagram()

    print("\n" + "=" * 60)
    print("Phase 5 complete.")
    print("=" * 60)
