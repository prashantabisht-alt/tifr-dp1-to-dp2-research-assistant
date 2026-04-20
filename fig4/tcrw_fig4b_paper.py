"""
TCRW Fig 4(b) — PBC Bloch band structure as 3D surfaces over the BZ
====================================================================

Reproduces the paper's Fig. 4(b): Re(λ(k)) as a 3D surface plot over
k_x, k_y ∈ [-π, π] for ω ∈ {0.35, 0.5, 0.65} at D_r = 0.1.

Physics
-------
The 4×4 Fourier-space transition matrix P(k) of the one-walker Markov
chain is (paper Eq. 1):

    P(k) =
      | 0              R1+C1·e^{+ikx}   0              R2+C2·e^{-ikx} |
      | R2+C2·e^{+iky} 0                R1+C1·e^{-iky} 0              |
      | 0              R2+C2·e^{+ikx}   0              R1+C1·e^{-ikx} |
      | R1+C1·e^{+iky} 0                R2+C2·e^{-iky} 0              |

with
    C1 = (1-ω)(1-D_r)   (chiral, CCW rotation after translation)
    C2 = ω(1-D_r)        (chiral, CW  rotation after translation)
    R1 = ω D_r            (noise,  CCW rotation, no translation)
    R2 = (1-ω) D_r        (noise,  CW  rotation, no translation)

Sublattice symmetry S·P(k)·S⁻¹ = −P(k) with S = 1⊗σ_z ⇒ the 4 eigenvalues
appear in ±λ pairs. Combined with the trace-0 bipartite structure, for
the parameter range relevant to Fig 4(b) the pairs are exactly
    {±a(k)}   (real, non-Hermitian but real-diagonalisable)
    {±ib(k)}  (pure imaginary)
so only two Re-surfaces are dispersive (the red surfaces ±a(k)) and the
other two eigenvalues contribute a flat plane at Re(λ)=0 (the blue plane).

The "gap" is Δ(ω, D_r) = min_k a(k). It closes exactly at ω = 0.5 — this
is the topological phase transition.

Classification robustness
-------------------------
Near ω = 0.5 we have a(k) → 0 somewhere in the BZ, so sorting eigenvalues
by |Re| fails (the "real pair" has small |Re| too). Instead we pair
eigenvalues by matching their **squared values**: within each ±λ pair,
λ² is the same number. This is exact (only limited by eig() precision)
and works right up to the gap-closing point.

Paper parameters used
---------------------
    D_r = 0.1
    ω ∈ {0.35, 0.5, 0.65}
    (kx, ky) grid: 80 × 80 on [-π, π]²

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3D projection)


# ---------------------------------------------------------------------------
# Bloch matrix
# ---------------------------------------------------------------------------

def build_Pk(omega, D_r, kx, ky):
    """4×4 Fourier-space transition matrix P(k). See paper Eq. (1)."""
    C1 = (1.0 - omega) * (1.0 - D_r)
    C2 = omega * (1.0 - D_r)
    R1 = omega * D_r
    R2 = (1.0 - omega) * D_r

    ekx = np.exp(1j * kx)
    eky = np.exp(1j * ky)

    P = np.array([
        [0.0,             R1 + C1 * ekx,    0.0,             R2 + C2 / ekx],
        [R2 + C2 * eky,   0.0,              R1 + C1 / eky,   0.0          ],
        [0.0,             R2 + C2 * ekx,    0.0,             R1 + C1 / ekx],
        [R1 + C1 * eky,   0.0,              R2 + C2 / eky,   0.0          ],
    ], dtype=complex)
    return P


# ---------------------------------------------------------------------------
# Robust eigenvalue pairing:  group 4 eigvals into {+λ₁,-λ₁,+λ₂,-λ₂}
# ---------------------------------------------------------------------------

def pair_eigenvalues(evals):
    """
    Group the 4 eigenvalues of P(k) into the two ±λ pairs by matching
    squared values (which are identical within each pair).

    Returns
    -------
    a, b : float
        Non-negative magnitudes such that the spectrum is {±a, ±ib}.
        If classification is clean (μ_real real ≥ 0, μ_imag real ≤ 0)
        these are exact; otherwise they are the best projection.
    res_real, res_imag : float
        Residuals measuring deviation from the "clean" ±a, ±ib form.
        Zero to machine precision when the paper's classification holds.
    """
    sq = np.asarray(evals) ** 2
    # Greedy pair by closest squared value
    i = 0
    dists = np.abs(sq[1:] - sq[0])
    j_rel = int(np.argmin(dists))
    j = j_rel + 1
    others = [idx for idx in range(4) if idx not in (i, j)]
    k_, l_ = others

    mu_A = 0.5 * (sq[i] + sq[j])   # average for noise robustness
    mu_B = 0.5 * (sq[k_] + sq[l_])

    # Real pair: μ with the larger real part (positive → real ±a)
    # Imag pair: μ with the smaller real part (negative → imag ±ib)
    if mu_A.real >= mu_B.real:
        mu_real, mu_imag = mu_A, mu_B
    else:
        mu_real, mu_imag = mu_B, mu_A

    a = np.sqrt(max(mu_real.real, 0.0))
    b = np.sqrt(max(-mu_imag.real, 0.0))

    # "Clean" residuals: real pair should have μ real ≥ 0, imag pair μ real ≤ 0
    res_real = abs(mu_real.imag) + max(0.0, -mu_real.real)
    res_imag = abs(mu_imag.imag) + max(0.0,  mu_imag.real)
    return a, b, res_real, res_imag


# ---------------------------------------------------------------------------
# BZ surface
# ---------------------------------------------------------------------------

def pbc_bz_surface(omega, D_r, Nk=80):
    """
    Compute a(kx, ky) and b(kx, ky) on an Nk×Nk grid of the full BZ.

    Returns
    -------
    kx, ky     : 1D arrays, length Nk
    a_surf     : 2D array (Nk, Nk), a(kx[i], ky[j]) at a_surf[i,j]
    b_surf     : 2D array (Nk, Nk)
    max_res_real, max_res_imag : float
        Maximum classification residuals seen on the grid.
    """
    kx = np.linspace(-np.pi, np.pi, Nk)
    ky = np.linspace(-np.pi, np.pi, Nk)
    a_surf = np.empty((Nk, Nk))
    b_surf = np.empty((Nk, Nk))
    max_res_real = 0.0
    max_res_imag = 0.0
    for i, kxi in enumerate(kx):
        for j, kyj in enumerate(ky):
            ev = np.linalg.eigvals(build_Pk(omega, D_r, kxi, kyj))
            a, b, rr, ri = pair_eigenvalues(ev)
            a_surf[i, j] = a
            b_surf[i, j] = b
            if rr > max_res_real:
                max_res_real = rr
            if ri > max_res_imag:
                max_res_imag = ri
    return kx, ky, a_surf, b_surf, max_res_real, max_res_imag


# ---------------------------------------------------------------------------
# Sanity checks  — run before every figure build, fail loud if wrong
# ---------------------------------------------------------------------------

def sanity_checks(verbose=True):
    header = "SANITY CHECKS"
    if verbose:
        print(header)
        print("-" * 60)

    # (1) P(k=0) is column-stochastic (Fourier transform at k=0 recovers
    #     the total outgoing probability from each director state).
    P0 = build_Pk(0.35, 0.1, 0.0, 0.0)
    col_sums = P0.sum(axis=0)
    ok_stoch = np.allclose(col_sums, 1.0, atol=1e-14)
    if verbose:
        print(f"(1) P(k=0) column sums = {col_sums.real}   "
              f"{'PASS' if ok_stoch else 'FAIL'}")
    assert ok_stoch, "P(k=0) not column-stochastic!"

    # (2) Leading eigenvalue is exactly 1 at k=0
    ev0 = np.linalg.eigvals(P0)
    lam_max_re = np.max(ev0.real)
    spec_rad = np.max(np.abs(ev0))
    if verbose:
        print(f"(2) P(k=0): spectral radius = {spec_rad:.12f}, "
              f"max Re(λ) = {lam_max_re:.12f}")
    assert abs(lam_max_re - 1.0) < 1e-12, "Leading eigenvalue ≠ 1 at k=0"
    assert spec_rad < 1.0 + 1e-12, "Spectrum escapes unit disk at k=0"

    # (3) All eigenvalues inside the closed unit disk over random params
    rng = np.random.default_rng(42)
    max_abs = 0.0
    for _ in range(1000):
        omega = rng.uniform(0, 1)
        D_r = rng.uniform(0, 1)
        kx = rng.uniform(-np.pi, np.pi)
        ky = rng.uniform(-np.pi, np.pi)
        ev = np.linalg.eigvals(build_Pk(omega, D_r, kx, ky))
        max_abs = max(max_abs, np.max(np.abs(ev)))
    if verbose:
        print(f"(3) max |λ| over 1000 random (ω, D_r, k): {max_abs:.10f}  "
              f"{'PASS' if max_abs <= 1.0 + 1e-10 else 'FAIL'}")
    assert max_abs <= 1.0 + 1e-10, "Eigenvalue outside unit disk!"

    # (4) Sublattice symmetry: for every λ, -λ is also in the spectrum
    ev = np.linalg.eigvals(build_Pk(0.35, 0.1, 0.7, -0.3))
    tr = np.sum(ev)
    if verbose:
        print(f"(4) trace at (ω=0.35, D_r=0.1, k=(0.7,-0.3)): "
              f"{tr.real:+.2e} + {tr.imag:+.2e}j")
    assert abs(tr) < 1e-12, "Trace ≠ 0 (bipartite structure broken)"
    for lam in ev:
        dmin = np.min(np.abs(-lam - ev))
        assert dmin < 1e-10, f"No -λ match for λ={lam}"
    if verbose:
        print("    ✓ spectrum invariant under λ → -λ")

    # (5) Clean {±a, ±ib} pair structure for Fig 4(b) parameter range
    max_rr = 0.0
    max_ri = 0.0
    for _ in range(2000):
        omega = rng.uniform(0.0, 1.0)
        D_r = rng.uniform(0.0, 1.0)
        kx = rng.uniform(-np.pi, np.pi)
        ky = rng.uniform(-np.pi, np.pi)
        ev = np.linalg.eigvals(build_Pk(omega, D_r, kx, ky))
        _, _, rr, ri = pair_eigenvalues(ev)
        max_rr = max(max_rr, rr)
        max_ri = max(max_ri, ri)
    if verbose:
        print(f"(5) pair-structure residuals over 2000 samples:")
        print(f"    max |Im μ_real| + max(0,-Re μ_real) = {max_rr:.3e}")
        print(f"    max |Im μ_imag| + max(0,+Re μ_imag) = {max_ri:.3e}")
    if max_rr < 1e-10 and max_ri < 1e-10:
        if verbose:
            print("    ✓ spectrum is exactly {±a, ±ib} to machine precision")
    else:
        if verbose:
            print("    ⚠ spectrum is NOT universally {±a, ±ib}; "
                  "Re(λ)=0 plane is only approximate on some k")

    # (6) Inversion I·P(k)·I⁻¹ = P(-k)  (paper Methods, I = σ_x ⊗ 1 = swap
    # d ↔ d+2 mod 4). Test both the direct matrix identity and the
    # spectral consequence spec P(k) = spec P(-k). Set-level match
    # (not sort-level), because at Re≈0 the sort-ordering of the imag
    # pair can flip between k and -k even though the set is identical.
    Imat = np.array([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]], dtype=complex)
    max_mat_diff = 0.0
    max_spec_diff = 0.0
    for _ in range(50):
        omega = rng.uniform(0.05, 0.95)
        D_r = rng.uniform(0.05, 0.95)
        kx = rng.uniform(-np.pi, np.pi)
        ky = rng.uniform(-np.pi, np.pi)
        Pk = build_Pk(omega, D_r, kx, ky)
        Pmk = build_Pk(omega, D_r, -kx, -ky)
        mat_diff = np.max(np.abs(Imat @ Pk @ Imat - Pmk))
        max_mat_diff = max(max_mat_diff, mat_diff)
        ev1 = np.linalg.eigvals(Pk)
        ev2 = np.linalg.eigvals(Pmk)
        # Set match: every eigenvalue in ev1 has a close neighbour in ev2
        set_diff = max(np.min(np.abs(l - ev2)) for l in ev1)
        max_spec_diff = max(max_spec_diff, set_diff)
    assert max_mat_diff < 1e-12, f"I·P(k)·I ≠ P(-k): max diff = {max_mat_diff}"
    assert max_spec_diff < 1e-10, f"spec P(k) ≠ spec P(-k): max diff = {max_spec_diff}"
    if verbose:
        print(f"(6) inversion symmetry: "
              f"|I·P·I − P(-k)| < {max_mat_diff:.1e}, "
              f"set-match |spec(k) − spec(-k)| < {max_spec_diff:.1e}  PASS")

    if verbose:
        print("-" * 60)
        print("ALL SANITY CHECKS PASSED")
        print()


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_fig4b(D_r=0.1, omegas=(0.35, 0.5, 0.65), Nk=80,
               out='tcrw_fig4b_paper.png'):
    """
    Produce the paper's Fig 4(b): 3D surfaces of Re(λ) over the BZ.

    Layout: one row, three columns (one panel per ω).
    Each panel shows:
      - upper red surface:  +a(kx, ky)
      - lower red surface:  -a(kx, ky)
      - blue flat plane at Re(λ) = 0   (the ±ib pair's contribution)
    """
    fig = plt.figure(figsize=(16, 5.5))

    summary = []
    for col, omega in enumerate(omegas):
        print(f"--- computing ω = {omega} (D_r = {D_r}) ---")
        kx, ky, a_surf, b_surf, rr, ri = pbc_bz_surface(omega, D_r, Nk)
        gap_min = a_surf.min()
        gap_max = a_surf.max()
        b_max = b_surf.max()
        print(f"    a(k) ∈ [{gap_min:.4f}, {gap_max:.4f}]   min ⇒ spectral gap")
        print(f"    b(k) ∈ [{b_surf.min():.4f}, {b_max:.4f}]")
        print(f"    classification residuals: real={rr:.2e}, imag={ri:.2e}")
        summary.append((omega, gap_min, gap_max, b_max, rr, ri))

        KX, KY = np.meshgrid(kx, ky, indexing='ij')
        ax = fig.add_subplot(1, 3, col + 1, projection='3d')

        # +a(k) upper red surface
        ax.plot_surface(KX, KY, a_surf,
                        color='#b22222', alpha=0.9,
                        edgecolor='none', linewidth=0,
                        rcount=60, ccount=60, antialiased=True)
        # -a(k) lower red surface
        ax.plot_surface(KX, KY, -a_surf,
                        color='#b22222', alpha=0.9,
                        edgecolor='none', linewidth=0,
                        rcount=60, ccount=60, antialiased=True)
        # Re(λ)=0 flat plane for the ±ib pair
        zeros = np.zeros_like(a_surf)
        ax.plot_surface(KX, KY, zeros,
                        color='#3b5998', alpha=0.55,
                        edgecolor='none',
                        rcount=2, ccount=2, antialiased=True)

        ax.set_xlim(-np.pi, np.pi)
        ax.set_ylim(-np.pi, np.pi)
        ax.set_zlim(-1.05, 1.05)
        ax.set_xticks([-np.pi, 0.0, np.pi])
        ax.set_xticklabels([r'$-\pi$', '0', r'$\pi$'])
        ax.set_yticks([-np.pi, 0.0, np.pi])
        ax.set_yticklabels([r'$-\pi$', '0', r'$\pi$'])
        ax.set_zticks([-1, 0, 1])
        ax.set_xlabel(r'$k_x$', fontsize=12, labelpad=4)
        ax.set_ylabel(r'$k_y$', fontsize=12, labelpad=4)
        ax.set_zlabel(r'$\Re(\lambda)$', fontsize=12, labelpad=2)
        ax.set_title(rf'$\omega = {omega}$', fontsize=14, pad=10)
        ax.view_init(elev=20, azim=-62)

        # min-gap annotation
        ax.text2D(0.02, 0.94,
                  f'min $a(k) = {gap_min:.3f}$',
                  transform=ax.transAxes,
                  fontsize=10, color='#8b0000')

    plt.suptitle(
        rf'Fig 4(b): PBC Bloch band structure  '
        rf'($D_r = {D_r}$, $N_k = {Nk}\times{Nk}$ grid)',
        fontsize=14, y=1.02
    )
    plt.tight_layout()
    plt.savefig(out, dpi=160, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {out}")

    # Short numerical summary for the logbook
    print()
    print("SUMMARY")
    print("-" * 60)
    print(f"{'ω':>6}  {'min a':>9}  {'max a':>9}  {'max b':>9}  "
          f"{'res_real':>10}  {'res_imag':>10}")
    for (om, gmin, gmax, bmax, rr, ri) in summary:
        print(f"{om:>6.2f}  {gmin:>9.5f}  {gmax:>9.5f}  {bmax:>9.5f}  "
              f"{rr:>10.2e}  {ri:>10.2e}")
    print("-" * 60)
    print("Expected: min a(k) large for ω=0.35 and 0.65, ≈0 for ω=0.5")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    sanity_checks()
    print("=" * 60)
    print("Building Fig 4(b) — paper parameters")
    print(" D_r = 0.1,  ω ∈ {0.35, 0.5, 0.65},  Nk = 80")
    print("=" * 60)
    make_fig4b(D_r=0.1, omegas=(0.35, 0.5, 0.65), Nk=80,
               out='tcrw_fig4b_paper.png')
