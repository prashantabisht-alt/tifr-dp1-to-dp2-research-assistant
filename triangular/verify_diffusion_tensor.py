"""Thorough cross-verification of the diffusion tensor computation.

Three independent methods:
1. Eigenvalue curvature (numerical finite differences)
2. Green-Kubo (pseudoinverse of rate matrix)
3. Direct eigendecomposition of the circulant R (analytical DFT)

Plus sanity checks:
- D_odd must vanish at b=0
- D_even must be isotropic (D_xx = D_yy)
- D_xy = -D_yx
- Known exact values at special points
- Eigenvalue method D_even must match Green-Kubo D_even
"""

import numpy as np
from triangular_chiral_rtw import (
    build_Mk_chiral_rtw, rates_from_gamma_bias, triangular_deltas,
    diffusion_tensor_green_kubo, numerical_diffusion_tensor,
    N_DIR
)

PASS = 0
FAIL = 0

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  PASS  {name}")
    else:
        FAIL += 1
        print(f"  FAIL  {name}  {detail}")


# ======================================================================
# CHECK 1: Probability conservation (column sums = 0) at k=0
# ======================================================================
print("=" * 70)
print("CHECK 1: Probability conservation at k=0")
print("=" * 70)
for v, gamma, b, gr in [(1,1,0,0), (1,0.8,0.35,0.2), (2,0.5,1.0,0.1), (1,1,1,0)]:
    gp, gm = rates_from_gamma_bias(gamma, b)
    M0 = build_Mk_chiral_rtw(0, 0, v=v, gamma_plus=gp, gamma_minus=gm, gamma_r=gr)
    col_sums = np.abs(M0.sum(axis=0))
    maxerr = np.max(col_sums)
    check(f"v={v} γ={gamma} b={b} γ_r={gr}", maxerr < 1e-14, f"max col sum = {maxerr:.2e}")


# ======================================================================
# CHECK 2: Eigenvalue method vs Green-Kubo — D_even match
# ======================================================================
print("\n" + "=" * 70)
print("CHECK 2: Eigenvalue D_even vs Green-Kubo D_even")
print("=" * 70)
test_params = [
    (1.0, 1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0, 0.2),
    (1.0, 0.8, 0.35, 0.2),
    (1.0, 0.8, 0.8, 0.2),
    (1.0, 1.0, 1.0, 0.0),
    (1.0, 1.0, 1.0, 0.5),
    (2.0, 0.5, 0.5, 0.1),
    (1.0, 0.3, 0.9, 0.0),
    (1.0, 2.0, 0.7, 0.3),
    (0.5, 0.5, 0.5, 0.5),
]

for v, gamma, b, gr in test_params:
    gp, gm = rates_from_gamma_bias(gamma, b)
    kw = dict(v=v, gamma_plus=gp, gamma_minus=gm, gamma_r=gr)
    De_eig = numerical_diffusion_tensor(**kw, dk=1e-5)['D_even']
    De_gk = diffusion_tensor_green_kubo(**kw)['D_even']
    reldiff = abs(De_eig - De_gk) / max(abs(De_gk), 1e-15)
    check(f"v={v} γ={gamma} b={b} γ_r={gr}: eig={De_eig:.8f} gk={De_gk:.8f}",
          reldiff < 1e-4, f"rel diff = {reldiff:.2e}")


# ======================================================================
# CHECK 3: Isotropy (D_xx = D_yy) and antisymmetry (D_xy = -D_yx)
# ======================================================================
print("\n" + "=" * 70)
print("CHECK 3: C6 isotropy and antisymmetry")
print("=" * 70)
for v, gamma, b, gr in test_params:
    gp, gm = rates_from_gamma_bias(gamma, b)
    kw = dict(v=v, gamma_plus=gp, gamma_minus=gm, gamma_r=gr)
    res = diffusion_tensor_green_kubo(**kw)
    iso_err = abs(res['D_xx'] - res['D_yy'])
    anti_err = abs(res['D_xy'] + res['D_yx'])
    check(f"v={v} γ={gamma} b={b} γ_r={gr}: D_xx-D_yy={iso_err:.2e}, D_xy+D_yx={anti_err:.2e}",
          iso_err < 1e-14 and anti_err < 1e-14)


# ======================================================================
# CHECK 4: D_odd = 0 at b=0 (no chirality → no odd diffusion)
# ======================================================================
print("\n" + "=" * 70)
print("CHECK 4: D_odd = 0 at b=0")
print("=" * 70)
for v, gamma, gr in [(1,1,0), (1,0.5,0.2), (2,1,0.5), (0.5,2,0.1)]:
    gp, gm = rates_from_gamma_bias(gamma, 0.0)
    res = diffusion_tensor_green_kubo(v=v, gamma_plus=gp, gamma_minus=gm, gamma_r=gr)
    check(f"v={v} γ={gamma} γ_r={gr}: D_odd={res['D_odd']:.2e}",
          abs(res['D_odd']) < 1e-14)


# ======================================================================
# CHECK 5: D_odd sign flips with b → -b
# ======================================================================
print("\n" + "=" * 70)
print("CHECK 5: D_odd(b) = -D_odd(-b)")
print("=" * 70)
for v, gamma, b, gr in [(1,1,0.3,0), (1,0.8,0.7,0.2), (2,0.5,0.5,0.1)]:
    gp_pos, gm_pos = rates_from_gamma_bias(gamma, b)
    gp_neg, gm_neg = rates_from_gamma_bias(gamma, -b)
    res_pos = diffusion_tensor_green_kubo(v=v, gamma_plus=gp_pos, gamma_minus=gm_pos, gamma_r=gr)
    res_neg = diffusion_tensor_green_kubo(v=v, gamma_plus=gp_neg, gamma_minus=gm_neg, gamma_r=gr)
    sum_odd = abs(res_pos['D_odd'] + res_neg['D_odd'])
    diff_even = abs(res_pos['D_even'] - res_neg['D_even'])
    check(f"v={v} γ={gamma} b=±{b} γ_r={gr}: D_odd sum={sum_odd:.2e}, D_even diff={diff_even:.2e}",
          sum_odd < 1e-14 and diff_even < 1e-14)


# ======================================================================
# CHECK 6: Known exact values at b=1, γ_r=0, v=γ=1
# ======================================================================
print("\n" + "=" * 70)
print("CHECK 6: Exact values at b=1, γ_r=0, v=γ=1")
print("=" * 70)
gp, gm = rates_from_gamma_bias(1.0, 1.0)
res = diffusion_tensor_green_kubo(v=1.0, gamma_plus=gp, gamma_minus=gm, gamma_r=0.0)
check(f"D_even = 1/2 ? got {res['D_even']:.15f}", abs(res['D_even'] - 0.5) < 1e-14)
check(f"D_odd = √3/4 ? got {res['D_odd']:.15f}", abs(res['D_odd'] - np.sqrt(3)/4) < 1e-14)
check(f"D_odd/D_even = √3/2 ? got {res['D_odd']/res['D_even']:.15f}",
      abs(res['D_odd']/res['D_even'] - np.sqrt(3)/2) < 1e-14)


# ======================================================================
# CHECK 7: Analytical DFT cross-check for the circulant case
# ======================================================================
print("\n" + "=" * 70)
print("CHECK 7: DFT-based analytical diffusion tensor (circulant case)")
print("=" * 70)

def diffusion_tensor_dft(*, v=1.0, gamma_plus=0.5, gamma_minus=0.5, gamma_r=0.0):
    """Compute D via explicit DFT diagonalization of the circulant rate matrix.

    R is circulant with first column:
      c = [-(γ+ + γ- + γ_r), γ+, 0, γ_r, 0, γ-]
    Eigenvalues: λ_k = Σ_j c_j ω^{jk}, ω = exp(2πi/6)

    The pseudoinverse in the DFT basis is:
      R^+_{m,m'} = (1/6) Σ_{k≠0} (1/λ_k) ω^{k(m-m')}

    Then D_ax = -(v²/6) Δ^T R^{+T} Δ + (v/12) Δ^T Δ
    """
    omega = np.exp(2j * np.pi / 6)
    c = np.array([
        -(gamma_plus + gamma_minus + gamma_r),  # m=0: diagonal
        gamma_plus,                                # m=1: CCW tumble
        0.0,                                       # m=2
        gamma_r,                                   # m=3: reversal
        0.0,                                       # m=4
        gamma_minus,                               # m=5: CW tumble
    ])

    # DFT eigenvalues
    lam = np.array([sum(c[j] * omega**(j*k) for j in range(6)) for k in range(6)])

    # Check λ_0 = 0
    assert abs(lam[0]) < 1e-14, f"λ_0 should be 0, got {lam[0]}"

    # Build R^+ in real space via inverse DFT of 1/λ_k (skip k=0)
    R_pinv = np.zeros((6, 6), dtype=complex)
    for m in range(6):
        for mp in range(6):
            R_pinv[m, mp] = (1.0/6) * sum(
                (1.0/lam[k]) * omega**(k*(mp - m)) for k in range(1, 6)
            )

    deltas = triangular_deltas()
    D_ax_persistent = -(v**2 / 6) * (deltas.T @ R_pinv.T @ deltas)
    D_ax_jump = (v / 2) * (deltas.T @ deltas) / 6
    D_ax = np.real(D_ax_persistent + D_ax_jump)

    sqrt3 = np.sqrt(3.0)
    J = np.array([[1.0, 0.5], [0.0, 0.5 * sqrt3]])
    D_cart = J @ D_ax @ J.T

    return {
        'D_even': 0.5 * (D_cart[0,0] + D_cart[1,1]),
        'D_odd': 0.5 * (D_cart[0,1] - D_cart[1,0]),
        'D_cart': D_cart,
    }

for v, gamma, b, gr in test_params:
    gp, gm = rates_from_gamma_bias(gamma, b)
    kw = dict(v=v, gamma_plus=gp, gamma_minus=gm, gamma_r=gr)
    res_gk = diffusion_tensor_green_kubo(**kw)
    res_dft = diffusion_tensor_dft(**kw)
    de_diff = abs(res_gk['D_even'] - res_dft['D_even'])
    do_diff = abs(res_gk['D_odd'] - res_dft['D_odd'])
    check(f"v={v} γ={gamma} b={b} γ_r={gr}: ΔD_even={de_diff:.2e} ΔD_odd={do_diff:.2e}",
          de_diff < 1e-12 and do_diff < 1e-12)


# ======================================================================
# CHECK 8: Eigenvalue method D_odd ≈ 0 (can't see antisymmetric part)
# ======================================================================
print("\n" + "=" * 70)
print("CHECK 8: Eigenvalue method gives D_odd ≈ 0")
print("=" * 70)
for v, gamma, b, gr in [(1,1,0.5,0), (1,0.8,0.8,0.2), (1,1,1,0)]:
    gp, gm = rates_from_gamma_bias(gamma, b)
    res = numerical_diffusion_tensor(v=v, gamma_plus=gp, gamma_minus=gm, gamma_r=gr)
    check(f"v={v} γ={gamma} b={b} γ_r={gr}: eig D_odd={res['D_odd']:.2e}",
          abs(res['D_odd']) < 1e-8)


# ======================================================================
# SUMMARY
# ======================================================================
print("\n" + "=" * 70)
total = PASS + FAIL
print(f"SUMMARY: {PASS}/{total} checks passed, {FAIL} failed")
print("=" * 70)
