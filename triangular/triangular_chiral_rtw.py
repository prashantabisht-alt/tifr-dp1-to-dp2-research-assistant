"""Starter code for the chiral run-and-tumble walker on a triangular lattice.

This is the Phase-1 model after the PI meeting:

- position is a triangular-lattice site in axial coordinates (n1, n2)
- director m has 6 states: 0, 1, 2, 3, 4, 5
- a run moves only along the current director, at rate v
- tumbles change the director:
    m -> m + 1 at rate gamma_plus
    m -> m - 1 at rate gamma_minus
    m -> m + 3 at rate gamma_r

Matrix convention:

    d P / dt = M P

Columns are source states and rows are destination states.
"""

from __future__ import annotations

import numpy as np


N_DIR = 6


def triangular_deltas() -> np.ndarray:
    """Return the six nearest-neighbor steps in axial coordinates."""
    return np.array(
        [
            [1, 0],    # m = 0
            [0, 1],    # m = 1
            [-1, 1],   # m = 2
            [-1, 0],   # m = 3
            [0, -1],   # m = 4
            [1, -1],   # m = 5
        ],
        dtype=float,
    )


def rates_from_gamma_bias(gamma: float, bias: float) -> tuple[float, float]:
    """Convert total adjacent-turn rate and chirality bias to +/- rates."""
    if not -1.0 <= bias <= 1.0:
        raise ValueError("bias must lie in [-1, 1]")

    gamma_plus = 0.5 * gamma * (1.0 + bias)
    gamma_minus = 0.5 * gamma * (1.0 - bias)
    return gamma_plus, gamma_minus


def build_Mk_chiral_rtw(
    k1: float,
    k2: float,
    *,
    v: float = 1.0,
    gamma_plus: float = 0.5,
    gamma_minus: float = 0.5,
    gamma_r: float = 0.0,
) -> np.ndarray:
    """Build the 6x6 Bloch generator M(k) for the triangular chiral RTW.

    The Fourier convention is

        P_tilde(k) = sum_n exp(i k.n) P(n).

    Therefore the incoming run term P_m(n - Delta_m) becomes

        exp(i k.Delta_m) P_tilde_m(k).
    """
    if min(v, gamma_plus, gamma_minus, gamma_r) < 0.0:
        raise ValueError("all rates must be nonnegative")

    deltas = triangular_deltas()
    phases = deltas @ np.array([k1, k2], dtype=float)
    total_out = v + gamma_plus + gamma_minus + gamma_r

    M = np.zeros((N_DIR, N_DIR), dtype=complex)

    for m in range(N_DIR):
        M[m, m] = v * np.exp(1j * phases[m]) - total_out

        # Column m is the source director. Rows are destination directors.
        M[(m + 1) % N_DIR, m] += gamma_plus
        M[(m - 1) % N_DIR, m] += gamma_minus
        M[(m + 3) % N_DIR, m] += gamma_r

    return M


def rotate_k_60_for_spectrum(k1: float, k2: float) -> tuple[float, float]:
    """Rotate a Fourier covector by one 60-degree step in axial coordinates.

    This convention is the one that shifts the Bloch phases as

        phi_m(k_rot) = phi_{m + 1}(k).

    Since the director labels are cyclic, M(k_rot) is similar to M(k), so the
    eigenvalue spectrum must be unchanged.
    """
    return k2, k2 - k1


def reflect_k_for_spectrum(k1: float, k2: float) -> tuple[float, float]:
    """Reflect a Fourier covector across the e0 axis in axial coordinates.

    This reflection maps director phases as

        phi_m(k_reflected) = phi_{5 - m}(k).

    For the chiral model this operation also swaps gamma_plus and gamma_minus.
    Therefore mirror symmetry holds only when gamma_plus == gamma_minus.
    """
    return k1 - k2, -k2


def reflection_permutation() -> np.ndarray:
    """Permutation matrix for director reflection m -> 5 - m."""
    perm = np.zeros((N_DIR, N_DIR), dtype=float)
    for m in range(N_DIR):
        perm[m, (5 - m) % N_DIR] = 1.0
    return perm


def sorted_eigenvalues(M: np.ndarray) -> np.ndarray:
    """Sort complex eigenvalues in a stable display/check order."""
    vals = np.linalg.eigvals(M)
    order = np.lexsort((vals.imag, vals.real))
    return vals[order]


def max_spectral_difference(a: np.ndarray, b: np.ndarray) -> float:
    """Maximum absolute difference between the sorted eigenvalue magnitudes.

    Using magnitudes avoids the lexsort ambiguity when complex-conjugate
    eigenvalues appear in different order for M(k) vs M(k_rot).
    """
    va = np.sort(np.abs(np.linalg.eigvals(a)))
    vb = np.sort(np.abs(np.linalg.eigvals(b)))
    return float(np.max(np.abs(va - vb)))


def top_eigenvalue(k1: float, k2: float, **rates) -> complex:
    """Return the eigenvalue of M(k) closest to zero (the hydrodynamic mode)."""
    M = build_Mk_chiral_rtw(k1, k2, **rates)
    vals = np.linalg.eigvals(M)
    return vals[np.argmin(np.abs(vals))]


def numerical_diffusion_tensor(
    *,
    v: float = 1.0,
    gamma_plus: float = 0.5,
    gamma_minus: float = 0.5,
    gamma_r: float = 0.0,
    dk: float = 1e-5,
) -> dict:
    """Extract the 2x2 diffusion tensor from the curvature of lambda_0(k).

    Uses centered finite differences:

        D_ij = -(1/2) d^2 lambda_0 / (dk_i dk_j)

    evaluated at k = 0.  Returns Cartesian D_ij after converting from axial
    coordinates via the triangular-lattice metric.

    The axial second derivatives give the tensor in the (k1, k2) basis.
    Cartesian coordinates use x = n1 + n2/2, y = sqrt(3) n2 / 2.
    The Jacobian from Cartesian k to axial k is:

        k1 = kx,  k2 = kx/2 + ky sqrt(3)/2

    so derivatives transform accordingly.
    """
    rates = dict(v=v, gamma_plus=gamma_plus, gamma_minus=gamma_minus, gamma_r=gamma_r)

    # --- second derivatives in axial coordinates via finite differences ---

    # d^2 lambda_0 / dk1^2
    lp = top_eigenvalue(dk, 0.0, **rates)
    lm = top_eigenvalue(-dk, 0.0, **rates)
    l0 = top_eigenvalue(0.0, 0.0, **rates)
    d2_dk1dk1 = (lp + lm - 2.0 * l0) / dk**2

    # d^2 lambda_0 / dk2^2
    lp = top_eigenvalue(0.0, dk, **rates)
    lm = top_eigenvalue(0.0, -dk, **rates)
    d2_dk2dk2 = (lp + lm - 2.0 * l0) / dk**2

    # d^2 lambda_0 / dk1 dk2  (mixed)
    lpp = top_eigenvalue(dk, dk, **rates)
    lpm = top_eigenvalue(dk, -dk, **rates)
    lmp = top_eigenvalue(-dk, dk, **rates)
    lmm = top_eigenvalue(-dk, -dk, **rates)
    d2_dk1dk2 = (lpp - lpm - lmp + lmm) / (4.0 * dk**2)

    # Axial symmetric diffusion tensor:
    # lambda_0(k) = -D_ij k_i k_j + O(k^3), so
    # D_ij = -(1/2) d^2 lambda_0 / dk_i dk_j.
    D_ax = -0.5 * np.array(
        [[d2_dk1dk1, d2_dk1dk2], [d2_dk1dk2, d2_dk2dk2]], dtype=complex
    )

    # --- convert to Cartesian ---
    # Axial basis vectors: a1 = (1, 0), a2 = (1/2, sqrt3/2)
    # Reciprocal relation: k_axial = J^T k_cartesian
    #   where J = [[1, 1/2], [0, sqrt3/2]]  (columns are a1, a2)
    # So d/dk_cart = J d/dk_ax, and D_cart = J D_ax J^T
    sqrt3 = np.sqrt(3.0)
    J = np.array([[1.0, 0.5], [0.0, 0.5 * sqrt3]])
    D_cart = J @ D_ax @ J.T

    D_even = 0.5 * np.real(D_cart[0, 0] + D_cart[1, 1])
    D_odd = 0.5 * np.real(D_cart[0, 1] - D_cart[1, 0])

    return {
        "D_axial": D_ax,
        "D_cartesian": D_cart,
        "D_even": float(D_even),
        "D_odd": float(D_odd),
        "D_xx": float(np.real(D_cart[0, 0])),
        "D_yy": float(np.real(D_cart[1, 1])),
        "D_xy": float(np.real(D_cart[0, 1])),
        "D_yx": float(np.real(D_cart[1, 0])),
    }


def diffusion_tensor_green_kubo(
    *,
    v: float = 1.0,
    gamma_plus: float = 0.5,
    gamma_minus: float = 0.5,
    gamma_r: float = 0.0,
) -> dict:
    """Full diffusion tensor (including odd part) via Green-Kubo.

    The velocity autocorrelation is

        C^v_ij(t) = v^2 sum_{m,m'} (1/6) Delta_m^i [exp(R t)]_{m',m} Delta_{m'}^j

    where R = M(k=0) is the pure-rotation generator.  The diffusion tensor is

        D_ij = integral_0^infty C^v_ij(t) dt
             = -v^2 (1/6) Delta^i . R^{-1}_reduced . Delta^j

    R has a zero eigenvalue (probability conservation), but the zero mode
    drops out because sum_m Delta_m = 0.  We use the pseudoinverse R^+.

    For Poisson lattice jumps there is also a local jump-noise contribution:

        D_jump = (v / 2) <Delta_i Delta_j>_uniform.

    This is symmetric and must be added to the integrated persistent
    velocity-correlation part. The odd part comes only from the persistent
    correlation contribution.

    The result is converted to Cartesian coordinates (x, y) where
    x = n1 + n2/2, y = sqrt(3) n2 / 2.
    """
    # build k=0 rotation matrix
    R = build_Mk_chiral_rtw(0.0, 0.0, v=v, gamma_plus=gamma_plus,
                             gamma_minus=gamma_minus, gamma_r=gamma_r)
    R = np.real(R)  # pure real at k=0

    # pseudoinverse (drops the zero mode)
    R_pinv = np.linalg.pinv(R)

    # axial displacement vectors: Delta_m as rows of (6, 2) array
    deltas = triangular_deltas()  # shape (6, 2), columns are (dn1, dn2)

    # Persistent part:
    # D^ax_ij = -v^2 (1/6) sum_{m,m'} Delta_m^i  R^+_{m',m}  Delta_{m'}^j
    #         = -v^2 / 6  *  deltas^T . R_pinv^T . deltas
    # Note: R_pinv_{m',m} means row m', col m. We sum over m (source) and m'.
    D_ax_persistent = -(v**2 / N_DIR) * (deltas.T @ R_pinv.T @ deltas)

    # Bare Poisson jump-noise part. It is present because jumps occur at random
    # times with finite step length, even before velocity persistence is counted.
    D_ax_jump = 0.5 * v * (deltas.T @ deltas) / N_DIR
    D_ax = D_ax_persistent + D_ax_jump

    # convert to Cartesian
    sqrt3 = np.sqrt(3.0)
    J = np.array([[1.0, 0.5], [0.0, 0.5 * sqrt3]])
    D_cart_persistent = J @ D_ax_persistent @ J.T
    D_cart_jump = J @ D_ax_jump @ J.T
    D_cart = J @ D_ax @ J.T

    D_even = 0.5 * (D_cart[0, 0] + D_cart[1, 1])
    D_odd = 0.5 * (D_cart[0, 1] - D_cart[1, 0])

    return {
        "D_axial": D_ax,
        "D_axial_persistent": D_ax_persistent,
        "D_axial_jump": D_ax_jump,
        "D_cartesian": D_cart,
        "D_cartesian_persistent": D_cart_persistent,
        "D_cartesian_jump": D_cart_jump,
        "D_even": float(D_even),
        "D_odd": float(D_odd),
        "D_xx": float(D_cart[0, 0]),
        "D_yy": float(D_cart[1, 1]),
        "D_xy": float(D_cart[0, 1]),
        "D_yx": float(D_cart[1, 0]),
    }


def diffusion_tensor_exact(
    *,
    v: float = 1.0,
    gamma: float = 1.0,
    b: float = 0.0,
    gamma_r: float = 0.0,
) -> dict:
    """Exact closed-form diffusion tensor for the chiral RTW.

    The rate matrix R = M(k=0) is a 6x6 circulant.  Its eigenvalues are

        lambda_k = -gamma - gamma_r
                   + gamma cos(pi k/3) + i gamma b sin(pi k/3)
                   + gamma_r (-1)^k,    k = 0, ..., 5.

    Only k = 1 and k = 5 = bar{1} contribute to the diffusion tensor
    (the hexagonal displacement DFT vanishes at k = 0, 2, 3, 4).

    Defining alpha = gamma/2 + 2 gamma_r,  beta = sqrt(3) gamma b / 2,
    and |lambda_1|^2 = alpha^2 + beta^2, the result is:

        D_even = v^2 alpha / (2 |lambda_1|^2) + v/4
        D_odd  = v^2 beta  / (2 |lambda_1|^2)

    The v/4 term is bare Poisson jump noise; the persistent velocity
    correlation contributes to both D_even and D_odd.  The odd part
    is proportional to Im(1/lambda_1) and vanishes when b = 0.
    """
    alpha = gamma / 2 + 2 * gamma_r
    beta = np.sqrt(3.0) * gamma * b / 2
    lam1_sq = alpha**2 + beta**2

    D_even = v**2 * alpha / (2 * lam1_sq) + v / 4
    D_odd = v**2 * beta / (2 * lam1_sq)

    return {
        "D_even": float(D_even),
        "D_odd": float(D_odd),
        "alpha": float(alpha),
        "beta": float(beta),
        "lam1_sq": float(lam1_sq),
    }


def run_self_checks() -> None:
    """Run the first sanity checks we trust before making any plots."""
    v = 1.0
    gamma = 0.8
    bias = 0.35
    gamma_r = 0.2
    gamma_plus, gamma_minus = rates_from_gamma_bias(gamma, bias)

    M0 = build_Mk_chiral_rtw(
        0.0,
        0.0,
        v=v,
        gamma_plus=gamma_plus,
        gamma_minus=gamma_minus,
        gamma_r=gamma_r,
    )

    print("parameters")
    print(f"v = {v}")
    print(f"gamma_plus = {gamma_plus}")
    print(f"gamma_minus = {gamma_minus}")
    print(f"gamma_r = {gamma_r}")
    print()

    print("M(0, 0) =")
    print(np.array2string(M0, precision=4, suppress_small=True))
    print()

    print("column sums at k=0 =")
    print(np.array2string(M0.sum(axis=0), precision=12, suppress_small=True))
    print()

    vals0 = sorted_eigenvalues(M0)
    print("eigenvalues at k=0 =")
    print(np.array2string(vals0, precision=12, suppress_small=True))
    print()

    print("smallest |lambda(k=0)| =")
    print(np.min(np.abs(vals0)))
    print()

    k1, k2 = 0.37, -0.21
    M = build_Mk_chiral_rtw(
        k1,
        k2,
        v=v,
        gamma_plus=gamma_plus,
        gamma_minus=gamma_minus,
        gamma_r=gamma_r,
    )
    rk1, rk2 = rotate_k_60_for_spectrum(k1, k2)
    M_rot = build_Mk_chiral_rtw(
        rk1,
        rk2,
        v=v,
        gamma_plus=gamma_plus,
        gamma_minus=gamma_minus,
        gamma_r=gamma_r,
    )

    print("nonzero-k test point =")
    print((k1, k2))
    print("rotated k =")
    print((rk1, rk2))
    print()

    print("max spectral difference under C6 rotation =")
    print(max_spectral_difference(M, M_rot))
    print()

    mk1, mk2 = reflect_k_for_spectrum(k1, k2)
    M_ref_same_rates = build_Mk_chiral_rtw(
        mk1,
        mk2,
        v=v,
        gamma_plus=gamma_plus,
        gamma_minus=gamma_minus,
        gamma_r=gamma_r,
    )
    M_ref_swapped_rates = build_Mk_chiral_rtw(
        mk1,
        mk2,
        v=v,
        gamma_plus=gamma_minus,
        gamma_minus=gamma_plus,
        gamma_r=gamma_r,
    )

    print("reflected k =")
    print((mk1, mk2))
    print()

    reflection = reflection_permutation()
    reflected_M = reflection @ M @ reflection.T

    print("max matrix-covariance error under mirror, same rates =")
    print(np.max(np.abs(M_ref_same_rates - reflected_M)))
    print()

    print("max matrix-covariance error under mirror, swapped +/- rates =")
    print(np.max(np.abs(M_ref_swapped_rates - reflected_M)))

    gamma_achiral_plus, gamma_achiral_minus = rates_from_gamma_bias(gamma, 0.0)
    M_achiral = build_Mk_chiral_rtw(
        k1,
        k2,
        v=v,
        gamma_plus=gamma_achiral_plus,
        gamma_minus=gamma_achiral_minus,
        gamma_r=gamma_r,
    )
    M_achiral_ref = build_Mk_chiral_rtw(
        mk1,
        mk2,
        v=v,
        gamma_plus=gamma_achiral_plus,
        gamma_minus=gamma_achiral_minus,
        gamma_r=gamma_r,
    )
    reflected_M_achiral = reflection @ M_achiral @ reflection.T
    print()
    print("max matrix-covariance error under mirror at b=0 =")
    print(np.max(np.abs(M_achiral_ref - reflected_M_achiral)))

    # --- diffusion tensor ---
    print()
    print("=" * 60)
    print("DIFFUSION TENSOR")
    print("=" * 60)

    for b_test in [0.0, 0.35, 0.8]:
        gp_t, gm_t = rates_from_gamma_bias(gamma, b_test)
        rates_kw = dict(v=v, gamma_plus=gp_t, gamma_minus=gm_t, gamma_r=gamma_r)

        D_eig = numerical_diffusion_tensor(**rates_kw)
        D_gk = diffusion_tensor_green_kubo(**rates_kw)
        D_ex = diffusion_tensor_exact(v=v, gamma=gamma, b=b_test, gamma_r=gamma_r)

        print(f"\nbias = {b_test}")
        print(f"  --- eigenvalue method (symmetric part only) ---")
        print(f"  D_even = {D_eig['D_even']:.8f}")
        print(f"  D_odd  = {D_eig['D_odd']:.8f}  (must be ~0)")
        print(f"  --- Green-Kubo (full tensor including odd part) ---")
        print(f"  D_even = {D_gk['D_even']:.8f}")
        print(f"  D_odd  = {D_gk['D_odd']:.8f}")
        print(f"  D_xx == D_yy? diff = {abs(D_gk['D_xx'] - D_gk['D_yy']):.2e}")
        print(f"  D_xy == -D_yx? diff = {abs(D_gk['D_xy'] + D_gk['D_yx']):.2e}")
        print(f"  --- exact closed form ---")
        print(f"  D_even = {D_ex['D_even']:.8f}")
        print(f"  D_odd  = {D_ex['D_odd']:.8f}")
        print(f"  eig vs GK:    ΔD_even = {abs(D_eig['D_even'] - D_gk['D_even']):.2e}")
        print(f"  exact vs GK:  ΔD_even = {abs(D_ex['D_even'] - D_gk['D_even']):.2e}")
        print(f"  exact vs GK:  ΔD_odd  = {abs(D_ex['D_odd'] - D_gk['D_odd']):.2e}")


if __name__ == "__main__":
    run_self_checks()
