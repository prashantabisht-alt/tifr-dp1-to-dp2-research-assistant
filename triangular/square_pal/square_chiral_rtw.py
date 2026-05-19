"""Chiral run-and-tumble walker on the SQUARE lattice (C4, 4 directors).

Same Mallikarjun-Pal-style dynamics as our triangular model, but on the
square lattice for side-by-side comparison.

Important: this is a continuous-time square-lattice analogue.  It is a useful
bridge toward Wójcik-Kalz 2026, but it is not literally their model because
Wójcik-Kalz uses a discrete-time coin-step walk.

Model:
  - Position: square-lattice site (n1, n2) in Z^2.
  - Director: m in {0,1,2,3} -> right, up, left, down.
  - Events (continuous-time Gillespie):
      run:      hop along director m,   rate v
      CCW:      m -> (m+1) mod 4,       rate gamma_+ = gamma(1+b)/2
      CW:       m -> (m-1) mod 4,       rate gamma_- = gamma(1-b)/2
      reversal: m -> (m+2) mod 4,       rate gamma_r

Matrix convention: dP/dt = M P (columns = source states).

No axial-to-Cartesian conversion needed — the square lattice IS Cartesian.
"""

from __future__ import annotations
import numpy as np


N_DIR = 4


def square_deltas() -> np.ndarray:
    """Return the four NN steps in Cartesian coordinates."""
    return np.array(
        [
            [+1,  0],   # m=0: right
            [ 0, +1],   # m=1: up
            [-1,  0],   # m=2: left
            [ 0, -1],   # m=3: down
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
    """Build the 4×4 Bloch generator M(k) for the square chiral RTW.

    Fourier convention: P̃(k) = Σ_n exp(i k·n) P(n).
    """
    if min(v, gamma_plus, gamma_minus, gamma_r) < 0.0:
        raise ValueError("all rates must be nonneg")

    deltas = square_deltas()
    phases = deltas @ np.array([k1, k2], dtype=float)
    total_out = v + gamma_plus + gamma_minus + gamma_r

    M = np.zeros((N_DIR, N_DIR), dtype=complex)
    for m in range(N_DIR):
        M[m, m] = v * np.exp(1j * phases[m]) - total_out
        M[(m + 1) % N_DIR, m] += gamma_plus    # CCW tumble
        M[(m - 1) % N_DIR, m] += gamma_minus   # CW tumble
        M[(m + 2) % N_DIR, m] += gamma_r       # reversal

    return M


def numerical_diffusion_tensor(
    *,
    v: float = 1.0,
    gamma_plus: float = 0.5,
    gamma_minus: float = 0.5,
    gamma_r: float = 0.0,
    dk: float = 1e-5,
) -> dict:
    """Diffusion tensor from λ₀(k) curvature (finite differences).

    On the square lattice, k = (kx, ky) are already Cartesian — no
    Jacobian needed.

    NOTE: D_odd from this method is always zero on the square lattice.
    This is correct, not a bug.  The antisymmetric part of D drops out
    of the quadratic form k^T D k = D_xx kx² + (D_xy+D_yx) kx ky + ...,
    so ∂²λ₀/∂kx∂ky is symmetric and cannot capture D_odd.  On the
    triangular lattice the axial→Cartesian Jacobian mixes symmetric and
    antisymmetric parts, so D_odd shows up there.  To get D_odd on the
    square lattice, use Green-Kubo (which uses the full rate matrix R⁺).
    """
    rates = dict(v=v, gamma_plus=gamma_plus, gamma_minus=gamma_minus, gamma_r=gamma_r)

    def top_eig(k1, k2):
        M = build_Mk_chiral_rtw(k1, k2, **rates)
        vals = np.linalg.eigvals(M)
        return vals[np.argmin(np.abs(vals))]

    l0 = top_eig(0, 0)

    # d²λ/dk1²
    d2_11 = (top_eig(dk, 0) + top_eig(-dk, 0) - 2*l0) / dk**2
    # d²λ/dk2²
    d2_22 = (top_eig(0, dk) + top_eig(0, -dk) - 2*l0) / dk**2
    # d²λ/dk1dk2
    d2_12 = (top_eig(dk, dk) - top_eig(dk, -dk)
             - top_eig(-dk, dk) + top_eig(-dk, -dk)) / (4*dk**2)

    D = -0.5 * np.array([[d2_11, d2_12], [d2_12, d2_22]], dtype=complex)

    D_even = 0.5 * np.real(D[0, 0] + D[1, 1])
    D_odd = 0.5 * np.real(D[0, 1] - D[1, 0])

    return {
        "D_cartesian": D,
        "D_even": float(D_even),
        "D_odd": float(D_odd),
        "D_xx": float(np.real(D[0, 0])),
        "D_yy": float(np.real(D[1, 1])),
        "D_xy": float(np.real(D[0, 1])),
        "D_yx": float(np.real(D[1, 0])),
    }


def diffusion_tensor_green_kubo(
    *,
    v: float = 1.0,
    gamma_plus: float = 0.5,
    gamma_minus: float = 0.5,
    gamma_r: float = 0.0,
) -> dict:
    """Full diffusion tensor via Green-Kubo (pseudoinverse of R).

    D_ij = -(v²/N) Δ^T R^{+T} Δ  +  (v/2)(1/N) Δ^T Δ

    On square, Δ is already Cartesian — no Jacobian.
    """
    R = build_Mk_chiral_rtw(0, 0, v=v, gamma_plus=gamma_plus,
                             gamma_minus=gamma_minus, gamma_r=gamma_r)
    R = np.real(R)
    R_pinv = np.linalg.pinv(R)

    deltas = square_deltas()

    D_pers = -(v**2 / N_DIR) * (deltas.T @ R_pinv.T @ deltas)
    D_jump = (v / 2) * (deltas.T @ deltas) / N_DIR
    D = D_pers + D_jump

    D_even = 0.5 * (D[0, 0] + D[1, 1])
    D_odd = 0.5 * (D[0, 1] - D[1, 0])

    return {
        "D_cartesian": D,
        "D_persistent": D_pers,
        "D_jump": D_jump,
        "D_even": float(D_even),
        "D_odd": float(D_odd),
        "D_xx": float(D[0, 0]),
        "D_yy": float(D[1, 1]),
        "D_xy": float(D[0, 1]),
        "D_yx": float(D[1, 0]),
    }


def diffusion_tensor_exact(
    *,
    v: float = 1.0,
    gamma: float = 1.0,
    b: float = 0.0,
    gamma_r: float = 0.0,
) -> dict:
    """Exact closed-form diffusion tensor for the square chiral RTW.

    R = M(k=0) is a 4×4 circulant.  Eigenvalues:

        λ_k = -γ - γ_r + γ cos(πk/2) + i γ b sin(πk/2) + γ_r (-1)^k
        for k = 0, 1, 2, 3.

    Only k=1 and k=3 = conj(1) contribute to diffusion (k=2 has
    vanishing structure factor |g(2)|² = 0).

        λ₁ = -(γ + 2γ_r) + i γ b

    Define:
        α = γ + 2γ_r       (note: triangular has α = γ/2 + 2γ_r)
        β = γ b             (note: triangular has β = √3 γ b / 2)
        |λ₁|² = α² + β²

    Result (UNIVERSAL formula — same factor of 2 on any lattice):
        D_even = v² α / (2|λ₁|²) + v/4
        D_odd  = v² β / (2|λ₁|²)

    Derivation:  D^{pers}_{ij} = -v² Σ_{k≠0} g_i(k) g_j*(k) / λ_k
    where g_i(k) = (1/N) Σ_m δ^i_m ω^{km}.  On the square lattice
    only the k=1,3 pair contributes, giving the 2 Re(1/λ₁) structure
    that yields the factor 1/2.  The same factor appears on the
    triangular lattice (where k=1,5 is the contributing pair).
    Only α and β differ between lattices.
    """
    alpha = gamma + 2 * gamma_r
    beta = gamma * b
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
    """Quick sanity checks, same structure as the triangular module."""
    v = 1.0
    gamma = 0.8
    bias = 0.35
    gamma_r = 0.2
    gp, gm = rates_from_gamma_bias(gamma, bias)

    M0 = build_Mk_chiral_rtw(0, 0, v=v, gamma_plus=gp, gamma_minus=gm, gamma_r=gamma_r)
    print("M(0,0) =")
    print(np.array2string(np.real(M0), precision=4))
    print()

    col_sums = np.abs(M0.sum(axis=0))
    print(f"max |col sum| at k=0 = {np.max(col_sums):.2e}")

    vals = np.linalg.eigvals(M0)
    print(f"min |eigenvalue| at k=0 = {np.min(np.abs(vals)):.2e}")
    print()

    print("=" * 50)
    print("DIFFUSION TENSOR")
    print("=" * 50)

    for b_test in [0.0, 0.35, 0.8, 1.0]:
        gp_t, gm_t = rates_from_gamma_bias(gamma, b_test)
        kw = dict(v=v, gamma_plus=gp_t, gamma_minus=gm_t, gamma_r=gamma_r)

        D_eig = numerical_diffusion_tensor(**kw)
        D_gk = diffusion_tensor_green_kubo(**kw)
        D_ex = diffusion_tensor_exact(v=v, gamma=gamma, b=b_test, gamma_r=gamma_r)

        print(f"\nbias = {b_test}")
        print(f"  eigenvalue:  D_even={D_eig['D_even']:.8f}  D_odd={D_eig['D_odd']:.8f}")
        print(f"  Green-Kubo:  D_even={D_gk['D_even']:.8f}  D_odd={D_gk['D_odd']:.8f}")
        print(f"  exact:       D_even={D_ex['D_even']:.8f}  D_odd={D_ex['D_odd']:.8f}")
        print(f"  exact-GK:    ΔD_even={abs(D_ex['D_even']-D_gk['D_even']):.2e}  "
              f"ΔD_odd={abs(D_ex['D_odd']-D_gk['D_odd']):.2e}")


if __name__ == "__main__":
    run_self_checks()
