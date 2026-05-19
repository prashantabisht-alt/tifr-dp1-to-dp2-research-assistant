"""TCRW (Osat et al. / JMVR) model on the TRIANGULAR lattice.

Self-contained module for the JMVR-style triangular chiral random walker,
repackaged with the same API as triangular_chiral_rtw.py and
square_pal/square_chiral_rtw.py for side-by-side comparison.

Model (translation chirality):
  - Position: triangular-lattice site (n₁, n₂) ∈ Z².
  - Director: m ∈ {0,1,...,5} → the six NN directions.
  - Rates (continuous-time):
      hop to NN j:  depends on j and current director m
      CCW tumble:   m → (m+1) mod 6,   rate γ/2
      CW tumble:    m → (m-1) mod 6,   rate γ/2

  The Bloch matrix diagonal for director m is:
      M_{m,m}(k) = B(k) - 1 - γ ± 2iε sin(phase_m(k))
  where B(k) = (1/3)(cos(2k₁) + cos(k₁+k₂) + cos(k₁-k₂)).

  The (k₁,k₂) are the JMVR parameterisation, related to axial Fourier
  coordinates (q₁,q₂) by:  q₁ = 2k₁,  q₂ = k₁ + k₂  (with a=b=1).

KEY DIFFERENCE from our M-P style model:
  - JMVR: walker hops to ALL 6 neighbours (with bias), rotation symmetric
  - M-P:  walker hops ONLY along director m, rotation asymmetric

  JMVR chirality is in the translation–rotation coupling (ε),
  not in asymmetric rotation.  With symmetric rotation (γ/2 each way),
  D_odd = 0 in the bulk infinite-space limit.

Matrix convention:  dP/dt = M P  (columns = source states).
"""

from __future__ import annotations
import numpy as np


N_DIR = 6


def triangular_deltas() -> np.ndarray:
    """Return the six NN steps in AXIAL coordinates."""
    return np.array(
        [
            [+1,  0],   # m=0
            [ 0, +1],   # m=1
            [-1, +1],   # m=2
            [-1,  0],   # m=3
            [ 0, -1],   # m=4
            [+1, -1],   # m=5
        ],
        dtype=float,
    )


def axial_to_cartesian_jacobian() -> np.ndarray:
    """Jacobian J such that r_cart = J @ r_axial.

    J = [[1, 1/2], [0, √3/2]]  (columns are a₁, a₂).
    """
    sqrt3 = np.sqrt(3.0)
    return np.array([[1.0, 0.5], [0.0, 0.5 * sqrt3]])


def _jmvr_k_to_axial_jacobian() -> np.ndarray:
    """Jacobian J_kq = dk/dq for the JMVR k ↔ axial q conversion.

    q₁ = 2k₁,  q₂ = k₁ + k₂   (with a=b=1)

    So dq/dk = [[2, 0], [1, 1]], and dk/dq = [[1/2, 0], [-1/2, 1]].
    """
    return np.array([[0.5, 0.0], [-0.5, 1.0]])


def _composite_jacobian() -> np.ndarray:
    """The product J @ J_kq^T maps from k-space D to Cartesian D:

        D_cart = S @ D_k @ S^T   where S = J @ J_kq^T

    For a=b=1 with our conventions, S = diag(1/2, √3/2) — just a
    diagonal rescaling.
    """
    J = axial_to_cartesian_jacobian()
    J_kq = _jmvr_k_to_axial_jacobian()
    return J @ J_kq.T


def build_Mk_tcrw(
    k1: float,
    k2: float,
    *,
    gamma: float = 0.01,
    epsilon: float = 0.15,
    a: float = 1.0,
    b: float = 1.0,
) -> np.ndarray:
    """Build the 6×6 Bloch generator M(k) for the JMVR/TCRW model.

    Parameters
    ----------
    k1, k2 : float
        Wavevector components in the JMVR convention.
    gamma : float
        Director rotation rate (symmetric: γ/2 each way).
    epsilon : float
        Translation chirality parameter.  |ε| ≤ 1/6.
    a, b : float
        Lattice parameters (default 1.0).
    """
    if abs(epsilon) > 1.0 / 6.0 + 1e-10:
        raise ValueError("|epsilon| must be <= 1/6")

    bulk = (1.0 / 3.0) * (
        np.cos(2.0 * a * k1)
        + np.cos(a * k1 + b * k2)
        + np.cos(a * k1 - b * k2)
    )

    # Corrected c₃ sign (negative)
    c1 = bulk - 1.0 - gamma + 2.0j * epsilon * np.sin(2.0 * a * k1)
    c2 = bulk - 1.0 - gamma + 2.0j * epsilon * np.sin(a * k1 + b * k2)
    c3 = bulk - 1.0 - gamma - 2.0j * epsilon * np.sin(a * k1 - b * k2)

    diag = [c1, c2, c3, np.conj(c1), np.conj(c2), np.conj(c3)]

    M = np.zeros((N_DIR, N_DIR), dtype=complex)
    for m in range(N_DIR):
        M[m, m] = diag[m]
        M[(m + 1) % N_DIR, m] += gamma / 2.0   # CCW tumble
        M[(m - 1) % N_DIR, m] += gamma / 2.0   # CW tumble

    return M


def _top_eig(k1, k2, gamma, epsilon):
    """Eigenvalue of M(k) closest to zero (the hydrodynamic mode)."""
    M = build_Mk_tcrw(k1, k2, gamma=gamma, epsilon=epsilon)
    vals = np.linalg.eigvals(M)
    return vals[np.argmin(np.abs(vals))]


def numerical_diffusion_tensor(
    *,
    gamma: float = 0.01,
    epsilon: float = 0.15,
    dk: float = 1e-5,
) -> dict:
    """Diffusion tensor from λ₀(k) curvature (finite differences).

    Computes d²λ/dk_i dk_j in JMVR (k₁,k₂) space, then converts:
      1. k → axial:  D_axial = J_kq^T D_k J_kq
      2. axial → Cartesian:  D_cart = J D_axial J^T

    Equivalently  D_cart = S D_k S^T  where S = J J_kq^T = diag(1/2, √3/2).
    """
    l0 = _top_eig(0, 0, gamma, epsilon)

    d2_11 = (_top_eig(dk, 0, gamma, epsilon)
             + _top_eig(-dk, 0, gamma, epsilon) - 2 * l0) / dk**2
    d2_22 = (_top_eig(0, dk, gamma, epsilon)
             + _top_eig(0, -dk, gamma, epsilon) - 2 * l0) / dk**2
    d2_12 = (_top_eig(dk, dk, gamma, epsilon)
             - _top_eig(dk, -dk, gamma, epsilon)
             - _top_eig(-dk, dk, gamma, epsilon)
             + _top_eig(-dk, -dk, gamma, epsilon)) / (4 * dk**2)

    D_k = -0.5 * np.array([[d2_11, d2_12], [d2_12, d2_22]], dtype=complex)

    # Convert to axial, then to Cartesian
    J_kq = _jmvr_k_to_axial_jacobian()
    D_axial = J_kq.T @ D_k @ J_kq

    J = axial_to_cartesian_jacobian()
    D_cart = J @ D_axial @ J.T

    D_even = 0.5 * np.real(D_cart[0, 0] + D_cart[1, 1])
    D_odd = 0.5 * np.real(D_cart[0, 1] - D_cart[1, 0])

    return {
        "D_cartesian": D_cart,
        "D_axial": D_axial,
        "D_k": D_k,
        "D_even": float(D_even),
        "D_odd": float(D_odd),
        "D_xx": float(np.real(D_cart[0, 0])),
        "D_yy": float(np.real(D_cart[1, 1])),
        "D_xy": float(np.real(D_cart[0, 1])),
        "D_yx": float(np.real(D_cart[1, 0])),
    }


def diffusion_tensor_green_kubo(
    *,
    gamma: float = 0.01,
    epsilon: float = 0.15,
) -> dict:
    """Full diffusion tensor via Green-Kubo in the JMVR k-space.

    The JMVR model has the walker hop to all 6 neighbours, so the
    velocity v_i(m) and the jump-noise second moment are different from
    the M-P model (where only the director neighbour is visited).

    Jump noise:  D^{jump}_{k,ij} = -(1/(2N)) Σ_m d²M_{mm}/dk_i dk_j |_{k=0}
      This is the rate-weighted second moment of individual hops.

    Persistent part: D^{pers}_{k,ij} = (1/N) v_i^T (-R⁺) v_j
      where v_i(m) = -Im(dM_{mm}/dk_i |_{k=0}) is the mean displacement
      rate in the k_i direction for director state m.

    Then D_cart = S D_k S^T with S = J J_kq^T.
    """
    dk = 1e-7

    # R = M(k=0) — the pure rotation generator
    R = np.real(build_Mk_tcrw(0, 0, gamma=gamma, epsilon=epsilon))
    R_pinv = np.linalg.pinv(R)

    # ---- Jump noise: diagonal second derivatives of M ----
    # At k=0, d²M_{mm}/dk_i dk_j is the same for all m when ε=0,
    # and varies with m when ε≠0 (but the ε-dependent part vanishes
    # at k=0, so jump noise is ε-independent).
    #
    # We compute this numerically from the diagonal of M(k).
    M00 = build_Mk_tcrw(0, 0, gamma=gamma, epsilon=epsilon)

    def diag_at(kk1, kk2):
        return np.diag(build_Mk_tcrw(kk1, kk2, gamma=gamma, epsilon=epsilon))

    d0 = np.diag(M00)

    # d²M_{mm}/dk₁²
    d2_diag_11 = (diag_at(dk, 0) + diag_at(-dk, 0) - 2 * d0) / dk**2
    # d²M_{mm}/dk₂²
    d2_diag_22 = (diag_at(0, dk) + diag_at(0, -dk) - 2 * d0) / dk**2
    # d²M_{mm}/dk₁dk₂
    d2_diag_12 = (diag_at(dk, dk) - diag_at(dk, -dk)
                  - diag_at(-dk, dk) + diag_at(-dk, -dk)) / (4 * dk**2)

    # Jump noise: average the diagonal second derivatives
    D_k_jump = np.zeros((2, 2))
    D_k_jump[0, 0] = -0.5 * np.real(np.mean(d2_diag_11))
    D_k_jump[1, 1] = -0.5 * np.real(np.mean(d2_diag_22))
    D_k_jump[0, 1] = -0.5 * np.real(np.mean(d2_diag_12))
    D_k_jump[1, 0] = D_k_jump[0, 1]

    # ---- Persistent part: velocity vectors from diagonal first derivatives ----
    dM_dk1 = (diag_at(dk, 0) - diag_at(-dk, 0)) / (2 * dk)
    dM_dk2 = (diag_at(0, dk) - diag_at(0, -dk)) / (2 * dk)

    # Mean velocity in k_i direction for each director state
    # v_i(m) = Im(dM_{mm}/dk_i) / 1  (since dM/dk is pure imaginary at k=0
    # for the translation terms)
    v_k1 = np.imag(dM_dk1)  # shape (6,)
    v_k2 = np.imag(dM_dk2)

    # Persistent part: D^{pers}_{k,ij} = (1/N) v_i^T (-R⁺) v_j
    neg_R_pinv = -R_pinv
    v_vecs = [v_k1, v_k2]
    D_k_pers = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            D_k_pers[i, j] = (1.0 / N_DIR) * v_vecs[i] @ neg_R_pinv @ v_vecs[j]

    D_k = D_k_jump + D_k_pers

    # ---- Convert k → axial → Cartesian ----
    J_kq = _jmvr_k_to_axial_jacobian()
    J = axial_to_cartesian_jacobian()

    D_axial = J_kq.T @ D_k @ J_kq
    D_cart = J @ D_axial @ J.T

    D_axial_jump = J_kq.T @ D_k_jump @ J_kq
    D_axial_pers = J_kq.T @ D_k_pers @ J_kq
    D_cart_jump = J @ D_axial_jump @ J.T
    D_cart_pers = J @ D_axial_pers @ J.T

    D_even = 0.5 * (np.real(D_cart[0, 0]) + np.real(D_cart[1, 1]))
    D_odd = 0.5 * (np.real(D_cart[0, 1]) - np.real(D_cart[1, 0]))

    return {
        "D_cartesian": np.real(D_cart),
        "D_persistent": np.real(D_cart_pers),
        "D_jump": np.real(D_cart_jump),
        "D_axial": D_axial,
        "D_k": D_k,
        "D_k_jump": D_k_jump,
        "D_k_persistent": D_k_pers,
        "D_even": float(D_even),
        "D_odd": float(D_odd),
        "D_xx": float(np.real(D_cart[0, 0])),
        "D_yy": float(np.real(D_cart[1, 1])),
        "D_xy": float(np.real(D_cart[0, 1])),
        "D_yx": float(np.real(D_cart[1, 0])),
    }


def run_self_checks() -> None:
    """Quick sanity checks for the TCRW model."""
    gamma = 0.01
    epsilon = 0.15

    print("=" * 60)
    print("TCRW (JMVR) model on triangular lattice")
    print("=" * 60)

    M0 = build_Mk_tcrw(0, 0, gamma=gamma, epsilon=epsilon)
    print(f"\nM(0,0) for gamma={gamma}, epsilon={epsilon}:")
    print(np.array2string(np.real(M0), precision=4))

    col_sums = np.abs(M0.sum(axis=0))
    print(f"\nmax |col sum| at k=0 = {np.max(col_sums):.2e}")

    vals = np.linalg.eigvals(M0)
    print(f"min |eigenvalue| at k=0 = {np.min(np.abs(vals)):.2e}")

    # C₆ symmetry check
    from itertools import permutations as _perms  # noqa: F811
    k1_test, k2_test = 0.3, 0.7
    q1, q2 = 2 * k1_test, k1_test + k2_test
    q1r, q2r = q1 - q2, q1
    k1r, k2r = q1r / 2.0, q2r - q1r / 2.0

    M_orig = build_Mk_tcrw(k1_test, k2_test, gamma=gamma, epsilon=epsilon)
    M_rot = build_Mk_tcrw(k1r, k2r, gamma=gamma, epsilon=epsilon)
    vals_orig = np.sort(np.linalg.eigvals(M_orig).real)
    vals_rot = np.sort(np.linalg.eigvals(M_rot).real)
    sym_err = np.max(np.abs(vals_orig - vals_rot))
    print(f"\nC₆ symmetry error (spectral) = {sym_err:.2e}")

    print("\n" + "=" * 60)
    print("DIFFUSION TENSOR")
    print("=" * 60)

    for eps_test in [0.0, 0.05, 0.10, 0.15, 1.0 / 6.0]:
        D_eig = numerical_diffusion_tensor(gamma=gamma, epsilon=eps_test)
        D_gk = diffusion_tensor_green_kubo(gamma=gamma, epsilon=eps_test)

        print(f"\nepsilon = {eps_test:.4f}")
        print(f"  eigenvalue:  D_even={D_eig['D_even']:.8f}  D_odd={D_eig['D_odd']:.8f}")
        print(f"  Green-Kubo:  D_even={D_gk['D_even']:.8f}  D_odd={D_gk['D_odd']:.8f}")
        print(f"  eig-GK:      ΔD_even={abs(D_eig['D_even'] - D_gk['D_even']):.2e}  "
              f"ΔD_odd={abs(D_eig['D_odd'] - D_gk['D_odd']):.2e}")

    # ε=0 sanity: should give the simple random walk result
    print("\n--- ε=0 check (symmetric random walk) ---")
    D0 = diffusion_tensor_green_kubo(gamma=gamma, epsilon=0.0)
    print(f"  D_even = {D0['D_even']:.8f}")
    print(f"  D_odd  = {D0['D_odd']:.8f}")
    print(f"  D_jump (Cartesian):")
    print(f"    {D0['D_jump']}")
    print(f"  D_pers (Cartesian):")
    print(f"    {D0['D_persistent']}")
    print(f"  Expected D_even for simple RW on triangular = 1/4 = {0.25:.8f}")


if __name__ == "__main__":
    run_self_checks()
