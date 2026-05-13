"""
Canonical triangular JMVR generator after the c3 sign diagnosis.

This file is the source of truth for future triangular JMVR calculations.
Use build_Mk_corrected(...) by default.  build_Mk_dipanjan(...) is kept only
to reproduce the old notebook/draft sign bug for comparisons.

Model:
    - 6 internal director states on a triangular lattice.
    - Hopping to all 6 nearest neighbours with baseline rate 1/6.
    - Translation bias epsilon along the current director axis.
    - Director switching d -> d +/- 1 at rate gamma/2.
"""
from __future__ import annotations

from itertools import permutations
from pathlib import Path

import numpy as np


L_DEFAULT = 30
T_DEFAULT = 50.0


def _validate_rates(epsilon: float) -> None:
    if abs(epsilon) > 1.0 / 6.0:
        raise ValueError("epsilon must satisfy |epsilon| <= 1/6")


def _bulk_factor(k1: float, k2: float, a: float, b: float) -> float:
    return (1.0 / 3.0) * (
        np.cos(2.0 * a * k1)
        + np.cos(a * k1 + b * k2)
        + np.cos(a * k1 - b * k2)
    )


def _build_Mk(
    gamma: float,
    epsilon: float,
    k1: float,
    k2: float,
    *,
    a: float = 1.0,
    b: float = 1.0,
    c3_sign: float,
) -> np.ndarray:
    """Build the 6x6 continuous-time Bloch generator M(k)."""
    _validate_rates(epsilon)
    bulk = _bulk_factor(k1, k2, a, b)

    c1 = bulk - 1.0 - gamma + 2.0j * epsilon * np.sin(2.0 * a * k1)
    c2 = bulk - 1.0 - gamma + 2.0j * epsilon * np.sin(a * k1 + b * k2)
    c3 = bulk - 1.0 - gamma + c3_sign * 2.0j * epsilon * np.sin(
        a * k1 - b * k2
    )

    diag = [c1, c2, c3, np.conj(c1), np.conj(c2), np.conj(c3)]
    M = np.zeros((6, 6), dtype=complex)
    for src, value in enumerate(diag):
        M[src, src] = value
        M[(src + 1) % 6, src] += gamma / 2.0
        M[(src - 1) % 6, src] += gamma / 2.0
    return M


def build_Mk_corrected(
    gamma: float,
    epsilon: float,
    k1: float,
    k2: float,
    a: float = 1.0,
    b: float = 1.0,
) -> np.ndarray:
    """
    Corrected triangular JMVR generator.

    The state-2/state-c3 chirality term is

        c3 = B - 1 - gamma - 2 i epsilon sin(a k1 - b k2).
    """
    return _build_Mk(gamma, epsilon, k1, k2, a=a, b=b, c3_sign=-1.0)


def build_Mk_dipanjan(
    gamma: float,
    epsilon: float,
    k1: float,
    k2: float,
    a: float = 1.0,
    b: float = 1.0,
) -> np.ndarray:
    """
    Original Dipanjan/draft generator with the old c3 sign.

    This is intentionally not the default.  Use only for reproducing the
    old theory-vs-KMC mismatch.
    """
    return _build_Mk(gamma, epsilon, k1, k2, a=a, b=b, c3_sign=+1.0)


def build_Mk(
    gamma: float,
    epsilon: float,
    k1: float,
    k2: float,
    a: float = 1.0,
    b: float = 1.0,
) -> np.ndarray:
    """Default public builder: the corrected matrix."""
    return build_Mk_corrected(gamma, epsilon, k1, k2, a=a, b=b)


def sorted_eigenvalues(M: np.ndarray) -> np.ndarray:
    vals = np.linalg.eigvals(M)
    return vals[np.lexsort((vals.imag, vals.real))]


def best_spectrum_distance(vals_a: np.ndarray, vals_b: np.ndarray) -> float:
    """Minimum max-pair distance between two 6-eigenvalue spectra."""
    best = np.inf
    for perm in permutations(range(len(vals_b))):
        dist = float(np.max(np.abs(vals_a - vals_b[list(perm)])))
        if dist < best:
            best = dist
    return best


def k_from_q(q1: float, q2: float, a: float = 1.0, b: float = 1.0) -> tuple[float, float]:
    """
    Convert lattice reciprocal coordinates q to Dipanjan's (k1,k2).

        q1 = 2 a k1
        q2 = a k1 + b k2
    """
    k1 = q1 / (2.0 * a)
    k2 = (q2 - 0.5 * q1) / b
    return k1, k2


def sixfold_symmetry_error(
    builder,
    gamma: float = 0.01,
    epsilon: float = 0.15,
    *,
    a: float = 1.0,
    b: float = 1.0,
    q_points: tuple[tuple[float, float], ...] | None = None,
) -> float:
    """
    Check the 60-degree spectral symmetry in q coordinates.

    The rotation used here is (q1, q2) -> (q1 - q2, q1).
    Returns the largest spectrum mismatch over the test points.
    """
    if q_points is None:
        q_points = (
            (0.41, 1.17),
            (1.30, -0.73),
            (2.20, 0.51),
            (-1.10, 2.40),
        )

    max_err = 0.0
    for q1, q2 in q_points:
        k1, k2 = k_from_q(q1, q2, a=a, b=b)
        rk1, rk2 = k_from_q(q1 - q2, q1, a=a, b=b)
        vals = np.linalg.eigvals(builder(gamma, epsilon, k1, k2, a, b))
        vals_rot = np.linalg.eigvals(builder(gamma, epsilon, rk1, rk2, a, b))
        max_err = max(max_err, best_spectrum_distance(vals, vals_rot))
    return float(max_err)


def probability_conservation_error(
    builder=build_Mk_corrected,
    gamma: float = 0.01,
    epsilon: float = 0.15,
    *,
    a: float = 1.0,
    b: float = 1.0,
) -> float:
    """Return max absolute column sum at k=0."""
    M0 = builder(gamma, epsilon, 0.0, 0.0, a, b)
    return float(np.max(np.abs(M0.sum(axis=0))))


def zero_eigenvalue_error(
    builder=build_Mk_corrected,
    gamma: float = 0.01,
    epsilon: float = 0.15,
    *,
    a: float = 1.0,
    b: float = 1.0,
) -> float:
    """Return distance of the closest k=0 eigenvalue to zero."""
    M0 = builder(gamma, epsilon, 0.0, 0.0, a, b)
    return float(np.min(np.abs(np.linalg.eigvals(M0))))


def load_kmc_counts(path: Path, L: int = L_DEFAULT) -> tuple[np.ndarray, int]:
    """Load n1 n2 count rows into counts[n2,n1]."""
    counts = np.zeros((L, L), dtype=np.int64)
    n_walkers = None
    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                if "N=" in line:
                    try:
                        n_walkers = int(line.split("N=")[1].strip())
                    except ValueError:
                        n_walkers = None
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            n1, n2, count = int(parts[0]), int(parts[1]), int(parts[2])
            counts[n2, n1] = count
    if n_walkers is None:
        n_walkers = int(counts.sum())
    return counts, n_walkers


def theory_probability_on_lattice(
    builder,
    gamma: float = 0.01,
    epsilon: float = 0.15,
    *,
    L: int = L_DEFAULT,
    t_final: float = T_DEFAULT,
    a: float = 1.0,
    b: float = 1.0,
) -> np.ndarray:
    """Exact P[n2,n1,t] by finite Fourier inversion on the LxL torus."""
    eig_data: list[tuple[float, float, complex]] = []
    initial = (1.0 / 6.0) * np.ones(6, dtype=complex)

    for m1 in range(L):
        for m2 in range(L):
            k1 = np.pi * m1 / (a * L)
            k2 = np.pi * (2.0 * m2 - m1) / (b * L)
            M = builder(gamma, epsilon, k1, k2, a, b)
            evals, evecs = np.linalg.eig(M)
            prefactor = np.linalg.solve(evecs, initial)
            ptilde = np.sum(evecs @ (prefactor * np.exp(evals * t_final)))
            eig_data.append((k1, k2, ptilde))

    P = np.zeros((L, L), dtype=float)
    for n1 in range(L):
        for n2 in range(L):
            x = 2.0 * a * n1 + a * n2
            y = b * n2
            total = 0.0
            for k1, k2, ptilde in eig_data:
                total += np.real(ptilde * np.exp(-1j * (k1 * x + k2 * y)))
            P[n2, n1] = total / (L * L)
    return P


def compare_saved_kmc(
    counts_path: Path | None = None,
    gamma: float = 0.01,
    epsilon: float = 0.15,
    *,
    L: int = L_DEFAULT,
    t_final: float = T_DEFAULT,
) -> dict[str, float]:
    """Compare saved KMC counts against Dipanjan and corrected exact theory."""
    if counts_path is None:
        counts_path = Path(__file__).resolve().parent / "kmc_triangular_counts.txt"

    counts, n_walkers = load_kmc_counts(counts_path, L=L)
    P_kmc = counts / n_walkers
    P_bug = theory_probability_on_lattice(
        build_Mk_dipanjan, gamma, epsilon, L=L, t_final=t_final
    )
    P_fix = theory_probability_on_lattice(
        build_Mk_corrected, gamma, epsilon, L=L, t_final=t_final
    )

    return {
        "n_walkers": float(n_walkers),
        "rms_buggy": float(np.sqrt(np.mean((P_kmc - P_bug) ** 2))),
        "rms_corrected": float(np.sqrt(np.mean((P_kmc - P_fix) ** 2))),
        "mc_noise_floor": float(np.sqrt(P_kmc.max() / n_walkers)),
    }


def run_self_checks() -> None:
    """Run the lightweight checks that protect the corrected sign."""
    col_err = probability_conservation_error(build_Mk_corrected)
    zero_err = zero_eigenvalue_error(build_Mk_corrected)
    sym_fix = sixfold_symmetry_error(build_Mk_corrected)
    sym_bug = sixfold_symmetry_error(build_Mk_dipanjan)

    print("Triangular JMVR corrected self-checks")
    print(f"  column-sum error at k=0       = {col_err:.3e}")
    print(f"  zero-eigenvalue error at k=0  = {zero_err:.3e}")
    print(f"  sixfold error, corrected      = {sym_fix:.3e}")
    print(f"  sixfold error, Dipanjan       = {sym_bug:.3e}")

    assert col_err < 1e-12
    assert zero_err < 1e-12
    assert sym_fix < 1e-12
    assert sym_bug > 1e-6


if __name__ == "__main__":
    run_self_checks()
