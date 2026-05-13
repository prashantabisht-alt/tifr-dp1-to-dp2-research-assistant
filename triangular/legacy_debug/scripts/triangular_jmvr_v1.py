"""
Day-1 starter: triangular JMVR 6×6 Bloch matrix.

Goal of this file:
  1. Build the corrected 6x6 generator M(k) for JMVR-style chirality
  2. Keep Dipanjan's old c3 sign only as a comparison in
     triangular_jmvr_corrected.py
  3. Diagonalize across the Γ-M-K-Γ path of the triangular BZ
  4. Plot all 6 bands at ε=0 (achiral) and ε=0.15 (Dipanjan's value)

Convention used (matches Confinement 2021 draft §II + Dipanjan's rtp_tl_2.nb):
  - 6 directors, three "axes" : ê₁ = (2a, 0), ê₂ = (a, b), ê₃ = (a, -b)
                                ê₄ = -ê₁,    ê₅ = -ê₂,    ê₆ = -ê₃
  - Per-state translation rate: 1/6 to each of 6 NN, plus ±ε bias along director's axis
  - Per-state rotation rate: γ/2 to each of the two adjacent directors (d ± 1 mod 6)
  - Total outgoing rate per state: 1 + γ

Dipanjan's coeff_i (rtp_tl_2.nb), for reference — Eq. (5) of the Confinement draft:

    coeff_i(k1, k2) = (1/3) [cos(2 a k1) + cos(a k1 + b k2) + cos(a k1 - b k2)]
                      - 1 - γ
                      + 2 i ε sin(k · ê_i)

The corrected calculation flips the third sign:

    coeff_3 = bulk - 1 - gamma - 2 i epsilon sin(a k1 - b k2)

where the sin argument for i = 1, 2 is k · (2a, 0), k · (a, b), and the
state-2/c3 axis is fixed from the Appendix-B master equation.
Conjugates fill in i = 4, 5, 6 (the opposite directors).

Off-diagonal entries: γ/2 at positions (i, i+1 mod 6) and (i, i-1 mod 6).

Author: Prashant Bisht, TIFR Hyderabad — Day 1, 2026-05-12
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from triangular_jmvr_corrected import build_Mk_corrected


# ---------------------------------------------------------------------------
# 1.  6×6 Bloch matrix M(k)
# ---------------------------------------------------------------------------
def build_Mk(gamma: float, epsilon: float, k1: float, k2: float,
             a: float = 1.0, b: float = 1.0) -> np.ndarray:
    """
    Corrected triangular JMVR 6x6 generator at Bloch momentum (k1, k2).

    Dipanjan's original c3 sign is available only through
    triangular_jmvr_corrected.build_Mk_dipanjan.
    """
    return build_Mk_corrected(gamma, epsilon, k1, k2, a, b)


# ---------------------------------------------------------------------------
# 2.  Sanity checks — RUN THESE FIRST before any plotting
# ---------------------------------------------------------------------------
def sanity_check(gamma: float = 0.01, epsilon: float = 0.15):
    """
    Three sanity tests that MUST pass before trusting the matrix.
    """
    print("=" * 60)
    print("Sanity check at γ = {}, ε = {}".format(gamma, epsilon))
    print("=" * 60)

    # Test 1: At k=0, the generator should have column sums = 0
    # (probability conservation for a CTRW).
    M0 = build_Mk(gamma, epsilon, 0.0, 0.0)
    col_sums = M0.sum(axis=0)
    max_col = np.max(np.abs(col_sums))
    print(f"  Test 1 — column sums at k=0: max|col sum| = {max_col:.2e}")
    print(f"          (should be 0 to machine precision)")
    assert max_col < 1e-12, "Column sums not zero — kernel violates probability conservation"

    # Test 2: At k=0, one eigenvalue should be exactly 0 (Perron-like, stationary state)
    # and the others should have negative real parts.
    evals0 = np.linalg.eigvals(M0)
    sorted_evals = sorted(evals0, key=lambda x: -x.real)
    print(f"  Test 2 — eigenvalues at k=0 (sorted by Re):")
    for i, e in enumerate(sorted_evals):
        print(f"          λ_{i} = {e:.6f}")
    print(f"          (λ_0 should be 0; others should have Re < 0)")
    assert abs(sorted_evals[0]) < 1e-10, "No zero eigenvalue at k=0"

    # Test 3: Term-by-term match for an unaffected coefficient at random (k1, k2).
    # The corrected-vs-Dipanjan difference is only in c3.
    k1_test, k2_test = 0.7, -1.3
    M = build_Mk(gamma, epsilon, k1_test, k2_test)
    expected_coeff_1 = (
        (1.0 / 3.0) * (
            np.cos(2 * k1_test)
            + np.cos(k1_test + k2_test)
            + np.cos(k1_test - k2_test)
        )
        - 1.0 - gamma
        + 2j * epsilon * np.sin(2 * k1_test)
    )
    diff = abs(M[0, 0] - expected_coeff_1)
    print(f"  Test 3 — Dipanjan coeff_1 match at (k1, k2) = ({k1_test}, {k2_test}):")
    print(f"          |M[0,0] - expected| = {diff:.2e}")
    print(f"          (should be < 1e-14)")
    assert diff < 1e-14, "Coefficient does not match Dipanjan's formula"

    print()
    print("  ALL SANITY CHECKS PASSED ✓")
    print()


# ---------------------------------------------------------------------------
# 3.  Band structure along Γ-M-K-Γ path of the triangular BZ
# ---------------------------------------------------------------------------
def bz_path(n_per_segment: int = 60):
    """
    High-symmetry path in the triangular BZ, in (k1, k2) coordinates
    appropriate for the (2a, 0), (a, b) lattice convention with a = b = 1.

    Approximate high-symmetry points:
      Γ = (0, 0)
      M = (π/2, 0)               (midpoint of first BZ along k1)
      K = (π/3, π/(3√3))         (corner of BZ — symbolic; we use a placeholder)
      back to Γ

    Note: these k-coords are convention-dependent. Tomorrow we can refine the
    exact K point; for Day 1 a generic Γ-M-?-Γ path showing the gap is enough.
    """
    Gamma = np.array([0.0, 0.0])
    M = np.array([np.pi / 2.0, 0.0])
    K = np.array([np.pi / 3.0, np.pi / (3 * np.sqrt(3))])
    legs = [(Gamma, M, "Γ-M"), (M, K, "M-K"), (K, Gamma, "K-Γ")]
    k_path = []
    seg_labels = [0]
    for (start, end, _) in legs:
        for i in range(n_per_segment):
            t = i / n_per_segment
            k_path.append(start + t * (end - start))
        seg_labels.append(len(k_path))
    k_path.append(Gamma)
    return np.array(k_path), seg_labels


def band_structure(gamma: float, epsilon: float, n_per_segment: int = 60):
    """Diagonalize M(k) along the BZ path; return eigenvalues (sorted)."""
    k_path, seg_labels = bz_path(n_per_segment)
    n_k = len(k_path)
    evals = np.zeros((n_k, 6), dtype=complex)
    for i, (k1, k2) in enumerate(k_path):
        M = build_Mk(gamma, epsilon, k1, k2)
        w = np.linalg.eigvals(M)
        # Sort by real part (descending — λ=0 first, then negatives)
        evals[i] = sorted(w, key=lambda x: -x.real)
    return k_path, seg_labels, evals


def plot_bands(gamma: float = 0.01, eps_list=(0.0, 0.15)):
    """Two-row plot: Re(λ) and Im(λ) along the BZ path, one column per ε."""
    n_eps = len(eps_list)
    fig, axes = plt.subplots(2, n_eps, figsize=(5 * n_eps, 7), sharex=True)
    if n_eps == 1:
        axes = axes[:, None]

    for col, eps in enumerate(eps_list):
        k_path, seg_labels, evals = band_structure(gamma, eps)
        x = np.arange(len(k_path))

        # Top row: Re(λ)
        ax_re = axes[0, col]
        for n in range(6):
            ax_re.plot(x, evals[:, n].real, "-", lw=1.2, label=f"band {n}")
        ax_re.set_ylabel(r"Re(λ)")
        ax_re.set_title(rf"γ = {gamma}, ε = {eps}")
        ax_re.axhline(0, color="k", lw=0.5, alpha=0.3)
        for s in seg_labels:
            ax_re.axvline(s, color="gray", lw=0.5, alpha=0.5)
        ax_re.set_xticks(seg_labels)
        ax_re.set_xticklabels(["Γ", "M", "K", "Γ"])
        ax_re.grid(True, alpha=0.2)

        # Bottom row: Im(λ)
        ax_im = axes[1, col]
        for n in range(6):
            ax_im.plot(x, evals[:, n].imag, "-", lw=1.2)
        ax_im.set_ylabel(r"Im(λ)")
        ax_im.set_xlabel("k path")
        ax_im.axhline(0, color="k", lw=0.5, alpha=0.3)
        for s in seg_labels:
            ax_im.axvline(s, color="gray", lw=0.5, alpha=0.5)
        ax_im.set_xticks(seg_labels)
        ax_im.set_xticklabels(["Γ", "M", "K", "Γ"])
        ax_im.grid(True, alpha=0.2)

    axes[0, 0].legend(fontsize=8, loc="best")
    fig.suptitle("Triangular JMVR 6×6 generator — bands along BZ path",
                 fontsize=12)
    plt.tight_layout()
    out = "triangular_jmvr_bands_day1.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# 4.  Main — run sanity checks then plot bands
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Triangular JMVR — Day 1 starter\n")
    sanity_check(gamma=0.01, epsilon=0.15)
    plot_bands(gamma=0.01, eps_list=(0.0, 0.15))
    print("\nDay 1 deliverable produced: triangular_jmvr_bands_day1.png")
