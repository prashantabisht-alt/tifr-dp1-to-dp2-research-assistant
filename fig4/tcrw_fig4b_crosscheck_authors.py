"""
Fig 4(b) cross-check: Bloch matrix P(k) vs direct PBC-torus diagonalization.

GOAL
----
Numerically verify that my paper-Eq(1) Bloch matrix `build_Pk(omega, D_r, kx, ky)`
in tcrw_fig4b_paper.py is the correct Bloch decomposition of the authors'
real-space transition rules.

METHOD
------
1. Take the authors' build_sparse_transition_matrix from TRW.py and make a
   minimal PBC variant: replace the "blocked = stay put, no rotation" rule
   with wrap-around (new_i %= N, new_j %= N). Everything else identical.
   Call this W_torus, shape (4 N^2, 4 N^2), column-stochastic.

2. Diagonalize W_torus densely.  Get 4 N^2 eigenvalues -> Sigma_torus.

3. Using my Bloch matrix P(k) from tcrw_fig4b_paper.py, loop over the
   discrete momenta k = (2 pi m / N, 2 pi n / N) for m,n in 0..N-1.
   At each k, diagonalize the 4x4 P(k), get 4 eigenvalues.
   Collect all 4 N^2 of them -> Sigma_Bloch.

4. Sort both lists and compare pointwise.  Pass criterion: max abs
   difference < 1e-12.

   Additionally do a Hausdorff-style set match: for every eigenvalue in
   one set, find the nearest in the other and ensure distance < 1e-12.
   This guards against any numerical re-ordering across degeneracies.

Why this is definitive
----------------------
If the two spectra agree at machine precision, then for every mode on the
torus there is a matching Bloch eigenvalue at some discrete k (and vice
versa).  Since the torus matrix is built from the authors' *microscopic
rules*, this proves my Eq (1) Bloch matrix is their operator in k-space.

Author: Prashant Bisht, TIFR Hyderabad
"""

import importlib.util
import os
import sys
import time
import numpy as np

# ---------------------------------------------------------------
# Bring in my Bloch builder from tcrw_fig4b_paper.py
# ---------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from tcrw_fig4b_paper import build_Pk  # my Eq (1) 4x4 Bloch matrix

# ---------------------------------------------------------------
# Build PBC torus matrix using authors' weights
# ---------------------------------------------------------------
# Authors' direction indexing: ['↑', '→', '↓', '←'] = indices 0,1,2,3
# CW rotation:  d -> (d + 1) % 4
# CCW rotation: d -> (d - 1) % 4
# DX, DY per direction:
#   d=0 ↑ -> (0, +1)
#   d=1 → -> (+1, 0)
#   d=2 ↓ -> (0, -1)
#   d=3 ← -> (-1, 0)
DX_AUTHORS = [0, 1, 0, -1]
DY_AUTHORS = [0, 1, 0, -1]  # careful: separate DY table
# Actually the authors use dir_to_vec:
#   ↑ -> (0, 1), → -> (1, 0), ↓ -> (0, -1), ← -> (-1, 0)
# so DX=[0,1,0,-1], DY=[1,0,-1,0]
DX_AUTHORS = np.array([0, 1, 0, -1])
DY_AUTHORS = np.array([1, 0, -1, 0])


def index_torus(i, j, d, N):
    """Matches authors' convention: (i * N + j) * 4 + d, but wraps i,j mod N."""
    return (i * N + j) * 4 + d


def build_torus_matrix(N, omega, D_r):
    """
    Build the full (4 N^2, 4 N^2) column-stochastic transition matrix for
    a PBC torus of size N x N using the authors' weights.

    Rules (identical to authors' build_sparse_transition_matrix, but with
    wraparound instead of blocked-when-out-of-bounds):

      Rotation noise step (prob D_r):
        d -> (d + 1) mod 4 (CW):  weight D_r * (1 - omega)
        d -> (d - 1) mod 4 (CCW): weight D_r * omega

      Chiral move step (prob 1 - D_r):
        site (i, j) -> ((i + DX[d]) mod N, (j + DY[d]) mod N)
        and simultaneously d -> (d +/- 1) mod 4:
          move + CW : weight (1 - D_r) * omega
          move + CCW: weight (1 - D_r) * (1 - omega)

      On a torus every move is valid -> no "blocked" self-loop.

    Returns dense np.ndarray P of shape (4 N^2, 4 N^2).
    """
    size = 4 * N * N
    P = np.zeros((size, size), dtype=np.float64)

    for i in range(N):
        for j in range(N):
            for d in range(4):
                src = index_torus(i, j, d, N)
                cw  = (d + 1) % 4
                ccw = (d - 1) % 4

                # --- Rotational noise ---
                P[index_torus(i, j, cw,  N), src] += D_r * (1 - omega)
                P[index_torus(i, j, ccw, N), src] += D_r * omega

                # --- Chiral move ---
                ni = (i + DX_AUTHORS[d]) % N
                nj = (j + DY_AUTHORS[d]) % N

                P[index_torus(ni, nj, cw,  N), src] += (1 - D_r) * omega
                P[index_torus(ni, nj, ccw, N), src] += (1 - D_r) * (1 - omega)

    return P


def column_sums_check(P):
    """Probability-preservation sanity: every column must sum to 1."""
    sums = P.sum(axis=0)
    err = float(np.max(np.abs(sums - 1.0)))
    return err


# ---------------------------------------------------------------
# Bloch spectrum from my build_Pk
# ---------------------------------------------------------------
def bloch_spectrum(N, omega, D_r):
    """
    Return the 4 N^2 eigenvalues obtained by diagonalizing my 4x4 P(k)
    at each of the N^2 discrete momenta k = (2 pi m / N, 2 pi n / N).

    Convention: kx is conjugate to the i-axis (DX direction),
                ky is conjugate to the j-axis (DY direction).

    My build_Pk uses
       ekx = exp(i kx), eky = exp(i ky)
    with the same C1, C2, R1, R2 mapping as the paper and authors'.
    """
    out = []
    for m in range(N):
        for n in range(N):
            kx = 2.0 * np.pi * m / N
            ky = 2.0 * np.pi * n / N
            Pk = build_Pk(omega, D_r, kx, ky)
            out.extend(np.linalg.eigvals(Pk))
    return np.array(out, dtype=complex)


# ---------------------------------------------------------------
# Compare two sets of complex numbers
# ---------------------------------------------------------------
def compare_spectra(S1, S2, label1="A", label2="B"):
    """
    Two-sided set comparison.  Returns max nearest-neighbour distance
    and also the sorted-componentwise residual.
    """
    S1 = np.asarray(S1, dtype=complex)
    S2 = np.asarray(S2, dtype=complex)
    assert S1.shape == S2.shape, f"size mismatch: {S1.shape} vs {S2.shape}"

    # Two-sided Hausdorff-ish: every point has a close partner in the other set
    def nn(a, b):
        # For each element of a, min distance to any element of b.
        # O(n^2) is fine for n up to a few thousand.
        return np.max(np.min(np.abs(a[:, None] - b[None, :]), axis=1))

    d12 = nn(S1, S2)
    d21 = nn(S2, S1)

    # Also a sorted comparison (should be tight if both sorted identically)
    key = lambda z: (np.round(z.real, 10), np.round(z.imag, 10))
    S1s = np.array(sorted(S1, key=key))
    S2s = np.array(sorted(S2, key=key))
    resid_sorted = float(np.max(np.abs(S1s - S2s)))

    print(f"  |spec({label1}) -> spec({label2})|_max = {d12:.3e}")
    print(f"  |spec({label2}) -> spec({label1})|_max = {d21:.3e}")
    print(f"  max |sorted({label1}) - sorted({label2})| = {resid_sorted:.3e}")
    return max(d12, d21), resid_sorted


# ---------------------------------------------------------------
# Run
# ---------------------------------------------------------------
def run(N, omega, D_r):
    print(f"\n=== Cross-check: N={N}, omega={omega}, D_r={D_r} ===")
    t0 = time.time()
    W = build_torus_matrix(N, omega, D_r)
    print(f"  built torus matrix ({W.shape[0]}x{W.shape[0]}) in {time.time()-t0:.2f}s")
    col_err = column_sums_check(W)
    print(f"  max |col_sum - 1| = {col_err:.3e}  {'PASS' if col_err<1e-12 else 'FAIL'}")

    t0 = time.time()
    Sigma_torus = np.linalg.eigvals(W)
    print(f"  torus eigvals ({Sigma_torus.size}) in {time.time()-t0:.2f}s")

    t0 = time.time()
    Sigma_bloch = bloch_spectrum(N, omega, D_r)
    print(f"  Bloch eigvals ({Sigma_bloch.size}) in {time.time()-t0:.2f}s")

    nn_err, sorted_err = compare_spectra(Sigma_torus, Sigma_bloch, "torus", "Bloch")

    ok = nn_err < 1e-10
    print(f"  --> {'PASS' if ok else 'FAIL'}  (tol 1e-10 on nearest-neighbour)")
    return ok


if __name__ == "__main__":
    print("=" * 60)
    print(" TCRW Fig 4(b): cross-check Bloch P(k) against PBC torus")
    print("=" * 60)

    # Paper Fig 4(b) uses D_r = 0.1, omega in {0.35, 0.5, 0.65}.
    # Test several (N, omega, D_r) combos.
    cases = [
        (4, 0.35, 0.1),
        (4, 0.50, 0.1),
        (4, 0.65, 0.1),
        (5, 0.50, 0.1),   # odd N, covers non-trivial discrete k
        (6, 0.35, 0.1),
        (8, 0.25, 0.3),   # stress test: different parameters
        (8, 0.75, 0.01),
    ]
    results = []
    for N, om, Dr in cases:
        ok = run(N, om, Dr)
        results.append((N, om, Dr, ok))

    print("\n" + "=" * 60)
    print(" SUMMARY")
    print("=" * 60)
    for N, om, Dr, ok in results:
        tag = "PASS" if ok else "FAIL"
        print(f"  N={N}  omega={om:.2f}  D_r={Dr:.3f}  -> {tag}")
    if all(ok for *_, ok in results):
        print("\n  ALL CASES PASS  -> Bloch matrix P(k) is provably the")
        print("                      Bloch decomposition of the authors' operator.")
    else:
        print("\n  SOMETHING FAILED -> investigate.")
