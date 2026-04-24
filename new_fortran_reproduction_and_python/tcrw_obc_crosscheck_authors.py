"""
Head-to-head OBC cross-check: user's tcrw_obc.build_transition_matrix vs
authors' TRW.build_sparse_transition_matrix, plus a direct comparison of
the resulting P(X,Y) steady states and the J_Dr, J_omega current
decomposition.

Convention notes (important for interpreting the diff):
  - Authors index is site-major:  idx_A = (i*(L_A+1) + j)*4 + d_A
      with d_A ordering ['↑','→','↓','←'] = 0,1,2,3
      grid sites 0..L_A (so (L_A+1) sites per axis)
  - User index is direction-major: idx_U = d_U*L_U^2 + y*L_U + x
      with d_U ordering  ↑=0, →=1, ↓=2, ←=3  (same direction meaning)
      grid sites 0..L_U-1 (so L_U sites per axis)

For equivalent physics the user must set L_U = L_A + 1.
"""

import os, sys, numpy as np, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
# Work inside either macOS client or Linux sandbox
for candidate in (
    "/Users/prashantbisht/Documents/Claude/Projects/TIFR DP1 to DP2 Research Assistant",
    "/sessions/elegant-wizardly-einstein/mnt/TIFR DP1 to DP2 Research Assistant",
):
    if os.path.isdir(candidate):
        ROOT = candidate
        break
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "fortran_reproduction"))

# The authors' module has an unusual filename.  Load it as TRW.
import importlib.util
spec = importlib.util.spec_from_file_location(
    "TRW_authors",
    os.path.join(ROOT, "fortran_reproduction",
                 "TRW._original_code_by_paperauthors.py"))
TRW = importlib.util.module_from_spec(spec)
spec.loader.exec_module(TRW)

import tcrw_obc
import tcrw_currents

# ------------------------------------------------------------------
def compare_steady_states(L_A, omega, D_r):
    """
    L_A is the authors' 'L' (grid is (L_A+1)x(L_A+1)).
    User runs with L_U = L_A + 1 sites.
    """
    L_U = L_A + 1
    print(f"\n=== (L_A={L_A}, L_U={L_U}, omega={omega}, D_r={D_r}) ===")

    # Authors
    W_A = TRW.build_sparse_transition_matrix(L_A, omega, D_r)
    pi_A, lam_A = TRW.solve_steady_state_sparse(W_A, L_A)

    # User
    P_U = tcrw_obc.build_transition_matrix(omega, D_r, L_U)

    # Column-stochasticity check on both
    col_A = np.asarray(W_A.sum(axis=0)).ravel()
    col_U = np.asarray(P_U.sum(axis=0)).ravel()
    print(f"  max|col_sum(W_A)-1| = {np.max(np.abs(col_A-1)):.2e}")
    print(f"  max|col_sum(P_U)-1| = {np.max(np.abs(col_U-1)):.2e}")

    # User steady state.  np.linalg.eig returns an arbitrary scaling / sign for
    # the eigenvector; force it positive before normalising.
    Pd = P_U.toarray()
    vals, vecs = np.linalg.eig(Pd)
    idx1 = np.argmin(np.abs(vals - 1.0))
    pi_U = vecs[:, idx1].real
    # Fix sign so that the sum is positive (probability, not its negative)
    if pi_U.sum() < 0:
        pi_U = -pi_U
    pi_U[pi_U < 0] = 0.0
    total = pi_U.sum()
    if total <= 0:
        raise RuntimeError("user steady state zero; inspect eigenvector")
    pi_U = pi_U / total

    # ------------ Reduce both to P(X,Y) on a common grid ------------
    # Authors: pi_A index = (i*(L_A+1)+j)*4 + d
    P_XY_A = np.zeros((L_U, L_U))
    for i in range(L_U):
        for j in range(L_U):
            for d_name in TRW.directions:
                P_XY_A[i, j] += pi_A[TRW.index(i, j, d_name, L_A)]

    # User: pi_U index = d*L_U^2 + y*L_U + x
    P_XY_U = np.zeros((L_U, L_U))
    for x in range(L_U):
        for y in range(L_U):
            for d in range(4):
                P_XY_U[x, y] += pi_U[tcrw_obc.state_index(x, y, d, L_U)]

    # Both should already be normalised; renormalise defensively
    P_XY_A /= P_XY_A.sum()
    P_XY_U /= P_XY_U.sum()

    diff = np.max(np.abs(P_XY_A - P_XY_U))
    rel  = diff / np.max(P_XY_A)
    print(f"  max|P_XY_A - P_XY_U| = {diff:.3e}   (relative {rel:.3e})")

    # ------------ Currents decomposition ------------
    # Authors decomposition (flux-flavour split)
    J1_A, J2_A = TRW.calculate_J1_J2_with_boundaries(
        pi_A, W_A, L_A, D_r, omega)

    # User decomposition (exact_currents in tcrw_currents.py).  The user
    # function uses L = number of sites per axis = L_U.  It calls
    # tcrw_obc.exact_steady_state internally, which is consistent since
    # that is the same build_transition_matrix.
    out = tcrw_currents.exact_currents(omega, D_r, L_U)
    Jx_Dr_U = out["Jx_Dr"]; Jy_Dr_U = out["Jy_Dr"]
    Jx_om_U = out["Jx_om"]; Jy_om_U = out["Jy_om"]

    # Authors J1 shape (L_A+1, L_A+1, 2) = (L_U, L_U, 2)
    Jx_Dr_A = J1_A[..., 0]
    Jy_Dr_A = J1_A[..., 1]
    Jx_om_A = J2_A[..., 0]
    Jy_om_A = J2_A[..., 1]

    d_Dr  = max(np.max(np.abs(Jx_Dr_A - Jx_Dr_U)),
                np.max(np.abs(Jy_Dr_A - Jy_Dr_U)))
    d_om  = max(np.max(np.abs(Jx_om_A - Jx_om_U)),
                np.max(np.abs(Jy_om_A - Jy_om_U)))
    print(f"  max|J_Dr_A - J_Dr_U|    = {d_Dr:.3e}")
    print(f"  max|J_omega_A - J_om_U| = {d_om:.3e}")

    return dict(P_diff=diff, J_Dr_diff=d_Dr, J_om_diff=d_om)


if __name__ == "__main__":
    cases = [
        (5,  1.0, 0.01),
        (5,  0.5, 0.1),
        (5,  0.0, 0.001),
        (9,  1.0, 1e-3),
        (9,  0.7, 1e-2),
    ]
    for L_A, om, Dr in cases:
        compare_steady_states(L_A, om, Dr)
