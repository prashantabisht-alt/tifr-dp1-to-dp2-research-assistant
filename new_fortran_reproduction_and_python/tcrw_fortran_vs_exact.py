"""
Cross-check the new Fortran Fig 2 occupancy/currents against the
exact-Python steady state solved from the transition matrix.

The new Fortran driver uses a wall-ring convention: grid is L x L but
the outer ring (x in {0, L-1}, y in {0, L-1}) is a WALL (unreachable).
So the physical playground is a (L-2) x (L-2) OBC box that lives at
lattice coordinates x,y in {1..L-2}.

Exact-Python equivalence:
    user Fortran at L_F                (playground (L_F-2)^2)
  <==>
    authors' TRW with L_A = L_F - 3    (grid 0..L_A, so (L_F-2)^2 sites)

For L_F = 10 -> authors L_A = 7 (grid 0..7 -> 8x8 playground).

We load the Fortran occupancy file, extract P on the inner (L_F-2)^2
playground, build the exact steady-state P on 0..7, and compare.
"""

import os, sys, numpy as np

for c in (
    "/sessions/elegant-wizardly-einstein/mnt/TIFR DP1 to DP2 Research Assistant",
    "/Users/prashantbisht/Documents/Claude/Projects/TIFR DP1 to DP2 Research Assistant"):
    if os.path.isdir(c):
        ROOT = c
        break

NEW = os.path.join(ROOT, "new_fortran_reproduction_and_python")
sys.path.insert(0, ROOT)
sys.path.insert(0, NEW)

# Load authors' module
import importlib.util
spec = importlib.util.spec_from_file_location(
    "TRW_authors",
    os.path.join(NEW, "TRW._original_code_by_paperauthors.py"))
TRW = importlib.util.module_from_spec(spec)
spec.loader.exec_module(TRW)


def load_fortran_occ(path, L_F):
    """Read Fortran occupancy file: columns 'x y P'."""
    arr = np.loadtxt(path)
    P = np.zeros((L_F, L_F))
    for x, y, p in arr:
        P[int(x), int(y)] = p
    return P


def load_fortran_current(path):
    """Read Fortran current file: columns 'x y Jx Jy |J|'."""
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        arr = arr[None, :]
    return arr  # raw; we'll index by (x,y) directly


def compare_fig2(omega, D_r=1e-3, L_F=10):
    """
    Compare Fortran Fig 2 P(X,Y) and J_Dr/J_omega at given omega against
    the exact-Python steady state on the 8x8 playground.
    """
    print(f"\n=== Fig 2, omega={omega}, D_r={D_r}, L_F={L_F} ===")

    # --- Fortran side ---
    occ_path = os.path.join(NEW, f"tcrw_fig2_occ_w{omega:.1f}.txt")
    P_F_full = load_fortran_occ(occ_path, L_F)
    # Extract inner playground (x,y in 1..L_F-2 -> indices 0..L_F-3)
    P_F = P_F_full[1:L_F-1, 1:L_F-1]
    P_F = P_F / P_F.sum()   # renormalize on playground

    # --- Exact-Python side: authors' L_A = L_F - 3 ---
    L_A = L_F - 3
    W = TRW.build_sparse_transition_matrix(L_A, omega, D_r)
    pi, lam = TRW.solve_steady_state_sparse(W, L_A)
    P_E = np.zeros((L_A + 1, L_A + 1))
    for i in range(L_A + 1):
        for j in range(L_A + 1):
            for dn in TRW.directions:
                P_E[i, j] += pi[TRW.index(i, j, dn, L_A)]
    P_E = P_E / P_E.sum()

    # --- Compare ---
    assert P_F.shape == P_E.shape, f"shape mismatch {P_F.shape} vs {P_E.shape}"
    diff = np.abs(P_F - P_E)
    print(f"  playground shape: {P_F.shape}")
    print(f"  max|P_F - P_E| = {diff.max():.3e}   (rel {diff.max()/P_E.max():.2%})")
    print(f"  RMS (P_F - P_E) = {np.sqrt(np.mean((P_F-P_E)**2)):.3e}")
    print(f"  sum P_F = {P_F.sum():.6f}   sum P_E = {P_E.sum():.6f}")
    print(f"  max P_F = {P_F.max():.3e}   max P_E = {P_E.max():.3e}")

    # --- Currents (if files exist) ---
    for kind_F, label in [("JDr", "J_Dr"), ("Jomega", "J_omega")]:
        cpath = os.path.join(NEW, f"tcrw_fig2_{kind_F}_w{omega:.1f}.txt")
        if not os.path.exists(cpath):
            continue
        arr = load_fortran_current(cpath)
        # build (L_F-2, L_F-2) Jx, Jy from columns
        Jx_F = np.zeros((L_A + 1, L_A + 1))
        Jy_F = np.zeros((L_A + 1, L_A + 1))
        for row in arr:
            x, y, Jx, Jy, Jmag = row
            ix = int(x) - 1; iy = int(y) - 1   # shift into playground
            if 0 <= ix < L_A + 1 and 0 <= iy < L_A + 1:
                Jx_F[ix, iy] = Jx
                Jy_F[ix, iy] = Jy

    # exact currents
    J1_E, J2_E = TRW.calculate_J1_J2_with_boundaries(pi, W, L_A, D_r, omega)
    # J1_E = J_Dr, J2_E = J_omega

    # Re-load each current file to match
    for kind_F, J_exact in [("JDr", J1_E), ("Jomega", J2_E)]:
        cpath = os.path.join(NEW, f"tcrw_fig2_{kind_F}_w{omega:.1f}.txt")
        if not os.path.exists(cpath):
            print(f"  [skip] {cpath} missing")
            continue
        arr = load_fortran_current(cpath)
        Jx_F = np.zeros((L_A + 1, L_A + 1))
        Jy_F = np.zeros((L_A + 1, L_A + 1))
        for row in arr:
            x, y, Jx, Jy, _ = row
            ix = int(x) - 1; iy = int(y) - 1
            if 0 <= ix < L_A + 1 and 0 <= iy < L_A + 1:
                Jx_F[ix, iy] = Jx
                Jy_F[ix, iy] = Jy
        # Fortran accumulates 1 count per successful translation.  Exact
        # calculation produces current in the same units (per step).
        dJx = Jx_F - J_exact[..., 0]
        dJy = Jy_F - J_exact[..., 1]
        mag = np.hypot(Jx_F, Jy_F)
        ref = max(np.abs(J_exact).max(), 1e-30)
        print(f"  {kind_F:6s}: max|Fortran - exact| = "
              f"(Jx) {np.abs(dJx).max():.3e}  (Jy) {np.abs(dJy).max():.3e}"
              f"   (ref scale {ref:.2e})")


if __name__ == "__main__":
    for om in (0.0, 0.5, 1.0):
        try:
            compare_fig2(omega=om)
        except FileNotFoundError as e:
            print(f"  [skip omega={om}] {e}")
