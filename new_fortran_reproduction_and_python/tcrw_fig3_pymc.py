"""
TCRW Fig 3 — self-contained Python reproduction (no TRW.py dependency)
========================================================================

Same multi-panel figure as `tcrw_fig3_authors.py`, but with the
transition matrix and current decomposition implemented inline (port of
authors' algorithm, bit-verified).  This makes the file forkable for
DP2 — change the step rule, change one line of `build_obc_matrix`, and
the rest of the pipeline still works.

What "exact" means here
-----------------------
At T = ∞ the MC walker's empirical histogram converges (by ergodicity)
to the Perron eigenvector π of the column-stochastic transition matrix.
We solve that eigenvalue problem directly:
  - dense `np.linalg.eig` for n ≤ 600 (L ≤ 11)
  - shift-invert `scipy.sparse.linalg.eigs(σ = 1)` for larger
Then we compute the current decomposition via the authors' flux-flavour
algorithm.  No randomness, no Monte Carlo.

Cross-check
-----------
At the bottom of this file: `crosscheck_authors` builds W and J at a
test point both via this module's `build_obc_matrix` and via authors'
`TRW.build_sparse_transition_matrix`, and verifies element-wise zero
diff.  This carries over from the bit-verification of `tcrw_fig4c` (the
same matrix layout) — repeating the check here keeps pymc honest in
isolation.

Convention
----------
Authors' L = max-index ⇒ grid 0..L inclusive, (L+1)² sites.
State index: s = (i*(L+1)+j)*4 + d  with directors d = 0(↑) 1(→) 2(↓) 3(←).
Noise (prob D_r):    P(CCW) = ω,   P(CW) = 1-ω
Chiral (prob 1-D_r): P(CW)  = ω,   P(CCW) = 1-ω    (translate THEN rotate)
Blocked chiral: self-loop weight (1-D_r), no rotation.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import sys
import time
import importlib.util
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _fig3_plot import make_figure


# ---------------------------------------------------------------------------
# 0. Convention constants
# ---------------------------------------------------------------------------
DX = np.array([0, 1, 0, -1], dtype=int)
DY = np.array([1, 0, -1, 0], dtype=int)


def state_index(i: int, j: int, d: int, L: int) -> int:
    return (i * (L + 1) + j) * 4 + d


# ---------------------------------------------------------------------------
# 1. OBC transition matrix (port of authors' build_sparse_transition_matrix)
# ---------------------------------------------------------------------------
def build_obc_matrix(omega: float, D_r: float, L: int) -> sp.csc_matrix:
    n = 4 * (L + 1) ** 2
    rows, cols, vals = [], [], []
    for i in range(L + 1):
        for j in range(L + 1):
            for d in range(4):
                src = state_index(i, j, d, L)
                d_ccw = (d - 1) % 4
                d_cw  = (d + 1) % 4
                # noise step
                rows.append(state_index(i, j, d_ccw, L)); cols.append(src)
                vals.append(D_r * omega)
                rows.append(state_index(i, j, d_cw,  L)); cols.append(src)
                vals.append(D_r * (1.0 - omega))
                # chiral step
                ni, nj = i + int(DX[d]), j + int(DY[d])
                if 0 <= ni <= L and 0 <= nj <= L:
                    rows.append(state_index(ni, nj, d_cw,  L)); cols.append(src)
                    vals.append((1.0 - D_r) * omega)
                    rows.append(state_index(ni, nj, d_ccw, L)); cols.append(src)
                    vals.append((1.0 - D_r) * (1.0 - omega))
                else:
                    rows.append(src); cols.append(src)
                    vals.append(1.0 - D_r)
    return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()


# ---------------------------------------------------------------------------
# 2. Steady state π via sparse Perron
# ---------------------------------------------------------------------------
def steady_state(W: sp.csc_matrix) -> np.ndarray:
    n = W.shape[0]
    if n <= 600:
        Wd = W.toarray()
        evals, evecs = np.linalg.eig(Wd)
        idx = int(np.argmin(np.abs(evals - 1.0)))
        pi = evecs[:, idx].real
    else:
        try:
            evals, evecs = spla.eigs(W, k=1, sigma=1.0 + 1e-8, which="LM",
                                     maxiter=2000, tol=1e-12)
            pi = evecs[:, 0].real
        except Exception:
            Wd = W.toarray()
            evals, evecs = np.linalg.eig(Wd)
            idx = int(np.argmin(np.abs(evals - 1.0)))
            pi = evecs[:, idx].real
    if pi.sum() < 0:
        pi = -pi
    pi[pi < 0] = 0.0
    return pi / pi.sum()


# ---------------------------------------------------------------------------
# 3. Current decomposition (port of authors' calculate_J1_J2_with_boundaries)
# ---------------------------------------------------------------------------
def calculate_J1_J2(pi: np.ndarray, W: sp.csc_matrix, L: int):
    """Returns J_Dr, J_omega each shape (L+1, L+1, 2)."""
    n = pi.size
    P_rot = np.zeros(n)
    P_chiral = np.zeros(n)
    Wc = W.tocoo()
    rows_a, cols_a, data_a = Wc.row, Wc.col, Wc.data

    # Step 1: classify incoming flux to each state
    for idx_dest, idx_src, prob in zip(rows_a, cols_a, data_a):
        if idx_src == idx_dest:
            continue
        i_d = (idx_dest // 4) // (L + 1); j_d = (idx_dest // 4) % (L + 1)
        i_s = (idx_src  // 4) // (L + 1); j_s = (idx_src  // 4) % (L + 1)
        flux = pi[idx_src] * prob
        if i_s == i_d and j_s == j_d:
            P_rot[idx_dest]    += flux
        else:
            P_chiral[idx_dest] += flux
    total = P_rot + P_chiral
    mask = total > 1e-12
    P_rot[mask]    = pi[mask] * (P_rot[mask]    / total[mask])
    P_chiral[mask] = pi[mask] * (P_chiral[mask] / total[mask])
    P_rot[~mask] = 0.0
    P_chiral[~mask] = 0.0

    # Step 2: outgoing currents
    J1 = np.zeros((L + 1, L + 1, 2))
    J2 = np.zeros((L + 1, L + 1, 2))
    for idx_dest, idx_src, prob in zip(rows_a, cols_a, data_a):
        if idx_src == idx_dest:
            continue
        i_d = (idx_dest // 4) // (L + 1); j_d = (idx_dest // 4) % (L + 1)
        i_s = (idx_src  // 4) // (L + 1); j_s = (idx_src  // 4) % (L + 1)
        if i_s == i_d and j_s == j_d:
            continue
        dx, dy = i_d - i_s, j_d - j_s
        J1[i_s, j_s, 0] += P_rot[idx_src]    * prob * dx
        J1[i_s, j_s, 1] += P_rot[idx_src]    * prob * dy
        J2[i_s, j_s, 0] += P_chiral[idx_src] * prob * dx
        J2[i_s, j_s, 1] += P_chiral[idx_src] * prob * dy
    return J1, J2


# ---------------------------------------------------------------------------
# 4. Per-panel observables
# ---------------------------------------------------------------------------
def _solve_one(omega: float, D_r: float, L: int):
    W = build_obc_matrix(omega, D_r, L)
    pi = steady_state(W)
    J_Dr, J_om = calculate_J1_J2(pi, W, L)
    return pi, J_Dr, J_om


def _edge_bulk_per_site(pi: np.ndarray, L: int):
    P = np.zeros((L + 1, L + 1))
    for i in range(L + 1):
        for j in range(L + 1):
            for d in range(4):
                P[i, j] += pi[state_index(i, j, d, L)]
    edge = np.zeros_like(P, dtype=bool)
    edge[0, :] = True; edge[L, :] = True
    edge[:, 0] = True; edge[:, L] = True
    bulk = ~edge
    n_e = int(edge.sum()); n_b = int(bulk.sum())
    return P[edge].sum() / n_e, P[bulk].sum() / n_b, n_e, n_b


def _wall_totals(J_Dr, J_om, L):
    """Full-wall sums (matches Fortran fig3b/3g/3hij; cross-check uses
    th_Dr_tot, th_om_tot) plus interior-only angles
    th_Dr_tot_inner/th_om_tot_inner (paper convention for plotting Fig
    3(e) and 3(j) — corner sites contaminate Jx and smooth the step)."""
    # Full-wall (y = 0..L) — matches Fortran
    Jx_Dr = float(J_Dr[0, :, 0].sum()); Jy_Dr = float(J_Dr[0, :, 1].sum())
    Jx_om = float(J_om[0, :, 0].sum()); Jy_om = float(J_om[0, :, 1].sum())
    abs_Dr = float(np.hypot(Jx_Dr, Jy_Dr))
    abs_om = float(np.hypot(Jx_om, Jy_om))
    # Interior-only (y = 1..L-1)
    Jx_Dr_in = float(J_Dr[0, 1:L, 0].sum()); Jy_Dr_in = float(J_Dr[0, 1:L, 1].sum())
    Jx_om_in = float(J_om[0, 1:L, 0].sum()); Jy_om_in = float(J_om[0, 1:L, 1].sum())
    return dict(
        abs_Dr=abs_Dr, abs_om=abs_om,
        ratio=abs_Dr / abs_om if abs_om > 0 else np.nan,
        # FULL-WALL angles (cross-check vs Fortran):
        th_Dr_tot=float(np.arctan2(Jy_Dr, Jx_Dr)),
        th_om_tot=float(np.arctan2(Jy_om, Jx_om)),
        # INTERIOR-ONLY angles (paper-faithful plotting):
        th_Dr_tot_inner=float(np.arctan2(Jy_Dr_in, Jx_Dr_in)),
        th_om_tot_inner=float(np.arctan2(Jy_om_in, Jx_om_in)),
    )


def _wall_per_y(J_Dr, J_om, L):
    Jx_Dr = J_Dr[0, :, 0]; Jy_Dr = J_Dr[0, :, 1]
    Jx_om = J_om[0, :, 0]; Jy_om = J_om[0, :, 1]
    return dict(
        y=np.arange(L + 1),
        Jx_Dr=Jx_Dr, Jy_Dr=Jy_Dr,
        Jx_om=Jx_om, Jy_om=Jy_om,
        absJ_Dr=np.hypot(Jx_Dr, Jy_Dr),
        absJ_om=np.hypot(Jx_om, Jy_om),
        th_Dr=np.arctan2(Jy_Dr, Jx_Dr),
        th_om=np.arctan2(Jy_om, Jx_om),
    )


# ---------------------------------------------------------------------------
# 5. Panel scans (matching tcrw_fig3_exact's API, but using local code)
# ---------------------------------------------------------------------------
def scan_panel_a(D_r_grid, L_list, omega: float = 1.0):
    out = {}
    for L in L_list:
        rec = []
        for D in D_r_grid:
            pi, *_ = _solve_one(omega, float(D), L)
            P_e, P_b, ne, nb = _edge_bulk_per_site(pi, L)
            rec.append((float(D), P_e, P_b, P_e / P_b if P_b > 0 else np.nan, ne, nb))
        out[L] = np.array(rec, dtype=float)
    return out


def scan_panel_b(D_r_grid, L_list, omega: float = 1.0):
    out = {}
    for L in L_list:
        rec = []
        for D in D_r_grid:
            _, J_Dr, J_om = _solve_one(omega, float(D), L)
            t = _wall_totals(J_Dr, J_om, L)
            rec.append((float(D), t["ratio"], t["abs_Dr"], t["abs_om"]))
        out[L] = np.array(rec, dtype=float)
    return out


def scan_panel_cde(D_r_grid, L_paper: int = 10, omega: float = 1.0):
    per_site, totals = [], []
    for D in D_r_grid:
        _, J_Dr, J_om = _solve_one(omega, float(D), L_paper)
        per_site.append(_wall_per_y(J_Dr, J_om, L_paper))
        totals.append(_wall_totals(J_Dr, J_om, L_paper))
    return per_site, totals


def scan_panel_f(omega_grid, L_list, D_r: float = 1e-3):
    out = {}
    for L in L_list:
        rec = []
        for w in omega_grid:
            pi, *_ = _solve_one(float(w), D_r, L)
            P_e, P_b, ne, nb = _edge_bulk_per_site(pi, L)
            rec.append((float(w), P_e, P_b, P_e / P_b if P_b > 0 else np.nan, ne, nb))
        out[L] = np.array(rec, dtype=float)
    return out


def scan_panel_g(omega_grid, L_list, D_r: float = 1e-3):
    out = {}
    for L in L_list:
        rec = []
        for w in omega_grid:
            _, J_Dr, J_om = _solve_one(float(w), D_r, L)
            t = _wall_totals(J_Dr, J_om, L)
            rec.append((float(w), t["ratio"], t["abs_Dr"], t["abs_om"]))
        out[L] = np.array(rec, dtype=float)
    return out


def scan_panel_hij(omega_grid, L_paper: int = 10, D_r: float = 1e-3):
    per_site, totals = [], []
    for w in omega_grid:
        _, J_Dr, J_om = _solve_one(float(w), D_r, L_paper)
        per_site.append(_wall_per_y(J_Dr, J_om, L_paper))
        totals.append(_wall_totals(J_Dr, J_om, L_paper))
    return per_site, totals


# ---------------------------------------------------------------------------
# 6. Author cross-check (bit-verification of build_obc_matrix)
# ---------------------------------------------------------------------------
def crosscheck_authors(omega: float, D_r: float, L: int = 2,
                       tol: float = 1e-12):
    p = os.path.join(HERE, "TRW._original_code_by_paperauthors.py")
    if not os.path.exists(p):
        print("  [warn] authors' TRW.py not found — cross-check skipped")
        return
    spec = importlib.util.spec_from_file_location("TRW", p)
    TRW = importlib.util.module_from_spec(spec); spec.loader.exec_module(TRW)
    W_A = TRW.build_sparse_transition_matrix(L, omega, D_r).toarray()
    W_U = build_obc_matrix(omega, D_r, L).toarray()
    diff = float(np.max(np.abs(W_A - W_U)))
    print(f"  ω={omega} D_r={D_r} L={L}: max|W_A − W_U| = {diff:.2e}")
    assert diff < tol, f"matrix mismatch {diff:.2e} > tol {tol:.2e}"


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 62)
    print(" TCRW Fig 3 — self-contained Python reproduction (pymc)")
    print("=" * 62)

    print("\n[0] Bit-check vs authors' TRW (build_obc_matrix)")
    crosscheck_authors(0.3, 0.1, L=2)
    crosscheck_authors(1.0, 1e-3, L=2)
    crosscheck_authors(0.5, 0.5, L=3)

    # Parameter grids — same as authors' version
    N_Dr  = 25
    # Stop just below 1: at D_r=1 the walker never translates, eigenspace
    # at λ=1 is degenerate, and currents → 0/0 → ratio NaN.
    D_r_grid    = np.logspace(-4.0, -0.01, N_Dr)
    N_omega   = 21
    omega_grid  = np.linspace(0.0, 1.0, N_omega)
    L_LIST_a    = (4, 9, 19, 49)
    L_LIST_f    = (10, 19, 49)
    L_HEATMAP   = 10
    D_r_FIXED   = 1e-3

    data = {}
    print("\n[1] Panel (a)")
    t0 = time.time()
    data["a"] = scan_panel_a(D_r_grid, L_LIST_a, omega=1.0)
    print(f"    cpu = {time.time()-t0:.1f}s")
    print("[2] Panel (b)")
    t0 = time.time()
    data["b"] = scan_panel_b(D_r_grid, L_LIST_a, omega=1.0)
    print(f"    cpu = {time.time()-t0:.1f}s")
    print("[3] Panel (c)(d)(e)")
    t0 = time.time()
    p_c, t_c = scan_panel_cde(D_r_grid, L_paper=L_HEATMAP, omega=1.0)
    data["cde_per"] = p_c; data["cde_tot"] = t_c
    print(f"    cpu = {time.time()-t0:.1f}s")
    print("[4] Panel (f)")
    t0 = time.time()
    data["f"] = scan_panel_f(omega_grid, L_LIST_f, D_r=D_r_FIXED)
    print(f"    cpu = {time.time()-t0:.1f}s")
    print("[5] Panel (g)")
    t0 = time.time()
    data["g"] = scan_panel_g(omega_grid, L_LIST_f, D_r=D_r_FIXED)
    print(f"    cpu = {time.time()-t0:.1f}s")
    print("[6] Panel (h)(i)(j)")
    t0 = time.time()
    p_w, t_w = scan_panel_hij(omega_grid, L_paper=L_HEATMAP, D_r=D_r_FIXED)
    data["hij_per"] = p_w; data["hij_tot"] = t_w
    print(f"    cpu = {time.time()-t0:.1f}s")

    print("\n[7] Plot")
    make_figure(data, D_r_grid=D_r_grid, omega_grid=omega_grid,
                L_LIST_a=L_LIST_a, L_LIST_f=L_LIST_f,
                savepath="tcrw_fig3_pymc.png",
                title_extra="(self-contained Python — no TRW.py dep)")
    print("\nDone.")
