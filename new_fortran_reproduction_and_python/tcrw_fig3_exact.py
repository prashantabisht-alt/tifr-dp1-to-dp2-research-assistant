"""
TCRW Fig 3 — exact-Python ("black-line") observables for every panel.
======================================================================

Single source of truth for the Fig 3 cross-checks.  Calls the authors'
reference module `TRW._original_code_by_paperauthors.py` for the
transition matrix, steady state, and J_Dr / J_omega decomposition; this
file just wraps the per-panel reductions (edge/bulk averages, wall
totals, per-y vector fields, angles) so each *_crosscheck.py can stay
small.

L convention
------------
We adopt **authors' L throughout this file** (`L_paper`).  At
`L_paper = 10` the lattice is 11 × 11 (sites 0..10 inclusive).
The Fortran summary files tag rows with `L_paper`, so the cross-check
is a 1:1 lookup.

Author convention summary (verified bit-exact against TRW.py):
  - directions ['↑','→','↓','←'] = 0..3
  - dir_to_vec ↑(0,1) →(1,0) ↓(0,-1) ←(-1,0)
  - state index = (i*(L+1)+j)*4 + d        (site-major)
  - noise (prob D_r):    P(CCW) = ω,   P(CW) = 1-ω
  - chiral (prob 1-D_r): P(CW) = ω,   P(CCW) = 1-ω
  - blocked chiral:      self-loop (1-D_r), no rotation

Author angle convention
-----------------------
θ = atan2(Jy, Jx) measured from +x axis, in radians ∈ [-π, π].
At ω=1, small D_r:  θ_J_ω = +π/4 (constant);  θ_J_Dr = -π/2.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations

import importlib.util
import os
from typing import Dict

import numpy as np
import scipy.sparse.linalg as spla


# ---------------------------------------------------------------------------
# 1.  Load authors' module
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_authors_trw():
    """Import authors' TRW.py.  The filename has dots, so use importlib."""
    candidates = [
        os.path.join(_HERE, "TRW._original_code_by_paperauthors.py"),
    ]
    for path in candidates:
        if os.path.exists(path):
            spec = importlib.util.spec_from_file_location("TRW_authors", path)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod
    raise FileNotFoundError(
        "TRW._original_code_by_paperauthors.py not found next to "
        f"{_HERE} — make sure tcrw_fig3_exact.py lives in the same folder."
    )


TRW = _load_authors_trw()


# ---------------------------------------------------------------------------
# 2.  Steady state π and current decomposition
# ---------------------------------------------------------------------------
def steady_state_and_currents(omega: float, D_r: float, L_paper: int):
    """
    Returns
    -------
    pi : (4 (L+1)^2,) probability vector (sums to 1, non-negative)
    J_Dr : (L+1, L+1, 2)   noise-attributed current per site
    J_om : (L+1, L+1, 2)   chiral-attributed current per site
    W : sparse W (kept for downstream uses)

    L = L_paper (authors' convention).  Sites are 0..L inclusive ⇒ (L+1)×(L+1).
    """
    L = L_paper
    W = TRW.build_sparse_transition_matrix(L, omega, D_r)
    n = (L + 1) * (L + 1) * 4

    if n <= 600:
        # tiny: dense eig (more reliable than ARPACK for L≤4)
        Wd = W.toarray()
        evals, evecs = np.linalg.eig(Wd)
        idx = int(np.argmin(np.abs(evals - 1.0)))
        pi = evecs[:, idx].real
    else:
        # large: shift-invert sparse around λ=1
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
    pi = pi / pi.sum()

    J_Dr, J_om = TRW.calculate_J1_J2_with_boundaries(pi, W, L, D_r, omega)
    return pi, J_Dr, J_om, W


# ---------------------------------------------------------------------------
# 3.  Per-site spatial probability map (sum over 4 directors)
# ---------------------------------------------------------------------------
def P_xy_from_pi(pi: np.ndarray, L_paper: int) -> np.ndarray:
    L = L_paper
    P = np.zeros((L + 1, L + 1))
    for i in range(L + 1):
        for j in range(L + 1):
            for dn in TRW.directions:
                P[i, j] += pi[TRW.index(i, j, dn, L)]
    return P


# ---------------------------------------------------------------------------
# 4.  Per-site averages — edge vs bulk (matches Fortran fig3a/3f convention)
# ---------------------------------------------------------------------------
def edge_bulk_per_site(pi: np.ndarray, L_paper: int):
    """
    Returns (P_edge_per_site, P_bulk_per_site, n_edge_sites, n_bulk_sites).
    Edge = any (i, j) on the outer frame (i ∈ {0, L} or j ∈ {0, L}).
    Bulk = strict interior.
    Direction-summed.  Authors' L convention.
    """
    L = L_paper
    P = P_xy_from_pi(pi, L)
    edge_mask = np.zeros((L + 1, L + 1), dtype=bool)
    edge_mask[0, :] = True
    edge_mask[L, :] = True
    edge_mask[:, 0] = True
    edge_mask[:, L] = True
    bulk_mask = ~edge_mask

    n_edge = int(edge_mask.sum())                 # = 4(L+1) - 4 = 4L
    n_bulk = int(bulk_mask.sum())                 # = (L-1)^2
    P_edge = float(P[edge_mask].sum() / n_edge)
    P_bulk = float(P[bulk_mask].sum() / n_bulk)
    return P_edge, P_bulk, n_edge, n_bulk


# ---------------------------------------------------------------------------
# 5.  Left-wall current totals — fig3b / fig3g observable
# ---------------------------------------------------------------------------
def left_wall_J_totals(J_Dr: np.ndarray, J_om: np.ndarray, L_paper: int):
    """
    J fields have shape (L+1, L+1, 2).  Left wall is i = 0, j = 0..L.

    Returns full-wall sums (matches Fortran fig3b/3g/3hij outputs:
    `th_Dr_tot`, `th_om_tot` are atan2(Σ_full Jy, Σ_full Jx)) AND
    interior-only sums (`th_Dr_tot_inner`, `th_om_tot_inner` over
    y = 1..L-1, the paper convention for plotting panels 3(e), 3(j)).

    Why two flavours
    ----------------
    Corner sites y = 0, y = L sit on TWO walls and have a constant Jx
    component (geometry, not chirality).  Including them in Σ J makes
    Jx_sum ≈ const ≠ 0, so atan2(Σ Jy, Σ Jx) sweeps SMOOTHLY from
    +π/2 → 0 → -π/2 through 0.  Excluding them gives Σ Jx ≈ 0 and
    a sharp ±π/2 step at ω = 0.5 — what paper Fig 3(j) shows.

    Cross-check files compare against Fortran's full-wall sums, so
    `th_*_tot` keeps the full-wall convention (no breaking change).
    Plotting code uses `th_*_tot_inner` for paper-faithful display.
    """
    L = L_paper
    # Full-wall sums — matches Fortran (cross-check uses these)
    Jx_Dr_tot = float(J_Dr[0, :, 0].sum())
    Jy_Dr_tot = float(J_Dr[0, :, 1].sum())
    Jx_om_tot = float(J_om[0, :, 0].sum())
    Jy_om_tot = float(J_om[0, :, 1].sum())
    abs_Dr = float(np.hypot(Jx_Dr_tot, Jy_Dr_tot))
    abs_om = float(np.hypot(Jx_om_tot, Jy_om_tot))
    ratio = abs_Dr / abs_om if abs_om > 0 else np.nan

    # Interior-only sums — paper convention for plotting panels 3(e), 3(j)
    Jx_Dr_in = float(J_Dr[0, 1:L, 0].sum())
    Jy_Dr_in = float(J_Dr[0, 1:L, 1].sum())
    Jx_om_in = float(J_om[0, 1:L, 0].sum())
    Jy_om_in = float(J_om[0, 1:L, 1].sum())

    return dict(
        Jx_Dr_tot=Jx_Dr_tot, Jy_Dr_tot=Jy_Dr_tot,
        Jx_om_tot=Jx_om_tot, Jy_om_tot=Jy_om_tot,
        abs_Dr=abs_Dr, abs_om=abs_om, ratio=ratio,
        # FULL-WALL angles (cross-check vs Fortran)
        th_Dr_tot=float(np.arctan2(Jy_Dr_tot, Jx_Dr_tot)),
        th_om_tot=float(np.arctan2(Jy_om_tot, Jx_om_tot)),
        # INTERIOR-ONLY angles (paper-faithful plotting)
        th_Dr_tot_inner=float(np.arctan2(Jy_Dr_in, Jx_Dr_in)),
        th_om_tot_inner=float(np.arctan2(Jy_om_in, Jx_om_in)),
    )


# ---------------------------------------------------------------------------
# 6.  Per-y left-wall current — fig3cde / fig3hij observable
# ---------------------------------------------------------------------------
def left_wall_J_per_y(J_Dr: np.ndarray, J_om: np.ndarray, L_paper: int):
    """
    Per-site current vectors on the left wall (i = 0).
    Returns dict of length-(L+1) arrays, indexed by y = 0..L.
    """
    Jx_Dr = J_Dr[0, :, 0].astype(float)
    Jy_Dr = J_Dr[0, :, 1].astype(float)
    Jx_om = J_om[0, :, 0].astype(float)
    Jy_om = J_om[0, :, 1].astype(float)
    return dict(
        y=np.arange(L_paper + 1, dtype=int),
        Jx_Dr=Jx_Dr, Jy_Dr=Jy_Dr,
        Jx_om=Jx_om, Jy_om=Jy_om,
        absJ_Dr=np.hypot(Jx_Dr, Jy_Dr),
        absJ_om=np.hypot(Jx_om, Jy_om),
        th_Dr=np.arctan2(Jy_Dr, Jx_Dr),
        th_om=np.arctan2(Jy_om, Jx_om),
    )


# ---------------------------------------------------------------------------
# 7.  Convenience: scan a parameter for each panel
# ---------------------------------------------------------------------------
def scan_panel_a(D_r_grid, L_paper_list, omega: float = 1.0):
    """Fig 3(a): per-site P_edge/P_bulk vs D_r at fixed ω, multiple L."""
    out = {}
    for L in L_paper_list:
        rec = []
        for D in D_r_grid:
            pi, *_ = steady_state_and_currents(omega, float(D), L)
            P_e, P_b, ne, nb = edge_bulk_per_site(pi, L)
            rec.append((float(D), P_e, P_b, P_e / P_b if P_b > 0 else np.nan, ne, nb))
        out[L] = np.array(rec, dtype=float)
    return out  # dict L -> array (n_Dr, 6)


def scan_panel_b(D_r_grid, L_paper_list, omega: float = 1.0):
    """Fig 3(b): |J_Dr|/|J_ω| at left wall vs D_r at fixed ω, multiple L."""
    out = {}
    for L in L_paper_list:
        rec = []
        for D in D_r_grid:
            _, J_Dr, J_om, _ = steady_state_and_currents(omega, float(D), L)
            tots = left_wall_J_totals(J_Dr, J_om, L)
            rec.append((float(D), tots["ratio"], tots["abs_Dr"], tots["abs_om"],
                        tots["Jx_Dr_tot"], tots["Jy_Dr_tot"],
                        tots["Jx_om_tot"], tots["Jy_om_tot"]))
        out[L] = np.array(rec, dtype=float)
    return out


def scan_panel_cde(D_r_grid, L_paper: int = 10, omega: float = 1.0):
    """Fig 3(c)(d)(e): per-y left-wall current and totals vs D_r."""
    per_site = []        # list of dicts (one per D_r)
    totals = []
    for D in D_r_grid:
        _, J_Dr, J_om, _ = steady_state_and_currents(omega, float(D), L_paper)
        per_site.append(left_wall_J_per_y(J_Dr, J_om, L_paper))
        totals.append(left_wall_J_totals(J_Dr, J_om, L_paper))
    return per_site, totals


def scan_panel_f(omega_grid, L_paper_list, D_r: float = 1e-3):
    """Fig 3(f): per-site P_edge/P_bulk vs ω at fixed D_r, multiple L."""
    out = {}
    for L in L_paper_list:
        rec = []
        for w in omega_grid:
            pi, *_ = steady_state_and_currents(float(w), D_r, L)
            P_e, P_b, ne, nb = edge_bulk_per_site(pi, L)
            rec.append((float(w), P_e, P_b, P_e / P_b if P_b > 0 else np.nan, ne, nb))
        out[L] = np.array(rec, dtype=float)
    return out


def scan_panel_g(omega_grid, L_paper_list, D_r: float = 1e-3):
    """Fig 3(g): |J_Dr|/|J_ω| at left wall vs ω at fixed D_r, multiple L."""
    out = {}
    for L in L_paper_list:
        rec = []
        for w in omega_grid:
            _, J_Dr, J_om, _ = steady_state_and_currents(float(w), D_r, L)
            tots = left_wall_J_totals(J_Dr, J_om, L)
            rec.append((float(w), tots["ratio"], tots["abs_Dr"], tots["abs_om"],
                        tots["Jx_Dr_tot"], tots["Jy_Dr_tot"],
                        tots["Jx_om_tot"], tots["Jy_om_tot"]))
        out[L] = np.array(rec, dtype=float)
    return out


def scan_panel_hij(omega_grid, L_paper: int = 10, D_r: float = 1e-3):
    """Fig 3(h)(i)(j): per-y left-wall current and totals vs ω."""
    per_site = []
    totals = []
    for w in omega_grid:
        _, J_Dr, J_om, _ = steady_state_and_currents(float(w), D_r, L_paper)
        per_site.append(left_wall_J_per_y(J_Dr, J_om, L_paper))
        totals.append(left_wall_J_totals(J_Dr, J_om, L_paper))
    return per_site, totals


# ---------------------------------------------------------------------------
# 8.  Self-test (smoke test on import or `python tcrw_fig3_exact.py`)
# ---------------------------------------------------------------------------
def _selftest():
    """Quick numerical sanity tests."""
    print("[fig3 exact] self-test")
    pi, J_Dr, J_om, W = steady_state_and_currents(1.0, 1e-3, 10)
    P_e, P_b, ne, nb = edge_bulk_per_site(pi, 10)
    tots = left_wall_J_totals(J_Dr, J_om, 10)
    print(f"  ω=1, D_r=1e-3, L=10:")
    print(f"    n_edge = {ne}  n_bulk = {nb}   (expect 40 / 81)")
    print(f"    P_edge_per_site = {P_e:.4e}")
    print(f"    P_bulk_per_site = {P_b:.4e}")
    print(f"    ratio = {P_e / P_b:.3f}")
    print(f"    |J_Dr|_wall = {tots['abs_Dr']:.4e}")
    print(f"    |J_om|_wall = {tots['abs_om']:.4e}")
    print(f"    ratio = {tots['ratio']:.3f}")
    print(f"    θ_J_Dr = {tots['th_Dr_tot']:.3f}   (expect ≈ -π/2 = -1.571)")
    print(f"    θ_J_om = {tots['th_om_tot']:.3f}   (expect ≈ +π/4 =  0.785)")


if __name__ == "__main__":
    _selftest()
