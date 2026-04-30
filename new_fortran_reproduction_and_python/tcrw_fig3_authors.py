"""
TCRW Fig 3 — exact reproduction using authors' TRW.py
======================================================

Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020, Fig 3.

10 panels (2 rows × 5 cols) computed via authors' TRW.py
(`build_sparse_transition_matrix`, `solve_steady_state_sparse`,
`calculate_J1_J2_with_boundaries`).

This script is the **TRW.py path**.  Its sister `tcrw_fig3_pymc.py`
computes the same panels with an inlined transition-matrix builder (no
TRW.py dependency); they are bit-verified against each other and
should produce pixel-identical plots.

Both scripts populate the same `data` dict shape and call
`_fig3_plot.make_figure`, so the plotting code is single-sourced.

Panel data shape (so DP2 forks know what to fill):
  data["a"][L]  (n_Dr , 6)  [D_r , P_e , P_b , ratio , n_e , n_b]
  data["b"][L]  (n_Dr , 4)  [D_r , ratio_wall , |J_Dr|_wall , |J_om|_wall]
  data["cde_per"]   list of dicts (per D_r)  per-y wall vectors and angles
  data["cde_tot"]   list of dicts (per D_r)  totals over wall + angles
  data["f"][L]  (n_omega, 6)  [ω , P_e , P_b , ratio , n_e , n_b]
  data["g"][L]  (n_omega, 4)  [ω , ratio_wall , |J_Dr|_wall , |J_om|_wall]
  data["hij_per"]   list of dicts (per ω)
  data["hij_tot"]   list of dicts (per ω)

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from tcrw_fig3_exact import (scan_panel_a, scan_panel_b, scan_panel_cde,
                             scan_panel_f, scan_panel_g, scan_panel_hij)
from _fig3_plot import make_figure


def main():
    print("=" * 64)
    print(" TCRW Fig 3 — exact via authors' TRW.py")
    print("=" * 64)

    # Parameter grids — must match tcrw_fig3_pymc.py exactly so the two
    # PNGs can be diffed pixel-by-pixel.
    N_Dr        = 25
    N_omega     = 21
    # Stop D_r grid just below 1.  At D_r = 1 the walker never translates,
    # the transition matrix has a degenerate eigenspace at λ = 1, and the
    # Perron solver returns an arbitrary stationary vector — currents are
    # 0/0 → ratio is NaN.  log10(0.99) ≈ -0.0044; we use -0.01 for safety.
    D_r_grid    = np.logspace(-4.0, -0.01, N_Dr)
    omega_grid  = np.linspace(0.0, 1.0, N_omega)
    L_LIST_a    = (4, 9, 19, 49)        # top row (D_r scan)
    L_LIST_f    = (10, 19, 49)          # bottom row (ω scan)
    L_HEATMAP   = 10                    # for panels c, d, e, h, i, j
    D_r_FIXED   = 1e-3

    data = {}
    print("\n[1] Panel (a) — P_edge/P_bulk vs D_r")
    t0 = time.time()
    data["a"] = scan_panel_a(D_r_grid, L_LIST_a, omega=1.0)
    print(f"    cpu = {time.time()-t0:.1f}s")

    print("[2] Panel (b) — |J_Dr|/|J_om| vs D_r")
    t0 = time.time()
    data["b"] = scan_panel_b(D_r_grid, L_LIST_a, omega=1.0)
    print(f"    cpu = {time.time()-t0:.1f}s")

    print("[3] Panel (c)(d)(e) — wall vectors vs D_r")
    t0 = time.time()
    p_c, t_c = scan_panel_cde(D_r_grid, L_paper=L_HEATMAP, omega=1.0)
    data["cde_per"] = p_c; data["cde_tot"] = t_c
    print(f"    cpu = {time.time()-t0:.1f}s")

    print("[4] Panel (f) — P_edge/P_bulk vs ω")
    t0 = time.time()
    data["f"] = scan_panel_f(omega_grid, L_LIST_f, D_r=D_r_FIXED)
    print(f"    cpu = {time.time()-t0:.1f}s")

    print("[5] Panel (g) — |J_Dr|/|J_om| vs ω")
    t0 = time.time()
    data["g"] = scan_panel_g(omega_grid, L_LIST_f, D_r=D_r_FIXED)
    print(f"    cpu = {time.time()-t0:.1f}s")

    print("[6] Panel (h)(i)(j) — wall vectors vs ω")
    t0 = time.time()
    p_w, t_w = scan_panel_hij(omega_grid, L_paper=L_HEATMAP, D_r=D_r_FIXED)
    data["hij_per"] = p_w; data["hij_tot"] = t_w
    print(f"    cpu = {time.time()-t0:.1f}s")

    print("\n[7] Plot")
    make_figure(data, D_r_grid=D_r_grid, omega_grid=omega_grid,
                L_LIST_a=L_LIST_a, L_LIST_f=L_LIST_f,
                savepath=os.path.join(HERE, "tcrw_fig3_authors.png"),
                title_extra="(authors' TRW.py via tcrw_fig3_exact)")
    print("\nDone.")


if __name__ == "__main__":
    main()
