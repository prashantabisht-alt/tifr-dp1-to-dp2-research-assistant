"""
Fig 3(g) cross-check: |J_Dr|_wall / |J_omega|_wall vs ω at fixed D_r=1e-3, L=10.

Loads tcrw_fig3g_summary.txt (raw counts), normalises by reconstructed
T_use, and overlays the exact-Python prediction.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os, sys, time
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tcrw_fig3_exact as ex
from _fig3_crosscheck_common import load_summary, reconstruct_T_use, rel_error


def crosscheck(summary_path: str = "tcrw_fig3g_summary.txt",
               plot_path: str = "tcrw_fig3g_crosscheck.png",
               D_r_fixed: float = 1e-3):
    print(f"loading {summary_path}")
    arr = load_summary(summary_path)
    # cols: L  ω  ratio  |J_Dr|  |J_om|  Jx_Dr  Jy_Dr  Jx_om  Jy_om
    L_list = sorted(int(x) for x in np.unique(arr[:, 0]))
    print(f"  L values: {L_list}, D_r = {D_r_fixed}")

    fig, (ax_ratio, ax_Jdr, ax_Jom) = plt.subplots(1, 3, figsize=(15, 4.5))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(L_list)))

    summary = []
    for color, L in zip(colors, L_list):
        m = arr[:, 0] == L
        rows = arr[m]
        ws = rows[:, 1]
        ratio_F = rows[:, 2]
        absDr_F_raw = rows[:, 3]
        absOm_F_raw = rows[:, 4]
        # fig3g uses K_meas = 100; τ_relax = max(L^2, 1/D_r)/D_r at fixed D_r.
        T_use = reconstruct_T_use(L, np.full_like(ws, D_r_fixed),
                                  K_meas=100.0, T_floor=1e8)
        absDr_F = absDr_F_raw / T_use
        absOm_F = absOm_F_raw / T_use

        # CONVENTION: current tcrw_fig3g.f90 writes L_paper directly.
        # Pass L straight through.  If you see a systematic offset, your
        # tcrw_fig3g_summary.txt is stale (older source wrote L_cur) —
        # rebuild and rerun `./tcrw_fig3g` first.
        L_authors = L

        ratio_x = np.zeros_like(ws)
        absDr_x = np.zeros_like(ws)
        absOm_x = np.zeros_like(ws)
        print(f"  L_file = L_authors = {L}: scanning {len(ws)} ω points...", end="", flush=True)
        t0 = time.time()
        for i, w in enumerate(ws):
            _, J_Dr, J_om, _ = ex.steady_state_and_currents(float(w), D_r_fixed, L_authors)
            t = ex.left_wall_J_totals(J_Dr, J_om, L_authors)
            ratio_x[i] = t["ratio"]
            absDr_x[i] = t["abs_Dr"]
            absOm_x[i] = t["abs_om"]
        print(f"  cpu = {time.time()-t0:.1f} s")

        r_max, _ = rel_error(ratio_F, ratio_x)
        d_max, _ = rel_error(absDr_F, absDr_x)
        o_max, _ = rel_error(absOm_F, absOm_x)
        summary.append((L, r_max, d_max, o_max))

        ax_ratio.scatter(ws, ratio_F, color=color, s=14, alpha=0.7, label=f"L={L} (MC)")
        ax_ratio.plot(ws, ratio_x, color=color, lw=1.0, alpha=0.9)
        ax_Jdr.scatter(ws, absDr_F, color=color, s=14, alpha=0.7)
        ax_Jdr.plot(ws, absDr_x, color=color, lw=1.0, alpha=0.9)
        ax_Jom.scatter(ws, absOm_F, color=color, s=14, alpha=0.7)
        ax_Jom.plot(ws, absOm_x, color=color, lw=1.0, alpha=0.9)

    for ax, ylab in [(ax_ratio, r"$|J_{Dr}|_{\rm wall}/|J_\omega|_{\rm wall}$"),
                     (ax_Jdr, r"$|J_{Dr}|_{\rm wall}/T$"),
                     (ax_Jom, r"$|J_\omega|_{\rm wall}/T$")]:
        ax.set_xlabel(r"$\omega$"); ax.set_ylabel(ylab)
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
    ax_ratio.legend(fontsize=9, loc="best")
    fig.suptitle(rf"Fig 3(g) cross-check — Fortran MC (markers) vs exact (lines), $D_r$ = {D_r_fixed}")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    print(f"  saved {plot_path}")

    print("\nrelative-error summary (max over ω per L):")
    print(f"  {'L':>3} {'ratio max':>10} {'|J_Dr| max':>11} {'|J_om| max':>11}")
    for L, r, d, o in summary:
        print(f"  {L:>3} {r:>10.2e} {d:>11.2e} {o:>11.2e}")


if __name__ == "__main__":
    crosscheck()
