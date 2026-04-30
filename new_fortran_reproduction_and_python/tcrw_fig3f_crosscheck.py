"""
Fig 3(f) cross-check: P_edge/P_bulk per-site vs ω at fixed D_r=1e-3, multiple L.

Loads tcrw_fig3f_summary.txt and overlays the exact-Python prediction.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os, sys, time
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tcrw_fig3_exact as ex
from _fig3_crosscheck_common import load_summary, rel_error


def crosscheck(summary_path: str = "tcrw_fig3f_summary.txt",
               plot_path: str = "tcrw_fig3f_crosscheck.png",
               D_r_fixed: float = 1e-3):
    print(f"loading {summary_path}")
    arr = load_summary(summary_path)
    # cols:  L  ω  ratio  P_edge  P_bulk  n_edge  n_bulk
    L_list = sorted(int(x) for x in np.unique(arr[:, 0]))
    print(f"  L values: {L_list}, D_r = {D_r_fixed}")

    fig, (ax_edge, ax_bulk, ax_ratio) = plt.subplots(1, 3, figsize=(15, 4.5))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(L_list)))

    summary = []
    for color, L in zip(colors, L_list):
        m = arr[:, 0] == L
        rows = arr[m]
        ws = rows[:, 1]
        ratio_F = rows[:, 2]
        Pe_F = rows[:, 3]
        Pb_F = rows[:, 4]

        # CONVENTION: current tcrw_fig3f.f90 writes L_paper directly
        # (matches authors' TRW.build_sparse_transition_matrix L = max-index).
        # Pass L straight through.  If you see a systematic ~10 % offset,
        # your tcrw_fig3f_summary.txt is stale (older source wrote L_cur =
        # L_paper + 1) — rebuild and rerun `./tcrw_fig3f` first.
        L_authors = L

        Pe_x = np.zeros_like(ws)
        Pb_x = np.zeros_like(ws)
        print(f"  L_file = L_authors = {L}: scanning {len(ws)} ω points...", end="", flush=True)
        t0 = time.time()
        for i, w in enumerate(ws):
            pi, *_ = ex.steady_state_and_currents(float(w), D_r_fixed, L_authors)
            Pe_x[i], Pb_x[i], _, _ = ex.edge_bulk_per_site(pi, L_authors)
        print(f"  cpu = {time.time()-t0:.1f} s")

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_x = Pe_x / Pb_x

        e_max, _ = rel_error(Pe_F, Pe_x)
        b_max, _ = rel_error(Pb_F, Pb_x)
        r_max, _ = rel_error(ratio_F, ratio_x)
        summary.append((L, e_max, b_max, r_max))

        ax_edge.scatter(ws, Pe_F, color=color, s=14, alpha=0.7, label=f"L={L} (MC)")
        ax_edge.plot(ws, Pe_x, color=color, lw=1.0, alpha=0.9)
        ax_bulk.scatter(ws, Pb_F, color=color, s=14, alpha=0.7)
        ax_bulk.plot(ws, Pb_x, color=color, lw=1.0, alpha=0.9)
        ax_ratio.scatter(ws, ratio_F, color=color, s=14, alpha=0.7)
        ax_ratio.plot(ws, ratio_x, color=color, lw=1.0, alpha=0.9)

    for ax, ylab in [(ax_edge, r"$\langle P\rangle_{\rm edge}$"),
                     (ax_bulk, r"$\langle P\rangle_{\rm bulk}$"),
                     (ax_ratio, r"ratio $= P_{\rm edge}/P_{\rm bulk}$")]:
        ax.set_xlabel(r"$\omega$"); ax.set_ylabel(ylab)
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
    ax_edge.legend(fontsize=9, loc="best")
    fig.suptitle(rf"Fig 3(f) cross-check — Fortran MC (markers) vs exact (lines), $D_r$ = {D_r_fixed}")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    print(f"  saved {plot_path}")

    print("\nrelative-error summary (max over ω per L):")
    print(f"  {'L':>3} {'edge max':>10} {'bulk max':>10} {'ratio max':>10}")
    for L, e, b, r in summary:
        print(f"  {L:>3} {e:>10.2e} {b:>10.2e} {r:>10.2e}")


if __name__ == "__main__":
    crosscheck()
