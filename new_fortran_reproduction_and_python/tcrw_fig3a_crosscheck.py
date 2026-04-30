"""
Fig 3(a) cross-check: P_edge/P_bulk per-site vs D_r at ω=1, multiple L.

Loads tcrw_fig3a_summary.txt, computes the same observable from authors'
exact transition matrix via tcrw_fig3_exact, and overlays.

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


def crosscheck(summary_path: str = "tcrw_fig3a_summary.txt",
               plot_path: str = "tcrw_fig3a_crosscheck.png"):
    print(f"loading {summary_path}")
    arr = load_summary(summary_path)
    # columns:  L  D_r  ratio  P_edge_norm  P_bulk_norm  n_edge  n_bulk
    L_list = sorted(int(x) for x in np.unique(arr[:, 0]))
    print(f"  L values: {L_list}")

    fig, (ax_edge, ax_bulk, ax_ratio) = plt.subplots(1, 3, figsize=(15, 4.5))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(L_list)))

    summary = []
    for color, L in zip(colors, L_list):
        m = arr[:, 0] == L
        rows = arr[m]
        Drs = rows[:, 1]
        ratio_F = rows[:, 2]
        Pe_F = rows[:, 3]
        Pb_F = rows[:, 4]

        Pe_x = np.zeros_like(Drs)
        Pb_x = np.zeros_like(Drs)
        print(f"  L = {L}: scanning {len(Drs)} D_r points...", end="", flush=True)
        t0 = time.time()
        for i, D in enumerate(Drs):
            pi, *_ = ex.steady_state_and_currents(1.0, float(D), L)
            Pe_x[i], Pb_x[i], _, _ = ex.edge_bulk_per_site(pi, L)
        print(f"  cpu = {time.time()-t0:.1f} s")

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_x = Pe_x / Pb_x

        # Skip the D_r = 1 row — pure-noise limit has a degenerate steady state
        keep = Drs < 0.95
        e_max, e_mean = rel_error(Pe_F, Pe_x, keep=keep)
        b_max, b_mean = rel_error(Pb_F, Pb_x, keep=keep)
        r_max, r_mean = rel_error(ratio_F, ratio_x, keep=keep)
        summary.append((L, e_max, b_max, r_max))

        ax_edge.scatter(Drs, Pe_F, color=color, s=14, alpha=0.7, label=f"L={L} (MC)")
        ax_edge.plot(Drs, Pe_x, color=color, lw=1.0, alpha=0.9)
        ax_bulk.scatter(Drs, Pb_F, color=color, s=14, alpha=0.7)
        ax_bulk.plot(Drs, Pb_x, color=color, lw=1.0, alpha=0.9)
        ax_ratio.scatter(Drs, ratio_F, color=color, s=14, alpha=0.7)
        ax_ratio.plot(Drs, ratio_x, color=color, lw=1.0, alpha=0.9)

    for ax, ylab in [(ax_edge, r"$\langle P\rangle_{\rm edge}$"),
                     (ax_bulk, r"$\langle P\rangle_{\rm bulk}$"),
                     (ax_ratio, r"ratio $= P_{\rm edge}/P_{\rm bulk}$")]:
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$D_r$"); ax.set_ylabel(ylab)
        ax.grid(True, which="both", alpha=0.3)
    ax_edge.legend(fontsize=9, loc="lower left")
    fig.suptitle("Fig 3(a) cross-check — Fortran MC (markers) vs exact Python (lines), ω = 1")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    print(f"  saved {plot_path}")

    print("\nrelative-error summary (max over D_r per L, D_r < 0.95):")
    print(f"  {'L':>3} {'edge max':>10} {'bulk max':>10} {'ratio max':>10}")
    for L, e, b, r in summary:
        print(f"  {L:>3} {e:>10.2e} {b:>10.2e} {r:>10.2e}")


if __name__ == "__main__":
    crosscheck()
