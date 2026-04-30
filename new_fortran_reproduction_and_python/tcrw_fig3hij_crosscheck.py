"""
Fig 3(h)(i)(j) cross-check
==========================

3(h) and 3(i) — quiver of left-wall current vs (ω, y) for J_Dr (3h) and
J_omega (3i) at fixed D_r = 1e-3, L = 10.
3(j) — angle of the TOTAL left-wall J_Dr and J_omega vs ω.

Loads:
  tcrw_fig3hij_summary.txt   (per-y per-ω, 11 sites × 21 ω = 231 rows)
  tcrw_fig3j_summary.txt     (totals per ω, 21 rows)

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os, sys, time
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tcrw_fig3_exact as ex
from _fig3_crosscheck_common import load_summary, reconstruct_T_use, rel_error, angle_diff


def crosscheck(per_site_path: str = "tcrw_fig3hij_summary.txt",
               totals_path:    str = "tcrw_fig3j_summary.txt",
               plot_path:      str = "tcrw_fig3hij_crosscheck.png",
               L_paper: int = 10, D_r_fixed: float = 1e-3):
    print(f"loading {per_site_path}")
    rec = load_summary(per_site_path)
    # cols:  iW  ω  y  Jx_Dr  Jy_Dr  Jx_om  Jy_om  |J_Dr|  th_Dr  |J_om|  th_om
    ws = np.unique(rec[:, 1])
    L_cur = L_paper + 1
    n_w = len(ws)
    print(f"  ω range: {ws.min():.3f} .. {ws.max():.3f}  ({n_w} pts)")
    print(f"  L_paper = {L_paper}, L_cur = {L_cur}, D_r = {D_r_fixed}")

    print(f"loading {totals_path}")
    tot = load_summary(totals_path)
    # cols:  iW  ω  Jx_Dr_tot  Jy_Dr_tot  Jx_om_tot  Jy_om_tot  th_Dr_tot  th_om_tot
    ws_t = tot[:, 1]
    th_Dr_F = tot[:, 6]
    th_om_F = tot[:, 7]

    K_meas = 300.0     # fig3hij uses K_meas = 300
    T_use = reconstruct_T_use(L_cur, np.full_like(ws_t, D_r_fixed),
                              K_meas=K_meas, T_floor=1e8)

    per_x = []
    tot_x = []
    print(f"  computing exact at {n_w} ω points...", end="", flush=True)
    t0 = time.time()
    for w in ws:
        _, J_Dr, J_om, _ = ex.steady_state_and_currents(float(w), D_r_fixed, L_paper)
        per_x.append(ex.left_wall_J_per_y(J_Dr, J_om, L_paper))
        tot_x.append(ex.left_wall_J_totals(J_Dr, J_om, L_paper))
    print(f"  cpu = {time.time()-t0:.1f} s")

    th_Dr_x = np.array([t["th_Dr_tot"] for t in tot_x])
    th_om_x = np.array([t["th_om_tot"] for t in tot_x])

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    (ax_qDr, ax_qOm), (ax_th, ax_resid) = axes

    ys = np.arange(L_cur)
    cmap = plt.cm.viridis
    for iy, y in enumerate(ys):
        absJ_Dr_F = np.zeros(n_w); absJ_om_F = np.zeros(n_w)
        absJ_Dr_x = np.zeros(n_w); absJ_om_x = np.zeros(n_w)
        for iw, w in enumerate(ws):
            mask = (rec[:, 1] == w) & (rec[:, 2] == y)
            if mask.any():
                row = rec[mask][0]
                Jx_Dr_raw, Jy_Dr_raw = row[3], row[4]
                Jx_om_raw, Jy_om_raw = row[5], row[6]
                T = reconstruct_T_use(L_cur, D_r_fixed, K_meas=K_meas, T_floor=1e8)
                absJ_Dr_F[iw] = np.hypot(Jx_Dr_raw, Jy_Dr_raw) / T
                absJ_om_F[iw] = np.hypot(Jx_om_raw, Jy_om_raw) / T
            absJ_Dr_x[iw] = per_x[iw]["absJ_Dr"][y]
            absJ_om_x[iw] = per_x[iw]["absJ_om"][y]
        c = cmap(iy / max(1, L_cur - 1))
        ax_qDr.plot(ws, absJ_Dr_x, color=c, lw=0.8, alpha=0.9)
        ax_qDr.scatter(ws, absJ_Dr_F, color=c, s=8, alpha=0.7,
                       label=(f"y={y}" if iy in (0, L_cur//2, L_cur-1) else None))
        ax_qOm.plot(ws, absJ_om_x, color=c, lw=0.8, alpha=0.9)
        ax_qOm.scatter(ws, absJ_om_F, color=c, s=8, alpha=0.7)

    for ax, ylab in [(ax_qDr, r"$|J_{Dr}(0,y)|/T$"),
                     (ax_qOm, r"$|J_\omega(0,y)|/T$")]:
        ax.set_yscale("log")
        ax.set_xlabel(r"$\omega$"); ax.set_ylabel(ylab)
        ax.grid(True, which="both", alpha=0.3)
    ax_qDr.legend(fontsize=8, loc="best", title="left-wall y")
    ax_qDr.set_title("Fig 3(h) — magnitude per y (Fortran scatter, exact lines)")
    ax_qOm.set_title("Fig 3(i) — magnitude per y")

    # Panel (j): θ_J_Dr_tot, θ_J_om_tot vs ω
    ax_th.scatter(ws_t, th_Dr_F, color="C0", s=20, alpha=0.7, label=r"$\theta_{J_{Dr}}$ (MC)")
    ax_th.plot   (ws_t, th_Dr_x, color="C0", lw=1.0, alpha=0.9)
    ax_th.scatter(ws_t, th_om_F, color="C3", s=20, alpha=0.7, label=r"$\theta_{J_\omega}$ (MC)")
    ax_th.plot   (ws_t, th_om_x, color="C3", lw=1.0, alpha=0.9)
    for h in (-np.pi/2, -np.pi/4, 0, np.pi/4, np.pi/2):
        ax_th.axhline(h, color="k", ls=":", lw=0.5, alpha=0.4)
    ax_th.set_xlabel(r"$\omega$"); ax_th.set_ylabel(r"$\theta$  (rad)")
    ax_th.set_yticks([-np.pi/2, -np.pi/4, 0, np.pi/4, np.pi/2])
    ax_th.set_yticklabels([r"$-\pi/2$", r"$-\pi/4$", "0", r"$\pi/4$", r"$\pi/2$"])
    ax_th.set_title("Fig 3(j) — total angles vs ω")
    ax_th.legend(fontsize=9, loc="best")
    ax_th.grid(True, which="both", alpha=0.3)

    dth_Dr = angle_diff(th_Dr_F, th_Dr_x)
    dth_om = angle_diff(th_om_F, th_om_x)
    ax_resid.scatter(ws_t, dth_Dr, color="C0", s=14, alpha=0.7, label=r"$\Delta\theta_{J_{Dr}}$")
    ax_resid.scatter(ws_t, dth_om, color="C3", s=14, alpha=0.7, label=r"$\Delta\theta_{J_\omega}$")
    ax_resid.axhline(0, color="k", lw=0.5, alpha=0.5)
    ax_resid.set_xlabel(r"$\omega$"); ax_resid.set_ylabel(r"MC $-$ exact  (rad)")
    ax_resid.set_title("Angle residuals (MC − exact)")
    ax_resid.legend(fontsize=9, loc="best")
    ax_resid.grid(True, which="both", alpha=0.3)

    fig.suptitle(rf"Fig 3(h)(i)(j) cross-check — Fortran (markers) vs exact (lines), L={L_paper}, $D_r$={D_r_fixed}")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    print(f"  saved {plot_path}")

    print("\nangle residual summary (max over ω):")
    print(f"  Δθ_J_Dr  max = {np.max(np.abs(dth_Dr)):.2e} rad")
    print(f"  Δθ_J_om  max = {np.max(np.abs(dth_om)):.2e} rad")


if __name__ == "__main__":
    crosscheck()
