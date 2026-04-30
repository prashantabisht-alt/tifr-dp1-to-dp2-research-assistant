"""
TCRW Fig 1 — Python MC simulator (vectorized numpy)
====================================================

Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020, Fig 1.

Self-contained reproduction of all three data panels of Fig 1 entirely
in Python.  Uses the same authors' convention for direction encoding and
step rules as `tcrw_fig4c.py` and `tcrw_step.f90`.

Three panels
------------
1(b) — sample UNWRAPPED trajectories at ω ∈ {0.0, 0.5, 0.75, 1.0},
       D_r = 10⁻³, T = 10⁶ steps.

1(c) — ensemble MSD ⟨|r|²⟩(t) at ω ∈ {0.5, 0.7, 1.0}, D_r = 10⁻³,
       T = 10⁶ steps, log-log axes.

1(d) — long-time diffusion coefficient D(ω, D_r) at
       ω ∈ {0.0, 0.1, ..., 1.0},  D_r ∈ {10⁻⁴, 10⁻³, 10⁻²}.
       Fit MSD ∝ 4 D t in the diffusive regime.  Cross-checked against
       the EXACT D from the small-k Bloch expansion in
       `tcrw_fig1_bloch_diffusion.py`.

Direction / step convention (verbatim from authors' TRW.py)
-----------------------------------------------------------
    d = 0 ↑   d = 1 →   d = 2 ↓   d = 3 ←
    DX = [0, 1, 0, -1]   DY = [1, 0, -1, 0]

    Per step (independent random draws):
       prob D_r       — NOISE  step (no translation, rotate only)
                          CCW  w.p. ω,    CW  w.p. 1 − ω
       prob 1 − D_r   — CHIRAL step (translate THEN rotate)
                          CW   w.p. ω,    CCW w.p. 1 − ω

Boundary
--------
For Fig 1 we use UNBOUNDED dynamics — single walkers on an infinite
lattice.  PBC with unwrapped positions is equivalent to unbounded for
single-walker MSD; this is the standard convention for diffusion plots.

Vectorisation
-------------
N_traj walkers are simulated in parallel as numpy arrays of shape
(N_traj,).  For each MC step we draw two random arrays and update
positions / directors with boolean masking.  Throughput scales linearly
with N_traj at fixed T (essentially free up to a few thousand walkers
on a modern laptop).

Cross-check
-----------
Fig 1(d): MC D(ω, D_r) is overlaid against the exact Bloch line from
`tcrw_fig1_bloch_diffusion.D_exact` (which is verified against the
authors' real-space operator at machine precision via fig4b's torus
crosscheck).  The MC must converge to the Bloch line as T → ∞.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
for cand in (HERE,
             "/sessions/elegant-wizardly-einstein/mnt/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python",
             "/Users/prashantbisht/Documents/Claude/Projects/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python"):
    if os.path.isfile(os.path.join(cand, "tcrw_fig1_bloch_diffusion.py")):
        sys.path.insert(0, cand)
        ROOT = cand
        break

from tcrw_fig1_bloch_diffusion import D_exact   # exact reference D(ω, D_r)


DX_T = np.array([0, 1, 0, -1], dtype=np.int64)
DY_T = np.array([1, 0, -1, 0], dtype=np.int64)


# ---------------------------------------------------------------------------
# 1. Vectorized unbounded MC walker
# ---------------------------------------------------------------------------
def simulate_unbounded(omega: float, D_r: float, T_steps: int,
                       N_traj: int = 1000, seed: int = 0,
                       record_traj: bool = False,
                       traj_save_pts: int = 1000,
                       msd_log_pts: int = 80):
    """
    Run N_traj independent unbounded TCRW walkers for T_steps.

    Returns
    -------
    out : dict with keys
        't_pts'   : (n_t,) log-spaced time points where MSD is sampled
        'msd'     : (n_t,) ⟨x² + y²⟩ at each t_pts
        'final_x', 'final_y', 'final_d'
        (if record_traj) 'traj_t', 'traj_x', 'traj_y' for the FIRST walker
    """
    rng = np.random.default_rng(seed)

    x = np.zeros(N_traj, dtype=np.int64)
    y = np.zeros(N_traj, dtype=np.int64)
    d = rng.integers(0, 4, size=N_traj)

    # log-spaced sampling grid for MSD
    t_pts = np.unique(
        np.logspace(0, np.log10(T_steps), msd_log_pts).astype(np.int64)
    )
    msd = np.zeros(t_pts.shape[0])

    if record_traj:
        traj_t = np.unique(
            np.concatenate([
                np.array([0]),
                np.logspace(0, np.log10(T_steps), traj_save_pts).astype(np.int64),
            ])
        )
        traj_x = np.zeros(traj_t.shape[0], dtype=np.int64)
        traj_y = np.zeros(traj_t.shape[0], dtype=np.int64)
        traj_x[0] = 0; traj_y[0] = 0
        traj_kk = 1
    else:
        traj_t = traj_x = traj_y = None
        traj_kk = -1

    k = 0
    for it in range(1, T_steps + 1):
        r_step = rng.random(N_traj)
        r_rot  = rng.random(N_traj)
        is_noise = r_step < D_r
        is_chiral = ~is_noise

        # Noise: rotate only.  CCW w.p. ω, CW w.p. 1 − ω.
        m_ccw_n = is_noise & (r_rot <  omega)
        m_cw_n  = is_noise & (r_rot >= omega)
        d[m_ccw_n] = (d[m_ccw_n] - 1) % 4
        d[m_cw_n]  = (d[m_cw_n]  + 1) % 4

        # Chiral: translate then rotate.  CW w.p. ω, CCW w.p. 1 − ω
        if is_chiral.any():
            dc = d[is_chiral]
            x[is_chiral] += DX_T[dc]
            y[is_chiral] += DY_T[dc]
            m_cw_c  = is_chiral & (r_rot <  omega)
            m_ccw_c = is_chiral & (r_rot >= omega)
            d[m_cw_c]  = (d[m_cw_c]  + 1) % 4
            d[m_ccw_c] = (d[m_ccw_c] - 1) % 4

        # Sample MSD at log-spaced points
        if k < t_pts.shape[0] and it == t_pts[k]:
            msd[k] = float(np.mean(x.astype(np.float64) ** 2
                                   + y.astype(np.float64) ** 2))
            k += 1

        # Sample first-walker trajectory
        if record_traj and traj_kk < traj_t.shape[0] and it == traj_t[traj_kk]:
            traj_x[traj_kk] = x[0]
            traj_y[traj_kk] = y[0]
            traj_kk += 1

    out = dict(
        t_pts=t_pts[:k], msd=msd[:k],
        final_x=x, final_y=y, final_d=d,
    )
    if record_traj:
        out.update(traj_t=traj_t[:traj_kk],
                   traj_x=traj_x[:traj_kk],
                   traj_y=traj_y[:traj_kk])
    return out


# ---------------------------------------------------------------------------
# 2. D(ω, D_r) via long-time MSD slope fit
# ---------------------------------------------------------------------------
def fit_diffusion(t_pts, msd, t_min: float, t_max: float | None = None):
    """
    Direct LINEAR fit MSD = a + b · t  (not log-log) over t ∈ [t_min, t_max].
    D_fit = b / 4   in 2D.

    Linear-MSD fit is more accurate than log-log for D extraction:
      - imposes the known long-time scaling MSD ∝ t (no anomalous slope);
      - log-log fit conflates D with the slope, biases D when the slope
        deviates from 1 due to finite-T transient.

    Returns (D_fit, slope_loglog_diagnostic, slope_err, n_pts_used).
    The "slope_loglog" is computed in log-log space too, as a diagnostic
    of whether the data is genuinely diffusive (slope≈1).
    """
    if t_max is None:
        t_max = t_pts[-1]
    keep = (t_pts >= t_min) & (t_pts <= t_max) & (msd > 0)
    if keep.sum() < 3:
        return float("nan"), float("nan"), float("nan"), 0
    t_keep = t_pts[keep].astype(np.float64)
    m_keep = msd[keep]
    n = len(t_keep)

    # ---- linear MSD = a + b·t fit ----
    b, a = np.polyfit(t_keep, m_keep, 1)        # slope, intercept
    yhat = b * t_keep + a
    s2 = float(np.sum((m_keep - yhat) ** 2) / (n - 2))
    var_t = float(np.sum((t_keep - t_keep.mean()) ** 2))
    b_err = (s2 / var_t) ** 0.5
    D_fit = b / 4.0
    D_err = b_err / 4.0

    # ---- diagnostic log-log slope (should be ≈ 1 if genuinely diffusive) ----
    x = np.log10(t_keep); y = np.log10(m_keep)
    slope_ll, _ = np.polyfit(x, y, 1)

    return D_fit, float(slope_ll), float(D_err), n


# ---------------------------------------------------------------------------
# 3. Fig 1(b)  trajectories at 4 ω values
# ---------------------------------------------------------------------------
def make_fig1b(D_r=1e-3, T_steps=10**6, omega_list=(0.0, 0.5, 0.75, 1.0),
               seed=123, savepath="tcrw_fig1b_pymc.png"):
    print(f"\n[fig1b] D_r={D_r}, T={T_steps}, ω={omega_list}")
    fig, axs = plt.subplots(2, 2, figsize=(9, 9))
    cmap = mpl.colormaps.get_cmap("viridis")
    for ax, w in zip(axs.ravel(), omega_list):
        t0 = time.time()
        out = simulate_unbounded(w, D_r, T_steps, N_traj=1, seed=seed,
                                  record_traj=True, traj_save_pts=2000,
                                  msd_log_pts=2)
        dt = time.time() - t0
        tt = out["traj_t"].astype(float)
        # color by log10(t) so the dominant time is visible at all scales
        log_t = np.log10(np.maximum(tt, 1.0))
        norm = mpl.colors.Normalize(vmin=0, vmax=np.log10(T_steps))
        # plot as connected segments coloured by time
        from matplotlib.collections import LineCollection
        pts = np.stack([out["traj_x"], out["traj_y"]], axis=1)
        segs = np.stack([pts[:-1], pts[1:]], axis=1)
        lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=0.6)
        lc.set_array(log_t[:-1])
        ax.add_collection(lc)
        ax.scatter([0], [0], c="green", s=30, zorder=5, label="start")
        ax.scatter([out["traj_x"][-1]], [out["traj_y"][-1]],
                   c="red", s=30, zorder=5, label="end")
        ax.set_xlim(pts[:, 0].min() - 5, pts[:, 0].max() + 5)
        ax.set_ylim(pts[:, 1].min() - 5, pts[:, 1].max() + 5)
        ax.set_aspect("equal")
        ax.set_title(rf"$\omega = {w}$  (cpu {dt:.1f} s)", fontsize=11)
        ax.tick_params(labelsize=8)
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize(0, np.log10(T_steps)))
    cbar = fig.colorbar(sm, ax=axs.ravel().tolist(),
                        fraction=0.04, pad=0.03)
    cbar.set_label(r"$\log_{10}\, t$"); cbar.ax.tick_params(labelsize=8)
    fig.suptitle(rf"Fig 1(b) — trajectories,  $D_r = {D_r}$,  $T = {T_steps}$",
                 fontsize=12)
    plt.savefig(savepath, dpi=200, bbox_inches="tight")
    print(f"  saved {savepath}")


# ---------------------------------------------------------------------------
# 4. Fig 1(c)  ensemble MSD
# ---------------------------------------------------------------------------
def make_fig1c(D_r=1e-3, T_steps=10**6, omega_list=(0.5, 0.7, 1.0),
               N_traj=1000, seed=7, savepath="tcrw_fig1c_pymc.png"):
    print(f"\n[fig1c] D_r={D_r}, T={T_steps}, N_traj={N_traj}, ω={omega_list}")
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = plt.cm.tab10(np.arange(len(omega_list)))
    for color, w in zip(colors, omega_list):
        t0 = time.time()
        out = simulate_unbounded(w, D_r, T_steps, N_traj=N_traj, seed=seed,
                                  record_traj=False)
        dt = time.time() - t0
        ax.loglog(out["t_pts"], out["msd"], "o-", color=color, ms=3, lw=1.0,
                  label=f"ω={w}  (cpu {dt:.1f}s)")
    # MSD ∝ t guide
    t_g = np.array([1.0, 1e6])
    ax.loglog(t_g, 0.05 * t_g, "k--", lw=1.0, alpha=0.5, label=r"MSD $\propto t$")
    ax.set_xlabel("t");  ax.set_ylabel(r"$\langle |r|^2 \rangle$")
    ax.set_xlim(1, T_steps); ax.set_ylim(0.5, None)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=9)
    ax.set_title(rf"Fig 1(c) — ensemble MSD,  $D_r = {D_r}$,  $N$ = {N_traj}",
                 fontsize=11)
    plt.tight_layout(); plt.savefig(savepath, dpi=200, bbox_inches="tight")
    print(f"  saved {savepath}")


# ---------------------------------------------------------------------------
# 5. Fig 1(d)  D(ω) for several D_r, with exact Bloch overlay
# ---------------------------------------------------------------------------
def make_fig1d(D_r_list=(1e-4, 1e-3, 1e-2),
               omega_list=tuple(np.round(np.linspace(0.0, 1.0, 11), 2)),
               T_steps=10**6, N_traj=1000, seed=42,
               savepath="tcrw_fig1d_pymc.png"):
    print(f"\n[fig1d] T={T_steps}, N_traj={N_traj}")
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = ["tab:purple", "tab:pink", "tab:orange"]
    omega_dense = np.linspace(0, 1, 51)

    for color, D_r in zip(colors, D_r_list):
        t_min = max(10.0 / D_r, 100.0)        # diffusive-regime fit floor
        D_mc = np.zeros(len(omega_list))
        D_err = np.zeros(len(omega_list))
        slopes = np.zeros(len(omega_list))
        print(f"  D_r = {D_r:.0e}  (fit window t > {t_min:.0e})")
        t0 = time.time()
        for i, w in enumerate(omega_list):
            out = simulate_unbounded(float(w), float(D_r), T_steps,
                                      N_traj=N_traj, seed=seed + i)
            D_fit, slope, slope_err, _ = fit_diffusion(
                out["t_pts"], out["msd"], t_min)
            D_mc[i] = D_fit
            slopes[i] = slope
            # propagate slope_err to D_err via δ(10^a)/10^a = ln10·δa
            D_err[i] = D_fit * np.log(10) * slope_err
        print(f"    elapsed {time.time()-t0:.1f}s   slope range [{slopes.min():.3f}, {slopes.max():.3f}]")

        # exact Bloch line (dense, smooth)
        D_dense = np.array([D_exact(float(w), float(D_r)) for w in omega_dense])

        ax.errorbar(omega_list, D_mc, yerr=D_err, fmt="o", color=color,
                    capsize=3, ms=6, label=f"MC  D_r={D_r:.0e}")
        ax.plot(omega_dense, D_dense, color=color, lw=1.6, alpha=0.7,
                label=f"Bloch exact  D_r={D_r:.0e}")
    ax.set_xlabel(r"$\omega$"); ax.set_ylabel(r"$D(\omega, D_r)$")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(0, None)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.legend(fontsize=8, ncol=2)
    ax.set_title(rf"Fig 1(d) — D(ω): MC (markers) vs exact Bloch (lines)",
                 fontsize=11)
    ax.grid(alpha=0.3)
    plt.tight_layout(); plt.savefig(savepath, dpi=200, bbox_inches="tight")
    print(f"  saved {savepath}")


# ---------------------------------------------------------------------------
# 6. Sanity test (very small T) before any production run
# ---------------------------------------------------------------------------
def sanity_test():
    """
    Tiny run that catches obvious convention bugs before paper-spec.
    Predictions:
      ω=0.5, D_r=0.5 (achiral, half-noise): D ≈ 0.125 (exact = (1-D_r)/4)
      ω=0.5, D_r=0.1 (achiral): exact D ≈ 0.225
      ω=1.0, D_r=0.5 (chiral but lots of noise): D ≈ 0.125 (mid-range)
      ω=1.0, D_r=0.001 (mostly closed loops): D ≈ 0.00025 (tiny)
    """
    print("Sanity tests (T = 50_000, N_traj = 500):")
    cases = [(0.5, 0.5), (0.5, 0.1), (1.0, 0.5), (1.0, 0.001)]
    print(f"  {'ω':>5} {'D_r':>7} {'D_exact':>10} {'D_MC':>10} {'rel diff':>10}")
    for w, D_r in cases:
        D_x = D_exact(w, D_r)
        out = simulate_unbounded(w, D_r, 50_000, N_traj=500, seed=11)
        t_min = max(10.0 / D_r, 100.0)
        D_mc, slope, _, _ = fit_diffusion(out["t_pts"], out["msd"], t_min)
        rel = abs(D_mc - D_x) / max(D_x, 1e-30)
        print(f"  {w:>5.2f} {D_r:>7.0e} {D_x:>10.5f} {D_mc:>10.5f} {rel:>10.2%}")


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 64)
    print(" TCRW Fig 1 — Python MC")
    print("=" * 64)

    print("\n[0] Quick sanity test")
    sanity_test()

    # production runs (paper parameters)
    make_fig1b(D_r=1e-3, T_steps=10**6,
               savepath=os.path.join(ROOT, "tcrw_fig1b_pymc.png"))
    make_fig1c(D_r=1e-3, T_steps=10**6, N_traj=1000,
               savepath=os.path.join(ROOT, "tcrw_fig1c_pymc.png"))
    make_fig1d(D_r_list=(1e-4, 1e-3, 1e-2),
               T_steps=10**6, N_traj=1000,
               savepath=os.path.join(ROOT, "tcrw_fig1d_pymc.png"))

    print("\nDone.")
