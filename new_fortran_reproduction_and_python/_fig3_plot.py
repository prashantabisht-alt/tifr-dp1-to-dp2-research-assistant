"""
Shared plotting helper for tcrw_fig3_authors.py and tcrw_fig3_pymc.py.

The plot only needs a `data` dict with the panel scan results; it does
not import or depend on authors' TRW.py.  Both the authors and pymc
front-ends populate the same dict shape and call `make_figure(data)`.

Panels (c)(d)(h)(i) are quiver plots matching paper Fig 3:
  - y axis labelled "Left Edge", showing y ∈ {1, …, 8}
  - x axis = D_r (panels c, d) or ω (panels h, i)
  - arrow direction = (Jx, Jy) at that (x, y) cell
  - arrow colour    = log10 |J| via the colorbar
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# y-range to display on the "Left Edge" axis (paper convention)
_Y_MIN = 1
_Y_MAX = 8


def _per_y_quiver(ax, x_axis, per_site_list, mag_key: str, dx_key: str,
                  dy_key: str, x_label: str, title: str,
                  xscale: str = "linear",
                  vmin_pow: float = -7, vmax_pow: float = -3,
                  length_mode: str = "fixed",
                  arrow_len_x: float = None,
                  arrow_len_y: float = 0.40):
    """
    Quiver plot matching paper Fig 3(c)(d)(h)(i).

    Paper convention (Osat et al. Fig 3):
      - arrow DIRECTION = (Jx, Jy) / |J|        (orientation glyph)
      - arrow COLOR     = log10|J|              (diverging blue-white-red)
      - arrow LENGTH    = mostly fixed, EXCEPT panel (h) where |J_Dr| genuinely
                          vanishes near ω = 0.5 → arrows shrink to dots there.
      - inner-wall sites y = 1..8 (skip both corner-adjacent sites)

    `length_mode` selects how arrow length is computed:
        "fixed"        : every nonzero J cell gets the full arrow_len_y / arrow_len_x.
                         Use for panels (c), (d), (i) — orientation-only plots.
        "linear_mag"   : length proportional to |J| / J_max within the panel.
                         Use for panel (h) — captures the genuine zero at ω = 0.5.
        "log_mag"      : length proportional to (log|J| - vmin_pow) / range.
                         Available but not paper-faithful for any panel.

    Geometry defaults (paper-faithful):
        arrow_len_y = 0.40 (in wall-site units)
        arrow_len_x = None ⇒ derived from grid spacing:
            log-x panels (c, d): 0.50 * dx_grid ≈ 0.20 decade-fractions
            linear-x panels (h, i): 0.50 * dx_grid ≈ 0.025 ω-units
    """
    n_x = len(x_axis)
    ys = np.arange(_Y_MIN, _Y_MAX + 1)         # 1..8 inclusive

    # Stack into 2-D arrays indexed [iy, ix]
    Jx = np.zeros((len(ys), n_x))
    Jy = np.zeros((len(ys), n_x))
    Jmag = np.zeros((len(ys), n_x))
    for ix, rec in enumerate(per_site_list):
        for iy, y in enumerate(ys):
            Jx[iy, ix]   = rec[dx_key][y]
            Jy[iy, ix]   = rec[dy_key][y]
            Jmag[iy, ix] = rec[mag_key][y]

    log_J = np.where(Jmag > 0, np.log10(Jmag), np.nan)
    norm = mcolors.Normalize(vmin=vmin_pow, vmax=vmax_pow)
    cmap = plt.cm.RdBu_r                # diverging palette like paper

    # ---- length-scaling factor s ∈ [s_floor, 1]  --------------------------
    if length_mode == "fixed":
        # All nonzero arrows full length; vanishing cells get a dot
        s_floor = 0.05
        s = np.where(Jmag > 0, 1.0, s_floor)
    elif length_mode == "linear_mag":
        # length ∝ |J| / J_max(panel)  (panel h convention)
        s_floor = 0.05
        Jmax = float(np.nanmax(Jmag)) if np.any(Jmag > 0) else 1.0
        s = s_floor + (1.0 - s_floor) * (Jmag / max(Jmax, 1e-30))
        s = np.where(Jmag > 0, s, s_floor)
    elif length_mode == "log_mag":
        s_floor = 0.05
        log_clamp = np.clip(np.log10(np.maximum(Jmag, 1e-30)),
                            vmin_pow, vmax_pow)
        s = s_floor + (1.0 - s_floor) * (log_clamp - vmin_pow) \
            / max(vmax_pow - vmin_pow, 1e-12)
        s = np.where(Jmag > 0, s, s_floor)
    else:
        raise ValueError(f"unknown length_mode {length_mode!r}; "
                         "expected 'fixed' / 'linear_mag' / 'log_mag'")

    eps = 1e-30
    Jx_dir = (Jx / (Jmag + eps)) * s
    Jy_dir = (Jy / (Jmag + eps)) * s

    # Place arrows on x-grid in display coordinates so log-axes are flat.
    if xscale == "log":
        x_pos = np.log10(x_axis)
    else:
        x_pos = x_axis
    dx_grid = x_pos[1] - x_pos[0]
    Xq, Yq = np.meshgrid(x_pos, ys)

    if arrow_len_x is None:
        # Log-x panels (c/d): dx_grid ≈ 0.166 over 25 D_r points; 0.50*dx_grid
        # is too compressed and the +π/4 diagonal in (d) looks too vertical.
        # Paper/old-gnuplot geometry uses arrow_len_x ≈ 0.20 in decade units
        # for c/d and ≈ 0.025 in ω units for h/i (= 0.50 * dx_grid for the
        # ω-grid).  Pick the right default per axis type.
        arrow_len_x = 0.20 if xscale == "log" else 0.50 * dx_grid

    # Slim arrow style (closer to paper Fig 3): thin shaft, small sharp head
    Q = ax.quiver(Xq, Yq,
                  Jx_dir * arrow_len_x, Jy_dir * arrow_len_y,
                  log_J, cmap=cmap, norm=norm,
                  angles="xy", scale_units="xy", scale=1.0,
                  pivot="middle", width=0.004,
                  headwidth=4.5, headlength=5.5, headaxislength=4.5)

    if xscale == "log":
        decades = np.arange(int(np.floor(x_pos.min())),
                            int(np.ceil(x_pos.max())) + 1)
        ax.set_xticks(decades)
        ax.set_xticklabels([rf"$10^{{{d}}}$" for d in decades])
    ax.set_xlim(x_pos[0] - 0.6 * dx_grid, x_pos[-1] + 0.6 * dx_grid)
    ax.set_ylim(_Y_MIN - 0.6, _Y_MAX + 0.6)
    ax.set_xlabel(x_label); ax.set_ylabel("Left Edge")
    ax.set_yticks([_Y_MIN, _Y_MAX])
    ax.set_title(title, fontsize=10)
    return Q


def make_figure(data: dict, D_r_grid, omega_grid,
                L_LIST_a, L_LIST_f,
                savepath: str, title_extra: str = ""):
    fig, axs = plt.subplots(2, 5, figsize=(22, 8.5),
                            gridspec_kw={"width_ratios": [1, 1, 1.05, 1.05, 1]})
    plt.subplots_adjust(wspace=0.35, hspace=0.35)
    cm_a = plt.cm.viridis(np.linspace(0.15, 0.85, len(L_LIST_a)))
    cm_f = plt.cm.viridis(np.linspace(0.15, 0.85, len(L_LIST_f)))

    # (a) P ratio vs D_r
    ax = axs[0, 0]
    for c, L in zip(cm_a, L_LIST_a):
        rec = data["a"][L]
        ax.plot(rec[:, 0], rec[:, 3], color=c, lw=1.2, label=f"L={L}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$D_r$"); ax.set_ylabel(r"$P_{\rm edge}/P_{\rm bulk}$")
    ax.set_title("(a)", fontsize=11, loc="left")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8, loc="best")

    # (b) |J_Dr|/|J_ω| vs D_r
    ax = axs[0, 1]
    for c, L in zip(cm_a, L_LIST_a):
        rec = data["b"][L]
        ax.plot(rec[:, 0], rec[:, 1], color=c, lw=1.2, label=f"L={L}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$D_r$"); ax.set_ylabel(r"$|J_{Dr}|/|J_\omega|$")
    ax.set_title("(b)", fontsize=11, loc="left")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8, loc="best")

    # (c) per-y J_Dr vs D_r — orientation glyphs, color = log|J|.
    # Paper Fig 3(c) cb [-5, -3].  Length is FIXED (arrows visible at all D_r).
    Q_c = _per_y_quiver(axs[0, 2], D_r_grid, data["cde_per"],
                        mag_key="absJ_Dr", dx_key="Jx_Dr", dy_key="Jy_Dr",
                        x_label=r"$D_r$",
                        title=r"(c) $\vec{J}_{Dr}$ on left edge",
                        xscale="log", vmin_pow=-5, vmax_pow=-3,
                        length_mode="fixed")
    fig.colorbar(Q_c, ax=axs[0, 2], fraction=0.045, pad=0.03,
                 label=r"$\log_{10}|J_{Dr}|$")

    # (d) per-y J_ω vs D_r — orientation glyphs (steady +π/4 across all D_r).
    # Paper cb [-5, -3].  Length FIXED — magnitude is in color only.
    Q_d = _per_y_quiver(axs[0, 3], D_r_grid, data["cde_per"],
                        mag_key="absJ_om", dx_key="Jx_om", dy_key="Jy_om",
                        x_label=r"$D_r$",
                        title=r"(d) $\vec{J}_\omega$ on left edge",
                        xscale="log", vmin_pow=-5, vmax_pow=-3,
                        length_mode="fixed")
    fig.colorbar(Q_d, ax=axs[0, 3], fraction=0.045, pad=0.03,
                 label=r"$\log_{10}|J_\omega|$")

    # (e) θ_J_Dr ONLY vs D_r (paper convention).  Paper y-range is the
    # tighter [−π/2, 0] since θ_J_Dr never goes positive in this scan.
    # Interior-only sum (y=1..L-1) drops corner Jx contamination.
    ax = axs[0, 4]
    th_Dr = np.array([t["th_Dr_tot_inner"] for t in data["cde_tot"]])
    ax.plot(D_r_grid, th_Dr, color="C0", lw=1.6, label=r"$\theta_{J_{Dr}}$")
    for h in (-np.pi/2, -np.pi/4, 0):
        ax.axhline(h, color="k", ls=":", lw=0.5, alpha=0.4)
    ax.set_xscale("log")
    ax.set_xlabel(r"$D_r$"); ax.set_ylabel(r"$\theta_{J_{Dr}}$  (rad)")
    ax.set_yticks([-np.pi/2, -np.pi/4, 0])
    ax.set_yticklabels([r"$-\pi/2$", r"$-\pi/4$", "0"])
    ax.set_ylim(-np.pi/2 - 0.05, 0.05)
    ax.set_title("(e)", fontsize=11, loc="left")
    ax.grid(True, which="both", alpha=0.3)

    # (f) P ratio vs ω — log-y on full [10¹, 10³] range so trivial L-shifts
    # don't get visually exaggerated by tight auto-scaling (paper convention)
    ax = axs[1, 0]
    for c, L in zip(cm_f, L_LIST_f):
        rec = data["f"][L]
        ax.plot(rec[:, 0], rec[:, 3], color=c, lw=1.2, label=f"L={L}")
    ax.set_yscale("log")
    ax.set_ylim(1e1, 1e3)
    ax.set_xlabel(r"$\omega$"); ax.set_ylabel(r"$P_{\rm edge}/P_{\rm bulk}$")
    ax.set_title("(f)", fontsize=11, loc="left")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8, loc="best")

    # (g) |J_Dr|/|J_ω| vs ω
    ax = axs[1, 1]
    for c, L in zip(cm_f, L_LIST_f):
        rec = data["g"][L]
        ax.plot(rec[:, 0], rec[:, 1], color=c, lw=1.2, label=f"L={L}")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\omega$"); ax.set_ylabel(r"$|J_{Dr}|/|J_\omega|$")
    ax.set_title("(g)", fontsize=11, loc="left")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8, loc="best")

    # (h) per-y J_Dr vs ω — SPECIAL CASE: |J_Dr| genuinely → 0 at ω=0.5,
    # so length really should shrink there (paper shows tiny dots near
    # achiral ω).  Use linear_mag scaling for length; cb [-7, -5].
    Q_h = _per_y_quiver(axs[1, 2], omega_grid, data["hij_per"],
                        mag_key="absJ_Dr", dx_key="Jx_Dr", dy_key="Jy_Dr",
                        x_label=r"$\omega$",
                        title=r"(h) $\vec{J}_{Dr}$ on left edge",
                        xscale="linear", vmin_pow=-7, vmax_pow=-5,
                        length_mode="linear_mag")
    fig.colorbar(Q_h, ax=axs[1, 2], fraction=0.045, pad=0.03,
                 label=r"$\log_{10}|J_{Dr}|$")

    # (i) per-y J_ω vs ω — orientation rotates -π/4 → +π/4 with magnitude
    # nearly constant (cb is the very narrow [-4.65, -4.50]).  Length FIXED
    # so the rotation is visible at every ω.
    Q_i = _per_y_quiver(axs[1, 3], omega_grid, data["hij_per"],
                        mag_key="absJ_om", dx_key="Jx_om", dy_key="Jy_om",
                        x_label=r"$\omega$",
                        title=r"(i) $\vec{J}_\omega$ on left edge",
                        xscale="linear", vmin_pow=-4.65, vmax_pow=-4.50,
                        length_mode="fixed")
    fig.colorbar(Q_i, ax=axs[1, 3], fraction=0.045, pad=0.03,
                 label=r"$\log_{10}|J_\omega|$")

    # (j) θ totals vs ω — paper convention: TWO y-axes:
    #   left  axis: θ_{J_{Dr}} ∈ [-π/2, π/2]   (sharp ±π/2 step at ω=0.5)
    #   right axis: θ_{J_ω}    ∈ [-π/4, π/4]   (linear sweep)
    # Both use INTERIOR-only sums (corner sites contaminate Jx, see
    # tcrw_fig3_exact.left_wall_J_totals docstring).
    ax = axs[1, 4]
    th_Dr = np.array([t["th_Dr_tot_inner"] for t in data["hij_tot"]])
    th_om = np.array([t["th_om_tot_inner"] for t in data["hij_tot"]])

    l1, = ax.plot(omega_grid, th_Dr, color="C0", lw=1.6, label=r"$\theta_{J_{Dr}}$")
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"$\theta_{J_{Dr}}$ (rad)", color="C0")
    ax.set_yticks([-np.pi/2, -np.pi/4, 0, np.pi/4, np.pi/2])
    ax.set_yticklabels([r"$-\pi/2$", r"$-\pi/4$", "0", r"$\pi/4$", r"$\pi/2$"])
    ax.set_ylim(-np.pi/2 - 0.05, np.pi/2 + 0.05)
    ax.tick_params(axis="y", labelcolor="C0")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_title("(j)", fontsize=11, loc="left")

    ax2 = ax.twinx()
    l2, = ax2.plot(omega_grid, th_om, color="C3", lw=1.6, label=r"$\theta_{J_\omega}$")
    ax2.set_ylabel(r"$\theta_{J_\omega}$ (rad)", color="C3")
    ax2.set_yticks([-np.pi/4, 0, np.pi/4])
    ax2.set_yticklabels([r"$-\pi/4$", "0", r"$\pi/4$"])
    ax2.set_ylim(-np.pi/4 - 0.05, np.pi/4 + 0.05)
    ax2.tick_params(axis="y", labelcolor="C3")
    ax.legend(handles=[l1, l2], fontsize=8, loc="center right")

    fig.suptitle(f"TCRW Fig 3  —  exact transition-matrix data  {title_extra}",
                 fontsize=13)
    plt.savefig(savepath, dpi=200, bbox_inches="tight")
    print(f"  saved {savepath}")
