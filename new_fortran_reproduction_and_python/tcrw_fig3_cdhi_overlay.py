"""
Fig. 3(c)(d)(h)(i) vector overlay diagnostic.

This is not a replacement for the paper-style reproduction plots.  It is a
sanity-check view that overlays:

  gray arrows    = exact steady-state field from the authors' transition matrix
  colored arrows = Fortran simulation field from the fresh summary files

The two fields are offset slightly in y so that both are visible at the same
(x, y) grid point.  The colored Fortran arrows use the same convention as the
main Python Fig. 3 plot: arrow direction from (Jx, Jy), color from log10|J|,
fixed arrow length for panels c/d/i, and magnitude-scaled length for panel h.
"""
from __future__ import annotations

import os
import sys
import time
from dataclasses import dataclass

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import tcrw_fig3_exact as ex
from _fig3_crosscheck_common import angle_diff, load_summary


L_PAPER = 10
Y_MIN = 1
Y_MAX = 8
YS = np.arange(Y_MIN, Y_MAX + 1, dtype=int)


@dataclass
class FieldScan:
    x: np.ndarray
    Jx_Dr: np.ndarray
    Jy_Dr: np.ndarray
    absJ_Dr: np.ndarray
    Jx_om: np.ndarray
    Jy_om: np.ndarray
    absJ_om: np.ndarray


def _empty_scan(x_axis: np.ndarray) -> FieldScan:
    shape = (len(YS), len(x_axis))
    return FieldScan(
        x=x_axis,
        Jx_Dr=np.zeros(shape), Jy_Dr=np.zeros(shape), absJ_Dr=np.zeros(shape),
        Jx_om=np.zeros(shape), Jy_om=np.zeros(shape), absJ_om=np.zeros(shape),
    )


def _load_fortran_scan(path: str) -> FieldScan:
    """Load Fortran summary columns into a y-by-x vector field."""
    rec = load_summary(path)
    x_axis = np.unique(rec[:, 1])
    out = _empty_scan(x_axis)

    for ix, x in enumerate(x_axis):
        block = rec[np.isclose(rec[:, 1], x, rtol=1e-12, atol=1e-15)]
        for row in block:
            y = int(round(row[2]))
            if y < Y_MIN or y > Y_MAX:
                continue
            iy = y - Y_MIN
            out.Jx_Dr[iy, ix] = row[3]
            out.Jy_Dr[iy, ix] = row[4]
            out.Jx_om[iy, ix] = row[5]
            out.Jy_om[iy, ix] = row[6]

            # New Fortran summaries include normalized |J| columns.  Fall back
            # to raw |J|/T_use if an older summary is ever inspected.
            if rec.shape[1] >= 14:
                out.absJ_Dr[iy, ix] = row[12]
                out.absJ_om[iy, ix] = row[13]
            elif rec.shape[1] >= 12:
                T_use = row[11]
                out.absJ_Dr[iy, ix] = np.hypot(row[3], row[4]) / T_use
                out.absJ_om[iy, ix] = np.hypot(row[5], row[6]) / T_use
            else:
                out.absJ_Dr[iy, ix] = row[7]
                out.absJ_om[iy, ix] = row[9]
    return out


def _exact_cde_scan(D_r_grid: np.ndarray, omega: float = 1.0) -> FieldScan:
    out = _empty_scan(D_r_grid)
    for ix, D_r in enumerate(D_r_grid):
        _, J_Dr, J_om, _ = ex.steady_state_and_currents(omega, float(D_r), L_PAPER)
        per = ex.left_wall_J_per_y(J_Dr, J_om, L_PAPER)
        for iy, y in enumerate(YS):
            out.Jx_Dr[iy, ix] = per["Jx_Dr"][y]
            out.Jy_Dr[iy, ix] = per["Jy_Dr"][y]
            out.absJ_Dr[iy, ix] = per["absJ_Dr"][y]
            out.Jx_om[iy, ix] = per["Jx_om"][y]
            out.Jy_om[iy, ix] = per["Jy_om"][y]
            out.absJ_om[iy, ix] = per["absJ_om"][y]
    return out


def _exact_hij_scan(omega_grid: np.ndarray, D_r_fixed: float = 1.0e-3) -> FieldScan:
    out = _empty_scan(omega_grid)
    for ix, omega in enumerate(omega_grid):
        _, J_Dr, J_om, _ = ex.steady_state_and_currents(float(omega), D_r_fixed, L_PAPER)
        per = ex.left_wall_J_per_y(J_Dr, J_om, L_PAPER)
        for iy, y in enumerate(YS):
            out.Jx_Dr[iy, ix] = per["Jx_Dr"][y]
            out.Jy_Dr[iy, ix] = per["Jy_Dr"][y]
            out.absJ_Dr[iy, ix] = per["absJ_Dr"][y]
            out.Jx_om[iy, ix] = per["Jx_om"][y]
            out.Jy_om[iy, ix] = per["Jy_om"][y]
            out.absJ_om[iy, ix] = per["absJ_om"][y]
    return out


def _unit_vectors(Jx: np.ndarray, Jy: np.ndarray, mag_for_zero: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mag = np.hypot(Jx, Jy)
    eps = 1e-300
    U = Jx / (mag + eps)
    V = Jy / (mag + eps)
    U = np.where(mag_for_zero > 0, U, 0.0)
    V = np.where(mag_for_zero > 0, V, 0.0)
    return U, V


def _length_factor(mag: np.ndarray, mode: str, panel_max: float | None = None) -> np.ndarray:
    floor = 0.055
    if mode == "fixed":
        return np.where(mag > 0, 1.0, floor)
    if mode == "linear_mag":
        scale = max(float(panel_max or np.nanmax(mag)), 1e-300)
        return np.where(mag > 0, floor + (1.0 - floor) * mag / scale, floor)
    raise ValueError(f"unknown length mode {mode!r}")


def _components(scan: FieldScan, kind: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if kind == "Dr":
        return scan.Jx_Dr, scan.Jy_Dr, scan.absJ_Dr
    if kind == "om":
        return scan.Jx_om, scan.Jy_om, scan.absJ_om
    raise ValueError(f"unknown current kind {kind!r}")


def _plot_overlay_panel(
    ax,
    exact: FieldScan,
    mc: FieldScan,
    kind: str,
    title: str,
    x_label: str,
    xscale: str,
    vmin_pow: float,
    vmax_pow: float,
    length_mode: str,
    y_offset: float = 0.12,
    exact_length_boost: float = 1.0,
    mc_length_boost: float = 1.0,
    exact_alpha: float = 0.62,
    mc_alpha: float = 0.88,
):
    x_axis = mc.x
    if xscale == "log":
        x_pos = np.log10(x_axis)
        arrow_len_x = 0.20
    else:
        x_pos = x_axis
        arrow_len_x = 0.50 * np.median(np.diff(x_pos))
    arrow_len_y = 0.40

    X, Y = np.meshgrid(x_pos, YS)
    Jx_e, Jy_e, mag_e = _components(exact, kind)
    Jx_m, Jy_m, mag_m = _components(mc, kind)

    panel_max = max(float(np.nanmax(mag_e)), float(np.nanmax(mag_m)), 1e-300)

    Ue, Ve = _unit_vectors(Jx_e, Jy_e, mag_e)
    Um, Vm = _unit_vectors(Jx_m, Jy_m, mag_m)
    se = _length_factor(mag_e, length_mode, panel_max=panel_max)
    sm = _length_factor(mag_m, length_mode, panel_max=panel_max)

    ax.quiver(
        X, Y + y_offset,
        Ue * se * arrow_len_x * exact_length_boost,
        Ve * se * arrow_len_y * exact_length_boost,
        angles="xy", scale_units="xy", scale=1.0,
        pivot="middle", color="0.18", alpha=exact_alpha,
        width=0.0032, headwidth=4.4, headlength=5.2, headaxislength=4.4,
    )

    log_mag_m = np.where(mag_m > 0, np.log10(mag_m), np.nan)
    norm = mcolors.Normalize(vmin=vmin_pow, vmax=vmax_pow)
    q_mc = ax.quiver(
        X, Y - y_offset,
        Um * sm * arrow_len_x * mc_length_boost,
        Vm * sm * arrow_len_y * mc_length_boost,
        log_mag_m, cmap=plt.cm.RdBu_r, norm=norm,
        angles="xy", scale_units="xy", scale=1.0,
        pivot="middle", alpha=mc_alpha,
        width=0.0038, headwidth=4.6, headlength=5.4, headaxislength=4.6,
    )

    dx_grid = np.median(np.diff(x_pos))
    ax.set_xlim(x_pos[0] - 0.65 * dx_grid, x_pos[-1] + 0.65 * dx_grid)
    ax.set_ylim(Y_MIN - 0.7, Y_MAX + 0.7)
    ax.set_yticks([Y_MIN, Y_MAX])
    ax.set_ylabel("Left edge")
    ax.set_xlabel(x_label)
    ax.set_title(title, loc="left", fontsize=11)
    ax.grid(True, alpha=0.18, lw=0.5)

    if xscale == "log":
        decades = np.arange(int(np.floor(x_pos.min())), int(np.ceil(x_pos.max())) + 1)
        ax.set_xticks(decades)
        ax.set_xticklabels([rf"$10^{{{d}}}$" for d in decades])

    return q_mc


def _angle_report(
    name: str,
    exact: FieldScan,
    mc: FieldScan,
    kind: str,
    min_exact_fraction: float = 1e-12,
) -> dict[str, float]:
    Jx_e, Jy_e, mag_e = _components(exact, kind)
    Jx_m, Jy_m, mag_m = _components(mc, kind)
    th_e = np.arctan2(Jy_e, Jx_e)
    th_m = np.arctan2(Jy_m, Jx_m)

    # Ignore cells where the exact field is essentially zero; the angle is not
    # meaningful there, especially panel h at omega = 0.5.
    mag_floor = max(float(np.nanmax(mag_e)) * min_exact_fraction, 1e-300)
    valid = (mag_e > mag_floor) & (mag_m > 0)
    dtheta = np.full_like(mag_e, np.nan, dtype=float)
    dtheta[valid] = angle_diff(th_m[valid], th_e[valid])
    abs_dtheta = np.abs(dtheta[valid])

    if abs_dtheta.size:
        flat_valid = np.argwhere(valid)
        worst_local = int(np.nanargmax(abs_dtheta))
        iy, ix = flat_valid[worst_local]
        worst = {
            "max": float(abs_dtheta[worst_local]),
            "median": float(np.nanmedian(abs_dtheta)),
            "p90": float(np.nanpercentile(abs_dtheta, 90)),
            "x": float(mc.x[ix]),
            "y": float(YS[iy]),
            "n": int(abs_dtheta.size),
            "skipped": int(valid.size - abs_dtheta.size),
        }
    else:
        worst = {"max": np.nan, "median": np.nan, "p90": np.nan,
                 "x": np.nan, "y": np.nan, "n": 0, "skipped": int(valid.size)}

    threshold_note = ""
    if min_exact_fraction > 1e-12:
        threshold_note = f", exact |J|>{min_exact_fraction:.0e} max"
    print(
        f"{name:8s}: |dtheta| median={worst['median']:.3e}, "
        f"p90={worst['p90']:.3e}, max={worst['max']:.3e} rad "
        f"at x={worst['x']:.5g}, y={worst['y']:.0f} "
        f"(used {worst['n']}, skipped {worst['skipped']}{threshold_note})"
    )
    return worst


def _small_Dr_panel_d_report(exact_cde: FieldScan, mc_cde: FieldScan):
    ix = 0
    th_e = np.arctan2(exact_cde.Jy_om[:, ix], exact_cde.Jx_om[:, ix])
    th_m = np.arctan2(mc_cde.Jy_om[:, ix], mc_cde.Jx_om[:, ix])
    dth = angle_diff(th_m, th_e)
    print("\npanel d at smallest D_r:")
    print(f"  D_r = {mc_cde.x[ix]:.5e}")
    print(f"  exact theta range = {th_e.min():.6f} .. {th_e.max():.6f} rad")
    print(f"  MC theta range    = {th_m.min():.6f} .. {th_m.max():.6f} rad")
    print(f"  max |MC-exact|   = {np.max(np.abs(dth)):.3e} rad")


def make_overlay(
    cde_path: str = "tcrw_fig3cde_summary.txt",
    hij_path: str = "tcrw_fig3hij_summary.txt",
    out_path: str = "tcrw_fig3_cdhi_overlay.png",
    out_path_noshift: str = "tcrw_fig3_cdhi_overlay_noshift.png",
):
    t0 = time.time()
    print("loading Fortran summaries...")
    mc_cde = _load_fortran_scan(cde_path)
    mc_hij = _load_fortran_scan(hij_path)

    print("computing exact steady-state fields...")
    exact_cde = _exact_cde_scan(mc_cde.x, omega=1.0)
    exact_hij = _exact_hij_scan(mc_hij.x, D_r_fixed=1.0e-3)
    print(f"  exact computation cpu = {time.time() - t0:.1f} s")

    print("\nangle mismatch summary, Fortran MC minus exact steady state:")
    _angle_report("panel c", exact_cde, mc_cde, "Dr")
    _angle_report("panel d", exact_cde, mc_cde, "om")
    _angle_report("panel h", exact_hij, mc_hij, "Dr", min_exact_fraction=5e-3)
    _angle_report("panel i", exact_hij, mc_hij, "om")
    _small_Dr_panel_d_report(exact_cde, mc_cde)

    def _draw_figure(
        save_path: str,
        suptitle: str,
        y_offset: float,
        exact_length_boost: float = 1.0,
        mc_length_boost: float = 1.0,
        exact_alpha: float = 0.62,
        mc_alpha: float = 0.88,
    ):
        fig, axs = plt.subplots(2, 2, figsize=(12.8, 9.2))
        common = dict(
            y_offset=y_offset,
            exact_length_boost=exact_length_boost,
            mc_length_boost=mc_length_boost,
            exact_alpha=exact_alpha,
            mc_alpha=mc_alpha,
        )
        q_c = _plot_overlay_panel(
            axs[0, 0], exact_cde, mc_cde, "Dr",
            r"(c) $\vec{J}_{D_r}$, overlay", r"$D_r$", "log",
            vmin_pow=-5, vmax_pow=-3, length_mode="fixed", **common,
        )
        q_d = _plot_overlay_panel(
            axs[0, 1], exact_cde, mc_cde, "om",
            r"(d) $\vec{J}_{\omega}$, overlay", r"$D_r$", "log",
            vmin_pow=-5, vmax_pow=-3, length_mode="fixed", **common,
        )
        q_h = _plot_overlay_panel(
            axs[1, 0], exact_hij, mc_hij, "Dr",
            r"(h) $\vec{J}_{D_r}$, overlay", r"$\omega$", "linear",
            vmin_pow=-7, vmax_pow=-5, length_mode="linear_mag", **common,
        )
        q_i = _plot_overlay_panel(
            axs[1, 1], exact_hij, mc_hij, "om",
            r"(i) $\vec{J}_{\omega}$, overlay", r"$\omega$", "linear",
            vmin_pow=-4.65, vmax_pow=-4.50, length_mode="fixed", **common,
        )

        for ax, q, label in [
            (axs[0, 0], q_c, r"MC $\log_{10}|J_{D_r}|$"),
            (axs[0, 1], q_d, r"MC $\log_{10}|J_{\omega}|$"),
            (axs[1, 0], q_h, r"MC $\log_{10}|J_{D_r}|$"),
            (axs[1, 1], q_i, r"MC $\log_{10}|J_{\omega}|$"),
        ]:
            fig.colorbar(q, ax=ax, fraction=0.045, pad=0.025, label=label)

        fig.suptitle(suptitle, fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.965])
        fig.savefig(save_path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        print(f"\nsaved {save_path}")

    _draw_figure(
        out_path,
        "Fig. 3 cdhi overlay: gray = exact steady state (+0.12 y), "
        "colored = Fortran simulation (-0.12 y)",
        y_offset=0.12,
    )
    _draw_figure(
        out_path_noshift,
        "Fig. 3 cdhi no-shift overlay: same grid; gray exact drawn longer, "
        "colored Fortran MC drawn shorter",
        y_offset=0.0,
        exact_length_boost=1.10,
        mc_length_boost=0.86,
        exact_alpha=0.58,
        mc_alpha=0.86,
    )


if __name__ == "__main__":
    make_overlay()
