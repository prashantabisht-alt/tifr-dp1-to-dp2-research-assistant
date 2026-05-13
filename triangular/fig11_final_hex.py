"""
Final polished Fig. 11 replacement.

This script keeps the corrected theory/KMC comparison, but fixes the visual
geometry. The probability distribution is displayed in minimum-image
triangular Cartesian coordinates, so the sixfold/hexagonal structure is visible
instead of being hidden by a rolled rectangular or skew-parallelogram view.

Outputs:
    fig11_final_hex.png
    fig11_final_hex.pdf
    fig11_final_hex_evidence.png
    fig11_final_hex_evidence.pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from triangular_jmvr_corrected import (
    L_DEFAULT,
    T_DEFAULT,
    build_Mk_corrected,
    build_Mk_dipanjan,
    load_kmc_counts,
    theory_probability_on_lattice,
)


HERE = Path(__file__).resolve().parent
L = L_DEFAULT
t_final = T_DEFAULT
gamma = 0.01
epsilon = 0.15


def load_fortran_kmc(path: Path) -> tuple[np.ndarray, int]:
    """Load the Fortran KMC count file into counts[n2,n1]."""
    return load_kmc_counts(path, L=L)


def theory_P_on_lattice(*, fix_c3: bool) -> np.ndarray:
    """Exact theory on the LxL torus using corrected or historical c3 sign."""
    builder = build_Mk_corrected if fix_c3 else build_Mk_dipanjan
    return theory_probability_on_lattice(
        builder,
        gamma,
        epsilon,
        L=L,
        t_final=t_final,
    )


def min_image_coordinates(L: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Return X[n2,n1], Y[n2,n1] using the closest periodic image to the origin.

    The triangular display basis is

        a1 = (2, 0)
        a2 = (1, sqrt(3))

    so all six nearest-neighbour jumps have equal Cartesian length 2.
    """
    X = np.zeros((L, L), dtype=float)
    Y = np.zeros((L, L), dtype=float)
    sqrt3 = np.sqrt(3.0)

    for n2 in range(L):
        for n1 in range(L):
            best_x = 0.0
            best_y = 0.0
            best_r2 = np.inf
            for dn1 in (-1, 0, 1):
                for dn2 in (-1, 0, 1):
                    nn1 = n1 + dn1 * L
                    nn2 = n2 + dn2 * L
                    x = 2.0 * nn1 + nn2
                    y = sqrt3 * nn2
                    r2 = x * x + y * y
                    if r2 < best_r2:
                        best_x = x
                        best_y = y
                        best_r2 = r2
            X[n2, n1] = best_x
            Y[n2, n1] = best_y
    return X, Y


def periodic_tiled_lattice(
    L: int,
    x_range: tuple[float, float] = (-32.0, 32.0),
    y_range: tuple[float, float] = (-32.0, 32.0),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    For every lattice site (n1, n2) in {0,...,L-1}^2, return ALL its
    periodic-image Cartesian positions that fall inside the visible
    rectangle [x_range] x [y_range].

    The torus lattice vectors are:
        T1 = L * a1 = (2L, 0)
        T2 = L * a2 = (L, sqrt(3) L)

    Returns flat arrays (xs, ys, n1s, n2s) of all (image_position, lattice_index)
    pairs that fall in the window. This tiles the lattice across the entire
    rectangle so every pixel of the figure has a real lattice site — matching
    the rendering style of the Confinement-draft Fig 11.
    """
    sqrt3 = np.sqrt(3.0)
    T1x, T1y = 2.0 * L, 0.0
    T2x, T2y = float(L), sqrt3 * L

    # How many image-copies to consider in each direction. With visible
    # range of ~64 and |T| ~ 60, we need at most |i| <= 2 to cover edges.
    n_images = 3

    xs, ys, n1s, n2s = [], [], [], []
    for n1 in range(L):
        for n2 in range(L):
            x_base = 2.0 * n1 + n2
            y_base = sqrt3 * n2
            for i1 in range(-n_images, n_images + 1):
                for i2 in range(-n_images, n_images + 1):
                    x = x_base + i1 * T1x + i2 * T2x
                    y = y_base + i1 * T1y + i2 * T2y
                    if x_range[0] <= x <= x_range[1] and y_range[0] <= y <= y_range[1]:
                        xs.append(x)
                        ys.append(y)
                        n1s.append(n1)
                        n2s.append(n2)
    return (
        np.array(xs, dtype=float),
        np.array(ys, dtype=float),
        np.array(n1s, dtype=int),
        np.array(n2s, dtype=int),
    )


def unwrapped_coordinates(L: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Return X[n2,n1], Y[n2,n1] for the straight unwrapped torus parallelogram:

        X = 2 a n1 + a n2,   Y = b n2,   a = 1, b = sqrt(3)

    This is the same convention used in the draft's Fig 11: a single
    unwrapped copy of the L x L torus, displayed as a parallelogram of
    hex pixels. To centre the walker (which started at the lattice
    origin) inside this parallelogram, roll P by (L/2, L/2) before
    plotting.
    """
    sqrt3 = np.sqrt(3.0)
    n1_grid, n2_grid = np.meshgrid(np.arange(L), np.arange(L), indexing="xy")
    # n2_grid varies along axis 0, n1_grid along axis 1, matching P[n2, n1].
    X = 2.0 * n1_grid + n2_grid
    Y = sqrt3 * n2_grid
    return X, Y


def horizontal_cross_section(P: np.ndarray, X: np.ndarray, Y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return P along the minimum-image horizontal line y=0, sorted by x."""
    mask = np.isclose(Y, 0.0)
    x = X[mask]
    p = P[mask]
    order = np.argsort(x)
    return x[order], p[order]


def plot_hex_panel(ax, xs, ys, values, title, cmap, vmax, x_range, y_range):
    """Hex-tile heatmap with periodic-image tiling.

    xs, ys, values are flat arrays where each entry is one lattice-site image
    inside the visible rectangle. The rectangle is fully tiled — every pixel
    is a real (n1, n2) under PBC, so corners decay to the colormap's low
    end naturally (no white space).
    """
    sc = ax.scatter(
        xs,
        ys,
        c=values,
        s=110,                # tile size matched to lattice spacing of 2 in Cartesian
        marker="h",
        cmap=cmap,
        vmin=0.0,
        vmax=vmax,
        edgecolors="none",
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(x_range[0], x_range[1])
    ax.set_ylim(y_range[0], y_range[1])
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_title(title, fontsize=12)
    ax.tick_params(direction="in", top=True, right=True, colors="black")
    return sc


def main() -> None:
    counts, n_walkers = load_fortran_kmc(HERE / "kmc_triangular_counts.txt")
    P_kmc = counts / n_walkers

    print(f"Loaded KMC counts: N = {n_walkers:,}")
    print("Computing corrected theory...")
    P_fix = theory_P_on_lattice(fix_c3=True)
    print("Computing buggy theory...")
    P_bug = theory_P_on_lattice(fix_c3=False)

    rms_bug = float(np.sqrt(np.mean((P_kmc - P_bug) ** 2)))
    rms_fix = float(np.sqrt(np.mean((P_kmc - P_fix) ** 2)))
    mc_noise = float(np.sqrt(P_kmc.max() / n_walkers))

    # Tile the L x L lattice periodically across a square Cartesian window.
    # Every pixel of the rectangle is a real lattice site under PBC; corners
    # naturally fall to low P (rendered as the colormap's darkest colour).
    # Rectangle dimensions: must be smaller than half the fundamental-domain
    # extent in each direction, otherwise PBC wraps high-P "ring" values
    # back into the corners. Fundamental-domain half-extents are
    #   x_half ≈ |T1|/2 = L = 30,    y_half ≈ |T2_y|/2 = sqrt(3) L / 2 ≈ 26.
    # We use 28 in x, 25 in y → leaves a small margin so corners are firmly
    # in the depleted "far side" of the torus.
    x_range = (-28.0, 28.0)
    y_range = (-25.0, 25.0)
    xs_tile, ys_tile, n1s_tile, n2s_tile = periodic_tiled_lattice(L, x_range, y_range)
    P_fix_tile = P_fix[n2s_tile, n1s_tile]
    P_kmc_tile = P_kmc[n2s_tile, n1s_tile]
    P_bug_tile = P_bug[n2s_tile, n1s_tile]
    vmax = max(float(P_fix_tile.max()), float(P_kmc_tile.max()))

    # Cross-section: use the min-image lattice for a clean horizontal slice
    # at y = 0 through the walker (which is at origin).
    X_mi, Y_mi = min_image_coordinates(L)
    x_cross, kmc_cross = horizontal_cross_section(P_kmc, X_mi, Y_mi)
    _, bug_cross = horizontal_cross_section(P_bug, X_mi, Y_mi)
    _, fix_cross = horizontal_cross_section(P_fix, X_mi, Y_mi)

    draft_palette = LinearSegmentedColormap.from_list(
        "draft_black_purple_red_yellow",
        [
            (0.00, "#050505"),
            (0.22, "#32105f"),
            (0.42, "#b42bd6"),
            (0.60, "#ed1c24"),
            (0.78, "#ff8c00"),
            (1.00, "#fff200"),
        ],
    )

    # Residual P_fix - P_kmc, in units of one standard deviation of MC noise.
    # Per-site MC noise σ_i = sqrt(P_kmc * (1-P_kmc) / N) ≈ sqrt(P_kmc/N) for small P.
    # A z-score map |fix-kmc|/σ_i shows whether residuals are within statistical noise.
    sigma_local = np.sqrt(np.maximum(P_kmc, 1.0 / n_walkers) / n_walkers)
    residual_z = (P_fix - P_kmc) / sigma_local
    residual_z_tile = residual_z[n2s_tile, n1s_tile]
    z_max = float(np.percentile(np.abs(residual_z), 99))

    # Main PI-facing figure: clean, not crowded.
    fig = plt.figure(figsize=(13.5, 8.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.18, 0.86], hspace=0.34, wspace=0.23)

    ax_theory = fig.add_subplot(gs[0, 0])
    ax_kmc = fig.add_subplot(gs[0, 1])
    ax_cross = fig.add_subplot(gs[1, :])

    sc0 = plot_hex_panel(
        ax_theory,
        xs_tile,
        ys_tile,
        P_fix_tile,
        rf"(a) Corrected exact theory  ($\gamma={gamma}$, $\epsilon={epsilon}$, $t={t_final:g}$)",
        draft_palette,
        vmax,
        x_range,
        y_range,
    )
    sc1 = plot_hex_panel(
        ax_kmc,
        xs_tile,
        ys_tile,
        P_kmc_tile,
        rf"(b) Kinetic Monte Carlo  ($N={n_walkers:,}$)",
        draft_palette,
        vmax,
        x_range,
        y_range,
    )

    cb0 = fig.colorbar(sc0, ax=ax_theory, fraction=0.046, pad=0.03)
    cb0.set_label("Probability")
    cb1 = fig.colorbar(sc1, ax=ax_kmc, fraction=0.046, pad=0.03)
    cb1.set_label("Probability")

    ax_cross.plot(x_cross, kmc_cross, "ko", ms=5.5, label=f"KMC (N={n_walkers:,})", zorder=3)
    ax_cross.plot(
        x_cross,
        bug_cross,
        color="#d62728",
        lw=2.2,
        label=rf"buggy theory  (RMS $={rms_bug:.2e}$)",
    )
    ax_cross.plot(
        x_cross,
        fix_cross,
        color="#0057ff",
        lw=2.2,
        label=rf"corrected theory  (RMS $={rms_fix:.2e}$)",
    )
    ax_cross.set_title(r"(c) Horizontal cross-section through the walker", fontsize=12)
    ax_cross.set_xlabel(r"$x$ at $y=0$  (minimum-image coordinates)")
    ax_cross.set_ylabel(r"$P(x,y=0,t)$")
    ax_cross.grid(alpha=0.25)
    ax_cross.tick_params(direction="in", top=True, right=True)
    ax_cross.legend(loc="upper center", ncol=3, framealpha=0.95)

    fig.suptitle(
        "Triangular active walker: corrected theory matches KMC",
        fontsize=13,
        y=0.98,
    )
    fig.text(
        0.5,
        0.025,
        rf"Minimum-image triangular display.  MC noise floor $\approx {mc_noise:.2e}$; "
        rf"corrected theory is at the noise floor.  The old $c_3$ sign gives the red curve.",
        ha="center",
        fontsize=10.5,
    )

    fig.savefig(HERE / "fig11_final_hex.png", dpi=220, bbox_inches="tight")
    fig.savefig(HERE / "fig11_final_hex.pdf", bbox_inches="tight")
    plt.close(fig)

    # Evidence/backup figure: includes residual z-score panel.
    fig = plt.figure(figsize=(18.5, 9.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 0.9], hspace=0.34, wspace=0.28)

    ax_theory = fig.add_subplot(gs[0, 0])
    ax_kmc = fig.add_subplot(gs[0, 1])
    ax_res = fig.add_subplot(gs[0, 2])
    ax_cross = fig.add_subplot(gs[1, :])

    sc0 = plot_hex_panel(
        ax_theory,
        xs_tile,
        ys_tile,
        P_fix_tile,
        rf"(a) Corrected exact theory  ($\gamma={gamma}$, $\epsilon={epsilon}$, $t={t_final:g}$)",
        draft_palette,
        vmax,
        x_range,
        y_range,
    )
    sc1 = plot_hex_panel(
        ax_kmc,
        xs_tile,
        ys_tile,
        P_kmc_tile,
        rf"(b) Kinetic Monte Carlo  ($N={n_walkers:,}$)",
        draft_palette,
        vmax,
        x_range,
        y_range,
    )

    # Residual panel: signed z-score, diverging colormap centered at 0.
    # If theory matches KMC within MC noise, |z| should be O(1) everywhere
    # and the panel will look like pure noise — that IS the success criterion.
    sc2 = ax_res.scatter(
        xs_tile,
        ys_tile,
        c=residual_z_tile,
        s=110,
        marker="h",
        cmap="RdBu_r",
        vmin=-z_max,
        vmax=+z_max,
        edgecolors="none",
    )
    ax_res.set_aspect("equal", adjustable="box")
    ax_res.set_xlim(x_range[0], x_range[1])
    ax_res.set_ylim(y_range[0], y_range[1])
    ax_res.set_xlabel(r"$x$")
    ax_res.set_ylabel(r"$y$")
    ax_res.set_title(
        rf"(c) Residual $(P_{{\rm theory}}-P_{{\rm KMC}})/\sigma_{{\rm MC}}$",
        fontsize=12,
    )
    ax_res.tick_params(direction="in", top=True, right=True)

    cb0 = fig.colorbar(sc0, ax=ax_theory, fraction=0.046, pad=0.03)
    cb0.set_label("Probability")
    cb1 = fig.colorbar(sc1, ax=ax_kmc, fraction=0.046, pad=0.03)
    cb1.set_label("Probability")
    cb2 = fig.colorbar(sc2, ax=ax_res, fraction=0.046, pad=0.03)
    cb2.set_label(r"z-score  $\Delta P / \sigma_{\rm MC}$")

    ax_cross.plot(x_cross, kmc_cross, "ko", ms=5.5, label=f"KMC (N={n_walkers:,})", zorder=3)
    ax_cross.plot(
        x_cross,
        bug_cross,
        color="#d62728",
        lw=2.2,
        label=rf"buggy theory  (RMS $={rms_bug:.2e}$)",
    )
    ax_cross.plot(
        x_cross,
        fix_cross,
        color="#0057ff",
        lw=2.2,
        label=rf"corrected theory  (RMS $={rms_fix:.2e}$)",
    )
    ax_cross.set_title(r"(d) Horizontal cross-section through the walker", fontsize=12)
    ax_cross.set_xlabel(r"$x$ at $y=0$  (minimum-image coordinates)")
    ax_cross.set_ylabel(r"$P(x,y=0,t)$")
    ax_cross.grid(alpha=0.25)
    ax_cross.tick_params(direction="in", top=True, right=True)
    ax_cross.legend(loc="upper center", ncol=3, framealpha=0.95)

    fig.suptitle(
        "Triangular active walker: corrected theory matches KMC; residuals are noise-like",
        fontsize=13,
        y=0.98,
    )
    fig.text(
        0.5,
        0.025,
        rf"MC noise floor $\approx {mc_noise:.2e}$; corrected theory is at the noise floor. "
        rf"The old $c_3$ sign gives the red curve.",
        ha="center",
        fontsize=10.5,
    )

    fig.savefig(HERE / "fig11_final_hex_evidence.png", dpi=220, bbox_inches="tight")
    fig.savefig(HERE / "fig11_final_hex_evidence.pdf", bbox_inches="tight")
    plt.close(fig)

    print("Saved:")
    print(f"  {HERE / 'fig11_final_hex.png'}")
    print(f"  {HERE / 'fig11_final_hex.pdf'}")
    print(f"  {HERE / 'fig11_final_hex_evidence.png'}")
    print(f"  {HERE / 'fig11_final_hex_evidence.pdf'}")
    print(f"RMS buggy  vs KMC = {rms_bug:.3e}")
    print(f"RMS fixed  vs KMC = {rms_fix:.3e}")
    print(f"MC noise floor    = {mc_noise:.3e}")


if __name__ == "__main__":
    main()
