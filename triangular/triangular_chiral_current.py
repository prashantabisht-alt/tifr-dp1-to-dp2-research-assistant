"""Current-field plots for the triangular chiral RTW exact propagator.

The scalar probability density can hide chirality, especially for a uniform
initial orientation.  This script plots the local run-current field

    J(n, t) = v sum_m P_m(n, t) Delta_m,

using the orientation-resolved exact propagator.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from triangular_chiral_propagator import (
    axial_to_cartesian,
    centered_axis,
    centered_probability,
    exact_orientation_probabilities,
    parse_biases,
    total_probability,
)
from triangular_chiral_rtw import triangular_deltas


def axial_vectors_to_cartesian(j1: np.ndarray, j2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Map vector components in the axial basis to Cartesian components."""
    jx = j1 + 0.5 * j2
    jy = 0.5 * np.sqrt(3.0) * j2
    return jx, jy


def current_field(probs_by_dir: np.ndarray, *, v: float) -> tuple[np.ndarray, np.ndarray]:
    """Return centered Cartesian current components from orientation probabilities."""
    deltas = triangular_deltas()
    j1 = v * np.tensordot(deltas[:, 0], probs_by_dir, axes=(0, 0))
    j2 = v * np.tensordot(deltas[:, 1], probs_by_dir, axes=(0, 0))
    jx, jy = axial_vectors_to_cartesian(j1, j2)
    return centered_probability(jx), centered_probability(jy)


def current_diagnostics(prob: np.ndarray, jx: np.ndarray, jy: np.ndarray) -> dict[str, float]:
    """Compute simple scalars that summarize the current field."""
    L = prob.shape[0]
    axis = centered_axis(L)
    n1, n2 = np.meshgrid(axis, axis, indexing="ij")
    x, y = axial_to_cartesian(n1, n2)
    mag = np.sqrt(jx**2 + jy**2)
    angular_current = x * jy - y * jx

    return {
        "norm": float(prob.sum()),
        "net_jx": float(jx.sum()),
        "net_jy": float(jy.sum()),
        "mean_current_magnitude": float(mag.mean()),
        "max_current_magnitude": float(mag.max()),
        "angular_current": float(angular_current.sum()),
        "density_weighted_angular_current": float((prob * angular_current).sum()),
    }


def save_current_table(path: Path, prob: np.ndarray, jx: np.ndarray, jy: np.ndarray) -> None:
    """Save a gnuplot-friendly table with density and current components."""
    L = prob.shape[0]
    axis = centered_axis(L)
    n1, n2 = np.meshgrid(axis, axis, indexing="ij")
    x, y = axial_to_cartesian(n1, n2)
    data = np.column_stack(
        [
            n1.ravel(),
            n2.ravel(),
            x.ravel(),
            y.ravel(),
            prob.ravel(),
            jx.ravel(),
            jy.ravel(),
            np.sqrt(jx.ravel() ** 2 + jy.ravel() ** 2),
        ]
    )
    header = "n1 n2 x y probability Jx Jy Jmag"
    np.savetxt(path, data, header=header)


def quiver_components(
    x: np.ndarray,
    y: np.ndarray,
    jx: np.ndarray,
    jy: np.ndarray,
    prob: np.ndarray,
    *,
    stride: int,
    density_cut: float,
    arrow_length: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Subsample and rescale current arrows for a readable quiver plot."""
    xs = x[::stride, ::stride]
    ys = y[::stride, ::stride]
    us = jx[::stride, ::stride]
    vs = jy[::stride, ::stride]
    ps = prob[::stride, ::stride]

    mag = np.sqrt(us**2 + vs**2)
    cutoff = density_cut * float(prob.max())
    mask = (ps >= cutoff) & (mag > 0.0)

    if not np.any(mask):
        return xs[mask], ys[mask], us[mask], vs[mask], mag[mask]

    ref = float(np.quantile(mag[mask], 0.90))
    if ref <= 0.0:
        ref = float(mag[mask].max())
    scale = arrow_length * np.minimum(mag / ref, 1.0) / np.maximum(mag, 1e-300)

    return xs[mask], ys[mask], us[mask] * scale[mask], vs[mask] * scale[mask], mag[mask]


def plot_current_panels(
    *,
    panels: list[tuple[float, np.ndarray, np.ndarray, np.ndarray]],
    path: Path,
    time: float,
    v: float,
    gamma: float,
    gamma_r: float,
    initial: str,
    stride: int,
    density_cut: float,
    arrow_length: float,
) -> None:
    """Plot density background plus current arrows for each bias."""
    n_panels = len(panels)
    fig, axes = plt.subplots(
        1,
        n_panels,
        figsize=(4.8 * n_panels, 4.7),
        constrained_layout=True,
        squeeze=False,
    )
    axes = axes[0]
    vmax = max(float(prob.max()) for _, prob, _, _ in panels)
    last_scatter = None

    for ax, (bias, prob, jx, jy) in zip(axes, panels):
        L = prob.shape[0]
        axis = centered_axis(L)
        n1, n2 = np.meshgrid(axis, axis, indexing="ij")
        x, y = axial_to_cartesian(n1, n2)

        last_scatter = ax.scatter(
            x.ravel(),
            y.ravel(),
            c=prob.ravel(),
            s=6,
            cmap="viridis",
            vmin=0.0,
            vmax=vmax,
            linewidths=0.0,
        )
        qx, qy, qu, qv, qmag = quiver_components(
            x,
            y,
            jx,
            jy,
            prob,
            stride=stride,
            density_cut=density_cut,
            arrow_length=arrow_length,
        )
        if len(qx) > 0:
            ax.quiver(
                qx,
                qy,
                qu,
                qv,
                qmag,
                cmap="magma",
                angles="xy",
                scale_units="xy",
                scale=1.0,
                width=0.006,
                headwidth=3.2,
                headlength=4.0,
            )

        ax.set_aspect("equal")
        ax.set_title(rf"$b={bias:g}$")
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        ax.grid(False)

    fig.suptitle(
        rf"Current field: $t={time:g}$, $v={v:g}$, $\gamma={gamma:g}$, "
        rf"$\gamma_r={gamma_r:g}$, initial={initial}",
        fontsize=12,
    )
    fig.colorbar(last_scatter, ax=axes, shrink=0.82, label="probability")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot exact RTW probability currents.")
    parser.add_argument("--L", type=int, default=91, help="finite axial torus size")
    parser.add_argument("--time", type=float, default=24.0, help="propagation time")
    parser.add_argument("--v", type=float, default=1.0, help="run/hop rate")
    parser.add_argument("--gamma", type=float, default=0.6, help="total +/- turn rate")
    parser.add_argument("--gamma-r", type=float, default=0.0, help="reversal rate")
    parser.add_argument("--biases", type=str, default="0,0.6,1.0")
    parser.add_argument("--initial", type=str, default="uniform")
    parser.add_argument("--stride", type=int, default=4, help="quiver subsampling stride")
    parser.add_argument(
        "--density-cut",
        type=float,
        default=0.08,
        help="only draw arrows where density is at least this fraction of max density",
    )
    parser.add_argument(
        "--arrow-length",
        type=float,
        default=3.0,
        help="maximum displayed arrow length in plot coordinates",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default="chiral_rtw_current",
        help="base name for outputs inside triangular/outputs",
    )
    args = parser.parse_args()

    if args.stride <= 0:
        raise ValueError("stride must be positive")

    biases = parse_biases(args.biases)
    out_dir = Path(__file__).resolve().parent / "outputs"
    out_dir.mkdir(exist_ok=True)

    panels: list[tuple[float, np.ndarray, np.ndarray, np.ndarray]] = []

    for bias in biases:
        probs_by_dir = exact_orientation_probabilities(
            L=args.L,
            time=args.time,
            v=args.v,
            gamma=args.gamma,
            bias=bias,
            gamma_r=args.gamma_r,
            initial=args.initial,
        )
        prob = centered_probability(total_probability(probs_by_dir))
        jx, jy = current_field(probs_by_dir, v=args.v)
        panels.append((bias, prob, jx, jy))

        suffix = f"b{bias:+.3f}".replace("+", "p").replace("-", "m").replace(".", "p")
        table_path = out_dir / f"{args.output_prefix}_{suffix}.txt"
        save_current_table(table_path, prob, jx, jy)
        diag = current_diagnostics(prob, jx, jy)

        print(f"bias = {bias:g}")
        print(f"  wrote {table_path}")
        print(
            "  "
            f"norm={diag['norm']:.12f}, "
            f"netJ=({diag['net_jx']:.6e}, {diag['net_jy']:.6e}), "
            f"max|J|={diag['max_current_magnitude']:.6e}, "
            f"sum(x Jy - y Jx)={diag['angular_current']:.6e}"
        )

    fig_path = out_dir / f"{args.output_prefix}.png"
    plot_current_panels(
        panels=panels,
        path=fig_path,
        time=args.time,
        v=args.v,
        gamma=args.gamma,
        gamma_r=args.gamma_r,
        initial=args.initial,
        stride=args.stride,
        density_cut=args.density_cut,
        arrow_length=args.arrow_length,
    )
    print(f"wrote {fig_path}")


if __name__ == "__main__":
    main()
