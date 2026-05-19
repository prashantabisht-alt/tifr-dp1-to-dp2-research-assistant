"""Exact real-space propagator for the triangular chiral RTW.

This script takes the verified Bloch generator M(k) from
``triangular_chiral_rtw.py`` and performs the inverse Fourier transform

    P(n, t) = L^{-2} sum_k exp(-i k.n) 1^T exp[M(k) t] p0.

It is the next layer after the diffusion-tensor checks: instead of only
extracting the long-time coefficients, we plot the full probability
distribution in real space.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

from triangular_chiral_rtw import (
    N_DIR,
    build_Mk_chiral_rtw,
    rates_from_gamma_bias,
)


def parse_biases(text: str) -> list[float]:
    """Parse a comma-separated list such as ``0,0.6,1``."""
    biases = [float(part.strip()) for part in text.split(",") if part.strip()]
    if not biases:
        raise ValueError("at least one bias value is required")
    for bias in biases:
        if not -1.0 <= bias <= 1.0:
            raise ValueError("all bias values must lie in [-1, 1]")
    return biases


def initial_director_distribution(initial: str) -> np.ndarray:
    """Return the initial internal-state probability vector p0."""
    if initial == "uniform":
        return np.ones(N_DIR, dtype=complex) / N_DIR

    try:
        director = int(initial)
    except ValueError as exc:
        raise ValueError("initial must be 'uniform' or an integer director 0..5") from exc

    if not 0 <= director < N_DIR:
        raise ValueError("initial director must be in 0..5")

    p0 = np.zeros(N_DIR, dtype=complex)
    p0[director] = 1.0
    return p0


def centered_axis(length: int) -> np.ndarray:
    """Coordinate labels after numpy.fft.fftshift for a periodic axis."""
    return np.arange(length) - length // 2


def axial_to_cartesian(n1: np.ndarray, n2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Map axial triangular-lattice coordinates to Cartesian plot coordinates."""
    x = n1 + 0.5 * n2
    y = 0.5 * np.sqrt(3.0) * n2
    return x, y


def exact_orientation_probabilities(
    *,
    L: int,
    time: float,
    v: float,
    gamma: float,
    bias: float,
    gamma_r: float,
    initial: str,
) -> np.ndarray:
    """Compute orientation-resolved probabilities on an L x L axial torus.

    Returns an array of shape ``(6, L, L)`` with axes
    ``(director, n1, n2)``.  The origin is at array index ``(0, 0)`` before
    centering.
    """
    if L <= 0:
        raise ValueError("L must be positive")
    if time < 0.0:
        raise ValueError("time must be nonnegative")

    gamma_plus, gamma_minus = rates_from_gamma_bias(gamma, bias)
    p0 = initial_director_distribution(initial)

    tilde = np.zeros((N_DIR, L, L), dtype=complex)
    two_pi_over_L = 2.0 * np.pi / L

    for q1 in range(L):
        k1 = two_pi_over_L * q1
        for q2 in range(L):
            k2 = two_pi_over_L * q2
            M = build_Mk_chiral_rtw(
                k1,
                k2,
                v=v,
                gamma_plus=gamma_plus,
                gamma_minus=gamma_minus,
                gamma_r=gamma_r,
            )
            tilde[:, q1, q2] = expm(M * time) @ p0

    probs = np.fft.fft2(tilde, axes=(1, 2)).real / (L * L)
    probs[np.abs(probs) < 1e-15] = 0.0
    return probs


def total_probability(probs_by_dir: np.ndarray) -> np.ndarray:
    """Sum over directors to get scalar spatial probability."""
    return probs_by_dir.sum(axis=0)


def centered_probability(prob: np.ndarray) -> np.ndarray:
    """Move the origin to the center of the plotting array."""
    return np.fft.fftshift(prob)


def probability_moments(prob_centered: np.ndarray) -> dict[str, float]:
    """Compute simple Cartesian moments for a centered probability grid."""
    L = prob_centered.shape[0]
    axis = centered_axis(L)
    n1, n2 = np.meshgrid(axis, axis, indexing="ij")
    x, y = axial_to_cartesian(n1, n2)

    norm = float(prob_centered.sum())
    mean_x = float((prob_centered * x).sum() / norm)
    mean_y = float((prob_centered * y).sum() / norm)
    mean_r2 = float((prob_centered * (x**2 + y**2)).sum() / norm)

    return {
        "norm": norm,
        "mean_x": mean_x,
        "mean_y": mean_y,
        "mean_r2": mean_r2,
        "min_probability": float(prob_centered.min()),
        "max_probability": float(prob_centered.max()),
    }


def save_probability_table(path: Path, prob_centered: np.ndarray) -> None:
    """Save a gnuplot-friendly table with n1, n2, x, y, probability."""
    L = prob_centered.shape[0]
    axis = centered_axis(L)
    n1, n2 = np.meshgrid(axis, axis, indexing="ij")
    x, y = axial_to_cartesian(n1, n2)
    data = np.column_stack(
        [
            n1.ravel(),
            n2.ravel(),
            x.ravel(),
            y.ravel(),
            prob_centered.ravel(),
        ]
    )
    header = "n1 n2 x y probability"
    np.savetxt(path, data, header=header)


def plot_probability_panels(
    *,
    panels: list[tuple[float, np.ndarray]],
    path: Path,
    time: float,
    v: float,
    gamma: float,
    gamma_r: float,
    initial: str,
) -> None:
    """Plot scalar probability panels for several chirality biases."""
    n_panels = len(panels)
    fig, axes = plt.subplots(
        1,
        n_panels,
        figsize=(4.25 * n_panels, 4.4),
        constrained_layout=True,
        squeeze=False,
    )
    axes = axes[0]

    vmax = max(float(prob.max()) for _, prob in panels)
    last_scatter = None

    for ax, (bias, prob) in zip(axes, panels):
        L = prob.shape[0]
        axis = centered_axis(L)
        n1, n2 = np.meshgrid(axis, axis, indexing="ij")
        x, y = axial_to_cartesian(n1, n2)

        last_scatter = ax.scatter(
            x.ravel(),
            y.ravel(),
            c=prob.ravel(),
            s=7,
            cmap="viridis",
            vmin=0.0,
            vmax=vmax,
            linewidths=0.0,
        )
        ax.set_aspect("equal")
        ax.set_title(rf"$b={bias:g}$")
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        ax.grid(False)

    title = (
        rf"Triangular chiral RTW exact propagator: "
        rf"$t={time:g}$, $v={v:g}$, $\gamma={gamma:g}$, "
        rf"$\gamma_r={gamma_r:g}$, initial={initial}"
    )
    fig.suptitle(title, fontsize=12)
    fig.colorbar(last_scatter, ax=axes, shrink=0.82, label="probability")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run_self_checks() -> None:
    """Small exact-propagator checks that should never fail."""
    probs0 = exact_orientation_probabilities(
        L=17,
        time=0.0,
        v=1.0,
        gamma=0.8,
        bias=0.4,
        gamma_r=0.2,
        initial="0",
    )
    total0 = total_probability(probs0)
    assert abs(total0.sum() - 1.0) < 1e-12
    assert abs(total0[0, 0] - 1.0) < 1e-12
    assert np.count_nonzero(np.abs(total0) > 1e-12) == 1

    probs = exact_orientation_probabilities(
        L=31,
        time=5.0,
        v=1.0,
        gamma=0.8,
        bias=0.4,
        gamma_r=0.2,
        initial="uniform",
    )
    total = total_probability(probs)
    assert abs(total.sum() - 1.0) < 1e-11
    assert total.min() > -1e-11


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Exact Fourier propagator for triangular chiral RTW."
    )
    parser.add_argument("--L", type=int, default=121, help="finite axial torus size")
    parser.add_argument("--time", type=float, default=30.0, help="propagation time")
    parser.add_argument("--v", type=float, default=1.0, help="run/hop rate")
    parser.add_argument("--gamma", type=float, default=0.6, help="total +/- turn rate")
    parser.add_argument("--gamma-r", type=float, default=0.0, help="reversal rate")
    parser.add_argument(
        "--biases",
        type=str,
        default="0,0.6,1.0",
        help="comma-separated chirality biases",
    )
    parser.add_argument(
        "--initial",
        type=str,
        default="0",
        help="'uniform' or an initial director index 0..5",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default="chiral_rtw_exact_propagator",
        help="base name for output files inside triangular/outputs",
    )
    parser.add_argument(
        "--self-check",
        action="store_true",
        help="run exact-propagator self-checks before generating outputs",
    )
    args = parser.parse_args()

    if args.self_check:
        run_self_checks()
        print("self-checks passed")

    biases = parse_biases(args.biases)
    out_dir = Path(__file__).resolve().parent / "outputs"
    out_dir.mkdir(exist_ok=True)

    panels: list[tuple[float, np.ndarray]] = []

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
        panels.append((bias, prob))

        suffix = f"b{bias:+.3f}".replace("+", "p").replace("-", "m").replace(".", "p")
        table_path = out_dir / f"{args.output_prefix}_{suffix}.txt"
        save_probability_table(table_path, prob)
        moments = probability_moments(prob)

        print(f"bias = {bias:g}")
        print(f"  wrote {table_path}")
        print(
            "  "
            f"norm={moments['norm']:.12f}, "
            f"min={moments['min_probability']:.3e}, "
            f"max={moments['max_probability']:.3e}, "
            f"<x>={moments['mean_x']:.6f}, "
            f"<y>={moments['mean_y']:.6f}, "
            f"<r^2>={moments['mean_r2']:.6f}"
        )

    fig_path = out_dir / f"{args.output_prefix}.png"
    plot_probability_panels(
        panels=panels,
        path=fig_path,
        time=args.time,
        v=args.v,
        gamma=args.gamma,
        gamma_r=args.gamma_r,
        initial=args.initial,
    )
    print(f"wrote {fig_path}")


if __name__ == "__main__":
    main()
