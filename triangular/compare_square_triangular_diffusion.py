"""Compare square and triangular chiral run-and-tumble walkers.

This is the first "what is new on triangular?" diagnostic.

We compare two lattice RTW models on equal footing:

- square: 4 directors, nearest-neighbor square steps
- triangular: 6 directors, nearest-neighbor triangular steps

Both use Poisson run events at rate v and rotation rates

    gamma_plus  = gamma (1 + b) / 2
    gamma_minus = gamma (1 - b) / 2
    gamma_r     = reversal rate

The plotted quantities are:

- D_even: ordinary isotropic diffusion coefficient
- D_odd: antisymmetric / chiral diffusion response

This is not yet the full paper result, but it is the first clean transport
comparison between C4 and C6 orientation spaces.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from triangular_chiral_rtw import rates_from_gamma_bias


def square_deltas() -> np.ndarray:
    """Four square-lattice directions: right, up, left, down."""
    return np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [-1.0, 0.0],
            [0.0, -1.0],
        ]
    )


def triangular_deltas_cartesian() -> np.ndarray:
    """Six triangular directions in Cartesian coordinates."""
    sqrt3 = np.sqrt(3.0)
    return np.array(
        [
            [1.0, 0.0],
            [0.5, 0.5 * sqrt3],
            [-0.5, 0.5 * sqrt3],
            [-1.0, 0.0],
            [-0.5, -0.5 * sqrt3],
            [0.5, -0.5 * sqrt3],
        ]
    )


def rotation_generator(
    n_dir: int,
    gamma_plus: float,
    gamma_minus: float,
    gamma_r: float,
) -> np.ndarray:
    """Pure director-tumbling generator with column-source convention."""
    R = np.zeros((n_dir, n_dir), dtype=float)
    total_out = gamma_plus + gamma_minus + gamma_r
    reversal_shift = n_dir // 2

    for m in range(n_dir):
        R[m, m] -= total_out
        R[(m + 1) % n_dir, m] += gamma_plus
        R[(m - 1) % n_dir, m] += gamma_minus
        R[(m + reversal_shift) % n_dir, m] += gamma_r

    return R


def diffusion_green_kubo_lattice(
    deltas_cart: np.ndarray,
    *,
    v: float,
    gamma_plus: float,
    gamma_minus: float,
    gamma_r: float,
) -> dict[str, float]:
    """Full diffusion tensor for a Poisson-jump chiral RTW.

    The total diffusion has two pieces:

        D_total = D_persistent + D_jump.

    D_persistent is the integrated velocity-correlation part caused by memory
    in the director. D_jump is the local Poisson shot-noise contribution from
    finite jumps.
    """
    n_dir = deltas_cart.shape[0]
    R = rotation_generator(n_dir, gamma_plus, gamma_minus, gamma_r)
    R_pinv = np.linalg.pinv(R)

    D_persistent = -(v**2 / n_dir) * (deltas_cart.T @ R_pinv.T @ deltas_cart)
    D_jump = 0.5 * v * (deltas_cart.T @ deltas_cart) / n_dir
    D = D_persistent + D_jump

    D_even = 0.5 * (D[0, 0] + D[1, 1])
    D_odd = 0.5 * (D[0, 1] - D[1, 0])

    return {
        "D_even": float(D_even),
        "D_odd": float(D_odd),
        "D_xx": float(D[0, 0]),
        "D_yy": float(D[1, 1]),
        "D_xy": float(D[0, 1]),
        "D_yx": float(D[1, 0]),
    }


def compute_curves(
    *,
    v: float = 1.0,
    gamma: float = 0.8,
    gamma_r: float = 0.2,
    n_bias: int = 101,
) -> dict[str, np.ndarray]:
    """Compute square/triangular diffusion curves versus chirality bias."""
    biases = np.linspace(-0.95, 0.95, n_bias)
    square = square_deltas()
    triangular = triangular_deltas_cartesian()

    out = {
        "bias": biases,
        "square_even": np.zeros_like(biases),
        "square_odd": np.zeros_like(biases),
        "tri_even": np.zeros_like(biases),
        "tri_odd": np.zeros_like(biases),
    }

    for i, bias in enumerate(biases):
        gamma_plus, gamma_minus = rates_from_gamma_bias(gamma, float(bias))
        common = dict(
            v=v,
            gamma_plus=gamma_plus,
            gamma_minus=gamma_minus,
            gamma_r=gamma_r,
        )
        d_square = diffusion_green_kubo_lattice(square, **common)
        d_tri = diffusion_green_kubo_lattice(triangular, **common)

        out["square_even"][i] = d_square["D_even"]
        out["square_odd"][i] = d_square["D_odd"]
        out["tri_even"][i] = d_tri["D_even"]
        out["tri_odd"][i] = d_tri["D_odd"]

    return out


def save_data(curves: dict[str, np.ndarray], path: Path) -> None:
    """Save comparison data as a gnuplot-friendly text table."""
    header = "bias square_even square_odd triangular_even triangular_odd"
    data = np.column_stack(
        [
            curves["bias"],
            curves["square_even"],
            curves["square_odd"],
            curves["tri_even"],
            curves["tri_odd"],
        ]
    )
    np.savetxt(path, data, header=header)


def plot_curves(curves: dict[str, np.ndarray], path: Path) -> None:
    """Save a compact two-panel comparison plot."""
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), constrained_layout=True)

    ax = axes[0]
    ax.plot(curves["bias"], curves["square_even"], label="square C4", lw=2)
    ax.plot(curves["bias"], curves["tri_even"], label="triangular C6", lw=2)
    ax.set_xlabel("chirality bias b")
    ax.set_ylabel(r"$D_{\rm even}$")
    ax.set_title("ordinary diffusion")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    ax = axes[1]
    ax.plot(curves["bias"], curves["square_odd"], label="square C4", lw=2)
    ax.plot(curves["bias"], curves["tri_odd"], label="triangular C6", lw=2)
    ax.axhline(0.0, color="0.2", lw=0.8)
    ax.set_xlabel("chirality bias b")
    ax.set_ylabel(r"$D_{\rm odd}$")
    ax.set_title("chiral response")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    curves = compute_curves()

    out_dir = Path(__file__).resolve().parent / "outputs"
    out_dir.mkdir(exist_ok=True)
    data_path = out_dir / "square_vs_triangular_diffusion.txt"
    plot_path = out_dir / "square_vs_triangular_diffusion.png"

    save_data(curves, data_path)
    plot_curves(curves, plot_path)

    for bias in [0.0, 0.35, 0.8]:
        idx = int(np.argmin(np.abs(curves["bias"] - bias)))
        print(f"bias = {curves['bias'][idx]:.3f}")
        print(
            "  square:     "
            f"D_even={curves['square_even'][idx]:.8f}, "
            f"D_odd={curves['square_odd'][idx]:.8f}"
        )
        print(
            "  triangular: "
            f"D_even={curves['tri_even'][idx]:.8f}, "
            f"D_odd={curves['tri_odd'][idx]:.8f}"
        )

    print()
    print(f"wrote {data_path}")
    print(f"wrote {plot_path}")


if __name__ == "__main__":
    main()
