"""
Reproduce Dipanjan's Mathematica calculation, then plot it.

This script is intentionally narrow and historical:
  1. Use the same parameters and old c3 sign as rtp_tl_2.nb by default.
  2. Compute P(x, y, t) from Dipanjan's Fourier-space generator.
  3. Write the y = constant line cut that Dipanjan printed.
  4. Plot that line cut and a full real-space probability heatmap.

Run from the triangular folder or from the repo root:
    python3 triangular/reproduce_dipanjan_exact.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from triangular_jmvr_corrected import build_Mk_dipanjan


def evolve_fourier_grid(
    gamma: float,
    epsilon: float,
    time: float,
    L: int,
    x0: float = 0.0,
    y0: float = 0.0,
    a: float = 1.0,
    b: float = 1.0,
) -> np.ndarray:
    """
    Compute director-summed P_tilde(k, t) on Dipanjan's rectangular k-grid.

    Dipanjan's notebook uses
        kx = 2 pi nx / (2 a L),  nx = 0, ..., 2L - 1
        ky = 2 pi ny / (b L),    ny = 0, ..., L - 1

    At each k:
        P_tilde(k, t) = sum_d [ exp(M(k) t) P_tilde(k, 0) ]_d
    with a delta initial condition at (x0, y0) and uniform initial director.
    """
    nx_count = 2 * L
    ny_count = L
    p_k = np.zeros((nx_count, ny_count), dtype=complex)

    for nx in range(nx_count):
        kx = 2.0 * np.pi * nx / (2.0 * a * L)
        for ny in range(ny_count):
            ky = 2.0 * np.pi * ny / (b * L)
            M = build_Mk_dipanjan(gamma, epsilon, kx, ky, a=a, b=b)
            evals, evecs = np.linalg.eig(M)

            init = np.full(6, 1.0 / 6.0, dtype=complex)
            init *= np.exp(1j * (kx * x0 + ky * y0))

            coeff = np.linalg.solve(evecs, init)
            p_vec_t = evecs @ (coeff * np.exp(evals * time))
            p_k[nx, ny] = np.sum(p_vec_t)

    return p_k


def inverse_fourier_cartesian(p_k: np.ndarray) -> np.ndarray:
    """
    Inverse transform with Dipanjan's sign convention.

    Mathematica code:
        P(x,y,t) = (1/Nk) sum_k P_tilde(k,t) exp[-i(kx x + ky y)]

    numpy.fft.fft2 has the same negative sign but no 1/N prefactor, so divide
    by the number of k-points.
    """
    norm = p_k.shape[0] * p_k.shape[1]
    return np.fft.fft2(p_k) / norm


def parity_masks(L: int) -> tuple[np.ndarray, np.ndarray]:
    """Return masks for physical and off-lattice parity sites."""
    xs = np.arange(2 * L)[:, None]
    ys = np.arange(L)[None, :]
    physical = ((xs + ys) % 2) == 0
    return physical, ~physical


def write_linecut(
    P: np.ndarray,
    y: int,
    out_path: Path,
    gamma: float,
    epsilon: float,
    time: float,
    L: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Write exactly Dipanjan's style of output: x = 0, 2, ..., 2L-2."""
    xs = np.arange(0, 2 * L, 2)
    vals = P[xs, y].real
    with out_path.open("w") as f:
        f.write("# Reproduction of Dipanjan rtp_tl_2.nb line cut\n")
        f.write(f"# gamma={gamma} epsilon={epsilon} L={L} time={time} y={y}\n")
        f.write("# x P(x,y,t)\n")
        for x, val in zip(xs, vals):
            f.write(f"{x:.8g} {val:.16e}\n")
    return xs, vals


def plot_linecut(xs: np.ndarray, vals: np.ndarray, y: int, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xs, vals, "o-", ms=3.5, lw=1.4, color="#2855a3")
    ax.set_xlabel("x")
    ax.set_ylabel(rf"$P(x,y={y},t)$")
    ax.set_title("Dipanjan exact theory line cut")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_physical_heatmap(P: np.ndarray, out_path: Path, title: str) -> None:
    """
    Plot only the physical triangular-lattice sites.

    In Dipanjan's Cartesian embedding, the walker started at (0, 0) only lives
    on sites with x + y even. The off-parity sites are not physical sites for
    this initial condition.
    """
    L = P.shape[1]
    physical, _ = parity_masks(L)
    xs, ys = np.where(physical)
    vals = P.real[physical]

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    sc = ax.scatter(xs, ys, c=vals, s=23, cmap="inferno", marker="s")
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(r"$P(x,y,t)$")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def plot_full_grid(P: np.ndarray, out_path: Path, title: str) -> None:
    """Plot the full 2L x L Cartesian grid, including off-parity zero sites."""
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    im = ax.imshow(P.real.T, origin="lower", aspect="equal", cmap="inferno")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$P(x,y,t)$")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma", type=float, default=0.01)
    parser.add_argument("--epsilon", type=float, default=0.15)
    parser.add_argument("--L", type=int, default=30)
    parser.add_argument("--time", type=float, default=10.0)
    parser.add_argument("--y", type=int, default=4)
    parser.add_argument("--outdir", type=Path, default=Path(__file__).parent / "outputs" / "dipanjan_reproduction")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    print("Reproducing Dipanjan's exact Fourier calculation")
    print(f"gamma={args.gamma}, epsilon={args.epsilon}, L={args.L}, time={args.time}, y={args.y}")

    p_k = evolve_fourier_grid(
        gamma=args.gamma,
        epsilon=args.epsilon,
        time=args.time,
        L=args.L,
    )
    P = inverse_fourier_cartesian(p_k)

    max_imag = np.max(np.abs(P.imag))
    physical, off_lattice = parity_masks(args.L)
    full_sum = P.real.sum()
    physical_sum = P.real[physical].sum()
    off_sum = P.real[off_lattice].sum()
    dipanjan_sampled_sum = P.real[0 : 2 * args.L : 2, :].sum()

    print("\nSanity checks")
    print(f"max |Im P|                         = {max_imag:.3e}")
    print(f"sum over full 2L x L grid           = {full_sum:.12f}")
    print(f"sum over physical x+y even sites    = {physical_sum:.12f}")
    print(f"sum over off-parity sites           = {off_sum:.12e}")
    print(f"sum over Dipanjan x=0,2,... grid    = {dipanjan_sampled_sum:.12f}")

    time_tag = f"t{args.time:g}".replace(".", "p")
    y_tag = f"y{args.y:g}".replace(".", "p")
    line_txt = args.outdir / f"dipanjan_linecut_{y_tag}_{time_tag}.txt"
    line_png = args.outdir / f"dipanjan_linecut_{y_tag}_{time_tag}.png"
    heat_png = args.outdir / f"dipanjan_physical_heatmap_{time_tag}.png"
    full_png = args.outdir / f"dipanjan_full_cartesian_grid_{time_tag}.png"

    xs, vals = write_linecut(P, args.y, line_txt, args.gamma, args.epsilon, args.time, args.L)
    plot_linecut(xs, vals, args.y, line_png)
    plot_physical_heatmap(P, heat_png, rf"Physical triangular sites, $t={args.time}$")
    plot_full_grid(P, full_png, rf"Full Cartesian embedding, $t={args.time}$")

    print("\nWrote")
    print(f"  {line_txt}")
    print(f"  {line_png}")
    print(f"  {heat_png}")
    print(f"  {full_png}")

    print("\nFirst five line-cut values")
    for x, val in zip(xs[:5], vals[:5]):
        print(f"x={x:2d}  P={val:.16f}")


if __name__ == "__main__":
    main()
