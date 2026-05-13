import numpy as np
import matplotlib.pyplot as plt
from itertools import permutations
from pathlib import Path

from triangular_jmvr_corrected import build_Mk_corrected, build_Mk_dipanjan


def build_Mk(gamma, epsilon, k1, k2, a=1.0, b=1.0):
    """Default JMVR triangular generator: corrected c3 sign."""
    return build_Mk_corrected(gamma, epsilon, k1, k2, a, b)


def match_eigenvalues(prev_vals, vals):
    """Order vals so they move continuously from prev_vals."""
    best_perm = None
    best_cost = np.inf
    for perm in permutations(range(len(vals))):
        ordered = vals[list(perm)]
        cost = np.sum(np.abs(ordered - prev_vals) ** 2)
        if cost < best_cost:
            best_cost = cost
            best_perm = perm
    return vals[list(best_perm)]


def plot_band_line(gamma=0.01, epsilon=0.15, npts=300, outfile="band_line.png", tracked=True):
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    k1s = np.linspace(0.0, 2.0 * np.pi, npts)
    eigvals = np.zeros((npts, 6), dtype=complex)

    for i, k1 in enumerate(k1s):
        M = build_Mk(gamma=gamma, epsilon=epsilon, k1=k1, k2=0.0)
        vals = np.linalg.eigvals(M)
        if i == 0:
            vals = vals[np.argsort(vals.real)]
        elif tracked:
            vals = match_eigenvalues(eigvals[i - 1, :], vals)
        else:
            vals = vals[np.argsort(vals.real)]
        eigvals[i, :] = vals

    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)

    for band in range(6):
        axes[0].plot(k1s, eigvals[:, band].real, lw=1)
        imag = eigvals[:, band].imag.copy()
        imag[np.abs(imag) < 1e-14] = 0.0
        axes[1].plot(k1s, imag, lw=1)

    axes[0].set_ylabel(r"Re $\lambda$")
    axes[1].set_ylabel(r"Im $\lambda$")
    axes[1].set_xlabel(r"$k_1$ at $k_2=0$")
    if np.max(np.abs(eigvals.imag)) < 1e-14:
        axes[1].set_ylim(-0.05, 0.05)

    fig.suptitle(rf"Triangular JMVR bands, $\gamma={gamma}$, $\epsilon={epsilon}$")
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)
    print(f"saved {outfile}")


def run_sanity_tests():
    gamma = 0.01
    epsilon = 0.15

    M = build_Mk(gamma=gamma, epsilon=epsilon, k1=0.0, k2=0.0)

    print("M(0,0) =")
    print(M)

    print("column sums =")
    print(M.sum(axis=0))

    print("eigenvalues =")
    print(np.linalg.eigvals(M))

    print("\nnonzero-k tests")

    M_eps0 = build_Mk(gamma=gamma, epsilon=0.0, k1=0.3, k2=0.4)
    print("max imaginary part at epsilon=0 =")
    print(np.max(np.abs(M_eps0.imag)))

    M_pos = build_Mk(gamma=gamma, epsilon=epsilon, k1=0.3, k2=0.4)
    M_neg = build_Mk(gamma=gamma, epsilon=epsilon, k1=-0.3, k2=-0.4)

    print("max |M(-k) - conj(M(k))| =")
    print(np.max(np.abs(M_neg - np.conj(M_pos))))

    print("column sums at nonzero k =")
    print(M_pos.sum(axis=0))

    print("\ncoefficient check at k=(0.3,0.4)")
    k1 = 0.3
    k2 = 0.4
    a = 1.0
    b = 1.0

    M = build_Mk(gamma, epsilon, k1, k2, a, b)

    B = (1.0 / 3.0) * (
        np.cos(2.0 * a * k1)
        + np.cos(a * k1 + b * k2)
        + np.cos(a * k1 - b * k2)
    )
    c1 = B - 1.0 - gamma + 2.0j * epsilon * np.sin(2.0 * a * k1)
    c2 = B - 1.0 - gamma + 2.0j * epsilon * np.sin(a * k1 + b * k2)
    c3 = B - 1.0 - gamma - 2.0j * epsilon * np.sin(a * k1 - b * k2)

    expected_diag = np.array([c1, c2, c3, np.conj(c1), np.conj(c2), np.conj(c3)])
    actual_diag = np.diag(M)

    print("max diagonal mismatch =")
    print(np.max(np.abs(actual_diag - expected_diag)))

    print("\nDipanjan comparison at k=(0.3,0.4)")
    M_old = build_Mk_dipanjan(gamma, epsilon, k1, k2, a, b)
    print("max |corrected - Dipanjan| =")
    print(np.max(np.abs(M - M_old)))


if __name__ == "__main__":
    run_sanity_tests()
    plot_band_line(gamma=0.01, epsilon=0.0, outfile="outputs/band_line_eps0_tracked.png")
    plot_band_line(gamma=0.01, epsilon=0.15, outfile="outputs/band_line_eps015_tracked.png")
