"""Compare KMC MSD data with the exact Bloch-matrix theory.

Reads kmc_chiral_msd.dat (output of kmc_chiral_rtw.f90) and overlays:
  1. The full time-dependent MSD from the matrix exponential
  2. The asymptotic diffusive line 4 D_even t
  3. The short-time Poisson-jump line v t

Also prints D_even extracted from the KMC late-time slope.

Usage:
    python plot_kmc_vs_theory.py [datafile]
    python plot_kmc_vs_theory.py run

If datafile is omitted, or if the argument is "run", defaults to
kmc_chiral_msd.dat.
"""

from __future__ import annotations
import re
import sys
import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from triangular_chiral_rtw import (
    build_Mk_chiral_rtw, rates_from_gamma_bias, triangular_deltas,
    diffusion_tensor_exact, N_DIR,
)


def msd_theory_full(t_arr, *, v, gamma, b, gamma_r):
    """Compute the exact MSD(t) from the matrix exponential of R.

    MSD(t) = 2 * integral_0^t (t - s) Tr[C^v_cart(s)] ds

    where the velocity autocorrelation trace in Cartesian is:
        Tr[C^v_cart(s)] = (v^2 / 6) sum_m sum_{m'} delta_m^cart . delta_{m'}^cart [e^{Rs}]_{m',m}

    We evaluate this numerically via the matrix exponential.
    The jump-noise contribution adds v * t to the MSD.
    """
    gp, gm = rates_from_gamma_bias(gamma, b)
    R = build_Mk_chiral_rtw(0, 0, v=v, gamma_plus=gp, gamma_minus=gm, gamma_r=gamma_r)
    R = np.real(R)

    # Cartesian displacements
    deltas_ax = triangular_deltas()
    sqrt3 = np.sqrt(3.0)
    J = np.array([[1.0, 0.5], [0.0, 0.5 * sqrt3]])
    deltas_cart = (J @ deltas_ax.T).T  # shape (6, 2)

    # Precompute delta_m . delta_{m'} (trace contribution)
    # Tr[C^v_cart(s)] = (v^2/6) sum_{m,m'} (delta_m . delta_{m'}) [e^{Rs}]_{m',m}
    #                 = (v^2/6) Tr( D^T e^{Rs} )
    # where D_{m,m'} = delta_m . delta_{m'}
    dot_matrix = deltas_cart @ deltas_cart.T  # shape (6,6), D_{m,m'} = delta_m . delta_{m'}

    # Numerical integration via trapezoidal rule
    msd = np.zeros_like(t_arr)
    for i, t in enumerate(t_arr):
        if t <= 0:
            msd[i] = 0.0
            continue
        # Adaptive number of integration steps
        n_steps = max(200, int(50 * t))
        n_steps = min(n_steps, 5000)
        s_arr = np.linspace(0, t, n_steps + 1)
        ds = s_arr[1] - s_arr[0]

        integrand = np.zeros(n_steps + 1)
        for j, s in enumerate(s_arr):
            eRs = expm(R * s)
            # Tr(D^T e^{Rs}) = sum_{m,m'} D_{m',m} [e^{Rs}]_{m',m}
            # = sum of element-wise product
            trace_Cv = (v**2 / N_DIR) * np.sum(dot_matrix * eRs)
            integrand[j] = (t - s) * trace_Cv

        # Persistent MSD from velocity correlation
        msd_pers = 2.0 * np.trapz(integrand, s_arr)

        # Jump-noise MSD: each Poisson jump contributes |delta|^2 = 1
        # Rate of jumps = v, so MSD_jump = v * t * <|delta|^2> = v * t
        msd_jump = v * t

        msd[i] = msd_pers + msd_jump

    return msd


def main():
    datafile = sys.argv[1] if len(sys.argv) > 1 else 'kmc_chiral_msd.dat'
    if datafile == 'run':
        datafile = 'kmc_chiral_msd.dat'

    # Parse header for parameters
    v, gamma, b, gamma_r = 1.0, 1.0, 0.5, 0.0
    num = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
    with open(datafile) as f:
        for line in f:
            if not line.startswith('#'):
                break
            for name in ('v', 'gamma', 'bias', 'gamma_r'):
                match = re.search(rf'\b{name}\s*=\s*({num})', line)
                if match is None:
                    continue
                value = float(match.group(1))
                if name == 'v':
                    v = value
                elif name == 'gamma':
                    gamma = value
                elif name == 'bias':
                    b = value
                elif name == 'gamma_r':
                    gamma_r = value

    print(f"Parameters: v={v}, gamma={gamma}, b={b}, gamma_r={gamma_r}")

    # Load KMC data
    data = np.loadtxt(datafile)
    t_kmc = data[:, 0]
    mean_x = data[:, 1]
    mean_y = data[:, 2]
    mean_x2 = data[:, 3]
    mean_y2 = data[:, 4]
    mean_xy = data[:, 5]
    msd_kmc = data[:, 6]

    # Exact diffusion tensor
    D_ex = diffusion_tensor_exact(v=v, gamma=gamma, b=b, gamma_r=gamma_r)
    D_even = D_ex['D_even']
    D_odd = D_ex['D_odd']
    print(f"D_even = {D_even:.8f},  D_odd = {D_odd:.8f}")

    # Full time-dependent MSD from theory
    print("Computing full MSD(t) from matrix exponential...")
    t_theory = np.logspace(np.log10(t_kmc[0]), np.log10(t_kmc[-1]), 200)
    msd_theory = msd_theory_full(t_theory, v=v, gamma=gamma, b=b, gamma_r=gamma_r)

    # Extract D_even from late-time KMC
    late = t_kmc > 0.5 * t_kmc[-1]
    if np.sum(late) >= 2:
        # Linear fit to MSD = 4 D t + const
        coeffs = np.polyfit(t_kmc[late], msd_kmc[late], 1)
        D_even_kmc = coeffs[0] / 4.0
        print(f"D_even (KMC fit, late times) = {D_even_kmc:.8f}")
        print(f"D_even (theory)              = {D_even:.8f}")
        print(f"relative error               = {abs(D_even_kmc - D_even)/D_even:.6f}")

    # ---- Plot 1: MSD(t) ----
    fig, axes = plt.subplots(2, 1, figsize=(8, 10))

    ax = axes[0]
    ax.loglog(t_kmc, msd_kmc, 'o', ms=4, color='C0', label='KMC', zorder=3)
    ax.loglog(t_theory, msd_theory, '-', color='C1', lw=2,
              label='Theory (matrix exp)')
    ax.loglog(t_theory, 4 * D_even * t_theory, '--', color='C2', lw=1.5,
              label=f'$4D_{{\\mathrm{{even}}}} t$ = {4*D_even:.4f} t')
    ax.loglog(t_theory, v * t_theory, ':', color='C3', lw=1.5,
              label=f'$v t$ (short-time jumps)')
    ax.set_xlabel('$t$')
    ax.set_ylabel('MSD  $\\langle |\\mathbf{r}(t)|^2 \\rangle$')
    ax.set_title(f'Chiral RTW on triangular lattice\n'
                 f'$v={v}$, $\\gamma={gamma}$, $b={b}$, $\\gamma_r={gamma_r}$')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # ---- Plot 2: D_even(t) = MSD / 4t ----
    ax = axes[1]
    D_kmc_t = msd_kmc / (4 * t_kmc)
    D_theory_t = msd_theory / (4 * t_theory)
    ax.semilogx(t_kmc, D_kmc_t, 'o', ms=4, color='C0', label='KMC')
    ax.semilogx(t_theory, D_theory_t, '-', color='C1', lw=2, label='Theory')
    ax.axhline(D_even, color='C2', ls='--', lw=1.5,
               label=f'$D_{{\\mathrm{{even}}}}$ = {D_even:.6f}')
    ax.axhline(v / 4, color='C4', ls=':', lw=1,
               label=f'$v/4$ = {v/4:.4f} (jump noise)')
    ax.set_xlabel('$t$')
    ax.set_ylabel('$D_{\\mathrm{eff}}(t) = \\mathrm{MSD}/(4t)$')
    ax.set_title('Effective diffusion coefficient vs time')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, max(D_kmc_t) * 1.3)

    plt.tight_layout()
    plt.savefig('kmc_vs_theory.png', dpi=150)
    print("Saved kmc_vs_theory.png")

    # ---- Plot 3: Isotropy check ----
    fig2, ax2 = plt.subplots(1, 1, figsize=(7, 5))
    ax2.semilogx(t_kmc, mean_x2, 's', ms=3, color='C0', label='$\\langle x^2 \\rangle$')
    ax2.semilogx(t_kmc, mean_y2, 'D', ms=3, color='C1', label='$\\langle y^2 \\rangle$')
    ax2.semilogx(t_kmc, 2 * D_even * t_kmc, '--', color='C2', lw=1.5,
                 label=f'$2 D_{{\\mathrm{{even}}}} t$')
    ax2.set_xlabel('$t$')
    ax2.set_ylabel('second moment')
    ax2.set_title('Isotropy check: $\\langle x^2 \\rangle \\approx \\langle y^2 \\rangle$')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('kmc_isotropy_check.png', dpi=150)
    print("Saved kmc_isotropy_check.png")


if __name__ == '__main__':
    main()
