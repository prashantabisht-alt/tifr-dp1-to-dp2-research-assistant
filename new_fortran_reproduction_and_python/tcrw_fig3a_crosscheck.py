"""
Fig 3(a) Fortran-vs-Python exact cross-check
=============================================

The Fortran driver `tcrw_fig3a.f90` runs MC at ω = 1, OBC, varying L and
D_r, and writes per-site averages
    P_edge_norm = (visits to any edge site)  / (n_edge × T_meas)
    P_bulk_norm = (visits to any bulk site)  / (n_bulk × T_meas)
    ratio       = P_edge_norm / P_bulk_norm
into `tcrw_fig3a_summary.txt`.

For each row of that file we:
  1. Build the OBC transition matrix at the same (ω = 1, D_r, L)
     using the same authors'-convention indexing as fig4c.
  2. Solve for the steady-state right eigenvector π (sparse Perron).
  3. Sum π over edge / bulk SITES (any director).
  4. Divide by n_edge / n_bulk to get the per-site averages.
  5. Overlay the exact line on the Fortran scatter.

Expectation: at MC T = 10^8–10^9 (which is what the Fortran header
implies), per-site P should agree with the exact value to relative
1/√(T·P) ≈ 10⁻³–10⁻⁴ on the bulk sites; the EDGE side is more
trapped (P_edge ~ 10⁻¹), so MC noise there is well below 10⁻⁴.

Author: Prashant Bisht, TIFR Hyderabad
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Make sure tcrw_fig4c is importable (gives us build_obc_matrix +
# state_index, all in authors' convention).
HERE = os.path.dirname(os.path.abspath(__file__))
for cand in (HERE,
             "/sessions/elegant-wizardly-einstein/mnt/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python",
             "/Users/prashantbisht/Documents/Claude/Projects/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python"):
    if os.path.isdir(cand):
        sys.path.insert(0, cand)
        ROOT = cand
        break

import tcrw_fig4c as fig4c  # build_obc_matrix, state_index


# ---------------------------------------------------------------------------
# 1. Load Fortran summary
# ---------------------------------------------------------------------------
def load_fortran_summary(path: str):
    """Returns dict: L -> (D_r_array, ratio_F, P_edge_F, P_bulk_F, n_edge, n_bulk).

    The Fortran writer emits "1.79769+308" (huge_number with no 'E' before
    the exponent sign) for divisions by zero.  np.loadtxt can't parse
    that, so we read line-by-line and patch the format.
    """
    rows = []
    with open(path) as f:
        for line in f:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            parts = line.split()
            fixed = []
            for p in parts:
                # Repair "1.79769+308" -> "1.79769E+308" (and similarly for -)
                if "E" not in p and "e" not in p:
                    # strip leading sign, then look for + or - in the middle
                    body = p.lstrip("+-")
                    sign = p[: len(p) - len(body)]
                    ipos = max(body.rfind("+"), body.rfind("-"))
                    if ipos > 0 and body[ipos - 1].isdigit():
                        body = body[:ipos] + "E" + body[ipos:]
                    p = sign + body
                fixed.append(p)
            rows.append([float(x) for x in fixed])
    arr = np.array(rows, dtype=float)
    cols = ("L", "D_r", "ratio", "P_edge", "P_bulk", "n_edge", "n_bulk")
    out = {}
    for L_val in np.unique(arr[:, 0].astype(int)):
        m = arr[:, 0] == L_val
        sub = arr[m]
        out[int(L_val)] = dict(
            D_r       = sub[:, 1],
            ratio     = sub[:, 2],
            P_edge    = sub[:, 3],
            P_bulk    = sub[:, 4],
            n_edge    = int(sub[0, 5]),
            n_bulk    = int(sub[0, 6]),
        )
    return out


# ---------------------------------------------------------------------------
# 2. Per-site edge/bulk masks at the spatial-only level
# ---------------------------------------------------------------------------
def site_masks(L: int):
    """
    Return two boolean masks of length 4(L+1)^2:
      - edge_states: True where (i, j) is on the OBC boundary (any director)
      - bulk_states: True where (i, j) is strictly interior
    These match the Fortran fig3a partition (which is per-site, any-director,
    NOT the Fig 4(a) directional partition used in fig4c).
    """
    n = 4 * (L + 1) ** 2
    edge_states = np.zeros(n, dtype=bool)
    bulk_states = np.zeros(n, dtype=bool)
    n_edge_sites = 0
    n_bulk_sites = 0
    for i in range(L + 1):
        for j in range(L + 1):
            is_boundary = (i == 0 or i == L or j == 0 or j == L)
            if is_boundary:
                n_edge_sites += 1
                for d in range(4):
                    edge_states[fig4c.state_index(i, j, d, L)] = True
            else:
                n_bulk_sites += 1
                for d in range(4):
                    bulk_states[fig4c.state_index(i, j, d, L)] = True
    return edge_states, bulk_states, n_edge_sites, n_bulk_sites


# ---------------------------------------------------------------------------
# 3. Exact steady state at (ω, D_r, L) — sparse Perron
# ---------------------------------------------------------------------------
def steady_state(omega: float, D_r: float, L: int):
    P = fig4c.build_obc_matrix(omega, D_r, L)            # CSC, column-stochastic
    n = P.shape[0]
    if n <= 600:
        # tiny: dense eig
        Pd = P.toarray()
        evals, evecs = np.linalg.eig(Pd)
        idx = np.argmin(np.abs(evals - 1.0))
        pi = evecs[:, idx].real
    else:
        # large: shift-invert sparse eigs around λ = 1
        try:
            evals, evecs = spla.eigs(P, k=1, sigma=1.0 + 1e-8, which="LM",
                                     maxiter=500, tol=1e-12)
            pi = evecs[:, 0].real
        except Exception:
            # fall back to dense if sparse fails
            Pd = P.toarray()
            evals, evecs = np.linalg.eig(Pd)
            idx = np.argmin(np.abs(evals - 1.0))
            pi = evecs[:, idx].real

    if pi.sum() < 0:
        pi = -pi
    pi[pi < 0] = 0.0
    pi /= pi.sum()
    return pi


def edge_bulk_per_site(omega: float, D_r: float, L: int):
    pi = steady_state(omega, D_r, L)
    edge, bulk, n_edge_sites, n_bulk_sites = site_masks(L)
    P_edge_norm = float(pi[edge].sum() / n_edge_sites)
    P_bulk_norm = float(pi[bulk].sum() / n_bulk_sites)
    return P_edge_norm, P_bulk_norm


# ---------------------------------------------------------------------------
# 4. Cross-check
# ---------------------------------------------------------------------------
def crosscheck(summary_path: str, plot_path: str = "tcrw_fig3a_crosscheck.png"):
    print(f"loading {summary_path}")
    data = load_fortran_summary(summary_path)
    print(f"  L values: {sorted(data.keys())}")

    fig, (ax_edge, ax_bulk, ax_ratio) = plt.subplots(1, 3, figsize=(15, 4.5))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(data)))

    summary = []
    for color, L in zip(colors, sorted(data.keys())):
        d = data[L]
        Drs = d["D_r"]
        # check site counts match expectation
        n_edge_exp = 4 * L
        n_bulk_exp = (L - 1) ** 2
        assert d["n_edge"] == n_edge_exp, f"L={L}: n_edge {d['n_edge']} != {n_edge_exp}"
        assert d["n_bulk"] == n_bulk_exp, f"L={L}: n_bulk {d['n_bulk']} != {n_bulk_exp}"

        Pe_x = np.zeros_like(Drs)
        Pb_x = np.zeros_like(Drs)
        print(f"  L = {L}: scanning {len(Drs)} D_r points...", end="", flush=True)
        t0 = time.time()
        for i, D in enumerate(Drs):
            Pe_x[i], Pb_x[i] = edge_bulk_per_site(1.0, float(D), L)
        print(f" cpu = {time.time()-t0:.1f} s")

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_x = Pe_x / Pb_x

        # Skip the D_r = 1 row — pure-noise limit has a degenerate
        # steady state (block-diagonal per site, λ = 1 is highly
        # degenerate).  Both Fortran and Python pick arbitrary
        # representatives; comparison is meaningless there.
        keep = Drs < 0.95
        rel_edge  = float(np.max(np.abs(d["P_edge"][keep] - Pe_x[keep]) / Pe_x[keep]))
        rel_bulk  = float(np.max(np.abs(d["P_bulk"][keep] - Pb_x[keep]) / Pb_x[keep]))
        rel_ratio = float(np.max(np.abs(d["ratio"][keep]  - ratio_x[keep]) / ratio_x[keep]))
        summary.append((L, rel_edge, rel_bulk, rel_ratio))

        ax_edge.scatter(Drs, d["P_edge"], color=color, s=14, alpha=0.7,
                        label=f"L={L} (MC)")
        ax_edge.plot   (Drs, Pe_x,        color=color, lw=1.0, alpha=0.9)
        ax_bulk.scatter(Drs, d["P_bulk"], color=color, s=14, alpha=0.7)
        ax_bulk.plot   (Drs, Pb_x,        color=color, lw=1.0, alpha=0.9)
        ax_ratio.scatter(Drs, d["ratio"],  color=color, s=14, alpha=0.7)
        ax_ratio.plot   (Drs, ratio_x,    color=color, lw=1.0, alpha=0.9)

    for ax, ylab in [(ax_edge, r"$\langle P\rangle_{\rm edge}$"),
                     (ax_bulk, r"$\langle P\rangle_{\rm bulk}$"),
                     (ax_ratio, r"ratio $= \langle P\rangle_{\rm edge}/\langle P\rangle_{\rm bulk}$")]:
        ax.set_xscale("log");  ax.set_yscale("log")
        ax.set_xlabel(r"$D_r$"); ax.set_ylabel(ylab)
        ax.grid(True, which="both", alpha=0.3)
    ax_edge.legend(fontsize=9, loc="lower left")
    fig.suptitle("Fig 3(a) cross-check — Fortran MC (markers) vs exact Python (lines)")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    print(f"  saved {plot_path}")

    print("\nrelative-error summary (max over D_r per L):")
    print(f"  {'L':>3} {'edge':>10} {'bulk':>10} {'ratio':>10}")
    for L, e, b, r in summary:
        print(f"  {L:>3} {e:>10.2e} {b:>10.2e} {r:>10.2e}")


if __name__ == "__main__":
    summary_file = os.path.join(ROOT, "tcrw_fig3a_summary.txt")
    out_file = os.path.join(ROOT, "tcrw_fig3a_crosscheck.png")
    crosscheck(summary_file, out_file)
