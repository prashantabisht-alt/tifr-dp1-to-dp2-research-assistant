"""
Cross-check Fortran fig3b MC vs exact Python (|J_Dr|/|J_omega| at left wall).

Loads tcrw_fig3b_summary.txt:
    columns: L  D_r  ratio  |J_Dr|_wall  |J_omega|_wall  Jx_w_Dr Jy_w_Dr Jx_w_om Jy_w_om

The Fortran J columns are integer COUNTS (sums over MC steps).  To get
current per step they are divided by T_meas (roughly N_traj × T per row).
For the RATIO |J_Dr|/|J_omega|, the T factor cancels — so we can compare
ratios directly without normalisation.

Compute Python exact ratio at the same (L, D_r), overlay.
"""
import os, sys, importlib.util
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
for cand in (HERE,
             "/sessions/elegant-wizardly-einstein/mnt/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python",
             "/Users/prashantbisht/Documents/Claude/Projects/TIFR DP1 to DP2 Research Assistant/new_fortran_reproduction_and_python"):
    p = os.path.join(cand, "TRW._original_code_by_paperauthors.py")
    if os.path.isfile(p):
        ROOT = cand; TRW_PATH = p
        break

spec = importlib.util.spec_from_file_location("TRW", TRW_PATH)
TRW = importlib.util.module_from_spec(spec); spec.loader.exec_module(TRW)


def load_fig3b(path):
    """Robust loader handling Fortran's '1.79+308' format."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            parts = line.split()
            fixed = []
            for p in parts:
                if "E" not in p and "e" not in p:
                    body = p.lstrip("+-")
                    sign = p[: len(p) - len(body)]
                    ipos = max(body.rfind("+"), body.rfind("-"))
                    if ipos > 0 and body[ipos - 1].isdigit():
                        body = body[:ipos] + "E" + body[ipos:]
                    p = sign + body
                fixed.append(p)
            rows.append([float(x) for x in fixed])
    arr = np.array(rows)
    out = {}
    for L in np.unique(arr[:, 0].astype(int)):
        m = arr[:, 0] == L
        out[int(L)] = dict(D_r=arr[m, 1], ratio=arr[m, 2])
    return out


def exact_J_ratio(omega, D_r, L):
    W = TRW.build_sparse_transition_matrix(L, omega, D_r)
    pi, _ = TRW.solve_steady_state_sparse(W, L)
    J_Dr, J_om = TRW.calculate_J1_J2_with_boundaries(pi, W, L, D_r, omega)
    Jdr_mag = np.sqrt(J_Dr[0, :, 0] ** 2 + J_Dr[0, :, 1] ** 2).mean()
    Jom_mag = np.sqrt(J_om[0, :, 0] ** 2 + J_om[0, :, 1] ** 2).mean()
    return Jdr_mag / max(Jom_mag, 1e-30)


def crosscheck():
    summary = os.path.join(ROOT, "tcrw_fig3b_summary.txt")
    data = load_fig3b(summary)
    print(f"loaded {summary}")
    print(f"  L values: {sorted(data.keys())}")

    fig, ax = plt.subplots(figsize=(7, 5))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(data)))

    summary_rows = []
    for color, L in zip(colors, sorted(data.keys())):
        d = data[L]
        Drs = d["D_r"]
        ratio_F = d["ratio"]
        ratio_E = np.array([exact_J_ratio(1.0, float(D), L) for D in Drs])

        keep = Drs < 0.95   # skip pure-noise degeneracy
        rel = np.abs(ratio_F[keep] - ratio_E[keep]) / np.maximum(ratio_E[keep], 1e-12)
        med = np.median(rel)
        mx = np.max(rel)
        summary_rows.append((L, med, mx))

        ax.scatter(Drs, ratio_F, color=color, s=14, alpha=0.7, label=f"L={L} MC")
        ax.plot(Drs, ratio_E, color=color, lw=1.0, alpha=0.9)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$D_r$"); ax.set_ylabel(r"$|J_{D_r}|/|J_\omega|$  (left wall)")
    ax.set_title("Fig 3(b) cross-check — MC dots vs exact lines")
    ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
    plt.tight_layout()
    out = os.path.join(ROOT, "tcrw_fig3b_crosscheck.png")
    plt.savefig(out, dpi=180, bbox_inches="tight")
    print(f"  saved {out}")

    print(f"\n  L  median rel  max rel")
    for L, med, mx in summary_rows:
        print(f"  {L:>2}  {med:>10.2e}  {mx:>10.2e}")


if __name__ == "__main__":
    crosscheck()
