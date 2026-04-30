"""
Cross-check the new Fortran Fig 2 occupancy/currents against the
exact-Python steady state solved from the transition matrix.

[Updated 2026-04-30] The Fortran fig2_clean.f90 / fig2_defects.f90
have been fixed to drop the outer wall-ring (mask = .true. now), so
the walker plays on the FULL L × L grid.  Convention bookkeeping:

    Fortran  L_F = 10    →  walker on 10×10 grid, sites (0..9)
    Authors' TRW         →  L_A = L_F - 1 = 9  (grid 0..L_A, L_A+1 = 10 sites)

These describe the same physical lattice.  No slicing, no coordinate
shift — the Fortran (x, y) and authors' (i, j) labels are identical.
"""

import os, sys, numpy as np
import matplotlib.pyplot as plt

for c in (
    "/sessions/elegant-wizardly-einstein/mnt/TIFR DP1 to DP2 Research Assistant",
    "/Users/prashantbisht/Documents/Claude/Projects/TIFR DP1 to DP2 Research Assistant"):
    if os.path.isdir(c):
        ROOT = c
        break

NEW = os.path.join(ROOT, "new_fortran_reproduction_and_python")
sys.path.insert(0, ROOT)
sys.path.insert(0, NEW)

# Load authors' module
import importlib.util
spec = importlib.util.spec_from_file_location(
    "TRW_authors",
    os.path.join(NEW, "TRW._original_code_by_paperauthors.py"))
TRW = importlib.util.module_from_spec(spec)
spec.loader.exec_module(TRW)


def load_fortran_occ(path, L_F):
    """Read Fortran occupancy file: columns 'x y P'."""
    arr = np.loadtxt(path)
    P = np.zeros((L_F, L_F))
    for x, y, p in arr:
        P[int(x), int(y)] = p
    return P


def load_fortran_current(path):
    """Read Fortran current file: columns 'x y Jx Jy |J|'."""
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        arr = arr[None, :]
    return arr  # raw; we'll index by (x,y) directly


def compare_fig2(omega, D_r=1e-3, L_F=10):
    """
    Compare Fortran Fig 2 P(X,Y) and J_Dr/J_omega at given omega against
    the exact-Python steady state on the same L_F × L_F grid.
    """
    print(f"\n=== Fig 2, omega={omega}, D_r={D_r}, L_F={L_F} ===")

    # --- Fortran side: full L_F × L_F grid (no wall ring after 2026-04-30 fix) ---
    occ_path = os.path.join(NEW, f"tcrw_fig2_occ_w{omega:.1f}.txt")
    P_F = load_fortran_occ(occ_path, L_F)
    P_F = P_F / P_F.sum()   # renormalize (already 1 by construction)

    # --- Exact-Python: authors' L_A = L_F - 1 (10 sites needs L_A=9) ---
    L_A = L_F - 1
    W = TRW.build_sparse_transition_matrix(L_A, omega, D_r)
    pi, lam = TRW.solve_steady_state_sparse(W, L_A)
    P_E = np.zeros((L_A + 1, L_A + 1))
    for i in range(L_A + 1):
        for j in range(L_A + 1):
            for dn in TRW.directions:
                P_E[i, j] += pi[TRW.index(i, j, dn, L_A)]
    P_E = P_E / P_E.sum()

    # --- Compare ---
    assert P_F.shape == P_E.shape, f"shape mismatch {P_F.shape} vs {P_E.shape}"
    diff = np.abs(P_F - P_E)
    print(f"  playground shape: {P_F.shape}")
    print(f"  max|P_F - P_E| = {diff.max():.3e}   (rel {diff.max()/P_E.max():.2%})")
    print(f"  RMS (P_F - P_E) = {np.sqrt(np.mean((P_F-P_E)**2)):.3e}")
    print(f"  sum P_F = {P_F.sum():.6f}   sum P_E = {P_E.sum():.6f}")
    print(f"  max P_F = {P_F.max():.3e}   max P_E = {P_E.max():.3e}")

    # --- Currents (if files exist) ---
    for kind_F, label in [("JDr", "J_Dr"), ("Jomega", "J_omega")]:
        cpath = os.path.join(NEW, f"tcrw_fig2_{kind_F}_w{omega:.1f}.txt")
        if not os.path.exists(cpath):
            continue
        arr = load_fortran_current(cpath)
        # build (L_F-2, L_F-2) Jx, Jy from columns
        Jx_F = np.zeros((L_F, L_F))
        Jy_F = np.zeros((L_F, L_F))
        for row in arr:
            x, y, Jx, Jy, Jmag = row
            ix = int(x); iy = int(y)
            if 0 <= ix < L_F and 0 <= iy < L_F:
                Jx_F[ix, iy] = Jx
                Jy_F[ix, iy] = Jy

    # exact currents
    J1_E, J2_E = TRW.calculate_J1_J2_with_boundaries(pi, W, L_A, D_r, omega)
    # J1_E = J_Dr, J2_E = J_omega

    # Re-load each current file to match
    for kind_F, J_exact in [("JDr", J1_E), ("Jomega", J2_E)]:
        cpath = os.path.join(NEW, f"tcrw_fig2_{kind_F}_w{omega:.1f}.txt")
        if not os.path.exists(cpath):
            print(f"  [skip] {cpath} missing")
            continue
        arr = load_fortran_current(cpath)
        Jx_F = np.zeros((L_F, L_F))
        Jy_F = np.zeros((L_F, L_F))
        for row in arr:
            x, y, Jx, Jy, _ = row
            ix = int(x); iy = int(y)
            if 0 <= ix < L_F and 0 <= iy < L_F:
                Jx_F[ix, iy] = Jx
                Jy_F[ix, iy] = Jy
        # Fortran accumulates 1 count per successful translation.  Exact
        # calculation produces current in the same units (per step).
        dJx = Jx_F - J_exact[..., 0]
        dJy = Jy_F - J_exact[..., 1]
        mag = np.hypot(Jx_F, Jy_F)
        ref = max(np.abs(J_exact).max(), 1e-30)
        print(f"  {kind_F:6s}: max|Fortran - exact| = "
              f"(Jx) {np.abs(dJx).max():.3e}  (Jy) {np.abs(dJy).max():.3e}"
              f"   (ref scale {ref:.2e})")

    # ---------- Save a side-by-side comparison figure ----------
    # Re-collect Fortran currents on the playground for plotting.
    J_F = {}
    for kind_F in ("JDr", "Jomega"):
        cpath = os.path.join(NEW, f"tcrw_fig2_{kind_F}_w{omega:.1f}.txt")
        if not os.path.exists(cpath):
            J_F[kind_F] = None
            continue
        arr = load_fortran_current(cpath)
        Jx_F = np.zeros((L_F, L_F))
        Jy_F = np.zeros((L_F, L_F))
        for row in arr:
            x, y, Jx, Jy, _ = row
            ix = int(x); iy = int(y)
            if 0 <= ix < L_F and 0 <= iy < L_F:
                Jx_F[ix, iy] = Jx
                Jy_F[ix, iy] = Jy
        J_F[kind_F] = (Jx_F, Jy_F)

    fig = plt.figure(figsize=(12, 8))
    fig.suptitle(
        rf"Fig 2 cross-check: Fortran MC vs exact   $\omega={omega}$, "
        rf"$D_r={D_r}$, playground ${P_F.shape[0]}\times{P_F.shape[1]}$",
        fontsize=12)

    pmax = max(P_F.max(), P_E.max())
    extent = [0, P_F.shape[0], 0, P_F.shape[1]]
    common = dict(origin="lower", extent=extent, cmap="viridis",
                  vmin=0, vmax=pmax, interpolation="nearest")

    ax1 = plt.subplot(2, 3, 1)
    im1 = ax1.imshow(P_F.T, **common)
    ax1.set_title(r"$P_{\mathrm{Fortran}}(x,y)$")
    plt.colorbar(im1, ax=ax1, fraction=0.046)

    ax2 = plt.subplot(2, 3, 2)
    im2 = ax2.imshow(P_E.T, **common)
    ax2.set_title(r"$P_{\mathrm{exact}}(x,y)$")
    plt.colorbar(im2, ax=ax2, fraction=0.046)

    ax3 = plt.subplot(2, 3, 3)
    im3 = ax3.imshow(np.abs(P_F - P_E).T, origin="lower", extent=extent,
                     cmap="magma", interpolation="nearest")
    ax3.set_title(r"$|P_F - P_E|$")
    plt.colorbar(im3, ax=ax3, fraction=0.046)

    # Currents: quiver overlays (Fortran red, exact black) — paper-style
    # fixed scale + magnitude threshold so bulk-noise arrows are hidden.
    xs = np.arange(P_F.shape[0]) + 0.5
    ys = np.arange(P_F.shape[1]) + 0.5
    Xg, Yg = np.meshgrid(xs, ys, indexing="ij")

    # Global max across both J types for shared scale + threshold
    g_max = max(np.sqrt(J1_E[..., 0]**2 + J1_E[..., 1]**2).max(),
                np.sqrt(J2_E[..., 0]**2 + J2_E[..., 1]**2).max())
    q_scale = g_max / 0.8 if g_max > 0 else 1.0
    q_thresh = 0.01 * g_max

    def _filter(Jx, Jy):
        m = np.sqrt(Jx**2 + Jy**2)
        keep = m > q_thresh
        return np.where(keep, Jx, 0.0), np.where(keep, Jy, 0.0)

    for k, (kind_F, J_exact, label) in enumerate(
        [("JDr", J1_E, r"$J_{D_r}$"), ("Jomega", J2_E, r"$J_\omega$")]):
        ax = plt.subplot(2, 3, 4 + k)
        if J_F.get(kind_F) is not None:
            Jx_F, Jy_F = J_F[kind_F]
            Jx_F2, Jy_F2 = _filter(Jx_F, Jy_F)
            ax.quiver(Xg, Yg, Jx_F2, Jy_F2, color="C3",
                      scale=q_scale, scale_units="xy",
                      angles="xy", pivot="mid", width=0.005,
                      label="Fortran")
        Jx_E2, Jy_E2 = _filter(J_exact[..., 0], J_exact[..., 1])
        ax.quiver(Xg, Yg, Jx_E2, Jy_E2, color="k",
                  scale=q_scale, scale_units="xy",
                  angles="xy", alpha=0.55, pivot="mid", width=0.005,
                  label="exact")
        ax.set_xlim(0, P_F.shape[0])
        ax.set_ylim(0, P_F.shape[1])
        ax.set_aspect("equal")
        ax.set_title(f"{label}: Fortran (red) vs exact (black)")
        ax.legend(loc="upper right", fontsize=8)

    # Stats panel (bottom-right)
    ax_stat = plt.subplot(2, 3, 6)
    ax_stat.axis("off")
    txt = (f"max|P_F - P_E| = {np.abs(P_F-P_E).max():.3e}\n"
           f"  rel to max P_E = {np.abs(P_F-P_E).max()/P_E.max():.2%}\n"
           f"RMS (P_F - P_E) = "
           f"{np.sqrt(np.mean((P_F-P_E)**2)):.3e}\n"
           f"sum P_F = {P_F.sum():.6f}\n"
           f"sum P_E = {P_E.sum():.6f}\n"
           f"max P_F = {P_F.max():.3e}\n"
           f"max P_E = {P_E.max():.3e}")
    ax_stat.text(0.02, 0.98, txt, va="top", ha="left",
                 family="monospace", fontsize=10)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(NEW, f"tcrw_fig2_fortran_vs_exact_w{omega:.1f}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


if __name__ == "__main__":
    for om in (0.0, 0.5, 1.0):
        try:
            compare_fig2(omega=om)
        except FileNotFoundError as e:
            print(f"  [skip omega={om}] {e}")
