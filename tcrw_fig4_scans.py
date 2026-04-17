"""
Fig 4(d),(e),(i): spectrum parameter scans.

- Fig 4(d): Re(lambda) vs D_r at omega=1, OBC, L=10. Shows how eigenvalues
            collapse/coalesce as D_r increases; edge modes separate from bulk.
- Fig 4(e): Re(lambda) vs omega at fixed D_r, OBC, L=10. Shows gap closing
            at omega=0.5.
- Fig 4(i): PBC band structure plotted in the (cos kx, cos ky) plane —
            the characteristic "circle/square" shape.

Uses existing obc_spectrum() and pbc_full_bz() from tcrw_spectrum.py.

Author: Prashant Bisht, TIFR Hyderabad
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from tcrw_spectrum import obc_spectrum, pbc_full_bz, build_Pk


# ============================================================
# Fig 4(d): Re(lambda) vs D_r for fixed omega
# ============================================================

def fig4d_real_vs_Dr(omega=1.0, L=10, D_r_values=None,
                      outfile='tcrw_fig4d_Re_vs_Dr.png'):
    if D_r_values is None:
        D_r_values = np.linspace(0.005, 0.5, 40)

    print(f"  Fig 4(d): omega={omega}, L={L}, {len(D_r_values)} D_r points")
    all_re = []
    all_edge = []
    all_Dr = []

    t0 = time.time()
    for i, D_r in enumerate(D_r_values):
        evals, edge_w = obc_spectrum(omega, D_r, L)
        all_re.append(evals.real)
        all_edge.append(edge_w)
        all_Dr.append(np.full_like(evals.real, D_r))
        if (i + 1) % 10 == 0:
            print(f"    {i+1}/{len(D_r_values)}: D_r={D_r:.3f} "
                  f"({time.time()-t0:.1f}s)")

    x = np.concatenate(all_Dr)
    y = np.concatenate(all_re)
    c = np.concatenate(all_edge)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    sc = ax.scatter(x, y, c=c, cmap='viridis', s=4, alpha=0.8,
                     norm=Normalize(vmin=0, vmax=c.max()))
    ax.axhline(1.0, color='k', lw=0.5, ls=':', alpha=0.5)
    ax.axhline(0.0, color='k', lw=0.5, ls=':', alpha=0.5)
    ax.set_xlabel('$D_r$', fontsize=13)
    ax.set_ylabel('$\\mathrm{Re}(\\lambda)$', fontsize=13)
    ax.set_title(f'Fig 4(d): OBC spectrum vs $D_r$ at $\\omega={omega}$, $L={L}$',
                 fontsize=12)
    fig.colorbar(sc, ax=ax, label='edge weight', pad=0.02)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outfile}")


# ============================================================
# Fig 4(e): Re(lambda) vs omega for fixed D_r
# ============================================================

def fig4e_real_vs_omega(D_r=0.1, L=10, omega_values=None,
                         outfile='tcrw_fig4e_Re_vs_omega.png'):
    if omega_values is None:
        omega_values = np.linspace(0.0, 1.0, 41)

    print(f"  Fig 4(e): D_r={D_r}, L={L}, {len(omega_values)} omega points")
    all_re = []
    all_edge = []
    all_om = []

    t0 = time.time()
    for i, om in enumerate(omega_values):
        evals, edge_w = obc_spectrum(om, D_r, L)
        all_re.append(evals.real)
        all_edge.append(edge_w)
        all_om.append(np.full_like(evals.real, om))
        if (i + 1) % 10 == 0:
            print(f"    {i+1}/{len(omega_values)}: omega={om:.3f} "
                  f"({time.time()-t0:.1f}s)")

    x = np.concatenate(all_om)
    y = np.concatenate(all_re)
    c = np.concatenate(all_edge)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    sc = ax.scatter(x, y, c=c, cmap='viridis', s=4, alpha=0.8,
                     norm=Normalize(vmin=0, vmax=c.max()))
    ax.axhline(1.0, color='k', lw=0.5, ls=':', alpha=0.5)
    ax.axhline(0.0, color='k', lw=0.5, ls=':', alpha=0.5)
    ax.axvline(0.5, color='tab:red', lw=1.0, ls='--', alpha=0.6,
                label='$\\omega=1/2$ (gap closes)')
    ax.set_xlabel('$\\omega$', fontsize=13)
    ax.set_ylabel('$\\mathrm{Re}(\\lambda)$', fontsize=13)
    ax.set_title(f'Fig 4(e): OBC spectrum vs $\\omega$ at $D_r={D_r}$, $L={L}$',
                 fontsize=12)
    ax.legend(fontsize=10, loc='lower right')
    fig.colorbar(sc, ax=ax, label='edge weight', pad=0.02)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outfile}")


# ============================================================
# Fig 4(i): PBC bands on (cos kx, cos ky) plane
# ============================================================

def fig4i_cos_plane(omega=1.0, D_r=0.1, Nk=80,
                     outfile='tcrw_fig4i_cos_plane.png'):
    """
    Plot PBC eigenvalues as a function of (cos kx, cos ky) over the
    full BZ. The 4-band spectrum traces out a characteristic surface.
    We scatter all 4 eigenvalues at each (kx, ky), coloured by Re(lambda).
    """
    print(f"  Fig 4(i): omega={omega}, D_r={D_r}, BZ {Nk}x{Nk}")
    kxs = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    kys = np.linspace(-np.pi, np.pi, Nk, endpoint=False)

    cx_list, cy_list, re_list, im_list = [], [], [], []
    for kx in kxs:
        for ky in kys:
            P = build_Pk(omega, D_r, kx, ky)
            evals = np.linalg.eigvals(P)
            for ev in evals:
                cx_list.append(np.cos(kx))
                cy_list.append(np.cos(ky))
                re_list.append(ev.real)
                im_list.append(ev.imag)

    cx = np.array(cx_list); cy = np.array(cy_list)
    re = np.array(re_list); im = np.array(im_list)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # (left) color by Re(lambda)
    sc0 = axes[0].scatter(cx, cy, c=re, s=3, cmap='RdBu_r',
                            vmin=-1, vmax=1, alpha=0.5)
    axes[0].set_xlabel('$\\cos k_x$', fontsize=13)
    axes[0].set_ylabel('$\\cos k_y$', fontsize=13)
    axes[0].set_title(f'Fig 4(i) left: bands coloured by $\\mathrm{{Re}}(\\lambda)$',
                       fontsize=11)
    axes[0].set_aspect('equal')
    fig.colorbar(sc0, ax=axes[0], shrink=0.85, pad=0.02)

    # (right) color by |lambda|
    mag = np.sqrt(re**2 + im**2)
    sc1 = axes[1].scatter(cx, cy, c=mag, s=3, cmap='viridis',
                            vmin=0, vmax=1, alpha=0.5)
    axes[1].set_xlabel('$\\cos k_x$', fontsize=13)
    axes[1].set_ylabel('$\\cos k_y$', fontsize=13)
    axes[1].set_title(f'Fig 4(i) right: bands coloured by $|\\lambda|$',
                       fontsize=11)
    axes[1].set_aspect('equal')
    fig.colorbar(sc1, ax=axes[1], shrink=0.85, pad=0.02)

    plt.suptitle(f'PBC spectrum on $(\\cos k_x, \\cos k_y)$ plane, '
                  f'$\\omega={omega}$, $D_r={D_r}$', fontsize=12, y=1.00)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outfile}")


if __name__ == '__main__':
    fig4d_real_vs_Dr()
    fig4e_real_vs_omega()
    fig4i_cos_plane()
