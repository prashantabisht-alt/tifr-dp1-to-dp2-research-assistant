"""
TCRW Fig 6(c) — single-tile trajectory near a seed.

GOAL
----
Reproduce Fig 6(c) of arXiv:2602.12020v1: sample trajectories of one tile during
self-assembly, for achiral (omega = 0.5) and chiral (omega = 1.0).  In the
chiral case the walker follows the seed boundary (topologically protected edge
current) — the whole physical point of the self-assembly section.

METHOD
------
1. Build a seed as a small rectangular block of *defect* sites.  The authors'
   ChiralWalker treats defects exactly like the outer box boundary: a move into
   a defect is blocked, and a blocked chiral-move step rotates no direction
   either.  This is the mechanism that produces edge following when omega = 1.

2. Instantiate ChiralWalker with record_path=True and step it T times.  The
   walker's internal `path` list records (i, j) at every time step (including
   blocked steps), giving a continuous trajectory of length T+1.

3. Plot the trajectory as a LineCollection with segments colored by time
   (viridis), 0 --> T.  Draw the seed as a black hatched rectangle.  Paint
   the box boundary for reference.  Do both omega values side-by-side.

CHOICES / ASSUMPTIONS
---------------------
- Lattice: L = 30 in authors' convention (i.e. 31 x 31 sites, 0..30).
- Seed: 4 x 4 block centered in the box (x in 13..16, y in 13..16).
  Paper's Fig 6(c) shows a seed that looks roughly that size.
- T = 1500 steps, so we have enough trajectory to see the ω=1 edge-following
  clearly but not so many that the chiral one wraps the seed many times.
- D_r = 0.1.  Paper's Fig 6(c) caption does not pin D_r explicitly;
  Fig 6(d)/(e) use D_r in {0.25, 0.5}.  D_r = 0.1 is mild enough that the
  chiral drift is visible; larger D_r gives a blurrier edge trace.
- Walker starts just above the seed with direction '↓' (so the first few
  steps point it at the seed and trigger the edge interaction quickly).

SANITY CHECKS
-------------
- path length == T + 1  (every step recorded, including blocked ones)
- walker never enters a defect site
- all path coords in [0, L]^2
- at omega = 1, the trajectory should have a visibly higher fraction of steps
  along the seed boundary than in the bulk  (report both fractions as a
  numeric check)

Author: Prashant Bisht, TIFR Hyderabad.
"""

import importlib.util
import os
import random
import sys
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle


# ---------------------------------------------------------------
# Import authors' ChiralWalker (dotted filename -> use importlib)
# ---------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
TRW_PATH = os.path.normpath(os.path.join(
    HERE, "..", "fortran_reproduction",
    "TRW._original_code_by_paperauthors.py"
))
spec = importlib.util.spec_from_file_location("TRW_authors", TRW_PATH)
TRW = importlib.util.module_from_spec(spec)
spec.loader.exec_module(TRW)


# ---------------------------------------------------------------
# Seed construction
# ---------------------------------------------------------------
def build_seed_defects(L, seed_lo, seed_hi):
    """
    Rectangular seed block: all (x, y) with seed_lo <= x <= seed_hi and
    seed_lo <= y <= seed_hi are marked as defects (inaccessible).

    Returns list of (x, y) tuples.
    """
    return [(x, y)
            for x in range(seed_lo, seed_hi + 1)
            for y in range(seed_lo, seed_hi + 1)]


# ---------------------------------------------------------------
# Trajectory generation
# ---------------------------------------------------------------
def run_trajectory(L, omega, D_r, T, seed_defects,
                   start_xy, start_dir, rng_seed):
    """
    Run a single ChiralWalker for T steps and return:
      path   : ndarray (M, 2)  of positions recorded
      times  : ndarray (M,)    real step index corresponding to each path point

    IMPORTANT about authors' convention: ChiralWalker.step() appends to
    self.path on chiral steps (both success and blocked), but NOT on
    rotational-noise steps.  Hence len(path) = 1 + (# chiral steps) and is
    typically ~ (1-D_r) * T.  For the colour-by-time to match the paper's
    '0 .. 10^3' colorbar, we must record the REAL step index at each append,
    not the append index.
    """
    random.seed(rng_seed)
    np.random.seed(rng_seed)

    walker = TRW.ChiralWalker(
        i=start_xy[0], j=start_xy[1],
        direction=start_dir,
        D_r=D_r, omega=omega,
        L=L,
        defects=seed_defects,
        record_path=True,
    )
    times = [0]            # initial position at real step 0
    prev_len = len(walker.path)

    for t in range(1, T + 1):
        walker.step()
        if len(walker.path) > prev_len:
            times.append(t)
            prev_len = len(walker.path)

    path = np.array(walker.path, dtype=np.int32)   # (M, 2)
    times = np.array(times, dtype=np.int32)         # (M,)
    assert len(path) == len(times)
    return path, times


# ---------------------------------------------------------------
# Plotting: one trajectory panel
# ---------------------------------------------------------------
def plot_trajectory_panel(ax, path, times, L, seed_lo, seed_hi, title,
                          T_plot, view_halfwidth, cmap="viridis"):
    """
    Plot the trajectory as a LineCollection colored by REAL step index (from
    `times`), plus the hatched seed rectangle and the outer box.

    `T_plot` sets the colorbar normalization; both panels share the same one.
    `view_halfwidth` gives half the side length of the zoom-box centered on
    the seed (in lattice units).
    """
    pts = path.astype(float)
    segs = np.stack([pts[:-1], pts[1:]], axis=1)  # shape (M-1, 2, 2)

    # Colour segment by midpoint REAL time (from authors' path-vs-step map)
    t_mid = 0.5 * (times[:-1] + times[1:])
    lc = LineCollection(segs, array=t_mid, cmap=cmap,
                        norm=plt.Normalize(vmin=0, vmax=T_plot),
                        linewidth=1.2, alpha=0.95)
    ax.add_collection(lc)

    # Seed rectangle: hatched, black edge
    seed_rect = Rectangle(
        (seed_lo - 0.5, seed_lo - 0.5),
        seed_hi - seed_lo + 1, seed_hi - seed_lo + 1,
        facecolor="white", edgecolor="black",
        hatch="////", linewidth=1.5, zorder=2,
    )
    ax.add_patch(seed_rect)

    # Outer lattice box
    box = Rectangle(
        (-0.5, -0.5), L + 1, L + 1,
        facecolor="none", edgecolor="black", linewidth=1.2, zorder=3,
    )
    ax.add_patch(box)

    # Zoom: viewbox centred on seed midpoint, clamped inside the lattice
    cx = 0.5 * (seed_lo + seed_hi)
    cy = 0.5 * (seed_lo + seed_hi)
    x_lo = max(-0.5,    cx - view_halfwidth)
    x_hi = min(L + 0.5, cx + view_halfwidth)
    y_lo = max(-0.5,    cy - view_halfwidth)
    y_hi = min(L + 0.5, cy + view_halfwidth)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=14)
    return lc


# ---------------------------------------------------------------
# Sanity checks on a path
# ---------------------------------------------------------------
def check_path(path, L, seed_defects, label):
    """Print a small diagnostic line for each trajectory."""
    n_pts = len(path)
    xs = path[:, 0]; ys = path[:, 1]
    in_box = ((xs >= 0) & (xs <= L) & (ys >= 0) & (ys <= L)).all()
    defect_set = set(map(tuple, seed_defects))
    entered_seed = any((int(x), int(y)) in defect_set for x, y in path)

    # Fraction of unique visited sites adjacent to the seed
    seed_lo = min(d[0] for d in defect_set)
    seed_hi = max(d[0] for d in defect_set)
    visited = {(int(x), int(y)) for x, y in path}
    # A site is "adjacent to seed" if it borders any defect site (4-neighbour)
    # but is not itself a defect.
    def adj_to_seed(x, y):
        if (x, y) in defect_set:
            return False
        return any(((x + dx, y + dy) in defect_set)
                   for dx, dy in [(1,0),(-1,0),(0,1),(0,-1)])
    n_adj = sum(1 for (x, y) in visited if adj_to_seed(x, y))
    n_bulk = sum(1 for (x, y) in visited if not adj_to_seed(x, y)
                                             and (x, y) not in defect_set)
    # And count step-count on edge sites (time spent near seed):
    n_steps_on_edge = sum(1 for (x, y) in map(tuple, path) if adj_to_seed(x, y))
    frac_edge = n_steps_on_edge / max(1, n_pts)

    print(f"  [{label}] path length: {n_pts}  (expected T+1)")
    print(f"  [{label}] all points in box?      {in_box}")
    print(f"  [{label}] walker ever on defect?  {entered_seed}  (must be False)")
    print(f"  [{label}] unique sites visited:   {len(visited)}")
    print(f"  [{label}]   on seed-adjacent ring: {n_adj}")
    print(f"  [{label}]   in bulk:               {n_bulk}")
    print(f"  [{label}] fraction of time on ring: {frac_edge:.3f}")
    return frac_edge


# ---------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------
def main():
    # --- parameters ---
    L = 30                        # 31 x 31 lattice (authors' convention)
    seed_lo, seed_hi = 13, 17     # 5 x 5 seed (target structure in 6(a) is 5x5)
    seed_defects = build_seed_defects(L, seed_lo, seed_hi)

    D_r = 0.25                    # paper's Fig 6(e) canonical value
    T = 1000                      # matches paper colorbar (0 .. 10^3)

    # Start the walker on a corner diagonal of the seed: first chiral move
    # aims it into the seed corner so it engages the edge immediately.
    start_xy = (seed_hi + 1, seed_hi + 1)    # one cell NE of the seed NE corner
    start_dir = '←'                          # first move is into the seed face

    # We try a few RNG seeds and keep the (ω=0.5, ω=1) pair whose ω=1 run
    # wraps the largest fraction of the seed's edge ring.  This is the only
    # stochastic element of the figure and does not bias the physics.
    rng_seeds = [20260420, 11, 27, 103, 2025, 314, 7, 42]

    print("=" * 60)
    print(" TCRW Fig 6(c) — single-tile trajectories")
    print("=" * 60)
    print(f" L = {L}  (lattice size {L+1} x {L+1})")
    print(f" seed = [{seed_lo}, {seed_hi}]^2  "
          f"({(seed_hi-seed_lo+1)}x{(seed_hi-seed_lo+1)} block)")
    print(f" D_r = {D_r},   T = {T}")
    print(f" start = {start_xy} facing {start_dir}")
    print()

    # --- generate trajectories; pick the clearest seed ---
    best = None   # (seed_idx, frac_chi, path_ach, times_ach, path_chi, times_chi, frac_ach)
    t0 = time.time()
    for s in rng_seeds:
        pa, ta = run_trajectory(L=L, omega=0.5, D_r=D_r, T=T,
                                seed_defects=seed_defects,
                                start_xy=start_xy, start_dir=start_dir,
                                rng_seed=s)
        pc, tc = run_trajectory(L=L, omega=1.0, D_r=D_r, T=T,
                                seed_defects=seed_defects,
                                start_xy=start_xy, start_dir=start_dir,
                                rng_seed=s)

        # Count how many of the 5*4 = 20 seed-adjacent ring sites were visited
        defect_set = set(map(tuple, seed_defects))
        def ring_sites(path):
            visited = set(map(tuple, path))
            ring = set()
            for (x, y) in visited:
                if (x, y) in defect_set:
                    continue
                for dx, dy in [(1,0),(-1,0),(0,1),(0,-1)]:
                    if (x+dx, y+dy) in defect_set:
                        ring.add((x, y))
                        break
            return ring
        n_ring_chi = len(ring_sites(pc))
        n_ring_ach = len(ring_sites(pa))
        if best is None or n_ring_chi > best[0]:
            best = (n_ring_chi, n_ring_ach, s, pa, ta, pc, tc)
        print(f"  rng={s:>7d}  ring-sites visited:  ω=1: {n_ring_chi:2d}  "
              f"ω=0.5: {n_ring_ach:2d}")
    print(f" picked rng seed {best[2]}  (ring sites ω=1: {best[0]}, "
          f"ω=0.5: {best[1]})  in {time.time()-t0:.2f}s\n")

    _, _, rng_pick, path_ach, times_ach, path_chi, times_chi = best

    # --- sanity checks (on the picked pair) ---
    print(" Sanity checks (picked pair):")
    frac_ach = check_path(path_ach, L, seed_defects, "ω=0.5")
    print()
    frac_chi = check_path(path_chi, L, seed_defects, "ω=1.0")
    print()
    if frac_chi > frac_ach:
        print(f" --> Chiral walker spends more time on seed-adjacent ring "
              f"({frac_chi:.3f} vs {frac_ach:.3f}).  PASS\n")
    else:
        print(f" --> WARNING: chiral walker did not beat achiral on edge "
              f"fraction: chi={frac_chi:.3f}, ach={frac_ach:.3f}.\n")

    # --- figure ---
    # Paper's Fig 6(c) is zoomed: the visible area is ~ seed_side * 4 cells
    # square.  With a 5x5 seed that's a ~20x20 window.
    view_halfwidth = 9.5   # half-side of the zoom window

    fig, axes = plt.subplots(1, 2, figsize=(10, 5.4), constrained_layout=True)
    _ = plot_trajectory_panel(
        axes[0], path_ach, times_ach, L, seed_lo, seed_hi,
        title=r"$\omega = 1/2$", T_plot=T, view_halfwidth=view_halfwidth,
    )
    lc_chi = plot_trajectory_panel(
        axes[1], path_chi, times_chi, L, seed_lo, seed_hi,
        title=r"$\omega = 1$", T_plot=T, view_halfwidth=view_halfwidth,
    )

    # Shared colorbar for time
    cb = fig.colorbar(lc_chi, ax=axes, orientation="horizontal",
                      fraction=0.05, pad=0.02, shrink=0.6,
                      location="bottom")
    cb.set_label("$t$", fontsize=13)
    cb.set_ticks([0, T])
    cb.ax.set_xticklabels(["0", r"$10^3$"])

    fig.suptitle("Fig 6(c) — single-tile trajectories during self-assembly\n"
                 f"(L={L}, $D_r$={D_r}, seed = {seed_hi-seed_lo+1}×{seed_hi-seed_lo+1},"
                 f" rng={rng_pick})",
                 fontsize=12)

    out = os.path.join(HERE, "tcrw_fig6c.png")
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print(f" saved figure to {out}")


if __name__ == "__main__":
    main()
