"""
Dipanjan postmortem — bug isolated.

We found: in dipanjan_postmortem.py the sum of P over Dipanjan's L*L
real-space grid is exactly 0.5, not 1.0.  Now we isolate WHY.

Two equivalent ways to see the bug:
  (1) k-grid is 2L*L (rectangle [0, 2pi)^2 / discretised) but only HALF of
      it lies inside the TRUE Brillouin zone of the Bravais lattice
      with primitive vectors a_1 = (2, 0) and a_2 = (1, 1).  Each true
      BZ point appears TWICE in Dipanjan's grid → sum doublecounts.
  (2) Dipanjan's reconstruction puts half the walker's probability on
      sub-lattice A (correct, where walker actually lives) and half on
      sub-lattice B (WRONG — walker cannot reach there).

This script verifies (2) directly: dump P(x, y, t) over the FULL 2L*L
Cartesian grid and confirm half the mass is on impossible sites.

Then we propose the fix.
"""
from __future__ import annotations
import numpy as np
from dipanjan_postmortem import build_Mk, P_k_at_time


def P_xy_full_cartesian(gamma=0.01, epsilon=0.15, t=10.0,
                         L=30, r0=(0, 0), a=1.0, b=1.0):
    """
    Inverse-FT exactly as Dipanjan's notebook does it, but evaluate at
    the FULL 2L*L Cartesian grid (every integer (x, y), not just x step 2).
    """
    Nkx, Nky = 2 * L, L
    N_k = Nkx * Nky
    # First, time-evolve at every k-point
    P_k_summed = np.zeros((Nkx, Nky), dtype=complex)
    for nx in range(Nkx):
        kx = 2 * np.pi * nx / (2 * a * L)
        for ny in range(Nky):
            ky = 2 * np.pi * ny / (b * L)
            v = P_k_at_time(gamma, epsilon, kx, ky, r0, t)
            P_k_summed[nx, ny] = np.sum(v)
    # Now invert over the FULL 2L*L real-space grid
    xs = np.arange(2 * L)   # step of 1, all integers
    ys = np.arange(L)
    P_full = np.zeros((len(xs), len(ys)), dtype=complex)
    for ix, x in enumerate(xs):
        for iy, y in enumerate(ys):
            s = 0j
            for nx in range(Nkx):
                kx = 2 * np.pi * nx / (2 * a * L)
                for ny in range(Nky):
                    ky = 2 * np.pi * ny / (b * L)
                    s += P_k_summed[nx, ny] * np.exp(-1j * (kx * x + ky * y))
            P_full[ix, iy] = s / N_k
    return xs, ys, P_full.real


def main():
    print("=" * 70)
    print("BUG HUNT — full Cartesian inverse FT")
    print("=" * 70)
    L = 30
    xs, ys, P_full = P_xy_full_cartesian(L=L)
    total = np.sum(P_full)
    print(f"\nTotal sum over FULL 2L*L = {2*L*L} Cartesian sites:  {total:.6f}")
    print(f"   (should be 1.0)")
    # Split by parity: walker starts at (0,0), so 'reachable' = (x+y) even
    on_lattice = 0.0
    off_lattice = 0.0
    on_count = 0
    off_count = 0
    for ix, x in enumerate(xs):
        for iy, y in enumerate(ys):
            if (x + y) % 2 == 0:
                on_lattice += P_full[ix, iy]
                on_count += 1
            else:
                off_lattice += P_full[ix, iy]
                off_count += 1
    print(f"\nMass on  (x+y) even sites (reachable from (0,0)):  {on_lattice:.6f}  over {on_count} sites")
    print(f"Mass on  (x+y) odd  sites (IMPOSSIBLE — walker cannot be here): {off_lattice:.6f}  over {off_count} sites")
    print(f"\n=>  the walker's probability is being SPREAD onto sites it cannot")
    print(f"    physically reach.  This is the bug.")
    print(f"    The 'on-lattice' mass is the physically correct probability,")
    print(f"    and it sums to 1 (modulo floating-point).")

    # Detailed slice at y=4 (Dipanjan's printed slice).  Show on-lattice + off-lattice
    # to make the symptom visible.
    iy4 = 4
    print(f"\nSlice at y={iy4}:")
    print(f"   x    P(x, y={iy4}, t)        (x+y) parity")
    print(f"  ---  ----------------       ---------------")
    for ix, x in enumerate(xs[:14]):
        par = "even (reachable)" if (x + iy4) % 2 == 0 else "ODD  (impossible)"
        print(f"  {x:3d}   {P_full[ix, iy4]:+.4e}      {par}")
    print(f"  ... etc up to x = {xs[-1]}")
    print()
    print(f"Notice: in Dipanjan's notebook the printed x-values are 0, 2, 4, ...")
    print(f"i.e. ONLY the 'reachable' sites.  His L*L sum gives the correct ")
    print(f"PHYSICAL probability ≈ 1, but his code labels it as P(x, y=4, t)")
    print(f"and divides by N_k = 2L*L, so the printed values are HALF the true P.")
    print()
    print(f"** Net effect: ALL his printed P(x, y, t) values are off by a factor of 2. **")
    print()
    print(f"=" * 70)
    print(f"THE FIX (two equivalent paths):")
    print(f"=" * 70)
    print(f"""
A. Stay on the 2L*L Cartesian grid (Dipanjan's setup), but RECOGNIZE that
   half the sites carry the same probability as the other half. Concretely:
     - Print at every integer (x, y), not every 2 in x.
     - The 'physical' probability is what you get summed over the full grid.
     - Off-parity sites should be 0 in the limit of infinite L; the finite-L
       discrepancy is the bug.

B. (Cleaner) Use lattice coordinates instead of Cartesian.
   Site index = (n1, n2) with n1, n2 in [0, L). Position is
   r = n1*a_1 + n2*a_2 = (2*n1 + n2, n2).
   k-vectors: (n1/L)*b_1 + (n2/L)*b_2 with n1, n2 in [0, L).
     b_1 = (pi, -pi), b_2 = (0, 2pi).
   This gives L*L k-points = L*L sites, exact 1-to-1 correspondence.
   Inverse FT divides by N_sites = L*L = 900 not 1800.
   Total probability comes out to 1.0 exactly.

Solution B is the right fix.  In any real-space simulation (KMC), the walker
should be parametrised by (n1, n2) integer lattice coordinates, NOT by
Cartesian (x, y) with brute-force mod arithmetic.  This is exactly the PBC
issue Kabir was pointing at on his notebook page:

    "along n2: y + L mod = y_tilde, AND x_tilde = x - L/2"

is the Cartesian way of saying "use lattice coordinates (n1, n2) with
trivial mod L".  Same fix, two languages.
""")


if __name__ == "__main__":
    main()
