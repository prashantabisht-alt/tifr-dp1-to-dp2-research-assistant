# Square-Lattice Chiral RTW Baseline

Purpose: keep the square-lattice version of our continuous-time chiral
run-and-tumble walker here.

This folder is for sanity checks against known square-lattice intuition:

- four directors instead of six,
- rates `v`, `gamma_plus`, `gamma_minus`, and `gamma_r`,
- even/odd diffusion formulas,
- comparison with the triangular six-director result.

Important distinction:

- Mallikarjun-Pal is the continuous-space chiral RTP reference.
- Wójcik-Kalz is the discrete-time square-lattice chiral-walk reference.
- `square_chiral_rtw.py` is our continuous-time square-lattice analogue, useful
  as a bridge, but not a literal implementation of Wójcik-Kalz.

