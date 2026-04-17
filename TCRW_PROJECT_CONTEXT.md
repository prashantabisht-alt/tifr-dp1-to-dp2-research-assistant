# TCRW Project: Full Context Summary

**Author:** Prashant Bisht, I.PhD Physics, TIFR Hyderabad (2024 intake)
**Supervisor:** Kabir Ramola
**Paper being reproduced:** Osat, Meyberg, Metson & Speck, "Topological Chiral Random Walker," arXiv:2602.12020
**Last updated:** 2026-04-07

---

## 1. Who I am and why this project exists

I'm a DP1-to-DP2 transition student. My DP1 ("Studies on single particle stochastic dynamics," Sep 2025) covered:
- **Ch 1:** Overdamped/underdamped Langevin, Fokker-Planck, large deviations
- **Ch 2:** Random walks (CSRW, CTRW, CS-CTRW, DSDT)
- **Ch 3:** Active Brownian Particles in 2D (exact MSDs, marginals)
- **Ch 4:** Undamped Jerky Active Particles (AOUP, MSD exponent sequences)
- **Ch 5:** Run-and-Tumble Particle

My coding toolkit: **Fortran 90** (Euler-Maruyama, MT19937 PRNG, Box-Muller), **gnuplot** for plotting, and now **Python** (NumPy/SciPy/Matplotlib).

Kabir sent my DP1 to **Stephy Jose** (Lowen group). She suggested five DP2 directions related to jerky active particles. Separately, Kabir gave me the TCRW paper to reproduce because it connects my random-walk + active-matter skills through non-Hermitian topology on a lattice, and is directly extendable toward DP2 (e.g., adding jerk/inertia to lattice models, which maps to Stephy's direction 4).

---

## 2. The TCRW model (paper rules, precisely)

2D square lattice walker at (x, y) with director d in {up, right, down, left} (encoded as 0,1,2,3).

**At each discrete time step:**
- With prob D_r: **NOISE STEP** -- walker stays put, director rotates CCW with prob omega, CW with prob (1-omega)
- With prob (1-D_r): **CHIRAL STEP** -- walker translates one step in direction d, THEN director rotates CW with prob omega, CCW with prob (1-omega) [OPPOSITE chirality to noise]
- **OBC blocking rule:** if chiral step would move walker off-grid, entire step is blocked (no move, no rotation)

**Parameters:**
- C1 = (1-omega)(1-D_r), C2 = omega(1-D_r), R1 = omega*D_r, R2 = (1-omega)*D_r
- omega = 0.5 is the achiral point (topological phase transition)
- D_r = rotational noise strength
- Topological phase: D_r < 0.5 AND omega != 0.5

---

## 3. What has been implemented (6 phases)

### Phase 1: PBC simulation engine (tcrw_core.py + tcrw_phase1_pbc.py)

**Core simulator** (`tcrw_core.py`, 363 lines):
- `simulate_tcrw_pbc(omega, D_r, L, T_steps, N_traj, ...)` -- vectorized over N_traj walkers, tracks unwrapped displacement for MSD, optional trajectory recording
- `simulate_tcrw_obc(omega, D_r, L, T_steps, N_traj, ...)` -- same with OBC blocking rule, optional visit histogram (`track_visits`), current tracking with J/J_Dr/J_omega decomposition (`track_currents`, `track_step_type`)
- `measure_diffusion_coeff(omega, D_r, ...)` -- linear fit to long-time MSD, returns D

**Figures produced:**
- `tcrw_fig1b_trajectories.png` -- Fig 1(b): sample trajectories for omega=0.5, 0.7, 1.0 at D_r=1e-3
- `tcrw_fig1c_MSD.png` -- Fig 1(c): MSD vs t confirming linear growth (normal diffusion)
- `tcrw_fig1d_D_vs_omega.png` -- Fig 1(d): D decreasing linearly with chirality

**Status:** Complete. Matches paper.

---

### Phase 2: OBC edge localization (tcrw_obc.py, 453 lines)

**Exact transition matrix approach:**
- `state_index(x, y, d, L)` -- maps (x,y,d) to linear index in 4L^2 space
- `build_transition_matrix(omega, D_r, L)` -- sparse column-stochastic (4L^2 x 4L^2) matrix
- `exact_steady_state(omega, D_r, L)` -- finds right eigenvector with eigenvalue 1 (dense for L<=30, sparse shift-invert for L>30); returns (Pxy, pi_vec) tuple
- `compute_edge_bulk_ratio(Pxy, L)` -- P_edge/P_bulk per-site ratio
- `mc_steady_state(omega, D_r, L, T_steps)` -- MC verification using long trajectory

**Figures produced:**
- `tcrw_fig2_Pxy_heatmaps.png` -- Fig 2(a),(f): P(X,Y) heatmaps for omega=1, 0.5, 0
- `tcrw_fig3a_Pedge_vs_Dr.png` -- Fig 3(a): P_edge/P_bulk vs D_r for various L (size-independent)
- `tcrw_fig3f_Pedge_vs_omega.png` -- Fig 3(f): P_edge/P_bulk vs omega (omega-independent)
- `tcrw_fig2_mc_verification.png` -- exact vs MC comparison

**Key result verified:** Edge localization depends on D_r but NOT on omega. This is non-obvious -- chirality affects the current, not the localization.

**Status:** Complete. Matches paper.

---

### Phase 3: Edge currents and decomposition (tcrw_currents.py, 633 lines)

**Exact current computation from steady state:**
- `build_split_matrices(omega, D_r, L)` -- builds P_noise and P_chiral separately
- `exact_currents(omega, D_r, L)` -- computes J, J_Dr, J_omega from pi and the split matrices
- `left_edge_currents(Jx, Jy, L)` -- extracts current vectors along x=0
- `current_angle(Jx, Jy)` -- angle from +x axis

**Figures produced:**
- `tcrw_fig2_currents.png` -- Fig 2(c)-(e),(f)-(h): J, J_omega, J_Dr vector fields for omega=1 and 0.5
- `tcrw_fig3b_Jratio_vs_Dr.png` -- Fig 3(b): |J_Dr|/|J_omega| vs D_r
- `tcrw_fig3cde_leftedge_vs_Dr.png` -- Fig 3(c)-(e): left-edge currents and angles vs D_r
- `tcrw_fig3ghij_leftedge_vs_omega.png` -- Fig 3(g)-(j): same vs omega

**Key physics:** Edge current runs OPPOSITE to walker chirality (like quantum Hall skipping orbits). For omega=0.5 (achiral), J_Dr = 0 exactly.

**Status:** Complete. Matches paper. Missing: Fig 2(i)-(o) defect/internal boundary panels (requires geometry modification).

---

### Phase 4: Transition matrix spectrum (tcrw_spectrum.py, 534 lines)

**PBC Fourier-space spectrum:**
- `build_Pk(omega, D_r, kx, ky)` -- 4x4 P(k) matrix (Eq 1 of paper)
- `pbc_band_structure(omega, D_r, Nk)` -- Re(lambda) along Gamma->X->M->Gamma
- `pbc_full_bz(omega, D_r, Nk)` -- all eigenvalues on Nk x Nk BZ grid
- `classify_eigenvalues(evals)` -- separates real pair (+/-a) and imaginary pair (+/-ib)

**OBC real-space spectrum:**
- `obc_spectrum(omega, D_r, L)` -- full diagonalization of (4L^2 x 4L^2) matrix, returns eigenvalues + edge-weight for each eigenvector

**Figures produced:**
- `tcrw_fig4b_pbc_bands.png` -- Fig 4(b): band structure for omega=0.35, 0.5, 0.65 showing gap closing
- `tcrw_fig4_obc_spectrum.png` -- Fig 4(c),(f)-(g): OBC spectrum in complex plane, colored by edge localization
- `tcrw_fig4_gap_closing.png` -- spectral gap vs omega at several D_r
- `tcrw_fig4_pbc_vs_obc.png` -- PBC bulk eigenvalues overlaid with OBC spectrum

**Key physics:** PBC eigenvalues are real (sublattice symmetry). OBC eigenvalues become complex -- imaginary parts encode edge current. Gap closes at omega=0.5 (topological phase transition). Edge modes appear inside the bulk gap for D_r < 0.5.

**Status:** Complete for core panels. Missing: Fig 4(a) state labeling diagram, Fig 4(d)-(e) real spectrum vs D_r/omega, Fig 4(h)-(i) hybrid BC + circle plot.

---

### Phase 5: Topological invariant -- Zak phase (tcrw_zak_phase.py, 511 lines)

**Biorthogonal Berry phase machinery:**
- `biorthogonal_eig(P)` -- left and right eigenvectors with biorthogonal normalization
- `wilson_loop_single_band(band_selector, omega, D_r, kperp, direction, Nk)` -- discrete Berry phase via Wilson loop
- `wilson_loop_multiband(band_selectors, omega, D_r, kperp, direction, Nk)` -- non-Abelian Wilson loop for band groups
- `select_top_band`, `select_bottom_band`, `select_imag_pos`, `select_imag_neg` -- band selectors
- `compute_zak_phases(omega, D_r, Nk)` -- Phi_x and Phi_y for top band
- `compute_zak_phases_multiband(omega, D_r, Nk)` -- 2-band Wilson loop version

**Figures produced:**
- `tcrw_fig8_zak_phase.png` -- Fig 8(a)-(c): Phi_x, Phi_y in (omega, D_r) plane with spectral gap boundary
- `tcrw_fig8_linecuts.png` -- linecuts through phase diagram

**Key result:** Phi_x = Phi_y = pi for D_r < 0.5 and omega != 0.5 (topological). Transitions exactly where spectral gap closes.

**Status:** Complete for main panels. Missing: Fig 8(d)-(l) nine individual spectrum panels at specific parameter points.

---

### Phase 6A: 1D effective edge model (tcrw_1d_edge.py, 723 lines)

**2-state Markov chain on edge (Eq 3 of paper):**
- States: left-arrow (into wall, d=3) and down-arrow (along wall, d=2) on left edge x=0
- `edge_rate_matrix(omega, D_r, k)` -- 2x2 A(k) matrix
- `edge_eigenvalues(omega, D_r, k)` -- analytical eigenvalues (Eq 4)
- `edge_spectrum(omega, D_r, Nk)` -- sweep k in [-pi, pi)
- `is_edge_state(xx, yy, dd)` -- checks boundary AND director state (wall-facing or edge-following)
- `measure_edge_residence_batch(omega, D_r, L, N_walkers, T_steps)` -- vectorized residence time measurement

**Rate matrix:**
```
A(k) = | 1-D_r         R1            |
       | R2 + C2*e^ik  0             |
```
where R1 = omega*D_r, R2 = (1-omega)*D_r, C2 = omega*(1-D_r).

**Loop splitting:** discriminant at k=pi is (1-D_r)(1-5*D_r). Vanishes at D_r = 1/5 = 0.2. Below 0.2: two loops in complex plane. Above 0.2: one loop.

**Figures produced:**
- `tcrw_fig10_pbc_obc_1d.png` -- Fig 10: PBC vs OBC vs 1D edge model spectra overlaid (3 panels: D_r=0.3, 0.2, 0.15)
- `tcrw_fig7e_edge_residence.png` -- Fig 7(e): P(tau_edge) distribution and exponential decay rate ~ D_r^2
- `tcrw_1d_edge_diagnostic.png` -- diagnostic varying D_r and omega

**Status:** Complete. 1D edge eigenvalues trace OBC edge band with mean distance < 0.05. Residence time scaling tau ~ 1/D_r^2 confirmed for D_r >= 0.01. Statistics insufficient for D_r < 5e-4 in Python (needs Fortran with 10^10 steps).

---

### Phase 6B: Maze solving (tcrw_maze.py, 532 lines)

- `generate_maze_prim(L, seed)` -- Prim's algorithm on (2L+1) x (2L+1) grid, entrance (1,0), exit (2L-1, 2L)
- `run_tcrw_maze(maze, entrance, exit_pos, omega, D_r, ...)` -- single walker through maze with OBC blocking
- `measure_mfpt(maze, entrance, exit_pos, omega, D_r, N_trials)` -- average MFPT
- Visit-frequency heatmap rendering with RGB composite image

**Figures produced:**
- `tcrw_fig5a_maze_traj.png` -- Fig 5(a): visit-frequency heatmaps for omega=0, 0.5, 1 (L=12)
- `tcrw_fig5c_mfpt_vs_omega.png` -- Fig 5(c): MFPT vs omega for D_r = 0.1, 0.01, 0.001
- `tcrw_fig5d_mfpt_vs_L.png` -- Fig 5(d): MFPT scaling (L^2 for chiral vs L^3 for achiral)

**Status:** Complete for core panels. Statistics limited by Python speed. Missing: Fig 5(b) MFPT heatmap by start position, Fig 5(e) disconnected mazes, Fig 5(f) wide-passage mazes.

---

## 4. Equation-level verification (10/10 passed)

| Check | Description | Result |
|-------|-------------|--------|
| A | Rate matrix A(k) matches paper Eq 3 | exact (err=0) |
| B | Eigenvalue formula vs numpy eig (7 test cases) | err < 1e-15 |
| C | Loop splitting at D_r=0.2: disc=(1-D_r)(1-5D_r)=0 | exact zero |
| D | All 1D edge eigenvalues sub-unitary (|lambda|<1) | confirmed |
| E | P(k=0) column-stochastic | exact |
| F | OBC transition matrix column-stochastic | exact |
| G | Steady state P@pi = pi | residual 2.5e-16 |
| H | Topology transition: two loops (D_r<0.2) -> one loop (D_r>0.2) | disc signs correct |
| I | Absorption rates conserve probability for all (omega, D_r) | exact |
| J | omega=1 special: left-arrow never absorbs, A(k=0) correct | exact |

---

## 5. Complete figure inventory: what's done vs what's missing

### Main text figures

**Fig 1 -- Chiral Random Walker:**
- (a) Schematic -- SKIP (diagram)
- (b) Trajectories -- DONE
- (c) MSD -- DONE
- (d) D vs omega -- DONE

**Fig 2 -- Chiral Edge Current:**
- (a) P(X,Y) heatmap -- DONE
- (b) Sample trajectory 10^6 steps -- MISSING
- (c)-(e) J, J_omega, J_Dr for omega=1 -- DONE
- (f)-(h) Same for omega=0.5 -- DONE
- (i)-(o) Defect/internal boundary panels -- MISSING (needs geometry modification to OBC builder)

**Fig 3 -- Impact of D_r and omega:**
- (a) P_edge/P_bulk vs D_r -- DONE
- (b) |J_Dr|/|J_omega| vs D_r -- DONE
- (c)-(e) Left-edge currents vs D_r -- DONE
- (f) P_edge/P_bulk vs omega -- DONE
- (g)-(j) Left-edge currents vs omega -- DONE

**Fig 4 -- Topological Spectrum:**
- (a) State labeling -- SKIP (diagram)
- (b) PBC band structure (gap closing) -- DONE
- (c) OBC spectrum complex plane (L=2) -- DONE
- (d) Real spectrum vs D_r for omega=1 -- MISSING
- (e) Real spectrum vs omega for fixed D_r -- MISSING
- (f)-(g) OBC complex plane L=10, varying D_r and omega -- DONE
- (h) Hybrid PBC/OBC boundary conditions -- MISSING
- (i) Band structure on circle (cos kx, cos ky) -- MISSING

**Fig 5 -- Maze Solving:**
- (a) Trajectories through maze -- DONE (visit heatmap, L=12)
- (b) MFPT heatmap by starting position -- MISSING
- (c) MFPT vs omega -- DONE
- (d) MFPT vs L scaling -- DONE
- (e) Disconnected mazes -- MISSING
- (f) Wide-passage mazes -- MISSING

**Fig 6 -- Self-Assembly:**
- (a)-(e) All panels -- NOT IMPLEMENTED (entirely separate multi-particle simulation with patchy tiles)

### Extended data figures

**Fig 7 -- Model Sketch:**
- (a)-(d) Schematics -- SKIP (diagrams)
- (e) P(tau_edge) -- DONE (limited stats at very small D_r)

**Fig 8 -- Zak Phase:**
- (a)-(c) Phase diagram + gap boundary -- DONE
- (d)-(l) Nine spectrum panels at specific points -- MISSING

**Fig 9 -- 1D Edge Model Sketch:**
- (a)-(b) Schematics -- SKIP (diagrams)

**Fig 10 -- Time-Scale Separation:**
- PBC vs OBC vs 1D edge overlay -- DONE

**Fig 11 -- Effect of Boundaries on Spectrum:**
- (a)-(d) All panels -- NOT IMPLEMENTED (needs internal defect geometry)

**Fig 12 -- Disconnected Boundaries and Nested Spectrum:**
- (a)-(d) All panels -- NOT IMPLEMENTED (needs multi-defect geometry)

### Summary count

Out of ~70+ data panels across the paper:
- **~35-40 DONE** (all core physics)
- **~15 MISSING but straightforward** (parameter scans of existing code, extra spectrum panels)
- **~10 MISSING, moderate effort** (defect geometries, hybrid BC, extra maze variants)
- **~5 MISSING, significant effort** (self-assembly = entirely new simulation)

---

## 6. Codebase structure

```
tcrw_core.py          (363 lines)  -- Walker engine: PBC + OBC simulators, diffusion measurement
tcrw_phase1_pbc.py    (177 lines)  -- Fig 1: trajectories, MSD, D vs omega
tcrw_obc.py           (453 lines)  -- Exact transition matrix, steady state, edge/bulk ratio
tcrw_currents.py      (633 lines)  -- Current decomposition J = J_Dr + J_omega
tcrw_spectrum.py      (534 lines)  -- Fourier-space P(k), PBC bands, OBC spectrum
tcrw_zak_phase.py     (511 lines)  -- Biorthogonal Berry phase, Wilson loop, phase diagram
tcrw_1d_edge.py       (723 lines)  -- 2-state edge model, edge residence time, Fig 10
tcrw_maze.py          (532 lines)  -- Prim's maze, TCRW solver, MFPT measurement
                     -----
Total:               3926 lines
```

All files in: `/TIFR DP1 to DP2 Research Assistant/` (the workspace folder)

**Import dependencies:**
```
tcrw_phase1_pbc.py  -->  tcrw_core.py
tcrw_obc.py         -->  (standalone, uses scipy.sparse)
tcrw_currents.py    -->  tcrw_obc.py
tcrw_spectrum.py    -->  tcrw_obc.py
tcrw_zak_phase.py   -->  tcrw_spectrum.py
tcrw_1d_edge.py     -->  tcrw_spectrum.py, tcrw_core.py
tcrw_maze.py        -->  (standalone)
```

**Direction encoding throughout:** d=0(up), 1(right), 2(down), 3(left). DX=[0,1,0,-1], DY=[1,0,-1,0].

---

## 7. Key technical decisions and conventions

**Matrix convention:** Paper uses row-vector convention p^T A. Our code computes eigenvalues of A directly. Eigenvalues of A and A^T are identical, so this is fine.

**P(k) is column-stochastic at k=0:** columns sum to 1. At general k, it's the Fourier-space transfer matrix and column sums != 1 (this is expected, not a bug).

**OBC blocking:** if chiral step would exit grid, ENTIRE step is blocked (no move AND no rotation). This is crucial for correct edge dynamics.

**Zak phase computation:** Uses biorthogonal left/right eigenvectors with discrete Wilson loop. Band selection by eigenvalue magnitude. The 2-band non-Abelian Wilson loop gives cleaner results than single-band near the gap closing.

**Edge-state criterion for residence time:** Walker is "on edge" only when on boundary AND in a wall-facing or edge-following director state. Not just being on the boundary -- the director must couple the walker to the edge dynamics.

```python
def is_edge_state(xx, yy, dd):
    left   = (xx == 0)   & ((dd == 3) | (dd == 2))   # left or down
    right  = (xx == L-1)  & ((dd == 1) | (dd == 0))   # right or up
    bottom = (yy == 0)    & ((dd == 2) | (dd == 1))   # down or right
    top    = (yy == L-1)  & ((dd == 0) | (dd == 3))   # up or left
    return left | right | bottom | top
```

---

## 8. Known limitations and what needs Fortran

1. **Fig 7(e) at small D_r:** P(tau_edge) for D_r = 5e-4 and 2e-4 needs ~10^10 steps to get enough edge-residence events. Python is too slow. Plan: port the OBC walker loop to Fortran.

2. **Fig 5(d) at large L:** MFPT for achiral walkers (omega=0.5, tau ~ L^3) at L=15,20 requires many maze trials. Fortran would help.

3. **Fig 5(b) MFPT heatmap:** Needs MFPT measured from every starting position x starting director -- a 4*L^2 sweep. Feasible in Fortran.

4. **All maze statistics** would benefit from Fortran for better averaging (1000+ mazes instead of 10).

---

## 9. Bugs encountered and fixed (reference for debugging)

- **OBC rotation convention:** Initially had CCW/CW swapped in the noise step. Fixed by carefully matching paper's convention: noise CCW with prob omega, chiral CW with prob omega.
- **Edge residence time too short:** Initially counted ALL boundary sites. Fixed with `is_edge_state()` checking director state.
- **Fig 10 colorbar overlap:** `tight_layout()` caused issues. Fixed with `fig.subplots_adjust(right=0.88)` + `fig.add_axes()`.
- **Fourier matrix column sums flagged as "not stochastic":** This is EXPECTED for P(k) at k != 0. Only the real-space OBC matrix is column-stochastic.
- **exact_steady_state return type:** Returns `(Pxy, pi_vec)` tuple, not a single array.
- **Achiral maze walker timeout:** omega=0.5 with L=15, D_r=0.001 exceeds 5M steps. Reduced L or increased max_steps.

---

## 10. DP2 connections

### Stephy Jose's five suggested directions (from email after reading my DP1):

1. **Jerky harmonic oscillator in 2D** -- extend Lowen's 1D to 2D, derive phase diagram. "Not particularly exciting but straightforward."
2. **Collective effects / jerk + MIPS** -- Stephy is actively working on this. "Tricky, highly unstable." Avoid competing.
3. **Snap-active particle (4th derivative)** -- MSD exponent vs highest-order derivative. Open question.
4. **Jerk + Tailleur's lattice model** -- exact calculations possible. Even inertia on a lattice is nontrivial.
5. **Active harmonic solids** -- jerk's effect on collective excitations ("entropons").

### How TCRW connects to DP2:

The TCRW is essentially a **lattice model with internal states and chirality** -- exactly the setting of Stephy's direction 4. Having reproduced this paper, I understand:
- How non-Hermitian topology arises from discrete stochastic dynamics
- Bulk-boundary correspondence in active systems
- Spectral analysis of non-Hermitian transfer matrices
- The role of chirality vs noise in creating edge states

**Most natural DP2 extension:** Add jerk/inertia/memory to the TCRW. This combines my DP1 jerky-particle expertise with the TCRW framework. The question "does higher-order dynamics change the topological phase?" is original and connects both halves of my work.

---

## 11. What to do next (prioritized)

**Quick wins (< 1 day each):**
- Fig 4(d)-(e): real spectrum vs D_r and omega (just parameter scans of existing `obc_spectrum`)
- Fig 8(d)-(l): spectrum panels at 9 specific (omega, D_r) points (same)
- Fig 2(b): sample OBC trajectory overlay (already have the simulator)

**Moderate effort (1-3 days each):**
- Port OBC walker to Fortran for 10^10 step statistics (Fig 7e, Fig 5 scaling)
- Add defect/internal boundary support to `build_transition_matrix` (for Fig 2(i)-(o), Fig 11, Fig 12)
- Fig 5(b) MFPT heatmap by starting position
- Fig 5(e)-(f) disconnected mazes, wide-passage mazes

**Significant effort (1+ weeks):**
- Fig 6 self-assembly (entirely new multi-particle simulation with patchy tile interactions)

**DP2 research directions:**
- Start formulating the "jerk on a lattice" model
- Literature review on inertial lattice models (Tailleur et al.)
- Prototype: add persistence/memory to TCRW and see how spectrum changes

---

## 12. Output files produced

```
tcrw_fig1b_trajectories.png       -- Fig 1(b): PBC trajectories
tcrw_fig1c_MSD.png                -- Fig 1(c): MSD vs t
tcrw_fig1d_D_vs_omega.png         -- Fig 1(d): D vs omega
tcrw_fig2_Pxy_heatmaps.png       -- Fig 2(a),(f): P(X,Y) heatmaps
tcrw_fig2_currents.png            -- Fig 2(c)-(e),(f)-(h): current vector fields
tcrw_fig2_mc_verification.png     -- MC vs exact verification
tcrw_fig3a_Pedge_vs_Dr.png        -- Fig 3(a): edge/bulk ratio vs D_r
tcrw_fig3b_Jratio_vs_Dr.png       -- Fig 3(b): current ratio
tcrw_fig3cde_leftedge_vs_Dr.png   -- Fig 3(c)-(e): left-edge vs D_r
tcrw_fig3f_Pedge_vs_omega.png     -- Fig 3(f): edge/bulk ratio vs omega
tcrw_fig3ghij_leftedge_vs_omega.png -- Fig 3(g)-(j): left-edge vs omega
tcrw_fig4b_pbc_bands.png          -- Fig 4(b): PBC band structure
tcrw_fig4_obc_spectrum.png        -- Fig 4(c),(f)-(g): OBC complex plane
tcrw_fig4_gap_closing.png         -- spectral gap vs omega
tcrw_fig4_pbc_vs_obc.png          -- PBC vs OBC overlay
tcrw_fig5a_maze_traj.png          -- Fig 5(a): maze visit heatmaps
tcrw_fig5c_mfpt_vs_omega.png      -- Fig 5(c): MFPT vs omega
tcrw_fig5d_mfpt_vs_L.png          -- Fig 5(d): MFPT scaling
tcrw_fig7e_edge_residence.png     -- Fig 7(e): P(tau_edge) distribution
tcrw_fig8_zak_phase.png           -- Fig 8(a)-(c): Zak phase diagram
tcrw_fig8_linecuts.png            -- phase diagram linecuts
tcrw_fig10_pbc_obc_1d.png         -- Fig 10: PBC vs OBC vs 1D edge
tcrw_1d_edge_diagnostic.png       -- 1D edge spectrum diagnostic
tcrw_phase3_verification.png      -- Phase 3 verification
```

---

## 13. How to run any figure

Each figure-generating function is self-contained. From the project directory:

```python
# Example: reproduce Fig 4(b)
import sys; sys.path.insert(0, '.')
from tcrw_spectrum import fig4b_band_structure
fig4b_band_structure()

# Example: reproduce Fig 8 phase diagram
from tcrw_zak_phase import fig8_phase_diagram
fig8_phase_diagram()

# Example: reproduce Fig 10
from tcrw_1d_edge import fig10_pbc_obc_1d_overlay
fig10_pbc_obc_1d_overlay()
```

Or run any file directly: `python tcrw_spectrum.py` (if it has `__main__` block).
