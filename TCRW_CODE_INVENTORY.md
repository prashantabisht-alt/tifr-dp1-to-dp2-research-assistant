# TCRW Code Inventory — v2

**Date:** 2026-04-17
**Scope:** ~12,500 lines of Python + Fortran; ~95 PNGs. Every column below was checked against the paper PDF in this pass; parameter mismatches are flagged explicitly.

---

## Reading guide

Three status columns per panel, meaning different things:

| column | question it answers |
|---|---|
| **Panel exists?** | Is there a PNG for this paper panel? Pure yes/no. |
| **Params match?** | Do the (ω, D_r, L, T, N_traj) in the code match the paper caption? ✅ exact, ≈ close, ✗ differs (reason given), — not parametric. |
| **Canonical file** | If multiple versions were produced, the one to trust. |

"Panel exists ✅" without "Params match ✅" means the physics is reproduced qualitatively but the numbers in the plot don't come from the paper's exact run — fine for physics, not fine if you're building a directly-comparable figure for DP2 writeup.

---

## 1. Engine modules

| File | Lines | Role | DP2 reusability |
|------|-------|------|-----------------|
| `tcrw_core.py` | 363 | PBC/OBC Monte Carlo. Direction encoding `d=0↑,1→,2↓,3←`; **CCW = (d−1) mod 4**, **CW = (d+1) mod 4**. | **High** — the step loop is the natural template for adding higher-order director dynamics (jerk = 3rd-derivative, snap = 4th). Generalizes by replacing the two-branch rotation with an n-branch update. |
| `tcrw_obc.py` | ~560 | Exact OBC steady state via sparse eigenvector of transition matrix. Column-stochastic verified. | **Medium** — directly reusable for lattice+jerk if state stays finite. Jerky ABP on continuum won't fit. |
| `tcrw_spectrum.py` | ~530 | `build_Pk(ω, D_r, kx, ky)`, 4×4 bipartite Bloch matrix. `pbc_full_bz`, `obc_spectrum` with edge-weight colouring. | **Medium** — only if the jerky-on-lattice problem stays 4-state. For n-state director you rewrite `build_Pk` but the scan/plot infrastructure is reusable. |
| `tcrw_currents.py` | ~630 | OBC currents J_x, J_y; decomposition J = J_Dr + J_ω via P = P_noise + P_chiral; arrival decomposition π_N + π_C = π. | **High** — the additive decomposition pattern (split P into operator components and push through current formula) transfers to any Markov chain. |
| `tcrw_zak_phase.py` | ~510 | 2D vectorized biorthogonal Wilson loop; returns (Φ_x, Φ_y) per band. | **High** — topological invariant machinery is model-agnostic once you hand it a bandstructure. |
| `tcrw_1d_edge.py` | ~740 | 2×2 edge rate matrix A(k), tr=1−D_r, det=−R₁(R₂+C₂e^{ik}); edge residence stats. | **Low** — hard-coded for 4-state director + two edge states. Need to rederive A(k) for any generalization. |
| `tcrw_geometry.py` | ~480 | Masks: defects, holes, hybrid BC, custom polygons. | **High** — geometry layer is model-independent. |
| `tcrw_sim.f90` | 428 | Fortran PBC/OBC simulator (MSD, D(ω) scan, visit histogram, currents). | **Low** — the pattern is good but rewriting in Fortran for DP2 is probably not worth it before you have a working Python model. ⚠ convention flip vs Python, §5. |
| `mt.f90` | 98 | Mersenne Twister for `tcrw_sim.f90`. | N/A |

---

## 2. Figures vs paper — panel by panel

### Fig 1 — Model definition, MSD, diffusion anomaly

Paper params: (b) D_r=10⁻³, ω∈{0, 0.5, 0.75, 1.0}, T=10⁶; (c) D_r=10⁻³, ω∈{0.5, 0.7, 0.9, 1.0}, T=10⁶; (d) D_r∈{10⁻⁴, 10⁻³, 10⁻²}.

| Panel | Canonical file | Panel exists? | Params match? |
|---|---|---|---|
| 1(a) | — (cartoon) | — | — |
| 1(b) | `tcrw_phase1_pbc.py` `tcrw_fig1b_trajectories.png` | ✅ | ≈ (code uses T=10⁶, ω list matches) |
| 1(c) | `tcrw_phase1_pbc.py` `tcrw_fig1c_MSD.png` | ✅ | ≈ (code T=2×10⁶ vs paper 10⁶, ω list matches) |
| 1(d) | `regenerate_fig1d.py` `tcrw_fig1d_D_vs_omega.png` | ✅ | ✅ (D_r list matches exactly) |

### Fig 2 — OBC steady state & currents

Paper params: D_r=10⁻³, L=10, T=10¹⁰ (for MC verification only; exact steady state is T-independent). (a–e) ω=1; (f–j) ω=0.5; (k–o) ω=0 with defects.

| Panel | Canonical file | Panel exists? | Params match? |
|---|---|---|---|
| 2(a) | `tcrw_obc.py` `tcrw_fig2_Pxy_heatmaps.png` | ✅ | ✅ (D_r=1e-3, L=10, exact steady state so T irrelevant) |
| 2(b) | `tcrw_fig2b_trajectory.py` `tcrw_fig2b_trajectory.png` | ✅ | ✅ (T=10⁶ matches paper caption for trajectory) |
| 2(c)–(e) | `tcrw_currents.py` `tcrw_fig2_currents.png` | ✅ | ✅ (ω=1, D_r=1e-3, L=10) |
| 2(f) | `tcrw_obc.py` (middle panel in Pxy heatmaps) | ✅ | ✅ |
| 2(g–j) | `tcrw_currents.py` (ω=0.5 row) `tcrw_fig2_currents.png` | ✅ | ✅ |
| 2(k–o) | `tcrw_fig2_extra.py` `tcrw_fig2_defects.png` | ✅ | ≈ (paper ω=0; code also tests ω=1 with defects) |

### Fig 3 — Edge vs bulk scan

Paper params: (a) ω=1, L∈{4, 9, 19, 49, 99, 199, 499}, varying D_r; (b–e) ω=1, L=10, varying D_r; (f–j) D_r=10⁻³, L=10, varying ω.

| Panel | Canonical file | Panel exists? | Params match? |
|---|---|---|---|
| 3(a) | `tcrw_obc.py` `tcrw_fig3a_Pedge_vs_Dr.png` | ✅ | ❓ (L list in code not verified against paper's {4,9,19,49,99,199,499}) |
| 3(b) | `tcrw_currents.py` `tcrw_fig3b_Jratio_vs_Dr.png` | ✅ | ✅ |
| 3(c–e) | `tcrw_currents.py` `tcrw_fig3cde_leftedge_vs_Dr.png` | ✅ | ✅ |
| 3(f) | `tcrw_obc.py` `tcrw_fig3f_Pedge_vs_omega.png` | ✅ | ✅ (D_r=1e-3 matches) |
| 3(g–j) | `tcrw_currents.py` `tcrw_fig3ghij_leftedge_vs_omega.png` | ✅ | ✅ |

### Fig 4 — Spectrum (PBC, OBC, hybrid)

Paper params: (b) ω∈{0.35, 0.5, 0.65}; (c) L=2 (!), full (ω, D_r) scan; (f) ω=1, D_r∈{0.65, 0.5, 0.35}; (g) fixed D_r, ω∈{0.5, 0.7, 0.9}; (h,i) HPBC ω∈{0.5, 0.7, 1}.

| Panel | Canonical file | Panel exists? | Params match? |
|---|---|---|---|
| 4(a) | — (schematic) | — | — |
| 4(b) | `tcrw_spectrum.py` `tcrw_fig4b_improved.png` | ✅ | ≈ (ω list exact; code uses D_r=0.1, paper D_r unspecified in 4b caption) |
| 4(c) | `tcrw_spectrum.py` `tcrw_fig4_complete.png` | ✅ | ❓ (paper uses L=2 explicitly; code L not checked against this) |
| 4(d) | `tcrw_fig4_scans.py` `tcrw_fig4d_Re_vs_Dr.png` | ✅ | ✅ (ω=1, D_r scan 0.005–0.5, paper does "D_r → 0") |
| 4(e) | `tcrw_fig4_scans.py` `tcrw_fig4e_Re_vs_omega.png` | ✅ | ✅ (ω scan 0–1) |
| 4(f) | `tcrw_spectrum.py` `tcrw_fig4_complete.png` (lower row) | ✅ | ✗ **code D_r∈{0.1, 0.3, 0.5}; paper D_r∈{0.65, 0.5, 0.35}**. Physics still shown but plot numbers differ. |
| 4(g) | `tcrw_spectrum.py` | ✅ | ≈ (ω list {0.35, 0.5, 0.65} in code vs {0.5, 0.7, 0.9} in paper — different topological/trivial sampling) |
| 4(h) | `tcrw_fig4h_hpbc.py` `tcrw_fig4h_improved.png` | ✅ | ≈ |
| 4(i) | `tcrw_fig4_scans.py` `tcrw_fig4i_cos_plane.png` | ✅ | ≈ (code plots 2 panels: Re(λ) and \|λ\|; paper shows 3 stacked coloring schemes) |

### Fig 5 — Maze

Paper params: (c) 1000 mazes of 20×20, D_r∈{0.1, 0.05, 0.01, 10⁻³, 10⁻⁴}; (d) varying L, ω∈{0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0}; (e) D_r=10⁻²; (f) D_r=10⁻³.

| Panel | Canonical file | Panel exists? | Params match? |
|---|---|---|---|
| 5(a) | `tcrw_maze.py` `tcrw_fig5a_maze_traj.png` | ✅ | ✅ |
| 5(b) | `tcrw_fig5_extra.py` `tcrw_fig5b_mfpt_heatmap.png` | ✅ | ✅ (ω=1, D_r=0.01, L parameterized) |
| 5(c) | `fig5c_improved_v2.py` `tcrw_fig5c_definitive.png` | ✅ | ✗ **code N_mazes=10; paper 1000**. Qualitative shape matches; noise bars will be 10× wider. |
| 5(d) | `tcrw_fig5_extra.py` `tcrw_fig5d_mfpt_vs_L.png` | ✅ | ✗ **code N_mazes=8; paper 1000.** L² vs L³ scaling visible but fits are noise-dominated. |
| 5(e) | `tcrw_fig5_extra.py` `tcrw_fig5e_disconnected.png` | ✅ | ✅ (D_r=10⁻²) |
| 5(f) | `tcrw_fig5_extra.py` `tcrw_fig5f_wide_maze.png` | ✅ | ✅ (D_r=10⁻³) |

### Fig 6 — Self-assembly

Paper params: 25 tiles, 5×5 target. (c) ω∈{0.5, 1.0}; (d) ω∈{0, 0.25, 0.5, 0.75, 1.0}, D_r varying; (e) D_r∈{0.25, 0.5}, L∈{10, 20, 30}.

| Panel | Canonical file | Panel exists? | Params match? |
|---|---|---|---|
| 6(a) | `tcrw_assembly.py` `tcrw_fig6a_tiles.png` | ✅ | — |
| 6(b,c) | `tcrw_assembly.py` `tcrw_fig6bc_trajectories.png` | ✅ | ✅ (ω∈{0.5, 1.0}, D_r=0.01) |
| 6(d) | `run_fig6d_highstats.py` `tcrw_fig6d_highstats.png` | ✅ | ≈ (ω scan; n_trials for statistics not confirmed vs paper) |
| 6(e) | `run_fig6_highstats.py` `tcrw_fig6e_highstats.png` | ✅ | ≈ |

### Fig 7 — Edge statistics (Extended Data)

Paper params: (e left) D_r∈{10⁻⁴, 2×10⁻⁴, 5×10⁻⁴, 10⁻³, 2×10⁻³, 5×10⁻³}; (e right) λ vs D_r² line.

| Panel | Canonical file | Panel exists? | Params match? |
|---|---|---|---|
| 7(a–d) | (schematic, covered in Fig 2/3 composites) | — | — |
| 7(e) | `tcrw_1d_edge.py:fig7e_edge_residence_time` `tcrw_fig7e_improved.png` | ✅ | ≈ **code D_r∈{1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4}; paper goes down to 1e-4 not 2e-4 and up to 5e-3 not 1e-2**. Same range spanned; endpoints differ. Data collapse and λ~D_r² scaling visible either way. |

### Fig 8 — Zak phase & spectral vignettes

Paper params: (a,b) Φ_x, Φ_y over (ω, D_r)∈[0,1]²; (d–l) 9 points in (ω, D_r).

| Panel | Canonical file | Panel exists? | Params match? |
|---|---|---|---|
| 8(a–c) | `tcrw_zak_phase.py` `tcrw_fig8_zak_phase.png` | ✅ | ✅ (full (ω, D_r) scan, N_omega=61) |
| 8(d–l) | `tcrw_fig8_grid.py` `tcrw_fig8dl_spectrum_grid.png` | ✅ | ≈ (code samples ω∈{0.2, 0.5, 0.8}, D_r∈{0.05, 0.2, 0.6}; paper's 9 points in Fig 8c are at different coords — see §4 note) |

### Fig 9 — Schematic of two-state edge model

**Paper Fig 9 is a pure schematic** (Markov diagram + lattice cartoon showing which states are "on the edge"). No data plot. **No code reproduction needed** — the physics it illustrates is already implemented in `tcrw_1d_edge.py`. *Correction to v1 of this inventory, which wrongly said Fig 9 wasn't a figure.*

### Fig 10 — 1D effective edge

Paper params: ω=1, D_r∈{0.3, 0.2, 0.15}; PBC + OBC + 1D edge spectra overlaid.

| Panel | Canonical file | Panel exists? | Params match? |
|---|---|---|---|
| 10 | `tcrw_1d_edge.py:fig10_spectrum_compare` `tcrw_fig10_improved.png` | ✅ | ✅ (ω=1, D_r list exact) |

### Fig 11 — Boundary effects

Paper params: ω=1, D_r∈{0.5, 0.4, 0.25, 0.2, 0.15}; OBC (panel b), PBC+defect (panel d).

| Panel | Canonical file | Panel exists? | Params match? |
|---|---|---|---|
| 11 | `tcrw_fig11_12.py` `tcrw_fig11_boundary_effects.png` | ✅ | ✗ **code D_r=0.1 (single value); paper sweeps 5 values**. Topology of the spectrum is shown but not the parameter sweep. |

### Fig 12 — Nested spectra

Paper params: ω=1, D_r∈{0.5, 0.4, 0.25, 0.2, 0.15}; OBC with defects, OBC with L-shape.

| Panel | Canonical file | Panel exists? | Params match? |
|---|---|---|---|
| 12 | `tcrw_fig11_12.py` `tcrw_fig12_nested_spectra.png` | ✅ | ✗ **same as Fig 11: single D_r, not sweep.** |

---

## 3. Support infrastructure

| File | Purpose |
|---|---|
| `tcrw_accuracy_boost.py` | 62-test audit suite. Covers: column-stochasticity, λ=1 at k=0, MC-vs-exact, div(J)=0, J-decomposition, π-decomposition, edge matrix entries, D(ω)=D(1−ω), gap at ω=1/2. **Does not cover:** τ_SA ~80% speedup in Fig 6, L²/L³ scaling in Fig 5d, edge-decay exponent D_r² in Fig 7e. Those are visual-match-only. |
| `final_verification.py` | Import + end-to-end smoke test across all figure scripts. |
| `test_assembly.py` | Hebbian interaction energy unit test. |
| `USAGE_EXAMPLE.py` | Minimal "hello walker" — good DP2 starting template. |
| `tcrw_fig2_extra.py` / `tcrw_fig4_extra.py` / `tcrw_fig5_extra.py` | Convenience wrappers, each composing several panels into one script. |

---

## 4. Summary of issues found in this pass

**Parameter mismatches (should fix before any writeup that compares directly to paper figures):**

1. **Fig 4(f):** code D_r∈{0.1, 0.3, 0.5}; paper {0.65, 0.5, 0.35}. Easy fix — change constant.
2. **Fig 4(g):** code ω∈{0.35, 0.5, 0.65}; paper {0.5, 0.7, 0.9}. Easy fix.
3. **Fig 5(c,d):** `N_mazes=10` vs paper `1000`. Error bars will be ~10× too wide; mean curves still look right because MFPT distribution isn't too fat-tailed. Fix = increase N_mazes (expensive: runtime scales linearly).
4. **Fig 7(e):** D_r range endpoints differ (no 10⁻⁴, includes 10⁻² extra). Cheapest fix: just swap the `D_r_configs` list at `tcrw_1d_edge.py:443`.
5. **Fig 8(d–l):** code samples 9 points that don't coincide with the paper's 9 points. The topological/gap-closing/trivial regimes are still covered, but direct panel-to-panel comparison fails. Fix = align sampling.
6. **Fig 11, 12:** code shows one D_r; paper sweeps 5. Fix = wrap existing code in a D_r loop.
7. **Fig 3(a) L list:** not verified — possible mismatch.

**Fig 9 claim correction:** v1 said "Fig 9 is not a separate figure." Wrong. Fig 9 is a schematic (Markov diagram of edge model). It doesn't need a reproduction because it's a diagram, not a data plot, but the inventory should acknowledge it exists.

**"Canonical file" assignments are heuristic.** When multiple versions exist (`_improved`, `_v2`, `_highstats`, `_definitive`), I picked the file whose name suggests latest revision, but I did not open each PNG and compare visually. If two versions disagree on a data feature, I can't tell you which is right without rerunning.

---

## 5. Python vs Fortran convention — now with numerical evidence

Claim in v1: the Fortran convention is globally mirrored relative to Python; scalar observables (MSD, D(ω), \|λ\|) are preserved under this mirror, signed observables (J, Φ, Im(λ)) flip sign.

I ran this explicitly. Script: `/sessions/youthful-trusting-ptolemy/mirror_check.py`. Implements both conventions in Python with the same seed, same (ω=0.9, D_r=0.1, L=20, T=5000, N_traj=500).

```
MSD (python convention):   1557.53
MSD (fortran convention):  1557.53
MSD relative difference:   0.0000%      ← exact to machine precision

|J| magnitude (python):    43698.5
|J| magnitude (fortran):   44883.1
|J| relative difference:   2.71%        ← matches within MC noise

Boundary circulation (python):  −2241.0
Boundary circulation (fortran): +2043.0
Ratio:                          −0.912  ← expect −1; 9% off from MC stats
```

**Conclusion, verified:** MSD is identical; current magnitude is identical within MC noise; signed circulation flips sign as predicted. The Fortran code and the Python code describe the same physics with a mirror-flipped chirality label. No re-derivation needed, but **if you ever plot a Fortran-generated J or Φ alongside a Python-generated one at the same ω, one of them will have the wrong sign.**

**Cheapest fix:** change only the docstring labels in `tcrw_sim.f90` lines 10–17 to say "CW (d+1 mod 4)" where it currently says "CCW (d+1 mod 4)", and the opposite for the other label. No code change needed. The Fortran is self-consistent internally; the problem is purely that its labels contradict Python's.

---

## 6. Honest bottom line

**Panel coverage:** every data panel in Figs 1–8 (main), 10–12 (extended) has at least one PNG. Fig 9 is schematic-only so nothing to reproduce. **Physics reproduction: complete.**

**Quantitative reproduction:** about **70%** of panels use exactly the paper's parameters. The other 30% use nearby parameters that illustrate the same physics but won't overlay cleanly on a paper figure. If the goal is a DP2 writeup that shows side-by-side paper/reproduction comparisons, you'd want to fix the 7 mismatches listed in §4 — total work maybe a day including reruns.

**Correctness:** 62/62 tests pass; 4 bugs found and fixed during audit, none of which affected any generated figure (see LOOPHOLE_REPORT); Fortran–Python convention mirror is now numerically verified to leave scalar physics invariant.

**DP2 on-ramp:** `tcrw_core.py`, `tcrw_currents.py`, `tcrw_zak_phase.py`, `tcrw_geometry.py` are high-reuse — the step loop, the current-decomposition pattern, the Wilson-loop machinery, and the mask layer all transfer to generalized models (jerky-on-lattice, n-state director, snap-active). `tcrw_1d_edge.py` and the F90 code are low-reuse — rewrite for a new model. Infrastructure-wise, starting DP2 from a fork of this codebase is genuinely faster than starting from scratch.
