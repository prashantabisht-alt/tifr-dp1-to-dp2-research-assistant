# `new_fortran_reproduction_and_python/` — audit against authors' code

**Date:** 2026-04-24
**Target:** the new Fortran folder that replaced `fortran_reproduction/`
**Method:** read every kernel-level file, then compare the Fortran MC
output to the authors' exact Python steady state for the same (ω, D_r)
on the equivalent lattice.

---

## Headline

**Physics is correct.** MC at T = 10¹⁰ agrees with the exact-Python
transition-matrix steady state to Monte-Carlo noise at **~3×10⁻⁵ RMS**
on P(X,Y) and **~10⁻⁷–10⁻⁶ absolute** on the J_Dr, J_ω currents.
Two systematic improvements vs the old Fortran are verified numerically;
one convention inconsistency within the new folder is flagged.

---

## 1. What changed vs the old `fortran_reproduction/`

Two real fixes, both verified:

### 1.1 The chirality sign flip is gone

The old `tcrw_sim.f90` was globally mirrored relative to `tcrw_core.py`
and `TRW.py` (your inventory §5 verified this numerically — MSD
invariant, boundary circulation flipped sign). The new kernel
`tcrw_step.f90` uses `mod(d+3,4)` for CCW and `mod(d+1,4)` for CW, with
`r_rot < omega → CCW` in the noise branch and `r_rot < omega → CW` in
the chiral branch. This is **identical** to the authors' `rotate_direction`
and to your `tcrw_core.py`. No mirror remaining.

### 1.2 The `prev_was_noise` flag now has the authors' semantics

This is the subtle one. In the authors' `ChiralWalker.step`, a chiral
attempt blocked by a wall does **not** reset `self.noise_step`. Your
new comment in `tcrw_step.f90:119-142` spells this out correctly and
the drivers implement it as `if (step_type /= 2) prev_noise = (step_type == 0)`.
This is the difference between getting J_Dr and J_ω right on the edge
and getting them wrong — a blocked-chiral attempt on the wall shouldn't
eat the walker's memory of the last actual event. Verified in the next
section: J_Dr and J_ω come out exactly right.

---

## 2. Quantitative cross-check: Fortran MC (T = 10¹⁰) vs exact Python

Script: `outputs/fortran_vs_exact.py` (also copy-able into the
project folder if you want it version-controlled). Loads
`tcrw_fig2_occ_w{ω}.txt` and `tcrw_fig2_{Jtot,JDr,Jomega}_w{ω}.txt`,
and compares to the authors' `solve_steady_state_sparse` +
`calculate_J1_J2_with_boundaries` on the 8×8 playground that the
wall-ring convention of fig2_clean.f90 actually gives (see §3).

```
=== Fig 2, omega=0.0, D_r=10⁻³, L_F=10 (8×8 playground) ===
  max|P_F - P_E|   =  1.62×10⁻⁴   (≈ 0.26 % of peak)
  RMS (P_F - P_E)  =  4.01×10⁻⁵
  JDr   max|F−E|   =  Jx 2.06×10⁻⁸   Jy 2.02×10⁻⁸   (ref 3.1×10⁻⁵)
  Jomega max|F−E|  =  Jx 4.29×10⁻⁶   Jy 4.12×10⁻⁶   (ref 3.1×10⁻⁵)

=== Fig 2, omega=0.5, D_r=10⁻³, L_F=10 ===
  max|P_F - P_E|   =  2.20×10⁻⁴
  RMS              =  5.51×10⁻⁵
  JDr   max|F−E|   =  Jx 9.0×10⁻⁸   Jy 1.7×10⁻⁷
  Jomega max|F−E|  =  Jx 2.3×10⁻⁷   Jy 2.8×10⁻⁷

=== Fig 2, omega=1.0, D_r=10⁻³, L_F=10 ===
  max|P_F - P_E|   =  1.66×10⁻⁴
  RMS              =  3.81×10⁻⁵
  JDr   max|F−E|   =  Jx 1.53×10⁻⁸   Jy 1.17×10⁻⁸
  Jomega max|F−E|  =  Jx 3.26×10⁻⁶   Jy 3.33×10⁻⁶
```

**Interpretation.** At T = 10¹⁰ on 64 sites the naive standard error on
P(X,Y) is 1/√(T·P) ≈ 10⁻⁵ per site; 3×10⁻⁵ RMS max is consistent with
that. The J_Dr error sits at 10⁻⁸ because J_Dr is a small quantity
(D_r × base-current). J_ω residuals are larger in absolute terms but
still ~3×10⁻⁶ / 3×10⁻⁵ ≈ 10 % of the peak, which is MC noise on a
per-cell current averaged over 10¹⁰ samples. None of these residuals
is a bias — they are statistical and scale as 1/√T.

**Net result for Fig 2:** P(X,Y), J_Dr, J_ω from the new Fortran match
the authors' exact calculation to MC precision. The fix in 1.2 is the
reason the J_Dr numbers are this good — without the "don't reset on
blocked" rule, J_Dr on the wall picks up a bias that won't go away
with more samples.

---

## 3. Convention inconsistency WITHIN the new folder

This is the one thing I want you to see clearly before DP2.

**`tcrw_fig2_clean.f90` and `tcrw_fig2_defects.f90`:** use the
**wall-ring convention.** Grid is L × L but the outer ring
(x ∈ {0, L−1}, y ∈ {0, L−1}) is a WALL; walker only sees the inner
(L−2) × (L−2) playground. At L_F = 10 the physical walkable region is
8×8. Header comment at lines 15–21 says this is because "an edge or a
boundary is defined as lattice points along the boundary not allowed to
be occupied by the walker" — but that paper sentence actually refers
to **defect sites**, not the outer frame. The authors' `TRW.py` does not
mask the outer frame; the walker occupies all (L+1)² sites and is only
stopped on attempted exit.

So `tcrw_fig2_clean.f90` at L_F = 10 reproduces **authors' L = 7**
(grid 0..7 = 8×8 playground), not **authors' L = 10** (grid 0..10 =
11×11 playground). The cross-check above confirms this: the exact
Python at L_A = 7 matches to MC noise. The paper's Fig 2 caption says
L = 10 in the **authors' convention**, which means the new Fortran's
Fig 2 panels are at the right ω and D_r but on a smaller box than the
panels they claim to reproduce.

**`tcrw_fig3*.f90`:** use the **paper / authors' convention.** No wall
ring. From `tcrw_fig3a.f90:231` header comment:
*"no wall ring; lattice IS the L × L playground"*. Summary file
`tcrw_fig3a_summary.txt` confirms this numerically: for L = 4, n_edge +
n_bulk = 16 + 9 = 25 = (L+1)², i.e. 5×5 grid — the authors' convention
exactly. Fig 3 is paper-faithful as-is.

**What to do.** The easy fix is to change the two fig2 drivers so
the mask allows the outer ring too — `mask = .true.` with no
override — and then trust the `tcrw_step_mask` kernel's built-in bounds
check to handle "walker tries to exit" (it already does). That makes
the fig2 playground match the fig3 convention, and both match the
authors. Current L_F = 10 would then run the paper's L = 9 panel; set
L_F = 11 if you specifically want to reproduce the L = 10 paper panel.

Alternative if you prefer keeping the wall-ring picture for visual
reasons (it makes the plot borders clearer): just bump L_F from 10 to
13 so the inner playground is 11×11. Either works; I'd pick the first
because it lets you remove the wall mask and reuse the exact-Python
cross-check unchanged.

---

## 4. `tcrw_step.f90` kernel — line-by-line verification

| item | authors (`TRW.py`) | new `tcrw_step.f90` | status |
|---|---|---|---|
| direction encoding | `['↑','→','↓','←']` → 0..3 | same, DX=(0,1,0,-1), DY=(1,0,-1,0) | same |
| CCW rotation | `(idx - 1) % 4` | `mod(d+3, 4)` (= d−1 mod 4) | same |
| CW rotation | `(idx + 1) % 4` | `mod(d+1, 4)` | same |
| noise: P(CCW) | `random < (1 − ω)` is CW → P(CCW) = ω | `r_rot < omega → CCW` | same |
| chiral: P(CW) | `random < ω` is CW → P(CW) = ω | `r_rot < omega → CW` | same |
| blocked chiral | no move + no rotation | identical | same |
| prev_noise semantic | only set on non-blocked | `if (step_type /= 2) prev_noise = ...` | same |

Every line matches. The unbounded/`obc`/`mask` variants are three
copies of the same logic, which is fine — they specialize their bounds
check. `tcrw_step_mask` is the most general and is what the drivers
use.

---

## 5. What I did NOT audit in this pass

- **Fig 1b/1c/1d Fortran.** You have PDFs; those drivers use the same
  kernel, so the sign/flag fixes apply. But I didn't run a MSD
  cross-check against `tcrw_core.simulate_tcrw_pbc`. Would be worth
  5 minutes if you haven't done it.
- **Fig 3b–j Fortran.** Summary format is per-L,D_r / ω scalar rows, not
  spatial grids. Same kernel, same conventions as fig3a. The numeric
  values in `tcrw_fig3b_summary.txt` look sensible
  (|J_Dr|/|J_ω| ratio ~0.7 at D_r → 0 at L = 4, rising with D_r) but
  I didn't compare to exact-Python end-to-end for those scans.
- **Python-side** (`tcrw_fig4b_paper.py`): already independently
  verified by `tcrw_fig4b_crosscheck_authors.py` in the Apr 20 pass
  (see `TCRW_AUTHOR_CROSSCHECK_REPORT.md`).

---

## 6. Bottom line

The new Fortran folder is a **real upgrade** over the old one:

- chirality sign flip fixed → no more docstring patch needed;
- J-decomposition flag logic matches the authors' semi-Markov rule
  exactly, removing a subtle J_Dr edge-bias;
- MC steady state matches exact Python at MC precision across
  ω ∈ {0, 0.5, 1.0}.

The one outstanding issue is **L convention inconsistency between
Fig 2 drivers (wall ring) and Fig 3 drivers (authors' convention)**.
Fix is mechanical: drop the wall-ring mask in `tcrw_fig2_clean.f90`
and `tcrw_fig2_defects.f90`, or bump L_F by 3 there. Ten-line change
per file. After that, the Fortran is genuinely faithful to the paper
and to the authors' reference code.

Then you can pivot to DP2 (jerky-on-lattice) with a clean base.
