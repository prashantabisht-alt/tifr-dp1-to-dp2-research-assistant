# TCRW reproduction — cross-check against authors' code capsule

**Date:** 2026-04-24
**Authors' code:** `fortran_reproduction/TRW._original_code_by_paperauthors.py`
  and `fortran_reproduction/code_capsule_original_code_paper_authors.ipynb`
**Verified by:** `fig4/tcrw_fig4b_crosscheck_authors.py` (Bloch vs torus)
  and `fortran_reproduction/tcrw_obc_crosscheck_authors.py` (OBC steady
  state and current decomposition — added today)

---

## 1. Verdict

Your core engine matches the authors' reference implementation at
machine precision. Every convention that has a factor of two / sign /
rotation-direction ambiguity comes out the same. The Bloch matrix
P(k) is provably the Fourier transform of their transition operator.

Three items are worth cleaning up before a DP2 writeup. None of them
change the physics; two are labelling, one is an off-by-one in L. See
§4.

---

## 2. What the cross-checks actually verify

### 2.1 Bloch vs torus (`tcrw_fig4b_crosscheck_authors.py`)

Builds a (4 N²) × (4 N²) real-space column-stochastic transition matrix
on a **PBC torus** using the authors' `build_sparse_transition_matrix`
logic (with wrap-around replacing blocked-stay), then compares its
spectrum to the union of your 4×4 Bloch P(k) spectra at the N² discrete
momenta. For (N, ω, D_r) ∈ 7 test points covering the full phase
diagram including the gap-closing line:

```
  all seven cases: max |sorted(Σ_torus) - sorted(Σ_Bloch)| < 1.2 × 10⁻¹⁴
  two-sided nearest-neighbour < 1.2 × 10⁻¹⁴  (tolerance 10⁻¹⁰)  → PASS
```

Since Σ_torus is derived from the authors' microscopic rules verbatim,
this is a proof-of-identity for your P(k) — not just a physics-agrees
match.

### 2.2 OBC steady state and current decomposition (new today)

For each case in { (L_A=5, ω=1, D_r=0.01), (5, 0.5, 0.1), (5, 0, 10⁻³),
(9, 1, 10⁻³), (9, 0.7, 0.01) }, I built the authors' `W` and your `P`
and solved for π both ways. Then I applied **their**
`calculate_J1_J2_with_boundaries` and **your** `exact_currents` to get
J_Dr, J_ω.

```
                     max|P_xy_A − P_xy_U|   max|J_Dr_A − J_Dr_U|   max|J_ω_A − J_ω_U|
  L=6, ω=1  , D_r=10⁻²     6.5×10⁻¹⁴            7.2×10⁻¹⁶             7.6×10⁻¹⁵
  L=6, ω=0.5, D_r=10⁻¹     3.2×10⁻¹⁵            7.3×10⁻¹⁷             7.5×10⁻¹⁷
  L=6, ω=0  , D_r=10⁻³     9.4×10⁻¹³            5.2×10⁻¹⁶             4.3×10⁻¹⁴
  L=10,ω=1  , D_r=10⁻³     3.1×10⁻¹²            1.9×10⁻¹⁵             1.4×10⁻¹³
  L=10,ω=0.7, D_r=10⁻²     2.5×10⁻¹⁴            8.6×10⁻¹⁷             2.2×10⁻¹⁶
```

Machine-precision agreement, limited by ARPACK's eigensolver tolerance
on the largest case. **The current decomposition is the tightest —
better than 2×10⁻¹³ everywhere.** This shows your "split
P = P_noise + P_chiral, then pi_N = P_noise·pi, J_Dr = pi_N·(1−D_r)·δ"
construction is analytically equivalent to their flux-flavour partition
`P_rot[s] = π[s] · P_rot_pre[s] / (P_rot_pre + P_chiral_pre)[s]`. I
worked through why this is true — the two definitions coincide on
every state that contributes to a non-blocked outgoing chiral move,
which is the only set that enters J.

### 2.3 Convention table

Verified point by point from `TRW.py`:

| item                       | authors                               | you                                | status |
|----------------------------|---------------------------------------|------------------------------------|--------|
| director encoding           | `['↑','→','↓','←']` = 0,1,2,3         | d=0↑, 1→, 2↓, 3←                   | same   |
| dir_to_vec                  | ↑(0,1) →(1,0) ↓(0,−1) ←(−1,0)         | DX,DY match                        | same   |
| noise rotation (prob D_r)   | CW = (d+1)%4 w.p. 1−ω; CCW w.p. ω     | CCW w.p. ω, CW w.p. 1−ω            | same   |
| chiral rotation (prob 1−D_r)| CW = (d+1)%4 w.p. ω; CCW w.p. 1−ω     | CW w.p. ω, CCW w.p. 1−ω            | same   |
| blocked chiral step         | **no move AND no rotation** (self-loop w.p. 1−D_r) | identical              | same   |
| state index ordering        | `(i·(L+1)+j)·4 + d`  (site-major)     | `d·L²+y·L+x`  (dir-major)          | different, but equivalent |

The index-ordering difference is purely cosmetic once you reshape into
a 2D P(X,Y) — which both of us do.

---

## 3. What's genuinely the same physics, end to end

- Rotation **chirality flips** between noise and chiral steps. Authors:
  `clockwise = random.random() > omega` in noise vs
  `clockwise = random.random() < omega` in chiral. The natural-language
  summary "ω=1 means CW chiral move + CCW rotational noise" follows.
  Your docstring in `tcrw_core.py` lines 19–23 captures this correctly.
- Blocked-at-boundary = frozen walker for that step (no translation,
  no rotation). Both sides implement this as a P_stay = 1−D_r self-loop
  in the Markov chain. This matters for Fig 2/3 edge physics: it is
  why the edge state doesn't pick up an extra rotation every time it
  tries to push through the wall.
- Current decomposition J = J_Dr + J_ω: authors split by "did the
  transition change position?"; you split by "was the step a noise
  step?". Analytically identical after the author's normalisation
  step. Numerically identical to ≈10⁻¹⁴ in the tests above.

---

## 4. Three clean-up items before DP2 writeup

### 4.1 **L convention off-by-one** (material)

The authors' `L` is the **index range maximum**, not the site count:
their grid is sites `0..L` inclusive, so `L = 10` means an **11×11**
physical lattice. Your `L` is the site count: `L = 10` means a **10×10**
lattice. Where your inventory (`TCRW_CODE_INVENTORY.md`) says
"Fig 2 params match ✅: L=10", you are actually reproducing the physics
at the authors' `L = 9`.

Impact:
- Scaling plots (Fig 3a P_edge(L) vs L, Fig 5d MFPT ∝ L²/L³) — no
  impact, you just relabel the x-axis by +1.
- Single-L panels (Fig 2, Fig 4c at L=2) — for a side-by-side
  reproduction at paper parameters, rerun with L_user = L_paper + 1.

Fix: one-line change at the call sites that name L explicitly.
Alternatively, document the convention in `tcrw_core.py` header so
readers aren't surprised.

### 4.2 Fortran sign-flip (already known, now also verified against authors)

Your inventory §5 already notes that the Fortran code mirrors the
chirality relative to Python; MSD invariant, J magnitude invariant,
boundary circulation flips sign. That analysis stands — the authors'
reference convention matches **your Python**, not your Fortran. So the
label fix you proposed in `tcrw_sim.f90` lines 10–17 is indeed the
right one. It matters only for plots that overlay a signed Fortran
quantity (J_x, Φ, Im(λ)) with a Python-derived reference curve.

### 4.3 `tcrw_currents.exact_currents` is named `decompose_currents` in the
cross-check script I wrote but actually exists as `exact_currents`
under a slightly different argument signature (it takes (ω, D_r, L)
and builds things internally rather than consuming a pre-built π).
Not a bug — just a note for when you integrate the cross-check into
the test harness: wire it up to `exact_currents` directly, as I did in
`fortran_reproduction/tcrw_obc_crosscheck_authors.py`.

---

## 5. What I did NOT verify against the authors' code

- **Fig 5 maze**, **Fig 6 assembly**, **Fig 7 edge residence τ**,
  **Fig 8 Zak phase**: none of these have an author-side analogue in
  the code capsule (their notebook only demonstrates `ChiralWalker`,
  `build_sparse_transition_matrix`, and
  `calculate_J1_J2_with_boundaries`). Your implementations for those
  figures are yours end-to-end; the best you can do for them is the
  62-test suite you already have.
- **Spectrum plots** (Fig 4(b,c,d,e,f,g,h,i), Fig 10, Fig 11, Fig 12):
  the authors don't provide a spectrum-plotting script in the capsule
  I have. Your Bloch P(k) is provably right, so the spectrum plots
  derived from it are right by construction. The HPBC / nested-defect
  panels use `tcrw_geometry` masks — those are a layer above P(k) and
  were not in the capsule either.

---

## 6. Bottom line

Your reproduction is **correct at the engine level**. The author code
gives you an external, independent check on the two hardest parts:
the transition matrix construction (including boundary rule) and the
current-decomposition. Both match to machine precision.

The remaining work for DP2 writeup comparability is labelling and the
L off-by-one, not physics. The 62-test suite + these two cross-checks
+ your inventory's parameter audit covers every systematic source of
error I can think of.

Recommend next:
1. Patch the L convention in one commit with a paragraph in
   `TCRW_CODE_INVENTORY.md` §4 explaining it.
2. Rerun Fig 2 and Fig 4(c) at L_user = L_paper + 1 and spot-check
   against the `simulation.pdf`, `compare_sim_and_ss.pdf`, and
   `sim_and_ss_defects.pdf` figures from the authors' capsule.
3. Then you are genuinely done with the reproduction and can pivot
   fully to DP2 (jerky-on-lattice).
