# Fig 4 — fresh Python reproduction plan

**Date:** 2026-04-24
**Location:** `new_fortran_reproduction_and_python/`
**Language:** Python only (largest matrix is ~484 × 484 dense complex).
**Keep:** `tcrw_fig4b_paper.py` (already verified at 10⁻¹⁴ vs authors).
**Rewrite fresh:** every other Fig 4 panel, one self-contained file each.
**Convention:** authors' L — grid 0..L inclusive, (L+1) sites per axis.
So when the code says `L = 10` that is the paper's `L = 10` directly,
on an 11×11 playground.
**Cross-check:** wherever possible, diagonalise the authors'
`build_sparse_transition_matrix` (densified) and compare eigenvalue
sets to the new code's OBC matrix. Target: machine precision (≤10⁻¹²).

---

## 1. Panel roadmap

Paper Fig 4(a) is a schematic → skip. Working files to create:

| # | file | panel | matrix size (paper L) | new math vs 4(b) |
|---|------|-------|------------------------|-------------------|
| 0 | `tcrw_fig4b_paper.py` | 4(b) | 4×4 × N_k² | — (already done) |
| 1 | `tcrw_fig4c.py` | 4(c) | 4(L+1)² = 36 at L=2 | build OBC matrix; eigenvalue scatter; (ω, D_r) scan |
| 2 | `tcrw_fig4d.py` | 4(d) | 4(L+1)² = 484 at L=10 | OBC eig-scan vs D_r at ω=1 |
| 3 | `tcrw_fig4e.py` | 4(e) | 484 | OBC eig-scan vs ω at D_r=0.1 |
| 4 | `tcrw_fig4f.py` | 4(f) | 484 | OBC complex plane + edge-weight coloring, 3 D_r |
| 5 | `tcrw_fig4g.py` | 4(g) | 484 | same but 3 ω |
| 6 | `tcrw_fig4i.py` | 4(i) | 4×4 × N_k² | PBC band circle, 3 coloring schemes |
| 7 | `tcrw_fig4h.py` | 4(h) | 4(L+1) = 44 at L=10, × N_ky | HPBC matrix + k_y scan |

Order chosen to climb difficulty smoothly. Items 4-6 reuse the
edge-weight machinery introduced in item 4. Item 7 (HPBC) is kept
last because it has no author cross-check.

Each file is self-contained: it carries its own `build_Pk`,
`build_obc_matrix`, or `build_hpbc_matrix`. The `tcrw_fig4b_paper.py`
`build_Pk` is reproduced verbatim in items 1 and 6 — deliberately,
per the "one file per panel, self-contained" choice.

---

## 2. Universal file skeleton

Every new `tcrw_fig4X.py` must have these four sections in this
order:

```
"""
TCRW Fig 4(X) — <one-line description>
=======================================
Paper: Osat, Meyberg, Metson, Speck, arXiv:2602.12020.

Physics           : ...
Convention        : authors' L (grid 0..L inclusive, (L+1)^2 sites).
Paper parameters  : ω = ..., D_r = ..., L = ..., k grid = ...
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp


# ==== 1. Matrix builder (self-contained) ====

def build_XXX_matrix(omega, D_r, L, ...):
    """4(L+1)² × 4(L+1)² column-stochastic OBC matrix."""
    ...


# ==== 2. Spectrum (and edge weights if needed) ====

def obc_spectrum(omega, D_r, L, return_eigvecs=False):
    ...


# ==== 3. Author cross-check ====

def crosscheck_authors(omega, D_r, L, tol=1e-12):
    """Compare eigenvalue sets: ours vs authors' TRW.py densified."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "TRW", "TRW._original_code_by_paperauthors.py")
    TRW = importlib.util.module_from_spec(spec); spec.loader.exec_module(TRW)
    W_A = TRW.build_sparse_transition_matrix(L, omega, D_r).toarray()
    vals_A = np.linalg.eigvals(W_A)
    vals_U = np.linalg.eigvals(build_XXX_matrix(omega, D_r, L))
    d = _hausdorff(vals_A, vals_U)
    assert d < tol, f"Author crosscheck failed: {d:.2e}"
    print(f"  [ok] author crosscheck (ω={omega}, D_r={D_r}, L={L}): {d:.2e}")


# ==== 4. Figure generation ====

def make_figure():
    """Reproduces Fig 4(X) at paper params."""
    ...
    fig.savefig("tcrw_fig4X_paper.png", dpi=220, bbox_inches="tight")


if __name__ == "__main__":
    crosscheck_authors(...)   # quick sanity, tiny L
    make_figure()             # full paper-param run
```

Non-negotiable: **every file must run `crosscheck_authors` at at least
one small-L point before producing the figure**, so a regression is
caught immediately. For 4(h) HPBC where authors have no counterpart,
substitute the self-consistency tests in §4.

---

## 3. Paper parameters (ground truth — read off the PDF, 2026-04-24)

Source: tcrw.pdf page 5 (Fig 4 itself) + caption + main-text excerpts.
All L values are in authors' convention (grid 0..L inclusive).

| panel | ω | D_r | L | axis / plot type |
|---|---|---|---|---|
| 4(b) | {0.35, 0.5, 0.65} | **0.1** | — (PBC Bloch) | 3D surface Re(λ) over (k_x, k_y), 80×80 BZ grid |
| 4(c) | [0, 1] (scan) | [0, 1] (scan) | **2** (3×3 playground) | 3D surface Re(λ) and Im(λ) over (D_r, ω); colour = edge weight |
| 4(d) | **1** | scan [small → ~0.5] | **10** | Re(λ) vs D_r, all 4(L+1)²=484 eigenvalues |
| 4(e) | scan [0, 1] | **10⁻¹ = 0.1** (panel label) | **10** | Re(λ) vs ω, all 484 eigenvalues |
| 4(f) | **1** | **{0.65, 0.5, 0.35}** (3 columns) | **10** | complex plane Re vs Im; **two coloring rows**: top=edge weight, bottom=CCW weight |
| 4(g) | **{0.5, 0.7, 0.9}** (3 columns) | **fixed — value NOT labelled on panel**; likely 0.1 (check SI) | **10** | same layout as (f) |
| 4(h) | **{0.5, 0.7, 1.0}** (3 columns) | **fixed — NOT labelled on panel**; likely 0.1 | **10** | Re(λ) (top row) + Im(λ) (bottom row) vs k_y; colour = edge weight |
| 4(i) | **{0.5, 0.7, 1.0}** (3 columns) | **fixed — NOT labelled on panel**; likely 0.1 | **10** | Scatter over (cos k_x, cos k_y) plane; **three coloring rows**: top=bulk, middle=CCW, bottom=to-edge |

### Ambiguities to resolve BEFORE writing the code

1. **4(g) and 4(h) D_r value.** Panel columns show only ω labels; D_r
   value is not on the panels themselves. Main text p.5 says
   "illustrated for a system of size L = 10 for fixed of ω and Dr, see
   Fig. 4(f)–(g)" — it says "fixed" but gives no number. Best guess:
   D_r = 0.1, matching 4(e) and 4(b). **Action:** search the SI /
   Methods for explicit number; if absent, default to 0.1 and flag in
   the docstring.

2. **4(i) "HPBC" in caption.** Caption reads:
   *"Band structure of the HPBC on a closed circle defined by
   (cos k_x, cos k_y)."*
   But HPBC (x open, y periodic) has only k_y as a good quantum
   number — there is no k_x. So the axes (cos k_x, cos k_y) require
   **two** independent momenta, which only full PBC provides. Most
   likely interpretation: **caption typo, 4(i) is full PBC**, matching
   your existing `fig4i_band_circle.png` convention. **Action:**
   proceed as PBC for 4(i) but leave a `# PAPER_CAPTION_AMBIGUITY`
   note in the docstring; ask the authors if a second pass is warranted.

3. **4(f) and 4(g) have two coloring rows** (edge, CCW) per
   parameter column. My earlier draft said three — the three-coloring
   scheme is panel 4(i), not 4(f/g).

### Caption text (verbatim, for reference)

> FIG. 4. Topological origin of the edge localization and edge current.
> (a) Labeling of states for ω = 1 where CCW edge current is expected.
> Edge states are union of "to edge" and "along edge" states.
> (b) Spectrum of the model in PBC. The gap closes at ω = 0.5, where a
> topological transition occurs.
> (c) Decomposition into real and imaginary part of the spectrum in
> the presence of OBC for a small lattice of L = 2. Note that the
> spectrum is colored based on the localization of corresponding
> eigenvectors on the edge states.
> (d) Real part of spectrum for a fully chiral case ω = 1. Note the
> coalescence of eigenvalues at Dr → 0.
> (e) Same as (d) but for fixed Dr.
> (f) Spectrum of the model on a complex plane for a fully chiral case.
> Coloring differentiates localization on edge and current along edge.
> (g) Same as (f) but for changing ω.
> (h) Spectrum of a model with hybrid boundary conditions (periodic
> along the y-direction and open along the x-direction). Note the
> bands localized on the edge.
> (i) Band structure of the HPBC on a closed circle defined by
> (cos k_x, cos k_y). Different coloring is used to show localization
> on different states.

---

## 4. HPBC self-consistency (for 4(h) where author cross-check is impossible)

`build_hpbc_matrix(ω, D_r, L, k_y)` is a 4(L+1) × 4(L+1) complex
matrix. To trust it:

1. **At k_y = 0 and k_y = π**, the matrix must be real (because all
   e^{±i k_y} phases become ±1). Check `np.max(|Im(H)|) < 1e-14`.
2. **Column-stochastic property** fails at k_y ≠ 0 (phases break
   column sums), but **at k_y = 0 it should hold exactly**.
3. **Large-L convergence**: at large L, the HPBC spectrum should
   densely fill the PBC spectrum integrated over k_x. Concretely,
   build the 2D PBC spectrum of `build_Pk(ω, D_r, k_x, k_y=k_y*)` at
   the same k_y*, sample over 200 k_x values — the resulting 4 × 200
   eigenvalues should overlap (set-wise) with the HPBC spectrum at
   the same k_y* up to edge-mode residuals.
4. **Edge modes**: at k_y = 0.5 (say) with ω = 1, D_r = 0.1, expect
   a small number of eigenvalues that are strongly localised on x = 0
   or x = L — isolated from the bulk cloud. Verify by inspecting the
   edge-weight histogram.

Run all four checks inside `tcrw_fig4h.py::crosscheck_self()`.

---

## 5. Implementation order with ETAs

Sequential, each item gated by the author cross-check passing:

| step | ETA | notes |
|------|-----|-------|
| 1. `tcrw_fig4c.py` (L=2 OBC) | 45 min | Smallest matrix first — verifies the OBC builder. Author cross-check at L=2 is instant. |
| 2. `tcrw_fig4d.py` | 30 min | Reuses 4(c) builder. D_r scan. Author cross-check at a few D_r points. |
| 3. `tcrw_fig4e.py` | 30 min | Reuses 4(c) builder. ω scan. |
| 4. `tcrw_fig4f.py` | 45 min | Adds edge-weight coloring. Verify weights sum to 1 per eigenvector. |
| 5. `tcrw_fig4g.py` | 30 min | Near-clone of 4(f) with different parameter axis. |
| 6. `tcrw_fig4i.py` | 45 min | Uses build_Pk only. The "three colorings" is a plotting task, not new matrix math. |
| 7. `tcrw_fig4h.py` | 90 min | HPBC — the only item with genuine risk. Allow the full hour for the self-consistency suite. |

**Total: ~5.5 hours** of focused work if nothing blows up. Worst-case
with an HPBC snag, budget a full day.

---

## 6. What "done" looks like

When the plan is complete:

1. Seven new files in `new_fortran_reproduction_and_python/` (plus
   the retained `tcrw_fig4b_paper.py`).
2. Each file runs standalone: `python tcrw_fig4X.py` prints
   `[ok] author crosscheck ...` or the self-consistency report, then
   writes `tcrw_fig4X_paper.png`.
3. A `tcrw_fig4_regression.py` harness that imports all seven and
   runs their cross-checks back-to-back. Takes <1 minute. This is the
   CI-style check that guards against later refactors.
4. Updated `TCRW_CODE_INVENTORY.md` §2.4 (Fig 4 row) with params =
   paper-faithful status ✅ across the board.

---

## 7. First concrete action

Create `new_fortran_reproduction_and_python/tcrw_fig4c.py`. It's the
smallest matrix (36×36 at L=2), so the author cross-check takes
milliseconds and forces you to get the matrix-build convention right
before the bigger panels.

Skeleton for `tcrw_fig4c.py`:

```python
"""
TCRW Fig 4(c) — OBC spectrum on L=2 lattice (complex plane).

Paper: arXiv:2602.12020 Fig 4(c).
Convention: authors' L — grid 0..L inclusive, (L+1)^2 sites.
Paper parameters: L = 2 (3x3 playground, 36 states).
                  Scan (ω, D_r) on [0,1]^2, say 15 x 15 points.
Output: Re(λ) vs Im(λ) scatter, coloured by ω (or D_r).
"""
```

Once 4(c) is in and the author cross-check prints `[ok]`, the next
six panels are mechanical extensions.
