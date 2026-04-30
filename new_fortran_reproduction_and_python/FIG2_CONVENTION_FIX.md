# Fig 2 Fortran convention fix — paper-faithful 10×10 playground

**Date:** 2026-04-30

## What was wrong

The earlier Fortran drivers `tcrw_fig2_clean.f90` and
`tcrw_fig2_defects.f90` used a **wall-ring** boundary convention:

```fortran
mask = .false.
do iy = 1, L - 2
   do ix = 1, L - 2
      mask(ix, iy) = .true.    ! only the inner (L-2)x(L-2) is allowed
   end do
end do
```

With `L = 10` this gave the walker an **8×8 playground** (sites
(1..8) × (1..8)), not the paper's 10×10 (sites (0..9) × (0..9)).

Earlier audit (`NEW_FORTRAN_AUDIT.md` §3) flagged this:
*"Fig 2 Fortran at L_F = 10 reproduces authors' L = 7 (8×8 playground),
not paper's L = 10."*

## What was changed

Replaced the wall-ring with `mask = .true.` everywhere.  The walker
now lives on the full L×L grid, and OBC behaviour comes from the
kernel's bounds check inside `tcrw_step_mask`:

```fortran
if (nx >= 0 .and. nx < Lx .and. ny >= 0 .and. ny < Ly) then
    if (mask(nx, ny)) then
        ! valid move: translate + rotate
    else
        step_type = 2          ! defect or out-of-bounds
    end if
else
    step_type = 2              ! out-of-bounds
end if
```

This is exactly what authors' `TRW.build_sparse_transition_matrix`
does in `TRW.py` (line `if not (0 <= new_i <= L and 0 <= new_j <= L)`).

### Specific edits

**`tcrw_fig2_clean.f90`:**
- `mask = .false.` + inner loop  →  `mask = .true.`
- Header comment block updated (lines 13-25) to describe the new
  no-wall-ring convention and reference authors' code.
- Sanity-print loop changed from `1..L-2` to `0..L-1` for edge / bulk
  averages.
- Initial position `x = y = L/2 = 5` unchanged — it's still the centre
  of the (now larger) playground.
- Output file comments updated: "every site is allowed" instead of
  "(inner) sites".

**`tcrw_fig2_defects.f90`:**
- Same `mask = .false. → mask = .true.` change.
- **Defect pattern updated** to match paper Fig 2(k) / authors'
  notebook L-shape:
  ```
  defects = (3,3), (3,4), (3,5), (3,6), (4,3), (5,3), (6,3)
  ```
  (7 cells; previous version had a 5-cell plus-sign at (4, 5) which
  was a custom test pattern, not the paper's exact one.)
- Header comment block updated with new ASCII-art diagram of the
  L-shape defect.
- Sanity-print loop ranges changed to `0..L-1`.
- Initial position `(2, 2)` unchanged — still in the bulk, away from
  defects.

## Result

After the fix, with `L = 10`:

| convention | grid | walkable | matches paper L=10? |
|---|---|---|---|
| **before (wall ring)** | 10×10 | 8×8 inner | no (1.25× too small) |
| **after (no ring)** | 10×10 | full 10×10 | **yes** |

This now matches the paper's Fig 2 axes (0..9) and matches the
authors' `TRW.build_sparse_transition_matrix(L=9, ω, D_r)`.

### Convention naming caveat

There are TWO L conventions in the literature for this model:

| name | meaning | example |
|---|---|---|
| **Paper L** | linear site count | "L = 10" → 10 sites, axes 0..9 |
| **Authors' code L** | max index | `L = 9` in TRW.py → grid 0..9 = 10 sites |

The relation is `paper L = authors L + 1`.  My Python files
(`tcrw_fig4c.py` etc.) use **authors' L**, so they take `L = 9` to
reproduce paper L = 10.  The Fortran files now use the **paper L**
naming: `L = 10` directly means 10 sites, which is what readers
expect.

Cross-checking these two conventions is straightforward — pass `L = 9`
to the Python and `L = 10` to the Fortran for the same physical
problem.

## How to verify

The repo doesn't have gfortran in the sandbox, so I couldn't compile
in-place.  On your Mac:

```bash
cd new_fortran_reproduction_and_python
gfortran -O2 -fno-range-check -ffree-line-length-none \
         tcrw_fig2_clean.f90 -o tcrw_fig2_clean
gfortran -O2 -fno-range-check -ffree-line-length-none \
         tcrw_fig2_defects.f90 -o tcrw_fig2_defects
```

Then run and verify against the exact Python:

```bash
./tcrw_fig2_clean       # ~15 minutes, T = 10^10
./tcrw_fig2_defects     # ~5 minutes, T = 10^10
python tcrw_fortran_vs_exact.py     # MC vs exact overlay
```

The cross-check should now pass at MC noise levels (~10⁻⁵ on P at
T = 10¹⁰), without the systematic L=7-vs-L=9 bias that the wall-ring
convention introduced.

## What this enables

- **Direct comparison** to paper Fig 2 panels at the same axes.
- **Direct cross-check** between `tcrw_fig2_pymc.py` (L = 9, authors'
  convention, exact) and the new `tcrw_fig2_clean.f90` /
  `tcrw_fig2_defects.f90` (L = 10, paper convention, MC).  Both now
  describe a 10×10 walkable playground.
- **Paper-faithful reproduction** of Fig 2 panels (a-e), (f-j),
  (k-o) directly from Fortran output.

## Files modified

- `tcrw_fig2_clean.f90`  — wall-ring removed
- `tcrw_fig2_defects.f90` — wall-ring removed + L-shape defect pattern

No changes needed to `tcrw_step.f90` (the kernel was always correct;
it was the driver-level mask that was wrong).

No changes needed to `tcrw_fortran_vs_exact.py` (the cross-check
already adjusted for the smaller-playground convention; it'll
auto-adapt to the new larger playground because it parses L from the
output files).
