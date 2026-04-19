# TIFR DP1 → DP2 Research Assistant

Code and notes from my DP1 year at TIFR Hyderabad (PhD student under Prof. Kabir Ramola), and the beginning of my DP2 research direction on **jerky chiral active particles**.

The bulk of the code is a figure-by-figure reproduction of

> **Topological chiral random walker** — Saeed Osat, Ellen Meyberg, Jakob Metson, Thomas Speck, [arXiv:2602.12020](https://arxiv.org/abs/2602.12020) (2026)

which I am using as the lattice-model foundation for a DP2 project on adding jerk (3rd-order director dynamics) to chiral active walkers.

---

## Scope of the TCRW reproduction

- **Python implementation** (~12,500 lines across `tcrw_*.py`) covering every **data** panel in Figs 1–8 (main) and Figs 10–12 (extended). Fig 9 is a schematic and is not reproduced.
- **Accuracy suite:** 62/62 tests pass (column-stochasticity of the OBC transition matrix, λ = 1 at k = 0 for PBC, MC-vs-exact steady state agreement, div J = 0 on allowed bonds, the J = J_ω + J_{D_r} decomposition, the D(ω) = D(1−ω) symmetry, the gap at ω = 1/2, etc.).
- **Parameter fidelity:** ~70% of panels use the paper's exact (ω, D_r, L, T, N) parameters. ~30% (Fig 4f, 4g, 5c, 5d, 7e, 8d–l, 11, 12) reproduce the physics qualitatively at nearby parameters; these are flagged in `TCRW_CODE_INVENTORY.md` and will be re-run at paper parameters before any DP2 writeup.
- **Fortran + gnuplot** implementation of Fig 1(b,c,d) and all panels of Fig 2 at the paper-spec `T = 10^10`. Lives in `fortran_reproduction/`. Used as the quantitative-check harness for the Python code.

## Repository map

| path | what's there |
|------|--------------|
| `tcrw_core.py` | PBC/OBC Monte Carlo step loop. Direction encoding and CCW/CW convention. |
| `tcrw_obc.py` | Exact OBC steady state via sparse eigenvector of the transition matrix. |
| `tcrw_spectrum.py` | 4×4 Bloch matrix `P(k)`; PBC/OBC band structure. |
| `tcrw_currents.py` | OBC currents and the additive decomposition `J = J_ω + J_{D_r}`. |
| `tcrw_zak_phase.py` | 2D vectorized biorthogonal Wilson loop, returns `(Φ_x, Φ_y)` per band. |
| `tcrw_1d_edge.py` | 2×2 effective edge rate matrix `A(k)`; edge residence statistics. |
| `tcrw_geometry.py` | Masks for internal defects, holes, hybrid boundaries, custom polygons. |
| `tcrw_assembly.py`, `tcrw_maze.py` | Fig 5 (maze) and Fig 6 (self-assembly) applications. |
| `fortran_reproduction/` | Fortran driver (`tcrw_sim.f90`, `tcrw_fig2_defects.f90`, `mt.f90`) + gnuplot scripts for Fig 1 and Fig 2 at `T = 10^10`. |
| `TCRW_CODE_INVENTORY.md` | Authoritative panel-by-panel inventory with parameter-match status. Start here if you want to know what corresponds to what. |

## How to run one figure

```bash
# Python, Fig 2 (OBC steady-state currents and defect comparison):
python tcrw_currents.py
python tcrw_fig2_extra.py

# Fortran, Fig 2 defects at T = 10^10 (takes ~hours):
cd fortran_reproduction
bash run.sh tcrw-fig2-defects        # simulate
bash run.sh tcrw-fig2-defects-all    # then produce all three gnuplot figures
```

## Status

Work in progress — this is my active research directory, not a release. The panel-by-panel map and the 62-test suite are current as of 2026-04-17; anything else may be mid-edit. Issues and comments welcome, especially from people reproducing the same paper.

## Contact

Prashant Bisht — DP1 student, TIFR Hyderabad — prashantabisht@gmail.com
