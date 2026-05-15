# Triangular JMVR Active Random Walker

This folder contains the cleaned triangular-lattice JMVR work after the PI
meeting, the triangular-torus k-grid check, and the Dipanjan \(c_3\) sign-error
diagnosis.

## Scientific Status

The active model is the triangular JMVR continuous-time random walker:

- six internal director states,
- hopping to all six nearest neighbours,
- forward/backward translation bias \(1/6 \pm \epsilon\),
- director switching \(d \to d\pm1\) at rate \(\gamma/2\).

There are two required corrections:

1. Use the sheared triangular-torus reciprocal grid,

\[
k_1=\frac{\pi m_1}{aL},
\qquad
k_2=\frac{\pi(2m_2-m_1)}{bL}.
\]

2. Use the corrected \(c_3\) chirality sign,

\[
c_3 = B - 1 - \gamma - 2 i \epsilon \sin(a k_1 - b k_2).
\]

Dipanjan's notebook used a rectangular Fourier grid and the opposite sign for
this \(c_3\) term. Fixing both brings the exact theory to the \(10^8\)-walker
Fortran KMC noise floor.

## Run These First

From this folder:

```bash
python3 triangular_jmvr_corrected.py
```

This checks:

- probability conservation at \(k=0\),
- the zero eigenvalue at \(k=0\),
- sixfold spectral symmetry for the corrected matrix,
- failure of sixfold symmetry for the old Dipanjan sign.

Then run the two-bug forensic comparison:

```bash
python3 forensic_two_bugs.py
```

Only the sheared-grid plus corrected-\(c_3\) case reaches the KMC noise floor.

For the direct finite-torus Bloch check:

```bash
python3 verify_realspace_bloch.py
```

This builds the full real-space generator on small triangular tori and checks
that its spectrum matches the sheared-grid corrected Bloch spectra to machine
precision.

## Final Fig. 11 Reproduction

Python version:

```bash
python3 fig11_final_hex.py
```

Outputs:

- `fig11_final_hex.png`
- `fig11_final_hex.pdf`
- `fig11_final_hex_evidence.png`
- `fig11_final_hex_evidence.pdf`

Gnuplot version:

```bash
python3 export_fig11_final_hex_gnuplot_data.py
gnuplot fig11_final_hex.gnu
```

Outputs:

- `fig11_final_hex_gnuplot.png`
- `fig11_final_hex_gnuplot.pdf`

Original-draft-style visual reproduction:

```bash
gfortran -O3 -fno-range-check -ffree-line-length-none kmc_triangular_jmvr_L60.f90 -o kmc_triangular_L60
./kmc_triangular_L60 > kmc_L60_run.log
python3 export_fig11_original_style_gnuplot_data.py
gnuplot fig11_original_style.gnu
```

Outputs:

- `fig11_original_style_gnuplot.png`
- `fig11_original_style_gnuplot.pdf`

This uses \(L=60\), not the \(L=30\) verification torus, so the active front
has not visibly wrapped around the periodic boundary. That is why it reproduces
the annular/hexagonal look of the original draft Fig. 11.

## Monte Carlo

The Monte Carlo is Fortran, not Python.

```bash
./run_kmc.sh
```

This builds and runs:

- `kmc_triangular_jmvr.f90`
- `mt.f90`

and writes:

- `kmc_triangular_counts.txt`
- `kmc_run.log`

The saved count file currently contains \(N=100,000,000\) walkers. Python only
reads this count file, normalizes counts by \(N\), and overlays it with the
exact matrix result.

## Active File Map

| File | Role |
|---|---|
| `triangular_jmvr_corrected.py` | Canonical corrected 6x6 Bloch generator and core self-checks |
| `forensic_two_bugs.py` | Four-way comparison showing both k-grid and \(c_3\) corrections are required |
| `verify_realspace_bloch.py` | Finite real-space generator versus Bloch spectrum check |
| `triangular_active_walker.py` | Lightweight sanity/band-line starter using the corrected matrix |
| `kmc_triangular_jmvr.f90` | Fortran KMC implementation of the real-space master equation |
| `kmc_triangular_jmvr_L60.f90` | L60 Fortran KMC used only for original-draft-style Fig. 11 visuals |
| `mt.f90` | Random-number generator used by the Fortran KMC |
| `run_kmc.sh` | Build/run helper for KMC |
| `kmc_triangular_counts.txt` | Saved \(10^8\)-walker KMC counts |
| `kmc_triangular_counts_L60.txt` | Saved \(2\times10^7\)-walker L60 KMC counts for original-style plotting |
| `fig11_final_hex.py` | Final Python Fig. 11 replacement |
| `export_fig11_final_hex_gnuplot_data.py` | Exports final `.txt` tables for gnuplot |
| `fig11_final_hex.gnu` | Final full gnuplot Fig. 11 |
| `export_fig11_original_style_gnuplot_data.py` | Exports L60 `.txt` tables for original-draft-style gnuplot |
| `fig11_original_style.gnu` | Gnuplot script matching the old Fig. 11 annular/hexagonal visual convention |
| `PI_NOTE_TWO_BUGS_AND_CHECKS.md` | PI-facing explanation of the k-grid/PBC and \(c_3\) sign corrections |
| `TRIANGULAR_LATTICE_PLAN.md` | Research plan and remaining tasks |

## Folders

| Folder | Role |
|---|---|
| `outputs/` | Final text data used by `fig11_final_hex.gnu` |
| `papers/` | Reference papers sent or used for this project |
| `gnuplot_panels/` | Optional one-panel gnuplot scripts |
| `legacy_debug/` | Archived debugging, forensic scripts, and old figure attempts |

## What Not To Use For New Work

Do not build future physics from files inside `legacy_debug/`. They are kept as
evidence of the bug hunt and old plotting attempts. Future calculations should
start from `triangular_jmvr_corrected.py`.

## Next Physics Step

The next scientific task is not topology yet. It is:

1. corrected \(\Gamma-M-K-\Gamma\) band structure using
   `triangular_jmvr_corrected.py`,
2. gap scan,
3. real-space/Bloch consistency check,
4. then later multi-particle or triangular TCRW/topology.
