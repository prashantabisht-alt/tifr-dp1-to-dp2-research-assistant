# Meeting Prep: Triangular JMVR Result

## One-Sentence Story

I checked Dipanjan's triangular active-walker calculation against the real-space
master equation and high-statistics Fortran KMC. The mismatch needs **two**
corrections: the triangular-torus Fourier grid must be sheared, and the
\(c_3\) chirality term has the wrong sign in the old matrix. Fixing both brings
the exact theory to the \(10^8\)-walker KMC noise floor.

## What Is Done

1. Built the triangular JMVR \(6\times6\) continuous-time generator.
2. Audited Dipanjan's Mathematica notebook `rtp_tl_2.nb`.
3. Checked the triangular-torus reciprocal grid.
4. Checked the Appendix-B real-space master equations term by term.
5. Found the \(c_3\) sign error.
6. Wrote independent Fortran KMC directly from the real-space rules.
7. Compared the four combinations:
   - rectangular notebook grid + old \(c_3\),
   - rectangular notebook grid + corrected \(c_3\),
   - sheared triangular grid + old \(c_3\),
   - sheared triangular grid + corrected \(c_3\).
8. Made final Python and gnuplot Fig. 11 replacements.

## The Model

There are six internal states \(d=0,\ldots,5\), one for each triangular-lattice
nearest-neighbour direction.

For a walker in director state \(d\):

\[
\text{rate to hop forward along }d = \frac16+\epsilon,
\]

\[
\text{rate to hop backward along }d = \frac16-\epsilon,
\]

and the other four nearest-neighbour hops have rate

\[
\frac16.
\]

The director switches symmetrically:

\[
d\to d+1 \quad \text{at rate } \frac{\gamma}{2},
\]

\[
d\to d-1 \quad \text{at rate } \frac{\gamma}{2}.
\]

So the total outgoing rate is

\[
1+\gamma.
\]

The constraint is

\[
|\epsilon|\leq \frac16,
\]

so that all hopping rates are non-negative.

## Correct Triangular-Torus k-Grid

The triangular real-space periods are

\[
T_1 = L\mathbf{a}_1 = (2aL,0),
\]

\[
T_2 = L\mathbf{a}_2 = (aL,bL).
\]

Allowed Fourier modes must satisfy

\[
\mathbf{k}\cdot T_1 = 2\pi m_1,
\qquad
\mathbf{k}\cdot T_2 = 2\pi m_2.
\]

Therefore

\[
k_1 = \frac{\pi m_1}{aL},
\]

\[
k_2 = \frac{\pi(2m_2-m_1)}{bL}.
\]

This is the sheared reciprocal grid used in `triangular_jmvr_corrected.py`.

Dipanjan's notebook used a rectangular grid,

\[
k_1 = \frac{2\pi n_x}{2aL},
\qquad
k_2 = \frac{2\pi n_y}{bL}.
\]

For the second triangular period,

\[
\mathbf{k}\cdot T_2 = \pi n_x + 2\pi n_y.
\]

For odd \(n_x\), this is not a multiple of \(2\pi\). So those modes are not
periodic on the triangular torus. This confirms Kabir's PBC suspicion.

## Correct Matrix Statement

The common bulk part is

\[
B(k_1,k_2)
=
\frac13\left[
\cos(2ak_1)+\cos(ak_1+bk_2)+\cos(ak_1-bk_2)
\right].
\]

The corrected diagonal terms are

\[
c_1 = B-1-\gamma+2i\epsilon\sin(2ak_1),
\]

\[
c_2 = B-1-\gamma+2i\epsilon\sin(ak_1+bk_2),
\]

\[
c_3 = B-1-\gamma-2i\epsilon\sin(ak_1-bk_2).
\]

Dipanjan's notebook used

\[
c_3^{\rm old}
=
B-1-\gamma+2i\epsilon\sin(ak_1-bk_2),
\]

which has the wrong sign.

The full matrix has these six diagonal entries:

\[
c_1,\ c_2,\ c_3,\ c_1^*,\ c_2^*,\ c_3^*,
\]

and off-diagonal director-switching entries \(\gamma/2\) between neighbouring
director states on the six-state ring.

## Why Basic Checks Did Not Catch The Bugs

The old calculation can still look plausible because simple checks are too
weak:

- At \(k=0\), column sums vanish.
- There is a zero eigenvalue at \(k=0\).
- Probability conservation is not obviously violated.

The reason is that the sign error is invisible at \(k=0\), because

\[
\sin(0)=0.
\]

The check that catches the sign error is sixfold lattice symmetry. The
corrected matrix is spectrally invariant under a \(60^\circ\)
reciprocal-space rotation; the old matrix is not.

Numerically:

\[
\text{sixfold error}_{\rm corrected}\approx 2.27\times10^{-15},
\]

\[
\text{sixfold error}_{\rm old}\approx 1.54\times10^{-3}.
\]

## KMC Verification

The KMC is not Python. It is Fortran.

Fortran code:

```text
kmc_triangular_jmvr.f90
```

Fortran output:

```text
kmc_triangular_counts.txt
```

The count file contains \(N=100,000,000\) walkers. Python reads this file and
normalizes:

\[
P_{\rm KMC}(n_1,n_2,t)
=
\frac{\text{count}(n_1,n_2)}{100000000}.
\]

The clean four-way comparison is:

| Case | RMS vs KMC | Ratio to MC noise | Meaning |
|---|---:|---:|---|
| rectangular notebook grid + old \(c_3\) | \(4.448\times10^{-4}\) | \(105.34\) | original-style theory, bad |
| rectangular notebook grid + corrected \(c_3\) | \(4.591\times10^{-4}\) | \(108.73\) | sign alone does not fix it |
| sheared triangular grid + old \(c_3\) | \(1.428\times10^{-4}\) | \(33.82\) | PBC alone improves but not enough |
| sheared triangular grid + corrected \(c_3\) | \(3.434\times10^{-6}\) | \(0.81\) | correct, at noise floor |

The MC noise floor is

\[
\text{MC noise floor}
\approx
4.223\times10^{-6}.
\]

So neither correction alone is enough. Fixing both gives agreement within
Monte Carlo noise.

The final Fig. 11 red curve isolates the sign error **after** the k-grid has
already been corrected:

\[
\mathrm{RMS}_{\rm sheared+old\,sign}
=
1.428\times10^{-4},
\]

\[
\mathrm{RMS}_{\rm sheared+corrected}
=
3.434\times10^{-6},
\]

## What To Show Tomorrow

Show in this order:

1. `README.md` so the folder structure is clear.
2. `triangular_jmvr_corrected.py` for the corrected matrix.
3. `forensic_two_bugs.py` for the four-way PBC/sign comparison.
4. `verify_realspace_bloch.py` if he asks how we know the sheared grid is the correct torus grid.
5. `PI_NOTE_TWO_BUGS_AND_CHECKS.md` for the written explanation.
6. `kmc_triangular_jmvr.f90` to show KMC is independent of the matrix.
7. `fig11_final_hex.pdf` or `fig11_final_hex_gnuplot.pdf`.

Do not start by showing every old debugging script. Those are archived in
`legacy_debug/`.

## How To Say It To Kabir

> I started with the triangular JMVR track, not triangular TCRW yet. I found
> two issues. First, Kabir's PBC suspicion was right: the old rectangular
> Fourier grid is not the correct grid for the triangular torus. Second, after
> fixing the grid, the \(c_3\) chirality term still has the wrong sign. Neither
> correction alone matches KMC. Fixing both gives RMS \(3.43\times10^{-6}\),
> which is at the Monte Carlo noise floor.

## What To Ask Kabir

Do not ask "Should I switch to TCRW now?" immediately.

Ask:

> Should I now finish the corrected triangular JMVR analysis first: proper
> \(\Gamma-M-K-\Gamma\) band structure, gap scan, and real-space/Bloch spectrum
> check?

Then add:

> After this JMVR track is clean, I can start Track B: triangular version of
> the TCRW/topology model.

This frames TCRW as the next research direction, not as an escape from the
current result.

## Likely PI Questions

### Is this exact or Monte Carlo?

Both. The blue/corrected theory curve is exact finite-torus Fourier theory.
The black dots / KMC heatmap are independent Fortran Monte Carlo.

### Why is the matrix \(6\times6\)?

Because triangular lattice has six nearest-neighbour directions, and the walker
has one internal director state for each direction.

### What does \(\epsilon\) do?

\(\epsilon\) biases forward versus backward translation along the current
director:

\[
\frac16+\epsilon
\quad \text{vs.} \quad
\frac16-\epsilon.
\]

### What does \(\gamma\) do?

\(\gamma\) controls orientation switching. The walker turns to each adjacent
director at rate \(\gamma/2\), so total turning rate is \(\gamma\).

### Why did the old result look plausible?

Because both issues are hidden from the simplest checks. The rectangular grid
still gives a normalized-looking answer, and the wrong \(c_3\) sign vanishes at
\(k=0\). You need a torus-periodicity check, KMC comparison, and sixfold
symmetry test to see the failures.

### Is the PBC bug still relevant?

Yes. The PBC/k-grid bug is one of the two necessary corrections. A real-space
Bloch/PBC spectrum check has now been done on small tori in
`verify_realspace_bloch.py`; a larger/final version can be included later if
needed for publication.

## Next Work After The Meeting

1. Corrected \(\Gamma-M-K-\Gamma\) band structure.
2. Gap scan in \((\gamma,\epsilon)\).
3. Real-space transition matrix versus Bloch spectrum check.
4. Only then Track B: triangular TCRW/topology.
