# Triangular active walker: two-bug diagnosis and why basic sanity checks missed it

## One-line result

The mismatch in Dipanjan's triangular active-walker calculation needs two corrections. First, the Fourier modes must live on the sheared triangular torus, not the rectangular notebook grid. Second, one chirality sign in the state-2 coefficient is inconsistent with the real-space master equation.

The correct triangular-torus grid is

$$
k_1=\frac{\pi m_1}{aL},
\qquad
k_2=\frac{\pi(2m_2-m_1)}{bL}.
$$

The old notebook grid was rectangular:

$$
k_1=\frac{2\pi n_x}{2aL},
\qquad
k_2=\frac{2\pi n_y}{bL}.
$$

For the triangular period \(T_2=(aL,bL)\), the old grid gives

$$
\mathbf{k}\cdot T_2=\pi n_x+2\pi n_y,
$$

which is not a multiple of \(2\pi\) for odd \(n_x\). This is the PBC/k-grid issue Kabir suspected.

The correction is

$$
c_3 = B - 1 - \gamma + 2 i \epsilon \sin(a k_1 - b k_2)
$$

to

$$
c_3 = B - 1 - \gamma - 2 i \epsilon \sin(a k_1 - b k_2).
$$

With both corrections, the exact theory agrees with high-statistics KMC at the Monte Carlo noise level.

## Model being checked

The triangular walker has six internal director states. In director state \(d\), the particle hops to all six nearest-neighbour directions. The forward/backward pair along the director axis has biased rates

$$
\frac{1}{6} + \epsilon,
\qquad
\frac{1}{6} - \epsilon,
$$

and the other four directions have rate \(1/6\). The director switches symmetrically,

$$
d \to d+1
\quad \text{and} \quad
d \to d-1,
$$

each at rate

$$
\frac{\gamma}{2}.
$$

So the Bloch generator is a \(6 \times 6\) continuous-time matrix \(M(k_1,k_2)\), not a discrete-time stochastic matrix.

## What Dipanjan's setup got right

Dipanjan's notebook has the correct overall structure:

1. The diagonal entries contain the Fourier transform of hopping.
2. The off-diagonal entries are \(\gamma/2\) between neighbouring director states on the six-state ring.
3. At \(k=0\), probability is conserved:

$$
\sum_i M_{ij}(0,0)=0.
$$

4. At \(\epsilon=0\), the sign issue disappears and the matrix becomes real.
5. Opposite director states are written as complex conjugate pairs.

These are all good signs, but they are not enough to prove the torus grid and every chirality sign are correct.

## The k-grid / PBC issue

The triangular lattice is naturally represented by

$$
\mathbf{r}=n_1\mathbf{a}_1+n_2\mathbf{a}_2,
\qquad
\mathbf{a}_1=(2a,0),
\qquad
\mathbf{a}_2=(a,b).
$$

On an \(L\times L\) triangular torus,

$$
n_1\equiv n_1+L,
\qquad
n_2\equiv n_2+L.
$$

Therefore the physical periods are

$$
T_1=L\mathbf{a}_1=(2aL,0),
\qquad
T_2=L\mathbf{a}_2=(aL,bL).
$$

Allowed Bloch modes must obey

$$
e^{i\mathbf{k}\cdot T_1}=1,
\qquad
e^{i\mathbf{k}\cdot T_2}=1.
$$

Solving these two constraints gives the sheared grid

$$
k_1=\frac{\pi m_1}{aL},
\qquad
k_2=\frac{\pi(2m_2-m_1)}{bL}.
$$

This is what the corrected code uses.

The old rectangular grid is a different torus convention. It gives a normalized-looking answer, but it does not match the triangular \(L\times L\) torus used by the KMC.

## The sign error

From Appendix B, equation B4 for state 2 contains the biased incoming terms

$$
\left(\frac{1}{6}-\epsilon\right)P_2(x-a,y+b)
$$

and

$$
\left(\frac{1}{6}+\epsilon\right)P_2(x+a,y-b).
$$

The \(+\epsilon\) term is incoming from source \((x+a,y-b)\) to destination \((x,y)\). Therefore the actual hop displacement associated with the forward-biased rate is

$$
(-a,b).
$$

Using the Fourier convention

$$
\tilde P(k_1,k_2)
=
\sum_{x,y}
e^{i(k_1 x+k_2 y)}
P(x,y),
$$

the biased part contributes

$$
\epsilon e^{i(-a k_1+b k_2)}
-
\epsilon e^{i(a k_1-b k_2)}
=
-2 i \epsilon \sin(a k_1-b k_2).
$$

So the correct coefficient is

$$
c_3
=
B - 1 - \gamma
- 2 i \epsilon \sin(a k_1-b k_2),
$$

not

$$
c_3
=
B - 1 - \gamma
+ 2 i \epsilon \sin(a k_1-b k_2).
$$

## Why the usual sanity checks did not catch it

### 1. Probability conservation at \(k=0\) cannot see the bug

At \(k_1=k_2=0\),

$$
\sin(a k_1-b k_2)=0.
$$

So both the wrong and corrected coefficients reduce to the same value:

$$
c_3(0,0)=B(0,0)-1-\gamma=-\gamma.
$$

Therefore the column-sum check

$$
\sum_i M_{ij}(0,0)=0
$$

passes even with the wrong sign.

### 2. The \(\epsilon=0\) check cannot see the bug

When

$$
\epsilon=0,
$$

the entire chirality term vanishes:

$$
\pm 2 i \epsilon \sin(a k_1-b k_2)=0.
$$

So the wrong-sign and correct-sign matrices become identical in the achiral limit.

### 3. The conjugation symmetry can still pass

Dipanjan wrote opposite director states as conjugates:

$$
c_6 = c_3^*.
$$

If \(c_3\) has the wrong sign, but \(c_6\) is still set to its conjugate, the matrix can still look internally symmetric. This catches broken conjugation, but not a wrong physical assignment of the biased direction.

### 4. Nonzero-\(k\) column sums are not required to vanish

For the real-space generator, probability conservation means the total probability is conserved. In Fourier space, that conservation condition is imposed at the zero mode:

$$
k=0.
$$

For \(k \neq 0\), the column sums of \(M(k)\) do not have to be zero. So a nonzero-\(k\) column-sum test is not a valid detector for this sign error.

### 5. The error is not "losing probability"; it is biasing the wrong diagonal direction

The matrix still represents a process with the same total outgoing rate. The problem is subtler: for one internal state, the active bias points along the wrong triangular-lattice diagonal. That changes the shape of the probability distribution, but does not obviously break normalization.

## How we verified it

We compared exact theory against independent kinetic Monte Carlo. The KMC was done in lattice coordinates and does not use Dipanjan's Fourier matrix.

The high-statistics run used

$$
N = 100,000,000
$$

walkers.

The four-way comparison is:

| Case | RMS vs KMC | Ratio to MC noise |
|---|---:|---:|
| rectangular notebook grid + old \(c_3\) | \(4.448\times10^{-4}\) | \(105.34\) |
| rectangular notebook grid + corrected \(c_3\) | \(4.591\times10^{-4}\) | \(108.73\) |
| sheared triangular grid + old \(c_3\) | \(1.428\times10^{-4}\) | \(33.82\) |
| sheared triangular grid + corrected \(c_3\) | \(3.434\times10^{-6}\) | \(0.81\) |

The Monte Carlo noise floor is

$$
4.223 \times 10^{-6}.
$$

So neither correction alone is enough. The k-grid/PBC correction reduces the error, but the wrong \(c_3\) sign remains visible. Fixing both reaches the Monte Carlo noise floor.

As an independent torus check, `verify_realspace_bloch.py` builds the full
finite real-space generator for small \(L\) and compares its spectrum with the
Bloch spectra. The corrected sheared grid matches the real-space generator to
machine precision:

$$
\max_\lambda \min_{\lambda'} |\lambda_{\rm real}-\lambda'_{\rm sheared}|
\sim 10^{-14}.
$$

The rectangular notebook grid contains the physical modes plus extra
off-torus/anti-periodic modes. That is why it can look normalized but still
produce the wrong inverse-Fourier probability distribution.

The final Fig. 11 cross-section isolates the sign error after using the corrected triangular grid:

$$
\mathrm{RMS}(\mathrm{KMC},\mathrm{sheared\ grid+old\ sign})
=
1.428 \times 10^{-4},
$$

while

$$
\mathrm{RMS}(\mathrm{KMC},\mathrm{sheared\ grid+corrected\ sign})
=
3.434 \times 10^{-6}.
$$

## Why they could have missed it

This is easy to miss because both bugs survive simple sanity checks.

The rectangular grid gives a normalized-looking inverse Fourier sum, but on the wrong torus. The sign error is local to one of the three independent axes. It does not break the matrix dimension, the \(k=0\) conservation law, the \(\epsilon=0\) limit, or the conjugate-pair structure.

It only appears when all three conditions are true:

1. \(\epsilon \neq 0\),
2. \(k_1,k_2\) are away from zero,
3. the real-space direction of the biased hop is checked against the Fourier phase.

In other words, this bug is invisible to broad matrix sanity checks. It is visible only through a term-by-term derivation from the master equation, or through an independent simulation.

## Retrospective checks that would have caught it

### 1. Sixfold lattice-symmetry check

The corrected triangular matrix should respect the sixfold rotation symmetry of the lattice. In the lattice-coordinate convention, define

$$
q_1 = 2 a k_1,
\qquad
q_2 = a k_1 + b k_2.
$$

A \(60^\circ\) rotation of the direct lattice induces the reciprocal-space transformation

$$
(q_1,q_2)
\to
(q_1-q_2,q_1),
$$

or equivalently the inverse rotation

$$
(q_1,q_2)
\to
(q_2,q_2-q_1).
$$

The spectrum should obey

$$
\mathrm{spec}\,M(q_1,q_2)
=
\mathrm{spec}\,M(q_1-q_2,q_1).
$$

The corrected matrix passes this check up to numerical precision. The buggy matrix fails it because one director's biased direction is inconsistent with the other five.

This check is powerful, but it is not usually the first check people run because it requires writing the lattice symmetry explicitly in the same coordinate convention as the code.

### 2. Direct master-equation integration

Another way to catch the bug is to take the real-space master equations literally and step them forward in time, without hand-transcribing the Fourier sine signs. For example,

$$
\frac{dP}{dt} = M_{\mathrm{real}}P.
$$

If \(M_{\mathrm{real}}\) is built directly from a jump table and then compared with the hand-written Bloch-matrix result, the discrepancy reveals that the Fourier transcription is wrong.

This is different from checking probability conservation. It checks whether the matrix encodes the same microscopic jumps as the master equation.

### 3. Independent kinetic Monte Carlo

The strongest practical check is what we did: simulate the stochastic process directly on the same triangular lattice, using the transition rates from the master equation, and compare to the exact solution.

This is also why Dipanjan could have seen the mismatch but not immediately found the cause. If the hand-written theory passes the simple conservation checks, then a mismatch with KMC can look like a boundary-condition or simulation bug. The KMC tells us that something is inconsistent. The torus-periodicity/spectrum check identifies the k-grid/PBC issue, and the term-by-term master-equation derivation identifies the \(c_3\) sign issue.

## Future-proof checks we should use

For future active-walker matrices, use these checks in order:

1. Check probability conservation at \(k=0\):

$$
\sum_i M_{ij}(0)=0.
$$

2. Check the achiral limit:

$$
\epsilon=0
$$

should remove all imaginary chirality terms.

3. Check triangular sixfold spectral symmetry:

$$
\mathrm{spec}\,M(q_1,q_2)
=
\mathrm{spec}\,M(q_1-q_2,q_1).
$$

4. Build the Bloch matrix automatically from a real-space jump table instead of hand-writing sine signs.

5. For each internal state, check the short-time drift:

$$
\left.\frac{d}{dt}\langle r\rangle\right|_{t=0}
=
\sum_{\Delta r}
(\Delta r)\,w(\Delta r).
$$

The drift must point along the intended director.

6. Compare against KMC at one or two parameter points before trusting the exact matrix.

7. For publication-quality results, keep the full pipeline:

$$
\text{master equation}
\to
\text{Bloch matrix}
\to
\text{exact theory}
\to
\text{KMC check}
\to
\text{gnuplot from .txt data}.
$$

## How to explain this to Kabir

The concise explanation is:

> Kabir's PBC suspicion was partly right. Dipanjan's rectangular Fourier grid is not the correct grid for the triangular \(L\times L\) torus. Fixing the grid improves the theory, but does not bring it to KMC. Then I checked Appendix B term by term and found that \(c_3\) also has the wrong chirality sign. Neither correction alone reaches the noise floor. Fixing both gives RMS \(3.43\times10^{-6}\), while the MC noise floor is \(4.22\times10^{-6}\).
