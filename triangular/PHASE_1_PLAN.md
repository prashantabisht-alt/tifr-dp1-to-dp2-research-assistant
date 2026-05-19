# Phase 1 plan — chiral random walker on the triangular lattice

**Project**: Chiral run-and-tumble walker on the triangular lattice (single
walker, transport + first-passage + edge modes).
**Goal**: Build the first lattice chiral RTW with a six-state director, on
the triangular lattice. Compute transport (including the odd-diffusivity
tensor $D^{\rm odd}_{ij}$), first-passage statistics, and the OBC edge
spectrum. Closest lattice prior art is Wójcik–Kalz 2026 (square, $C_4$,
discrete-time coin-step); the chirality-rate mechanism is borrowed from
Mallikarjun–Pal 2023 (continuum, 4-director). Write a paper.

**Date**: 15 May 2026 · **Author**: P. Bisht (TIFR Hyderabad)
**Companion** : the bug-fixed `triangular_jmvr_corrected.py` infrastructure
(JMVR translation-chirality version, kept for reference and reuse).

**Revision 1 (15 May 2026 evening)** : technical corrections after a careful
read-through:
- $C_6$ rotational symmetry of the spectrum holds for *all* $b$, not just
  $b = 0$. What chirality breaks is the *mirror* (left/right) symmetry.
- The Bloch matrix is intrinsically non-Hermitian at any $b$, because
  translation is deterministic-forward. "Achiral" means $\gamma_+ = \gamma_-$,
  not "real-symmetric matrix."
- Steady-state drift $\vec{v} = 0$ by symmetry (uniform director distribution
  $\Rightarrow$ no preferred direction). Chirality shows up in $D_{ij}^{\rm odd}$
  and in circulation, not in $\langle \vec{r}(t) \rangle$.
- Odd diffusivity $D_{ij}^{\rm odd}$ cannot be extracted from $\lambda_0(\vec{k})$
  alone (antisymmetric tensor contracted with $k_i k_j$ vanishes). Must use
  velocity correlator / response formalism.
- The 120° rotation channel is *stretch*, not default. Phase 1 first does
  $\{+1, -1, +3\}$ — the direct ref-44 analogue.
- Whitelam-Klymko-Mandal moved to Phase 2 reading; Phase 1 is single-walker
  only.

**Revision 2 (19 May 2026, after follow-up with Kabir)** : positioning sharpened
after verifying ref 44 directly from the authors' companion code.
- **Mallikarjun–Pal 2023 (TCRW ref 44) is a continuum $\mathbb{R}^2$ model**
  with a 4-state discrete director, *not* a square-lattice walker. Verified
  from `papers/RMStatMech_companion_code/Chiral_RTP.ipynb`: position update
  is `x += v*tau*cos(theta_state)` in continuous space, absorbing walls live
  at continuum $x = \pm L$. Kabir was also unaware of this until today.
- Consequence: refs 42–44 do *not* form a "square-lattice chiral RTW" trio.
  Ref 42 (Hargus–Epstein–Mandadapu) is the abstract $D^{\rm odd}$ concept,
  ref 43 (Sevilla) is continuum, ref 44 (M–P) is continuum-with-discrete-director.
  The only true lattice prior art with chirality + odd diffusion is
  **Wójcik–Kalz 2026 (square)**.
- Re-framing: our paper is *not* "the triangular version of ref 44." It is
  *the lattice analogue of Mallikarjun–Pal's continuum chirality mechanism,
  the natural $C_6$ counterpart of Wójcik–Kalz, and the first lattice CRTP
  with a six-state director.* See revised §1 below.
- Practical consequence: the right sanity-check baseline before going to
  triangular is **Wójcik–Kalz** (lattice, square), not Mallikarjun–Pal
  (continuum). M–P's MFPT plot cannot serve as a lattice bit-match: the
  two models only agree in the diffusive limit, not at finite $v$.

**Revision 3 (19 May 2026 night, after careful Wójcik–Kalz read)** :
Wójcik–Kalz is the right lattice reference, but it is **not** the same
stochastic process as ours.
- Wójcik–Kalz use a **discrete-time internal-degree-of-freedom walk**:
  $|\rho(t+1)\rangle = U|\rho(t)\rangle$, with
  $U = S(I\otimes C)$ and
  $C_{\rm CRW}=(1-p)C_{\rm rand}+pC_{\rm chir}$.
- Our model is a **continuous-time Poisson run-and-tumble process** with
  rates $v,\gamma_+,\gamma_-,\gamma_r$. Therefore Wójcik–Kalz is a
  structural/baseline paper, not a bit-match target for the same code.
- Their square model gives
  $D_{xx}=D_{yy}=\frac{a^2}{2\tau}\frac{1-p^2}{1+p^2}$ and
  $D_{xy}=-D_{yx}=\frac{a^2}{2\tau}\frac{2p}{1+p^2}$, so
  $D_{xy}/D_{xx}=2p/(1-p^2)$ diverges as $p\to1$. Our continuous-time
  triangular model has a finite jump-noise floor in $D_{\rm even}$; for
  $v=\gamma=1,\gamma_r=0,b=1$ it gives
  $D_{\rm even}=1/2$, $D_{\rm odd}=\sqrt3/4$, and
  $D_{\rm odd}/D_{\rm even}=\sqrt3/2$.
- New framing: Wójcik–Kalz develop a **discrete-time topological coin-walk
  picture on square**; we develop the **continuous-time chiral RTW picture
  on triangular**. Complementary, not duplicate.
- Action change: implement a minimal Wójcik–Kalz discrete-time reproducer
  separately when needed. Keep `square_pal/square_chiral_rtw.py` as our
  continuous-time $C_4$ analogue, not as a literal Wójcik–Kalz module.

---

## 1. Mission statement, in one paragraph

We have a verified $6\times 6$ Bloch matrix and KMC infrastructure for
the JMVR-style triangular walker (translation chirality $\epsilon$).
We now extend to the **rotation-chirality** variant: a chiral
run-and-tumble walker on the triangular lattice where the director
$m \in \{0,\ldots,5\}$ rotates clockwise and counter-clockwise at
unequal rates ($\gamma_+ \ne \gamma_-$), with optional reversal
$\gamma_r$ and deterministic translation in the current director.
This model has not been studied before: Mallikarjun–Pal 2023 used a
related left/right/reversal chirality mechanism but in continuum
$\mathbb{R}^2$ with only four director states, and Wójcik–Kalz 2026
developed a different discrete-time coin-step chiral walk on the square
($C_4$) lattice. We compute the symmetric and antisymmetric
parts of the diffusion tensor ($D^{\rm sym}_{ij}$, $D^{\rm odd}_{ij}$),
real-space propagators, first-passage statistics, search-time
optimisation, and the open-boundary spectrum (edge modes). The headline
$C_6$-only result is the closed form $D^{\rm odd}/D^{\rm even} = \sqrt 3/2$
at $b=1, \gamma_r=0$, $v=\gamma=1$ — larger than the corresponding
continuous-time $C_4$ analogue and qualitatively different from the
diverging Wójcik–Kalz discrete-time ratio (verified to machine precision;
see `triangular/verify_diffusion_tensor.py`).

## 2. The model, precisely

### 2.1 State space

Lattice indices $(n_1, n_2) \in \mathbb{Z}^2$, director $m \in \{0,\ldots,5\}$.
NN displacements as before:
$\hat{e}_0 = (+1, 0),\, \hat{e}_1 = (0, +1),\, \hat{e}_2 = (-1, +1),\,
\hat{e}_3 = (-1, 0),\, \hat{e}_4 = (0, -1),\, \hat{e}_5 = (+1, -1)$.

Implementation convention: use these integer lattice coordinates internally.
For plotting, map to Cartesian triangular coordinates
\[
\vec{r}=n_1\vec{a}_1+n_2\vec{a}_2,\qquad
\vec{a}_1=(1,0),\quad \vec{a}_2=\left(\tfrac12,\tfrac{\sqrt3}{2}\right).
\]
This avoids the old rectangular-grid/PBC mistake. In axial coordinates the
finite-torus Fourier labels are simple; in Cartesian coordinates the same
labels appear as a sheared reciprocal grid.

### 2.2 Rates (chiral RTW / Mallikarjun–Pal style)

Per state, three event types:

- **Translation**: at rate $v$, walker hops deterministically in current
  director $m$: $(n_1, n_2) \to (n_1, n_2) + \Delta(m)$. No bias to other
  NN, no backward, no sideways. (This is the "run" of run-and-tumble.)
- **Rotation by +1 (CCW one step)**: at rate $\gamma_+ = \gamma (1 + b)/2$.
- **Rotation by −1 (CW one step)**: at rate $\gamma_- = \gamma (1 - b)/2$.
- **Reversal $m \to m + 3$**: at rate $\gamma_r$. Mallikarjun–Pal include this.
  We keep it in the equations/code from the beginning, but can set
  $\gamma_r=0$ for the simplest first plots.

Here $\gamma$ is the total rotation rate and $b \in [-1, 1]$ is the
chirality bias. Achiral: $b = 0$.

Total event rate per walker: $v + \gamma + \gamma_r$.

### 2.3 Decisions still open (mark as we lock them in)

- [x] Include the reversal channel $\gamma_r$ in the model API.
- [ ] Whether to also allow $\pm 2$ rotation (turn by $120^\circ$, no
  square analogue — *new triangular physics*). Not Phase 1 default. If yes
  later, there are two independent chirality biases, $b_{60}$ and $b_{120}$.
- [ ] Normalisation convention: fix $v = 1$ and report $\gamma / v$,
  or fix $v + \gamma = 1$ and report $\gamma$?
  *Recommendation*: $v = 1$, $\gamma$ free, matches Mallikarjun–Pal.
- [x] Start with direction-only translation (run-and-tumble style), not the
  old JMVR all-six-neighbour translation-bias rule.

### 2.4 Master equation

\[
\partial_t P_m(\vec{n}, t) = v\, P_m(\vec{n} - \Delta(m), t)
+ \gamma_+\, P_{m-1}(\vec{n}, t) + \gamma_-\, P_{m+1}(\vec{n}, t)
+ \gamma_r\, P_{m+3}(\vec{n}, t)
- (v + \gamma + \gamma_r)\, P_m(\vec{n}, t)
\]

Note signs of rotation: $\gamma_+$ feeds inflow *from* $m-1$ into $m$
(walker was in $m-1$, rotated $+1$, landed in $m$); $\gamma_-$ feeds
from $m+1$.

### 2.5 Bloch matrix $M(\vec{k})$

Use axial Fourier coordinates internally. On a finite triangular torus these
are equivalent to the sheared Cartesian reciprocal grid from the JMVR bug fix.
Diagonal:
\[
M_{m,m}(\vec{k}) = v\, e^{i \vec{k} \cdot \hat{e}_m} - (v + \gamma + \gamma_r)
\]
Off-diagonal: rotation contributions
\[
M_{m, m+1}(\vec{k}) = \gamma_- \quad \text{(inflow from $m+1$)}, \qquad
M_{m, m-1}(\vec{k}) = \gamma_+ \quad \text{(inflow from $m-1$)}.
\]
And the reversal: $M_{m, m+3}(\vec{k}) = \gamma_r$.

**Symmetries (corrected).** Under a 60° lattice rotation paired with
the director relabel $m \to m+1$, we have $M(R_{60}\vec{k}) = P M(\vec{k}) P^{-1}$
with $P$ the cyclic-permutation matrix on director indices. So
$\text{spec}\,M(R_{60}\vec{k}) = \text{spec}\,M(\vec{k})$ holds for *all*
$b$ — including the chiral case. **The $C_6$ rotation symmetry is *not*
broken by chirality.**

What chirality breaks is the *mirror* (left/right) symmetry. Under a
reflection (e.g., $k_y \to -k_y$ paired with $m \to -m$), the operations
$+1$ and $-1$ on the director swap, so $\gamma_+ \leftrightarrow \gamma_-$.
For the spectrum to be mirror-invariant we need $\gamma_+ = \gamma_-$,
i.e. $b = 0$.

**The matrix is non-Hermitian even at $b = 0$.** Translation is
deterministic forward, so the diagonal $v e^{i\vec{k}\cdot\hat{e}_m}$
is intrinsically complex (no forward-backward symmetrisation as in
JMVR). At $b = 0$ the matrix gains $\gamma_+ = \gamma_-$ which makes the
*rotation block* Hermitian, but the full $M(\vec{k})$ stays
non-Hermitian whenever $\vec{k} \ne 0$.

## 3. Reading list, in priority order

| # | Paper | Why | What to extract |
|---|---|---|---|
| 1 | **Wójcik & Kalz, arXiv:2602.09920 (2026)** | **Closest lattice prior art, but discrete-time.** Square-lattice coin-step chiral random walk with odd diffusivity and edge modes. It is not our continuous-time RTW, but it sets the lattice/topology language | Their $D^{\rm odd}$ formula, edge-current/topological framing, fidelity-decay diagnostic, and their explicit call for other lattice geometries. Use as a structural benchmark, not a bit-match target |
| 2 | Mallikarjun & Pal, Physica A 622, 128821 (2023), arXiv:2209.05912 | **The chirality mechanism reference.** Continuous space with four discrete orientations, *not* a square-lattice spatial model. Rate convention $\Gamma_\pm$ for left/right tumbling + $\Gamma_2$ reversal | Their rate convention; their MFPT-via-Laplace-transform technique (`CRTP.nb` in our `papers/RMStatMech_companion_code/`); their optimal-bias observable. Note: cannot serve as a lattice bit-match baseline |
| 3 | Hargus, Epstein, Mandadapu, PRL 127, 178001 (2021) | Definition of the odd diffusivity tensor; this is what makes $D^{\rm odd}$ a publishable observable | The Green–Kubo / velocity-correlator definition; the chirality-signature framing |
| 4 | Osat, Meyberg, Metson, Speck arXiv:2602.12020 (TCRW) | Square-lattice edge modes + topological band structure; provides the recipe we copy in Step 6 | Edge-mode protocol; OBC spectrum visualisation; the C_4 baseline figure templates |
| 5 | Sevilla, PRE 94, 062120 (2016) | Continuum chiral active particle (no lattice) | Continuum-limit cross-check formulas; sanity check on small-$k$ expansion |
| 6 | Gilbert & Sanders, PRE 80, 041121 (2009) | Persistent random walk / triangular Lorentz-gas near miss | How persistence and triangular geometry enter diffusion; useful citation to avoid overclaiming "first triangular" |
| 7 | Marris, Sarvaharman & Giuggioli (2023); Marris & Giuggioli (2024) | Persistent/anti-persistent lattice walks and first passage in domains | FPT methods and boundary-domain results relevant to our MFPT section |
| 8 | Oropesa, de Castro, Löwen & Liarte, arXiv:2602.04732 (2026) | Triangular-lattice RTP-adjacent active matter on networks | Distinguish their trail-mediated triangular RTP from our clean single-particle chiral RTW |
| -- | Whitelam-Klymko-Mandal, arXiv:1709.03951 (2017) | *Phase 2 reading* — multi-walker hard-core lattice ABP | Hard-core lattice rules; MIPS observable (defer until Phase 1 done) |

**Action status (17 May 2026)**: the relevant PDFs are now collected in
`triangular/papers/`, including the dangerous near-miss papers and the
Panagiotopoulos JCP 2005 lattice hard-sphere reference for Phase 3.

**Read first** (revised priority): read **Wójcik–Kalz** *thoroughly*. It is
the closest lattice prior art and the paper our results will be most
forcefully compared against, but remember it is a discrete-time coin walk.
After that, read Mallikarjun–Pal for the
chirality mechanism + the MFPT Laplace-transform technique (use their
`CRTP.nb`). Hargus–Epstein–Mandadapu next for the $D^{\rm odd}$ definition.
Then do a quick novelty audit of Gilbert–Sanders and Marris/Giuggioli
before writing the first paragraph of the Phase-1 paper.

**Priority-claim guardrail**: do not write "chiral random walk on triangular
lattice has never been done" as a blanket statement. Papers exist with
nearby titles and triangular persistent/chiral walks. The sharper claim is:
we study the **active six-director rotation-chiral run-and-tumble walker on
the triangular lattice**, and compute the transport (including the odd
diffusivity tensor), first-passage, and edge/localisation observables that
Wójcik–Kalz studied on the square lattice and that Mallikarjun–Pal and
Hargus–Epstein–Mandadapu studied in the continuum.

Do not deep-read the hard-hexagon/equation-of-state papers during Phase 1.
They motivate Phase 2/3, but the immediate deliverable is still the
single-particle chiral triangular walker.

## 4. Infrastructure inventory (reuse / extend / write new)

| Module | What it does | Phase-1 status |
|---|---|---|
| `triangular_jmvr_corrected.py` | $6\times 6$ Bloch matrix (JMVR, translation chirality) | **Reference only** — do not mix the new model into this file |
| `triangular_chiral_rtw.py` | New canonical Phase-1 Python module | **Write new** — rotation-chiral run-and-tumble walker |
| `kmc_triangular_jmvr.f90` | 100M-walker Fortran KMC | **Use as template** — write a separate chiral-RTW KMC, do not overwrite the JMVR verifier |
| `forensic_two_bugs.py` | $\{$rect, shear$\} \times \{$buggy, fixed$\}$ vs KMC | **Reuse** as cross-check template — same structure, different model |
| `fig11_final_hex.py`, `fig11_original_style.gnu` | Paper-style figures | **Reuse** for new $P(\vec{r}, t)$ heatmaps |
| `verify_realspace_bloch.py` | Bloch ↔ real-space generator consistency | **Reuse** verbatim |
| `jmvr_single_walker.html`, `tcrw_single_walker.html` | Live walker viewer | **Reuse**; the TCRW one is already the rotation-chirality version |

What's genuinely new:
- Small-$\vec{k}$ expansion of $\lambda_0(\vec{k})$ → zero drift plus the
  symmetric diffusion tensor $D^{\rm sym}_{ij}$ analytically.
- Velocity-correlation / response calculation for the odd part
  $D^{\rm odd}_{ij}$.
- OBC strip generator + edge-localisation computation.
- MFPT solver (real-space generator with absorbing boundary).

## 5. Work sequence

No calendar. Step $n$ starts when step $n-1$'s result is solid.

### Step 1 · Build the chiral-RTW Bloch matrix

Write `triangular_chiral_rtw.py` as a **new canonical file**, separate
from `triangular_jmvr_corrected.py`, so the rotation-chirality model
does not get muddled with the old translation-chirality JMVR.

Module exposes:
- `build_Mk_chiral_rtw(v, gamma, b, gamma_r, k1, k2, a, b_lat)` returning
  the $6\times 6$ matrix.

Self-checks (corrected):
- $\sum_m M_{m,j}(\vec{k}=0) = 0$ for each column $j$ (probability
  conservation; holds for all $b$).
- $\lambda_0(\vec{k}=0) = 0$ exactly (Perron eigenvalue; holds for all $b$).
- **$C_6$ spectral symmetry holds for all $b$**:
  $\text{spec}\,M(R_{60}\vec{k}) = \text{spec}\,M(\vec{k})$. *Do not*
  expect this to fail at $b \ne 0$.
- **Mirror symmetry fails when $b \ne 0$**: under
  $\vec{k} \to (k_x, -k_y)$ paired with director reflection, the spectrum
  differs unless $\gamma_+ = \gamma_-$.
- Bit-check at a few random $\vec{k}$ against the master-equation form.

### Step 2 · Small-$\vec{k}$ expansion → drift, symmetric diffusion, and (separately) the odd part

The scalar dispersion expansion is
\[
\lambda_0(\vec{k}) = -i \vec{v} \cdot \vec{k} - \tfrac{1}{2} D^{\rm sym}_{ij} k_i k_j + O(k^3)
\]
where $D^{\rm sym}_{ij}$ is the *symmetric* part of the diffusion tensor.
The antisymmetric (odd) part $D^{\rm odd}_{ij}$ does **not** appear in
$\lambda_0(\vec{k})$ — it contracts with $k_i k_j$ to zero.

**Drift expectation (corrected).** For uniform initial director the
stationary director distribution is $P_m = 1/6$ regardless of $b$. The
mean drift is $\vec{v} = v\sum_m P_m \hat{e}_m = (v/6)\sum_m \hat{e}_m = 0$,
since the six NN vectors sum to zero. So $\vec{v} = 0$ identically for
any $(b, \gamma, \gamma_r)$. This is the correct statement of "chirality
does not produce permanent drift on triangular." Chirality manifests as
**circulation / transient angular motion**, not as a steady $\vec{v}$.

**How to extract $D^{\rm odd}_{ij}$ (corrected).** Two equivalent routes:

(a) *Off-diagonal second moment* of the propagator:
\[
D^{\rm odd}_{xy} \;=\; \tfrac{1}{2}\,\lim_{t\to\infty}\frac{d}{dt}\langle x(t)y(t) - y(t)x(t)\rangle.
\]
Compute $\langle x(t) y(t) \rangle$ from the Bloch matrix by taking
$\partial_{k_x} \partial_{k_y} \Ptilde(\vec{k}, t)|_{\vec{k}=0}$ (or
equivalently as the Kubo-style $t\to\infty$ value of the integrated
cross-correlation).

(b) *Hargus–Epstein–Mandadapu velocity-correlator* form:
\[
D^{\rm odd}_{ij} \;=\; \tfrac{1}{2}\int_0^\infty dt\,\langle v_i(t) v_j(0) - v_j(t) v_i(0)\rangle,
\]
which on the lattice is a sum over director states weighted by the
$e^{Mt}$ propagator.

Both should agree. We compute (a) analytically from the Bloch matrix
and (b) numerically from the KMC, then check.

**Triangular-specific check.** Does the symmetric diffusion tensor at
leading order have $D^{\rm sym}_{ij} \propto \delta_{ij}$ (full
isotropy)? At higher order ($k^4$), do the leading anisotropic
corrections fall as 4-fold or 6-fold? Kabir's prediction: triangular
isotropises *more* than square, removing the 4-fold lattice signature.
Compute and compare to the square-lattice result.

**Output**: a sympy notebook with the symbolic expansion + a plot of the
leading $D^{\rm sym}_\perp$, $D^{\rm sym}_\parallel$, $D^{\rm odd}_{xy}$
coefficients as functions of $(b, \gamma, \gamma_r)$.

### Step 3 · $P(\vec{r}, t)$ on the torus: KMC vs analytic

- Write a separate Fortran KMC using the JMVR KMC as a template:
  deterministic translation, $\gamma_\pm$ rotation, optional reversal.
- **Achiral sanity check ($b = 0$).** Verify $P(\vec n, t)$ is isotropic
  ($\langle x^2 \rangle = \langle y^2 \rangle$, $\langle xy \rangle = 0$),
  and that at $t \gg 1/\gamma$ the lattice propagator approaches a 2D
  isotropic Gaussian with the analytic $D^{\rm sym}$ from Step 2. *Do not*
  attempt to bit-match Mallikarjun–Pal's continuum $P(x,y,t)$: their
  propagator is supported on $\mathbb{R}^2$, ours on the lattice; the two
  only agree in the diffusive limit, not at finite $v$.
- **Two square baselines, kept separate.**
  1. Build a minimal **Wójcik–Kalz discrete-time** reproducer:
     coin-step update $U=S(I\otimes C)$, reproduce their MSD formula
     $D_0^p=\frac{a^2}{2\tau}\frac{1-p^2}{1+p^2}$ and odd tensor
     $D_{xy}=\frac{a^2}{2\tau}\frac{2p}{1+p^2}$. This checks that we
     understand the lattice/topology reference.
  2. Build our own **continuous-time $C_4$ reduction** of the triangular
     RTW. This lives in `square_pal/square_chiral_rtw.py` and is the true
     apples-to-apples comparison for our triangular rates. It should match
     our Green-Kubo/exact formulas, not Wójcik–Kalz's discrete-time
     numbers exactly.
- Run at $b \ne 0$: visualise the propagator at $t \sim 1/\gamma$ and
  $t \gg 1/\gamma$. Expect transient circulation / rotating anisotropy in
  the ballistic-to-diffusive crossover, but no permanent centre-of-mass drift
  for a uniform initial director ensemble.

**Output**: figure analogous to our `fig11_original_style_gnuplot.png`
but for the chiral-RTW model. Verifies that the Bloch + KMC tooling
ports without bugs.

### Step 4 · First-passage on a half-plane and a strip

Real-space generator $M_\text{real}$ on a finite lattice with an absorbing
boundary. Two natural geometries:

- **Half-plane**: PBC in $n_1$ (cylinder), absorbing at $n_2 = 0$. MFPT
  to the boundary as a function of starting position $(n_1, n_2)$.
- **Strip**: absorbing at $n_2 = 0$ and $n_2 = L_y$, PBC in $n_1$.
  Splitting probability (left vs right edge), survival probability.

MFPT obtained as $-\langle \mathbf{1}^T (M_\text{real})^{-1} \mathbf{e}_{\vec{r}_0} \rangle$
where the inverse is on the absorbing subspace.

**Output**: MFPT heatmap over starting position, for several $(b, \gamma)$.

### Step 5 · Search-time optimisation

For fixed reversal rate $\gamma_r$ (and fixed $\gamma$), scan $b$ and
find the value $b^*$ that **minimises** the MFPT averaged over a
uniform initial distribution.

- Test whether the *qualitative phenomenon* Mallikarjun–Pal found in the
  continuum (an optimal chirality minimises the search time) survives on
  the triangular lattice. We are not claiming numerical agreement with
  their values — different geometry — only that the qualitative
  optimisation curve exists.
- **New question**: is $b^*$ larger or smaller on the triangular lattice
  than on the square lattice (from our own $C_4$ reduction, the same
  one we benchmark against Wójcik–Kalz in Step 3)? The expected answer
  is "different", and *why* it's different becomes a paragraph of the
  paper.

**Output**: $b^*(\gamma, \gamma_r)$ surface plot for triangular, plus the
analogous square ($C_4$) plot from our own code for comparison.

### Step 6 · Edge modes (OBC spectrum)

Build the OBC generator on a strip (say PBC in $n_1$, OBC in $n_2$).
Diagonalise. Look for eigenmodes localised at the edges:
- For each eigenvalue $\lambda$, compute $w_\partial = \sum_{\partial} |v|^2$
  (fraction of weight on the boundary rows).
- Plot eigenvalues coloured by $w_\partial$.

Two edge geometries (triangular only):
- **Zigzag** edge
- **Armchair** edge

**Compare** the edge-mode spectrum on triangular to the known
TCRW square-lattice edge modes. Kabir explicitly predicted:
*"same edge modes will come out, it will just be slightly different."*
We verify, and quantify what "slightly different" means.

**Output**: edge-mode spectrum + edge-current visualisation per
geometry, in the style of TCRW Fig.~3.

### Step 7 · Stress-test predictions and write up

- Sweep $(\gamma, b, \gamma_r)$ across the physically relevant range.
- For each observable (zero drift check, $D^{\rm sym}_{ij}$,
  $D^{\rm odd}_{ij}$, MFPT, $b^*$, edge-mode count and localisation),
  record the qualitative + quantitative result.
- Cross-reference, in two tiers:
  - *Lattice* (apples-to-apples): how does the triangular result compare
    to **Wójcik–Kalz 2026 (square, discrete-time)** and to our own
    continuous-time $C_4$ reduction at matched rates? This is where the
    $C_6$ story has to be told.
  - *Continuum* (limit-only): how does the diffusive-limit triangular
    result compare to Sevilla and to the diffusive limit of Mallikarjun–Pal?
    Hargus–Epstein–Mandadapu is the conceptual touchstone for $D^{\rm odd}$
    but is not a numerical comparator.

**Output**: a paper draft. Aim for the Physica A / PRE level (same
journals as refs 42–44).

## 6. Observables checklist (the "by the time we're done" list)

- [ ] Zero-drift check $\vec{v}=0$ for uniform initial director, all
      $(b,\gamma,\gamma_r)$.
- [ ] Symmetric diffusion tensor $D^{\rm sym}_{ij}$ from $\lambda_0(\vec{k})$.
- [ ] Odd diffusivity / antisymmetric response $D^{\rm odd}_{ij}$ from
      velocity correlations or response, not from scalar dispersion alone.
- [ ] Anisotropy at higher order (4-fold vs 6-fold? Kabir's prediction).
- [ ] Real-space propagator $P(\vec{r}, t)$ at three regimes
      (ballistic, crossover, diffusive).
- [ ] First-passage MFPT on half-plane, function of starting position.
- [ ] Splitting probability on a strip.
- [ ] Survival probability tail (exponential decay rate vs $b$).
- [ ] Optimal bias $b^*$ for search-time minimisation.
- [ ] Edge-mode spectrum, OBC strip, zigzag and armchair edges.
- [ ] Edge-current localisation, defect robustness check.
- [ ] (Optional, stretch) Berry curvature / Chern number on the
      $\vec{k}$-grid if the spectrum has a clean gap.

## 7. Verification protocol

Every analytic observable gets a KMC check at the appropriate scale:

- Bloch matrix vs real-space generator at small $L$: machine precision
  (`verify_realspace_bloch.py`-style).
- $P(\vec{r}, t)$ vs KMC: $10^7$ walkers for slider-friendly checks,
  $10^8$ for the headline figure (same as the JMVR bug-fix protocol).
- MFPT vs KMC: track first-passage time for $10^6$ independent walkers
  starting at a fixed site; histogram and compare to the analytic
  $-(M_\text{real})^{-1}$.
- Edge modes vs KMC: simulate on the OBC strip, accumulate the
  steady-state current along the boundary, compare to the eigenmode
  prediction.

## 8. Deliverables

By the end of Phase 1, we have:

- [ ] `triangular_chiral_rtw.py` — analytic infrastructure (analogue of
      `triangular_jmvr_corrected.py`).
- [ ] `kmc_triangular_chiral_rtw.f90` — high-stats KMC.
- [ ] `mfpt_solver.py` — real-space MFPT and splitting-probability solver.
- [ ] `edge_modes_obc.py` — OBC spectrum and edge-localisation analysis.
- [ ] One paper-ready figure for each major observable.
- [ ] A LaTeX paper draft modelled on the same elegant style as the
      derivation/project-report PDFs.

## 9. Risks and contingencies

- **Triangular result is only a mild extension of square-lattice ideas.**
  Then the paper becomes "continuous-time chiral RTW and geometry
  dependence of odd diffusion" — still publishable but less dramatic.
  *Mitigation*: lean into the sharp comparison:
  Wójcik–Kalz's discrete-time ratio diverges as $p\to1$, our
  continuous-time triangular ratio is finite at the fully chiral point,
  and our continuous-time $C_6$ value differs from the $C_4$ reduction.
  Also keep the 120°-rotation channel as a possible Phase 1.5 extension;
  it has no square analogue and may contain the genuinely new triangular
  physics.
- **MFPT solver is numerically unstable** (real-space generator is
  $6 L^2 \times 6 L^2$; inversion can be ill-conditioned). *Mitigation*:
  use sparse linear solver `scipy.sparse.linalg.spsolve`; check against
  iterative power method.
- **Edge modes are not well-defined on triangular** (because the gap
  closes everywhere). *Mitigation*: this would itself be a finding;
  document and discuss.
- **Kabir asks for something we didn't plan.** *Mitigation*: keep the
  scope modular so any one observable can be skipped or postponed
  without breaking the rest.

## 10. Connection to Phase 2 (and beyond)

Phase 1 builds the single-walker chiral RTW on triangular. Phase 2
adds hard-core exclusion (multi-walker) and asks about MIPS / clustering,
à la Whitelam–Klymko–Mandal. The single-walker propagator we compute
here is the clean microscopic motion rule that Phase 2 can build on,
but Phase 2 will need its own many-body simulation and theory.

Phase 3 is the long-term DP2 / PhD chapter: pressure and equation of
state for active triangular hard-hexagon, where the equilibrium answer
is exactly known. The key question Kabir pointed to is a small-activity
correction to the exact equilibrium equation of state:
\[
\mu(\rho,\alpha)=\mu_{\rm eq}(\rho)+\alpha\,\mu_1(\rho)+O(\alpha^2),
\]
and similarly for compressibility. This is a later multi-particle
project, not part of the first single-walker code.

---

## Appendix · Today's decision points

Mark with [x] when locked in.

- [x] Include the reversal channel $\gamma_r$ in the equations/code API
      (default: yes, matches Mallikarjun–Pal; set $\gamma_r=0$ when we
      want the simpler adjacent-turn model).
- [ ] Include the 120° rotation channel? **(default: NO for initial
      Phase 1 — direct ref-44 analogue uses only $\{+1, -1, +3\}$. Add
      the $\pm 2$ channel as a Phase 1.5 stretch once the baseline
      result is in hand.)**
- [ ] Normalisation: $v = 1$ vs $v + \gamma = 1$? (recommend: $v = 1$)
- [ ] Lattice convention for paper: algebraic $a = b = 1$ or isotropic
      $b = \sqrt{3} a$? (recommend: stick with the bug-fix convention,
      algebraic for math, isotropic only for visualisation)
- [ ] First-passage geometry: half-plane (cylinder) first, strip second.
- [ ] Edge geometries: zigzag first (matches TCRW square strip
      orientation), armchair second.
- [ ] Edge-current observable: $\sum_\partial |v|^2$ (TCRW Fig 4 style) — yes.
