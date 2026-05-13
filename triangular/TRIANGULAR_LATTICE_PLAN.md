# Triangular Lattice Active Random Walker — research plan after PI meeting

**Meeting date:** 2026-05-04 (Mon)  ·  **Today:** 2026-05-13 (Wed)  ·  **Deadline:** Fri 2026-05-15
**PI:** Kabir Ramola
**Recording:** transcribed from Apple Watch
**Last revised:** 2026-05-13 evening, after two-bug k-grid/sign diagnosis + KMC verification

**Strategic note (from PI lesson learned):** show small results frequently, not one big dump. Each working-day target = one figure or one result PI can absorb in 2 minutes. Friday is the soft deadline for first band-structure result.

---

## 0. Current status as of 2026-05-13

**Main result so far:** Dipanjan's triangular calculation failed for two reasons: the rectangular notebook k-grid is not the correct sheared grid for the triangular \(L\times L\) torus, and the Fourier coefficient \(c_3\) has a local chirality-sign error. The corrected torus grid and corrected sign are now the default in `triangular_jmvr_corrected.py`.

The corrected triangular-torus k-grid is

$$
k_1=\frac{\pi m_1}{aL},
\qquad
k_2=\frac{\pi(2m_2-m_1)}{bL}.
$$

The sign correction is

$$
c_3 = B - 1 - \gamma + 2 i \epsilon \sin(a k_1-b k_2)
\quad\longrightarrow\quad
c_3 = B - 1 - \gamma - 2 i \epsilon \sin(a k_1-b k_2).
$$

**Done**
- 6-state triangular JMVR generator derived and coded.
- Dipanjan's `rtp_tl_2.nb` audited.
- Triangular-torus k-grid/PBC issue diagnosed.
- Appendix-B master equations checked term by term.
- \(c_3\) sign error diagnosed.
- Four-way k-grid/sign comparison made in `forensic_two_bugs.py`.
- Finite real-space generator versus Bloch spectrum check made in `verify_realspace_bloch.py`.
- Independent \(10^8\)-walker KMC run completed.
- Corrected sheared-grid theory matches KMC at the Monte Carlo noise floor:

$$
\mathrm{RMS}_{\rm rectangular+old}=4.448\times 10^{-4},
\qquad
\mathrm{RMS}_{\rm sheared+old}=1.428\times 10^{-4},
\qquad
\mathrm{RMS}_{\rm sheared+corrected}=3.434\times 10^{-6},
\qquad
\sigma_{\rm MC}\approx 4.223\times 10^{-6}.
$$

- Final Fig. 11 replacement made: `fig11_final_hex.png`.
- Evidence figure with residuals made: `fig11_final_hex_evidence.png`.
- PI-facing note made: `PI_NOTE_TWO_BUGS_AND_CHECKS.md`.

**Partly done**
- Band plots exist, including tracked line plots, but they still need to be regenerated from the canonical corrected matrix and a proper triangular BZ path.

**Remaining**
- Corrected \(\Gamma-M-K-\Gamma\) band structure.
- Gap scan versus \(\epsilon\) at fixed \(\gamma\).
- Larger/polished real-space/Bloch finite-torus verification if needed for a paper note.
- Multi-particle/clustering extension.
- Topological invariants.
- Track B triangular TCRW.

---

## 1. What PI actually asked for — two separate tracks

**Critical framing point (corrected 2026-05-11 after review):** PI mentioned both TCRW (the paper you reproduced) and his own JMVR / Confinement-draft work. **These are two different models with two different chirality conventions.** Don't merge them. Do JMVR first (it's his group's model, has an unfinished paper); only attempt TCRW-on-triangular separately, after the JMVR triangular pass works.

### Track A — Triangular JMVR (THIS WEEK; PRIMARY)
Reproduce/finish Kabir-Dipanjan's triangular active random walker. This is the **Confinement 2021 draft §IV.B** (currently empty). Specifically:
- 6 directors (one per NN lattice direction); CTRW dynamics.
- Hopping rate: $\frac{1}{6} \pm \epsilon$ along the director's forward/backward axis, $\frac{1}{6}$ to the other 4 NN. (Translation chirality, controlled by $\epsilon$.)
- Director switching at rate $\gamma/2$ to each of the two adjacent directors (no chirality in rotation).
- 6×6 Bloch matrix $M(\mathbf{k})$.
- **Match equations to the Confinement draft**, verify against Dipanjan's `rtp_tl_2.nb`.
- **Diagnose** the old theory-vs-MC mismatch Kabir mentioned ("I think he implemented the boundary conditions wrong").
- Single-walker first. Then maybe multi-walker non-monotonic clustering in 2D (the draft's §V, also empty).

### Track B — Triangular TCRW (SECONDARY, only if Track A finishes early)
Same Osat–Speck rules you reproduced on square, but on triangular. Distinct from Track A: chirality is in the **rotation** (CW vs CCW with bias $\omega$), not in translation.
- 6 directors, $\omega \in [0, 1]$ rotation chirality parameter, $D_r$ noise rate.
- 6×6 Bloch matrix $P(\mathbf{k})$ with phases on chiral entries.
- Does the topology survive 4-fold → 6-fold? Does the gap still close at $\omega = 1/2$?
- This is the **research-question side** of PI's ask: whether the topological story from TCRW generalizes.

### Task C — Topological invariant question (when Track A or B gives a non-trivial spectrum)
PI explicitly admitted he's not sure of the answer:
> "Yeh aapko padhna padega. Berry curvature, wahan se integral leke aapko milega. As soon as you have this flat bands, there is something known as topological invariant... This is a question I don't honestly know."

Compute the **Chern number** (when TR-symmetry is broken, i.e. $\epsilon \neq 0$ in JMVR or $\omega \neq 1/2$ in TCRW) **OR the 2D Zak phase** (when TR-symmetric). See §3.6 for the recipe.

### Task D — Multi-particle interactions (LATER, NOT NOW)
Connection to:
- **Surajit Chakraborty's PhD work** (just finishing) — multi-particle chiral walkers with friction
- **Odd elasticity (Scheibner–Cohen–Vitelli)** — chirality becomes physically interesting once particles interact
- **Topological mechanics** (Kane & Lubensky line of work) — exploratory, "I'm not 100% sure it will be of interest"

PI was explicit: *"abhi ke liye triangular lattice karo"* — for now, just do triangular. Multi-particle is the next chapter.

---

## 2. Context PI revealed (things I didn't know before)

### 2.1 The JMVR paper is the home reference, not TCRW

**Confirmed citation (web-searched 2026-05-11):**
**Jose S., Mandal D., Barma M., Ramola K.** — *"Active Random Walks in One and Two Dimensions"* — Phys. Rev. E **105**, 064103 (2022) — arXiv:[2202.02995](https://arxiv.org/abs/2202.02995).

> **Note:** Kabir said "Jose, Mandal, *Verma*, Ramola" during the meeting. The actual co-author is **Mustansir Barma** (TIFR Hyderabad senior physicist), not Verma. Easy mishearing; use the correct citation in any write-up.

The JMVR model:
- **Continuous-time random walker (CTRW)** on a 1D or 2D lattice.
- Walker carries an internal "state" (1D: $\pm$; 2D triangular: 6 directors aligned with NN directions).
- **Translation rates** (per direction) — taking 1D as the canonical case from the Confinement draft Eq. (1):
  - Walker in state +: hops in +x at rate $\frac{1}{2} + \epsilon$, in −x at rate $\frac{1}{2} - \epsilon$.
  - Walker in state −: rates swapped.
  - For 6-direction triangular (Confinement draft §IV.B): per-direction baseline rate is $\frac{1}{6}$, with $\pm \epsilon$ bias along the director's forward/backward axis. Sideways hops are unbiased at rate $\frac{1}{6}$.
- **State-switching rate** $\gamma$: this is the orientation-flip rate (1D: $\pm \to \mp$; 2D triangular: $d \to d \pm 1$ mod 6 at rate $\gamma/2$ to each). $\gamma$ is **NOT** the translation baseline — it's the rotation rate.
- **Output**: exact occupation probabilities $P(x, t)$ via inverse Fourier; large-deviation free energy in the continuum limit.

> **What JMVR did NOT do:** topological edge modes or non-Hermitian band topology. That's the TCRW paper's contribution. JMVR is exact stochastic statistical mechanics; TCRW is non-Hermitian topology. Don't conflate them.

> "We were the first people to define this model... in both directions equally, gamma, gamma, gamma/2 forward and backward. When you break it, you go in this direction more and you go in this direction less." — PI was using approximate shorthand; the actual rates are 1/2 ± ε per direction in 1D (or 1/6 ± ε for triangular forward/backward), not γ/2 ± ε.

**Important:** Kabir's group's mental model is JMVR, not TCRW. The TCRW reproduction was a test. **The triangular work should be framed as "JMVR-extended-to-triangular" first, with TCRW as a parallel topology-side calculation.**

### 2.2 They previously tried triangular (Dipanjan, 2019) — never published
PI's draft notes + Mathematica code from **Dipanjan Mandal** (Sep 2019):
> "Yeh paper humne kabhi likha nahi... I think he implemented the boundary conditions wrong. That's what I'm trying to figure out."

So this is **an actual open project in the group**, not just a textbook exercise. There's a real chance you can finish what Dipanjan started and turn it into a paper.

**The Mathematica file `rtp_tl_2.nb` — audited 2026-05-11:**

Parameters at top of file:
- `gamma = 0.01` — rotation rate (γ/2 between adjacent directors)
- `epsilon = 0.15` — **chirality parameter** (forward/backward translation bias)
- `a = 1, b = 1` — lattice constants
- `L = 30` — finite grid size for any real-space work

The 6×6 Bloch matrix `mat[k1, k2]`:
- **Off-diagonal entries = γ/2** at positions $(i, i\pm 1) \mod 6$ — these are director-switching rates (orientational diffusion, symmetric).
- **Diagonal entries** `coeff_i(k1, k2)` for i=1,2,3 and `Conjugate(coeff_i)` for the opposite directors:
$$
\text{coeff}_i = \frac{1}{3}\Big(\cos(2 a k_1) + \cos(a k_1 + b k_2) + \cos(a k_1 - b k_2)\Big) - 1 - \gamma + 2 i \epsilon \sin(\mathbf{k} \cdot \hat{e}_i)
$$
where $\hat{e}_1 = (2a, 0)$, $\hat{e}_2 = (a, b)$, $\hat{e}_3 = (a, -b)$.

Reading it row-by-row (corrected after review 2026-05-11):
- **The walker in any director hops to all 6 NN.** The $(1/3) \cdot (3 \text{ cosines})$ structure is exactly $(1/6) \cdot (6 \text{ cosines})$, with each $\cos(\mathbf{k} \cdot \hat{e}_d) + \cos(-\mathbf{k} \cdot \hat{e}_d) = 2 \cos(\mathbf{k} \cdot \hat{e}_d)$ paired together → 3 distinct cosine terms with prefactor $1/3$. So each director hops symmetrically to all 6 NN with baseline rate $1/6$.
- **Chirality** enters as an imaginary $i \epsilon \sin(\mathbf{k} \cdot \hat{e}_d)$ term — this is the translation bias along director $d$'s axis: $1/6 + \epsilon$ forward, $1/6 - \epsilon$ backward. The other 4 NN remain at rate $1/6$.
- Chirality is in **translation rates**, NOT in rotation rates. (Distinguishes from TCRW.)
- At $\mathbf{k} = 0$ the column sums vanish (probability conservation for a CTRW generator): $\text{coeff}|_{\mathbf{k}=0} = 1 - 1 - \gamma = -\gamma$ balanced by the two off-diagonal $\gamma/2$ entries summing to $\gamma$. **This is a generator condition, not column-stochasticity** — for continuous-time Markov chains, the condition is column sums = 0, not column sums = 1 (column sums = 1 is for discrete-time chains).

**Two possible issues I identified — now revised:**

1. **Lattice parametrization (likely not a bug).** With `a = b = 1`, the NN vectors $(2a, 0), (a, b), (a, -b)$ have lengths $2, \sqrt{2}, \sqrt{2}$ — **not isotropic in Cartesian**. But this is a **deliberate algebraic parametrization** (matches the Confinement draft Fig 1 explicitly), not a geometric mistake. For a true 60°-isotropic triangular metric, use $b = \sqrt{3}\, a$; for the Confinement-draft convention, use $a = b$. The two are related by an affine reparametrization that preserves all topological invariants. **The original "this is a bug" hypothesis was wrong; this is just convention.**

2. **PBC asymmetry in Fourier wavevectors.** Dipanjan defines $k_x[n] = 2\pi n/(2aL)$ and $k_y[n] = 2\pi n/(bL)$ — the x-period is doubled relative to y. **This is the correct BZ for the chosen parametrization** (the lattice vectors $(2a, 0)$ and $(a, b)$ span a parallelogram with that exact reciprocal). This is probably NOT the bug Kabir flagged.

**Update 2026-05-13:** the dominant mismatch is now diagnosed as a sign error in \(c_3\), not as a PBC error. PBC verification is still useful, but now as a robustness check: build the real-space finite-torus generator and confirm its block spectrum matches the corrected Bloch matrix at the allowed discrete momenta.

### 2.3 The PBC trap on triangular
PI drew the geometry explicitly and walked through the trap:

Using Cartesian coordinates $(x, y)$, with primitive vectors $\mathbf{a}_1 = (1, 0)$ and $\mathbf{a}_2 = (1/2, \sqrt{3}/2)$ (so $\mathbf{a}_2$ at 60° from $\mathbf{a}_1$):

- **Periodicity along $\mathbf{a}_1$ (n1-axis):** trivial. $x \to x \mod L$.
- **Periodicity along $\mathbf{a}_2$ (n2-axis):** NOT trivial. If the walker crosses the n2-boundary, you wrap $y$ by the period of the lattice along $\mathbf{a}_2$, AND you also shift $x$ by $-L/2$. Quoting PI from his notes:
$$
\text{along } \hat{n}_1: \quad (x + L) \mod L = \tilde{x}
$$
$$
\text{along } \hat{n}_2: \quad (y + L) \mod L = \tilde{y}, \quad \tilde{x} = x - L/2
$$

This is the **standard skew-PBC for triangular/hexagonal lattices** — well known but easy to get wrong if you naively copy the square-lattice modular arithmetic.

**Simpler implementation:** parametrize sites by integer lattice coordinates $(m, n)$ where site position is $m \mathbf{a}_1 + n \mathbf{a}_2$. Then PBC is trivially $(m, n) \equiv (m + L, n) \equiv (m, n + L)$. Convert to Cartesian only for plotting. **This is what I'd recommend** to avoid the trap entirely.

### 2.4 PI doesn't fully grok the topological argument
> "I am not sure ki kya ho raha hai. But that is our theoretical calculation."

Explicit admission. So when you compute the Berry curvature / Zak phase for triangular, **you'll be teaching him**, not the other way around. This is good — it's a real research contribution rather than running someone else's recipe.

### 2.5 Anshuman is a resource for triangular-lattice PBC
> "You can ask Anshuman. He does triangular lattice all the time."

Probably Anshuman Pandey or another group member. Has correct PBC code already running. **Use him**, especially when implementing the lattice mod arithmetic.

---

## 3. The model on triangular lattice (math setup)

### 3.1 Two distinct chirality conventions — DECIDE FIRST

This is the most important conceptual point. There are **two physically different chiral active walker models on triangular**, and they need to be kept separate:

| Model | Chirality lives in | Parameter | Bloch matrix structure |
|---|---|---|---|
| **JMVR** (Jose-Mandal-Barma-Ramola, Dipanjan's Mathematica) | **Translation rates** (forward ≠ backward) | $\epsilon$ | Off-diagonal γ/2 (symmetric rotations); diagonal has $i \epsilon \sin(\mathbf{k} \cdot \hat{e})$ chirality term |
| **TCRW** (Osat-Speck, the paper you reproduced) | **Rotation rates** (CW ≠ CCW) | $\omega$ | Off-diagonal phases $e^{\pm i \mathbf{k} \cdot \hat{e}}$ on chiral entries; entire entries weighted by ω vs 1−ω |

**These are NOT the same model.** Both have edge currents and non-Hermitian Bloch matrices, but the topology, gap-closure conditions, and edge-mode physics may differ.

**My recommendation:** start with **JMVR-style** (matches Dipanjan's setup, matches Kabir's group's mental model, gives a direct comparison to the old Mathematica). Add TCRW-style as a second variant only if there's time / interest.

### 3.2 JMVR on triangular — the lattice setup

6 nearest-neighbor directions:
$$
\hat{e}_d = (\cos(d\pi/3),\ \sin(d\pi/3)), \quad d = 0, 1, 2, 3, 4, 5
$$
(or, equivalently, Dipanjan's parametrization with $(2a, 0), (a, b), (a, -b)$ and conjugates — see §2.2; either gives equivalent topology).

**Director-switching (rotation) rules — symmetric, no chirality:**
- $d \to (d + 1) \mod 6$ at rate $\gamma/2$
- $d \to (d - 1) \mod 6$ at rate $\gamma/2$
- Total rotation rate per state: $\gamma$. $\gamma$ is the **orientation-switching rate**, nothing else.

**Translation rates — biased along director's axis:**
- Walker in director $d$ hops to all 6 NN.
- Rate to NN in direction $+\hat{e}_d$ (forward): $\tfrac{1}{6} + \epsilon$
- Rate to NN in direction $-\hat{e}_d$ (backward): $\tfrac{1}{6} - \epsilon$
- Rate to other 4 NN (sideways): $\tfrac{1}{6}$ each
- Total translation rate per state: $\tfrac{4}{6} + (\tfrac{1}{6} + \epsilon) + (\tfrac{1}{6} - \epsilon) = 1$
- Achiral case: $\epsilon = 0$ → uniform 1/6 to all 6 NN, walker is a symmetric random walker.

**This is the rate convention to verify against Confinement draft Eq. (1) and §IV.B before writing any code.** The 1D version of the draft uses $\frac{1}{2} \pm \epsilon$ (the 2-direction analog of the 6-direction $\frac{1}{6} \pm \epsilon$).

**Total outgoing rate per state** = (translation total) + (rotation total) = $1 + \gamma$. The generator's diagonal then has $-(1 + \gamma)$ at $\mathbf{k} = 0$ minus the rotation-rate contribution (already accounted for via the off-diagonal $\gamma/2$ terms).

### 3.3 Bloch matrix $M(\mathbf{k})$ — 6×6 generator

The CTRW generator at momentum $\mathbf{k}$:
- **Off-diagonal entries**: rotation rates $\gamma/2$ at positions $(d, d \pm 1) \mod 6$.
- **Diagonal entries**: $-[\text{total outgoing rate}] + [\text{hopping in Fourier}]$, with the chirality $i \epsilon \sin(\mathbf{k} \cdot \hat{e}_d)$ term.

This is exactly Dipanjan's structure but with the lattice geometry corrected (if his $a = b = 1$ was wrong).

### 3.4 Symmetries to expect

**Square (4-fold) → triangular (6-fold)** changes the symmetry group. Expect:
- **C_6 cyclic** symmetry: rotation of all directors by one step ($d \to d+1$) plus rotation of $\mathbf{k}$ by 60° leaves the matrix invariant.
- **Inversion** ($\mathbf{k} \to -\mathbf{k}$): combined with director reflection.
- The square TCRW had a **chiral / sublattice symmetry** giving $\pm \lambda$ eigenvalue pairs. Whether this survives on triangular is **an open question** — it depends on the bipartite-ness of the rotation graph, which is **6-cycle = bipartite!** (6 is even, so the cycle is 2-colorable: alternating directors form two sets of 3). So chiral symmetry may survive. **Worth checking analytically once you have the matrix.**

### 3.5 Gap closure: does it stay at the achiral limit?

In JMVR, the achiral limit is $\epsilon = 0$. At that point, the matrix is real (no imaginary terms in the diagonal), the spectrum becomes manifestly Hermitian, and any gap closure can only be a regular band crossing.

When you turn on $\epsilon > 0$, the imaginary $i \epsilon \sin(...)$ terms break Hermiticity and the spectrum can move into the complex plane. **The first concrete numerical question is: at $\epsilon = 0$ do bands cross, and how does the structure of those crossings change as $\epsilon$ grows?**

This is conceptually different from TCRW (where the gap closes at $\omega = 1/2$ because that's where rotation chirality vanishes). The JMVR "gap closure" condition will involve $\epsilon$ and $\gamma$ in some other combination.

---

### 3.6 Topological invariants — what to compute and how

PI explicitly raised this question at the meeting (*"Berry curvature... topological invariant... I don't honestly know"*). This subsection is the bulletproof reference for when you get to Day 5+.

#### 3.6.1 Three objects, NOT one — keep them straight

| Object | Symbol | Gauge-invariant? | Topological? |
|---|---|---|---|
| Berry **connection** | $\mathbf{A}_n(\mathbf{k}) = i \langle u_n(\mathbf{k}) \vert \nabla_\mathbf{k} \vert u_n(\mathbf{k})\rangle$ | NO | NO |
| Berry **curvature** | $\Omega_n(\mathbf{k}) = \nabla_\mathbf{k} \times \mathbf{A}_n(\mathbf{k})$ (z-component in 2D) | YES | NO — smooth field |
| Topological **invariant** | e.g. Chern number $C_n = \frac{1}{2\pi}\int_{\rm BZ} \Omega_n\, d^2k$ | YES | YES — integer |

**Common confusion (answer for PI):** "Berry curvature" alone is **not** a topological invariant. Its **integral over a closed manifold** (the BZ torus in 2D) IS. The invariant is what's quantized; the curvature is the integrand.

#### 3.6.2 Which invariant does our model have?

Depends on symmetries. Two candidate invariants:

| Invariant | When non-trivial | Computed via |
|---|---|---|
| **Chern number** $C_n$ | Time-reversal broken | FHS plaquette algorithm — see §3.6.4 |
| **2D Zak phase / Wilson loop** $\theta_n$ | Inversion / chiral / sublattice symmetry preserved | Path-ordered Berry connection along 1D loops — see §3.6.5 |

**JMVR on triangular (Track A):**
- At $\epsilon = 0$ (achiral): the generator is real and symmetric → Hermitian. C = 0 trivially. Zak phase may still be non-trivial (similar to SSH-like systems).
- At $\epsilon > 0$ (chiral): the $i \epsilon \sin(\mathbf{k} \cdot \hat{e})$ term makes $M$ non-Hermitian. Under $\mathbf{k} \to -\mathbf{k}$, the matrix satisfies $M(-\mathbf{k}) = M(\mathbf{k})^*$ — a TR$^+$-like anti-unitary symmetry. **Whether this constrains the Chern number is a representation-class question** (Bernard–LeClair / Kawabata classification of non-Hermitian topology). Don't assume C = 0; compute it.

**Track A is primarily about exact stat-mech ($P(\mathbf{r}, t)$, clustering), not topology.** The topology question may not have an interesting answer for JMVR specifically. If you want non-trivial topology, **Track B (TCRW-style chirality on triangular) is the more promising direction** — TCRW on square already has a non-trivial 2D Zak phase, and the triangular extension is where PI's research question lives.

**Compute Chern and Zak phase on both Track A and Track B** when you reach the topology phase of the project. Sequence: Track A finished cleanly → Track B matrix built → compute topology on Track B (TCRW-style) first, since that's where it's most likely interesting.

#### 3.6.3 Non-Hermitian biorthogonal subtlety

JMVR's generator $M(\mathbf{k})$ is non-Hermitian: $M^\dagger \neq M$ because of the $i \epsilon \sin$ term. So $|u_n^R\rangle$ (right eigenvector) and $|u_n^L\rangle$ (left eigenvector, eigenvector of $M^\dagger$) are different.

The correct gauge-invariant Berry connection uses the **biorthogonal pair**:
$$
\mathbf{A}_n(\mathbf{k}) = i \frac{\langle u_n^L(\mathbf{k}) \vert \nabla_\mathbf{k} \vert u_n^R(\mathbf{k}) \rangle}{\langle u_n^L(\mathbf{k}) \vert u_n^R(\mathbf{k}) \rangle}
$$

Normalize so $\langle u_n^L | u_n^R \rangle = 1$ (biorthonormal). All FHS plaquette formulas below use this biorthogonal inner product.

**Pitfall:** at exceptional points (where two eigenvectors coalesce), $\langle u^L | u^R \rangle = 0$ and the connection blows up. Avoid these points by perturbing the BZ grid slightly. Numerically detect them as plaquettes with anomalously large flux.

#### 3.6.4 Fukui-Hatsugai-Suzuki (FHS) algorithm for Chern number

The standard discretized algorithm — quantizes the answer automatically on any grid.

```python
import numpy as np
from scipy.linalg import eig

def fhs_chern(build_Mk, n_band, Nk1=40, Nk2=40, k1_max=2*np.pi, k2_max=2*np.pi):
    """
    Chern number of band `n_band` of a non-Hermitian Bloch matrix M(k).
    Uses biorthogonal link variables.

    build_Mk(k1, k2) -> (N, N) complex matrix (the generator).
    """
    # 1. Set up k-grid on the BZ torus.
    k1s = np.linspace(0, k1_max, Nk1, endpoint=False)
    k2s = np.linspace(0, k2_max, Nk2, endpoint=False)
    N = build_Mk(0.0, 0.0).shape[0]

    # 2. Precompute right and left eigenvectors at every k.
    #    Sort by some consistent labeling (real part of eigenvalue here;
    #    in production use Hungarian band-matching like in TCRW Fig 4b).
    vR = np.zeros((Nk1, Nk2, N, N), dtype=complex)
    vL = np.zeros((Nk1, Nk2, N, N), dtype=complex)
    for i, k1 in enumerate(k1s):
        for j, k2 in enumerate(k2s):
            M = build_Mk(k1, k2)
            wR, VR = eig(M)
            wL, VL = eig(M.conj().T)        # left eigvecs from M†
            # Sort both by real part of eigenvalue; assume bands don't cross
            iR = np.argsort(wR.real); iL = np.argsort(wL.conj().real)
            VR = VR[:, iR]; VL = VL[:, iL]
            # Biorthonormalize: scale so <vL | vR> = 1
            for n in range(N):
                norm = np.vdot(VL[:, n], VR[:, n])
                VR[:, n] /= np.sqrt(norm)
                VL[:, n] /= np.sqrt(norm).conj()
            vR[i, j] = VR
            vL[i, j] = VL

    # 3. Compute plaquette Berry flux via link variables, take imaginary log.
    flux = np.zeros((Nk1, Nk2))
    for i in range(Nk1):
        for j in range(Nk2):
            ip = (i + 1) % Nk1
            jp = (j + 1) % Nk2
            # 4 links around a plaquette
            U1 = np.vdot(vL[i, j, :, n_band],  vR[ip, j, :, n_band])
            U2 = np.vdot(vL[ip, j, :, n_band], vR[ip, jp, :, n_band])
            U3 = np.vdot(vL[ip, jp, :, n_band], vR[i, jp, :, n_band])
            U4 = np.vdot(vL[i, jp, :, n_band], vR[i, j, :, n_band])
            # Normalize each to unit modulus (gauge-invariant phase)
            phase = np.angle(U1 / np.abs(U1) * U2 / np.abs(U2)
                             * U3 / np.abs(U3) * U4 / np.abs(U4))
            flux[i, j] = phase

    # 4. Sum and quantize.
    C = np.sum(flux) / (2 * np.pi)
    return C, flux  # C should be an integer up to grid discretization error
```

**Sanity checks:**
- Sum over the full BZ → quantized integer (verify $|C - \text{round}(C)| < 0.01$ at $40 \times 40$ grid).
- $\sum_n C_n = 0$ over all bands (sum rule from $\det(M)$).
- At $\epsilon = 0$ (TR-symmetric): all $C_n = 0$ by symmetry. **If you get non-zero here, your code is buggy.**

#### 3.6.5 Wilson loop / 2D Zak phase

For the achiral case (or any TR-preserving regime), the relevant invariant is the Wilson loop (Berry phase along a 1D loop in the BZ).

For each $k_1$ on a 1D grid, compute the path-ordered product of overlap matrices along $k_2$:
$$
W(k_1) = \prod_{j=0}^{Nk_2 - 1} \langle u^L(k_1, k_2^{(j+1)}) \vert u^R(k_1, k_2^{(j)}) \rangle
$$
then the Wilson-loop phase $\theta(k_1) = \arg W(k_1)$. A non-trivial Zak phase means $\theta(k_1)$ winds by $2\pi$ as $k_1$ traverses the BZ.

```python
def wilson_phase(build_Mk, n_band, Nk1=40, Nk2=200, k1_max=2*np.pi, k2_max=2*np.pi):
    k1s = np.linspace(0, k1_max, Nk1, endpoint=False)
    k2s = np.linspace(0, k2_max, Nk2, endpoint=True)        # closed loop
    theta = np.zeros(Nk1)
    for i, k1 in enumerate(k1s):
        # Path-ordered product of <u^L(j+1) | u^R(j)> along k2
        W = 1.0 + 0j
        prev_vR = None; prev_vL = None
        for j, k2 in enumerate(k2s):
            M = build_Mk(k1, k2)
            wR, VR = eig(M); wL, VL = eig(M.conj().T)
            iR = np.argsort(wR.real); iL = np.argsort(wL.conj().real)
            vR_here = VR[:, iR[n_band]]; vL_here = VL[:, iL[n_band]]
            # Biorthonormalize
            norm = np.vdot(vL_here, vR_here)
            vR_here /= np.sqrt(norm); vL_here /= np.sqrt(norm).conj()
            if j > 0:
                W *= np.vdot(prev_vL, vR_here)
            prev_vR = vR_here; prev_vL = vL_here
        theta[i] = np.angle(W)
    return k1s, theta   # plot theta vs k1: count windings → Z_2 invariant
```

**Sanity check:** at $\epsilon = 0$ on triangular at the achiral fixed point, the Zak phase should be quantized to $0$ or $\pi$ (Z₂ topology). Non-zero ⇒ non-trivial.

#### 3.6.6 Algorithm pitfalls

1. **Band crossings inside the BZ** → can't track $n$-th band continuously by argsort. Solution: Hungarian band-matching (you already have the code in `tcrw_fig4b_paper.py`). Port it.
2. **Exceptional points** → biorthogonal norm vanishes. Detect by anomalously large plaquette flux (>> π) and avoid by grid perturbation.
3. **Choice of starting band index at degeneracy** → not unique. Solution: compute invariants over **groups of degenerate bands** (multi-band Chern), not individual bands.
4. **Grid resolution** → $40 \times 40$ usually enough for $|C| \leq 3$; bump to $80 \times 80$ if quantization is off by > 0.05.

#### 3.6.7 Reading order — pragmatic, not exhaustive

| Source | Why | Pages |
|---|---|---|
| **Asbóth, Oroszlány, Pályi, "A Short Course on Topological Insulators"** (lecture notes, 2016, free PDF) | Chapter 1 (SSH model) and 2 (Berry phase) — the cleanest intro to the Zak phase calculation | ~40 pages |
| **Vanderbilt, "Berry Phases in Electronic Structure Theory"** (CUP 2018), §3.1–§3.4 | Chern number, BZ as torus, FHS algorithm | ~50 pages |
| **Fukui, Hatsugai, Suzuki, JPSJ 74, 1674 (2005)** | Original FHS paper. 4 pages. Read once for the formula. | 4 pages |
| **Kawabata et al., PRX 9, 041015 (2019)** | Non-Hermitian topology classification. Skim. | skim only |

Skip Hasan & Kane (too long for this purpose). Skip the Bernevig book (too long).

#### 3.6.8 What "bulletproof" means for the topology calculation

You're done with the topology side when:
- [ ] Chern number computed at multiple $(\gamma, \epsilon)$ points on a $40 \times 40$ grid. Quantized to integers with error < 0.01.
- [ ] Wilson loop phase computed at $\epsilon = 0$ and shown to be 0 or π.
- [ ] At least one regime found where Chern or Zak changes — that's the topological phase transition. Plot its location in $(\gamma, \epsilon)$ space.
- [ ] Sanity checks pass: $\sum_n C_n = 0$ over all 6 bands.
- [ ] If transitions found, verify they coincide with gap closures from §3.5.

---

### 3.7 What we reuse from the month of TCRW reproduction work

The TCRW month was not wasted; large parts carry over directly. Here's the honest accounting:

#### Direct architectural reuse (port, don't rewrite)

| TCRW artifact | What it does | Triangular reuse |
|---|---|---|
| `tcrw_fig4b_paper.py::build_Pk` | 4×4 Bloch matrix at $\mathbf{k}$ for square TCRW | Pattern: same skeleton, swap to 6×6 with JMVR / TCRW-on-triangular coefficients. Saves ~half the boilerplate. |
| Hungarian band-matching in `tcrw_fig4b_paper.py::fig4b_band_structure` | Tracks band identity across $\mathbf{k}$ in the presence of crossings | Port verbatim to 6 bands — same algorithm, just $n=6$ instead of $n=4$. |
| `tcrw_fig4{c,d,e}.py::build_obc_matrix` | OBC matrix builder with director, site indexing, state mapping `(i,j,d) → s` | Same pattern. Triangular needs different lattice connectivity (6 NN), but the indexing logic ports. |
| `tcrw_fig4{f,g}.py::compute_edge_weight` | $\sum_{\partial}\|v\|^2$ over right eigenvectors | Reuse verbatim. The math is geometry-independent. |
| `_fig3_plot.py::_per_y_quiver` | Quiver visualization on a 1D edge | Port to triangular edge geometry. |
| Bit-identity test pattern (we used 6 test points for TCRW) | `max\|W_user − W_authors\| = 0` cross-check | **Critical for Track A.** We'll bit-check our triangular JMVR matrix against Dipanjan's `coeff_i` formulas the same way. |
| FHS Berry curvature algorithm (§3.6.4 of this plan) | Topological invariant on any 2D BZ | Universal — works on triangular BZ identically. |

#### What we do NOT carry over

- **`TRW._original_code_by_paperauthors.py`** — Osat-Speck's square-only Python. Useful as a reference for TCRW conventions but cannot be reused for triangular (their code has 4 hardcoded directors). Don't try to monkey-patch it.
- **TCRW Fortran (`tcrw_step.f90` etc.)** — TCRW is discrete-time; JMVR is continuous-time. Different walker rules. **No Fortran needed this week** anyway (eigensolve only, see §5).
- **Specific TCRW physics observables** (J_Dr / J_ω decomposition, θ angles) — JMVR has a different chirality convention; observables will be different. The decomposition machinery doesn't transfer.

#### "Help from authors" — what's available

- **Osat-Speck (TCRW authors):** they sent us `TRW.py` for square. **They did not write a triangular version.** No useful additional code to ask for; do not email again unless we hit a specific TCRW-square issue.
- **Dipanjan Mandal's `rtp_tl_2.nb`** — this IS our closest equivalent to "authors' triangular code". Already audited in §2.2. The `coeff_i` expressions are our ground-truth reference.
- **Kabir (PI)** — direct collaborator. He's an author of both JMVR and Confinement-2021 draft. Ask him in person, not by email.
- **Stephy Jose** (JMVR co-author, currently with Löwen) — copy-cc later if/when results are ready; she might be a co-author on the Track A paper given her involvement with the original draft.

#### Practical implication for week-1 plan

This reuse cuts the Day 1–2 work nearly in half. Instead of writing a 6×6 Bloch matrix from scratch, we can:

```python
# Pseudocode for Day 1 - structurally identical to tcrw_fig4b_paper.build_Pk
def build_Mk_jmvr_triangular(gamma, epsilon, k1, k2, a=1.0, b=1.0):
    """6x6 generator for triangular JMVR — matches Dipanjan's coeff_i."""
    coeff = [None] * 3
    for i, (vx, vy) in enumerate([(2*a, 0), (a, b), (a, -b)]):
        kd = vx * k1 + vy * k2
        bulk = (1/3) * (np.cos(2*a*k1) + np.cos(a*k1 + b*k2) + np.cos(a*k1 - b*k2))
        coeff[i] = bulk - 1 - gamma + 2j * epsilon * np.sin(kd)
    # Build 6x6 with γ/2 off-diagonals at (i, i±1 mod 6)
    M = np.zeros((6, 6), dtype=complex)
    M[0, 0] = coeff[0]; M[1, 1] = coeff[1]; M[2, 2] = coeff[2]
    M[3, 3] = np.conj(coeff[0]); M[4, 4] = np.conj(coeff[1]); M[5, 5] = np.conj(coeff[2])
    for i in range(6):
        M[i, (i + 1) % 6] += gamma / 2
        M[i, (i - 1) % 6] += gamma / 2
    return M
```

That's ~15 lines. Plus the diagonalization and plot infrastructure from `tcrw_fig4b_paper.py`. **Day 1 is realistically 2 hours not 3.**

---

## 4. Papers — reading order

### 4.1 PRIORITY 1: JMVR 2022 — Jose, Mandal, Barma, Ramola — arXiv:[2202.02995](https://arxiv.org/abs/2202.02995)
Confirmed citation: *"Active Random Walks in One and Two Dimensions"*, Phys. Rev. E **105**, 064103 (2022).
This is the home model. Read §II (model definition) and §III (2D square-lattice results) tonight or tomorrow morning. **Don't proceed without this read.**

### 4.2 PRIORITY 2: `rtp_tl_2.nb` (Dipanjan's Mathematica) — AUDITED 2026-05-11
Already audited (see §2.2). Structure is clear. **Two suspected bug spots: lattice anisotropy (a=b=1 not isotropic) and PBC asymmetry (kx period ≠ ky period).** Use this as a reference for the matrix structure, not as ground truth.

### 4.3 PRIORITY 3: TCRW paper (you already know it)
Re-read paper §IV (topology) — relevant when you compute Berry curvature on triangular.
**Remember: TCRW uses rotation-bias chirality, JMVR uses translation-bias chirality. Don't conflate them.**

### 4.4 BACKGROUND: Lubensky-Kane-Mao-Souslov-Sun 2015 review
*"Phonons and elasticity in critically coordinated lattices"* — Rep. Prog. Phys. **78**, 073901 (2015).
**This is NOT an active-matter review** (corrected from earlier draft of this plan). It's the foundational review on **topological mechanics and Maxwell rigidity** — isostatic lattices, floppy modes, topological zero modes in mechanical metamaterials. The connection to your work is via topology, not via active matter directly. Useful if/when you read about topological mechanics (Kane–Lubensky line); skip for the JMVR triangular calculation itself.

### 4.5 LATER: Dolai_single_file_active_2020
*"Universal scaling in active single-file dynamics"* — multi-particle, single-file. Reference for the multi-particle direction PI mentioned. Not now.

### 4.6 LATER: Confinement_enhanced_clustering
Multi-particle paper. Same — defer.

### 4.7 LATER: Odd elasticity (Scheibner-Cohen-Vitelli)
Not sent but mentioned. Find it: "Scheibner Cohen Vitelli odd elasticity Nature Physics 2020". For when you do multi-particle.

---

## 5. The work — driven by physics, not by deadlines

**Revised 2026-05-11 evening.** Previous version of this section had a day-by-day schedule. Dropped. Research isn't a sprint; it's a sequence of physics questions. Each step is "what does the physics ask next", not "what's tomorrow's deliverable".

### Where we are right now

- Corrected 6×6 generator $M(\mathbf{k})$ is now centralized in `triangular_jmvr_corrected.py`.
- `build_Mk_corrected(...)` is the default scientific matrix; `build_Mk_dipanjan(...)` is kept only to reproduce the old sign bug.
- The root folder is now organized. Active workflows are documented in `README.md`; old bug-hunt scripts and superseded plots live in `legacy_debug/`.
- The active starter script `triangular_active_walker.py` points to the corrected matrix. The older day-1 starter `triangular_jmvr_v1.py` is archived in `legacy_debug/scripts/`.
- Eigenvalues at $\mathbf{k}=0$ (for $\gamma=0.01$): $\{0,\ -\gamma/2,\ -\gamma/2,\ -3\gamma/2,\ -3\gamma/2,\ -2\gamma\}$ — Perron eigenvalue $\lambda_0 = 0$ confirms a stationary state exists.
- At $\epsilon=0$ the spectrum is real; at $\epsilon=0.15$ it picks up imaginary parts.
- Sixfold spectral symmetry catches the bug: corrected matrix passes to numerical precision; Dipanjan's old sign fails.
- Four-way comparison confirms neither correction alone is enough; sheared grid plus corrected \(c_3\) matches \(10^8\)-walker KMC at the MC noise floor.
- Fig. 11 exact theory now matches \(10^8\)-walker KMC at the MC noise floor.
- Band plots exist, but the next publishable band result must be regenerated from the corrected canonical matrix and a proper BZ path.

### What the physics demands next, in logical order

Move to the next step only when the previous one's answer is solid. No time budget on any of them.

**Step 1 — Make the corrected band structure actually visible.** Use `triangular_jmvr_corrected.py` as the only matrix source. Hungarian band-matching from `tcrw_fig4b_paper.py::fig4b_band_structure` ports directly to 6 bands. Until this is done, we cannot answer "is there a gap?" cleanly.

**Step 2 — Compute the BZ path properly.** Use the reciprocal lattice $\mathbf{b}_1, \mathbf{b}_2$ from $\mathbf{a}_1 = (2, 0), \mathbf{a}_2 = (1, 1)$ (Dipanjan's convention) and find the exact $\Gamma, M, K$. Until this is right, features at the corners might be sampling artefacts and aren't trustworthy.

**Step 3 — Look at the cleaned bands. Is there a gap? At what $(\gamma, \epsilon)$?** This is the first real physics result. Scan $\epsilon$ at fixed $\gamma$, plot minimum gap. Three outcomes are equally interesting:
- Gap closes at some $\epsilon^* > 0$ → phase transition. Map its location.
- Gap stays open for all $\epsilon \in [0, 1/6]$ → JMVR triangular is topologically trivial in this sense; topology question moves to Track B (TCRW-style).
- Gap closes at one or more isolated $\mathbf{k}$ points only (Dirac-like) → richer story, look for monopoles in Berry curvature.

**Step 4 — Compute $P(\mathbf{r}, t)$ for a single walker.** This is the missing §IV.B of the Confinement draft. Inverse Fourier of $\tilde{P}(\mathbf{k}, t) = \sum_n c_n \mathbf{v}_n^R e^{\lambda_n t}$ over a finite BZ grid. Plot at late times. Look for the chirality signature (petal-like asymmetry similar to draft's Fig 7 for square).

**Step 5 — Multi-walker $P(\mathbf{r}, t)$ with step initial condition.** The §V hole of the draft. The 1D analog (draft Fig 2d) shows non-monotonic clustering in number of walkers $W$. Does that survive on triangular? Answer that and §V is filled.

**Step 6 — Topological invariants (Berry curvature, Wilson loops).** Only meaningful once Step 3 gives a non-trivial gap structure. FHS algorithm from §3.6.4. If JMVR turns out trivial, this step pivots to Track B (TCRW-style chirality on triangular).

**Step 7 — Real-space/Bloch verification.** The main mismatch is already explained by the \(c_3\) sign. Still build the real-space transition matrix on a finite triangular lattice and compare its discrete spectrum to the corrected Bloch eigenvalues at quantized $\mathbf{k}$-points. This future-proofs the PBC implementation and is also the foundation for edge-current calculations in OBC.

### What's deliberately NOT in this sequence

- **New expensive Monte Carlo.** The existing \(10^8\)-walker KMC already verifies the sign fix. Do not rerun it unless a new result specifically needs it.
- **TCRW-style chirality on triangular** (Track B). Only after Steps 1–5 give a clean JMVR picture. If Step 6 shows JMVR is topologically trivial, that's when Track B becomes the active path.
- **Multi-particle interactions, odd elasticity, topological mechanics.** PI's longer-term direction. Reading material for after the JMVR triangular calculation is done.

### Sanity instinct — what to do at each step before moving on

Look at the result. Spend a few minutes. Ask:
- Does it make physical sense at the limits ($\epsilon = 0$ should give isotropic spreading; $\gamma \to 0$ should approach ballistic motion in each director)?
- Does it match what you'd guess from the square version?
- Does it match Dipanjan's Mathematica where it can?

If something looks off, stop and investigate. **30 minutes of "wait, this doesn't add up" investigation is worth 3 hours of building on a wrong foundation.**

### When to talk to PI

When there's something concrete to say, not on a calendar. Likely natural moments:
- After Step 3: "Sir, here's the band structure on triangular. The gap [closes at ε=? / stays open]. Compared to square TCRW where it closes at ω=1/2, this is [the same / different]."
- After Step 4: "Sir, here's the §IV.B figure for the Confinement draft — single walker on triangular. The chirality shows up as ___."
- After Step 6: "Sir, the Chern number / Wilson phase is ___ in the regime ε > ε*."
- After Step 7: "Sir, verified the real-space matrix and Bloch matrix agree on the triangular torus."

PI doesn't want weekly progress updates. He wants research answers. Don't artificially fragment them.

---

## 6. Open questions

### Resolved 2026-05-11

1. ~~JMVR paper citation~~ — **RESOLVED.** Jose, Mandal, **Barma**, Ramola, Phys. Rev. E **105**, 064103 (2022), arXiv:[2202.02995](https://arxiv.org/abs/2202.02995). Note: it's **Barma**, not Verma, as Kabir said.
2. ~~Mathematica access~~ — **RESOLVED.** Don't need Mathematica installed; the .nb file structure was extracted via sandbox `grep` (it's plain XML-ish text). Structure is documented in §2.2.
3. ~~Timeline~~ — **RESOLVED.** Soft Friday deadline (2026-05-15) for first band-structure plot.
4. ~~Edge geometry~~ — **PROVISIONALLY RESOLVED.** Start with rhombus (simpler). Hexagon if time permits.

### Resolved 2026-05-13

5. ~~Where is the old theory-vs-KMC mismatch?~~ — **RESOLVED.** It is a two-part issue: Dipanjan's rectangular Fourier grid is not the sheared triangular-torus grid, and the draft/Mathematica \(c_3\) coefficient has the wrong chirality sign. Neither correction alone reaches KMC; both corrections together agree with \(10^8\)-walker KMC at the MC noise floor. See `PI_NOTE_TWO_BUGS_AND_CHECKS.md` and `forensic_two_bugs.py`.
6. ~~Does the corrected matrix pass a stronger check than probability conservation?~~ — **RESOLVED.** The corrected matrix passes the sixfold spectral-symmetry check; Dipanjan's old sign fails it.

### Still open (ask Kabir when convenient, or proceed with default)

7. **"γ/2_A vs γ/2_B" in PI's handwritten notes** — best guess: PI was using approximate shorthand for the forward/backward translation rates ($\frac{1}{6} + \epsilon$ vs $\frac{1}{6} - \epsilon$ in the triangular case, or $\frac{1}{2} \pm \epsilon$ in 1D). The "γ/2" was probably mis-named — the actual rotation rate γ doesn't enter the translation rates. Default: use 1/6 ± ε for triangular forward/backward bias.
8. **Lattice parametrization** — Dipanjan used $a = b = 1$ which gives Cartesian NN lengths $(2, \sqrt{2}, \sqrt{2})$. **Corrected (2026-05-11): this is a deliberate parametrization choice, not a bug.** The Confinement draft Fig 1 uses these exact coordinates. Topologically equivalent to true 60°-isotropic ($\mathbf{a}_1 = (1, 0)$, $\mathbf{a}_2 = (1/2, \sqrt{3}/2)$). Either works; pick one and stick with it. Recommend: **use Dipanjan's $(2a, 0), (a, b), (a, -b)$ convention** for algebra and the isotropic triangular display only for plotting.
9. **Anshuman's contact** — Anshuman Pandey? Anshuman Kumar? Don't know. **Plan:** ask Kabir or the group's mailing list for contact. Useful for PBC code reference but not blocking.
10. **Which chirality convention to do FIRST** — JMVR-style (translation bias) is now the primary track (Track A). TCRW-style (rotation bias) is Track B, secondary.
11. **PBC verification** — no longer the leading bug hypothesis, but still important. Build a finite real-space generator on a small triangular torus and compare its block spectrum to the corrected Bloch spectrum at discrete momenta.

---

## 7. The publishable angle (if it works)

**Corrected framing (2026-05-11 review):** JMVR is **exact stochastic stat-mech** (occupation probabilities, large-deviation function), NOT topology. TCRW is the topological-edge-current paper. Don't conflate them.

Two distinct possible papers come out of this work — pursue whichever Track succeeds:

### Track A paper — Complete the Confinement 2021 draft
**Title:** *Active random walks on a triangular lattice* (or extend the Confinement 2021 working title).

**Abstract sketch:**
1. Mandal-Barma-Ramola 2021 (unpublished draft) studied confinement-enhanced clustering for active random walkers in 1D and on the 2D square lattice. The 2D triangular case was sketched but never completed.
2. We correct a sign error in the triangular Bloch coefficient \(c_3\) and verify the corrected exact solution against \(10^8\)-walker KMC.
3. We fill in §IV.B (single walker, triangular) and §V (multi-walker clustering, triangular) of that draft.
4. Closed-form (or numerically computed) occupation distribution $P(\mathbf{r}, t)$ on triangular for chirality $\epsilon \in [0, 1/6]$.
5. Verify the non-monotonic clustering signature (1D Fig 2(d) of the draft) survives 6-fold geometry.
6. Connect to JMVR 2022 (the published square-lattice version) and discuss the role of lattice symmetry.

This is the **natural completion** of an existing draft. Authors would likely include Mandal, Barma, Ramola, and you. **Probability of becoming a paper: high**, since it finishes work the group already started.

### Track B paper — Triangular topology
**Title:** *Topological edge currents of a chiral random walker on a triangular lattice* (only if you pursue TCRW-style chirality on triangular and find non-trivial topology).

**Abstract sketch:**
1. The Osat-Speck 2026 TCRW model has non-trivial 2D-Zak-phase topology on the square lattice; gap closes at $\omega = 1/2$.
2. We extend to 6-fold triangular geometry with rotation-bias chirality.
3. Result: [the topology survives / shifts / a new invariant appears].
4. Edge modes confirmed in OBC.

Probability of becoming a paper: **conditional** on the topology turning out non-trivial in a different way than square. If identical, it's just a robustness note.

**Track A is the safer bet.** It's an in-group continuation of established work. Track B is a riskier research bet but a more original contribution if the topology genuinely differs.

The Dipanjan-2019 work being incomplete + the JMVR paper being the home stat-mech reference + the Confinement 2021 draft being unfinished + PI's own admission that he doesn't know the topology question = **this is genuinely open territory, not a textbook exercise.**

---

## 8. Risks

| Risk | Mitigation |
|---|---|
| PBC implementation gets wrong | Use lattice coordinates everywhere; verify by checking the achiral limit gives the JMVR result | 
| 6×6 matrix has degeneracies that complicate band tracking | Use Hungarian band-matching algo (already implemented in TCRW Fig 4(b)) — port the same logic |
| Berry curvature is numerically ill-defined at band crossings | Track the lowest-band invariant only; use Fukui-Hatsugai discretization |
| Triangular case gives identical physics to square — less interesting | Pivot to multi-particle (Task 3) earlier than planned |
| Mathematica file has more bugs than expected — uses for reference only | Don't trust Dipanjan's code blindly; verify everything from scratch in Python |

---

## 9. Bottom line

This is **a real research project**, not a reproduction exercise. PI's framing is:
- *Test*: can you generalize a model competently?
- *Research*: does the topology survive 4-fold → 6-fold?
- *Group strategic*: pick up a 5-year-old abandoned thread and bring it to closure.

The triangular extension is the right next step from the TCRW reproduction. If it goes well, this can be a paper with you as first author within ~6 months. The multi-particle direction (Task 3) is a 1–2 year project after this.
