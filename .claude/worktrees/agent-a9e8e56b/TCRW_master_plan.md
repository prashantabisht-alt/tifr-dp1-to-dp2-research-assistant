# Master Plan: Working Through "Topological Chiral Random Walker" (Speck et al.)

## Why this paper matters for you

Your DP1 built three pillars: (1) discrete random walks with internal structure (Ch. 2), (2) active particles with rotational dynamics in 2D (Ch. 3, ABP), (3) higher-order stochastic dynamics (Ch. 4, jerky particles). The jcABP paper extended pillar 3 into chirality in the continuum. This paper closes the loop — it puts chirality onto a lattice and discovers that topology emerges. It connects your random-walk skills to your active-matter skills through a genuinely new physics: non-Hermitian topological phases in stochastic systems.

Kabir gave you this because it's tractable (the simulation is a discrete Markov chain — no SDEs, no numerical instabilities), theoretically rich (non-Hermitian band theory, Zak phase, bulk-boundary correspondence), and directly extendable (add jerk? add memory? go 3D? change lattice geometry?).

---

## The six phases

### PHASE 1: The bare simulation engine (Days 1–2)
**Goal:** Build a single-walker TCRW simulator from scratch. Understand the dynamics by watching it.

**The model rules (be precise):**
- 2D square lattice, walker at (i, j) with director d ∈ {↑, ↓, →, ←}
- At each discrete time step, choose:
  - With probability D_r: **noise step** — walker stays put, director rotates CCW with prob ω, CW with prob (1−ω)
  - With probability (1−D_r): **chiral step** — walker translates one step in direction d, THEN director rotates CW with prob ω, CCW with prob (1−ω) [OPPOSITE chirality to noise]
- Boundary conditions: PBC or OBC (hard walls — if translation would exit grid, it's blocked, and the rotation part of the chiral step also doesn't happen)

**Figures to reproduce:**
- Fig 1(b): Sample trajectories for different (ω, D_r) — PBC, color by time
- Fig 1(c): MSD vs t for different ω at D_r = 10⁻³ — confirm linear growth (normal diffusion)
- Fig 1(d): Diffusion coefficient D vs ω for different D_r — confirm D decreases linearly with chirality

**Code:** Python (NumPy vectorized over trajectories). This is a discrete chain — no Fortran needed for speed. A single trajectory is ~10 lines of logic.

**Sanity checks:**
- ω = 0.5: ordinary random walk, D should match known value
- ω = 0, D_r = 0: deterministic closed orbit of period 4 on PBC
- ω = 1, D_r = 1: pure spinor (localized, no translation)
- MSD should be exactly linear at long times for all (ω, D_r) with D_r > 0

**What you learn:** The basic dynamics. How chirality suppresses diffusion. Why ω = 0.5 is special (achiral point).

---

### PHASE 2: Edge localization and probability distributions — OBC (Days 3–4)
**Goal:** See the edge states with your own eyes. Measure them quantitatively.

**What to compute:**
- P(x,y): steady-state probability of finding walker at each lattice point. Run one long trajectory (T = 10^10 steps or so) on an L×L grid with OBC, histogram the visits.
- Compare ω = 0 (achiral), ω = 0.5 (achiral), ω = 1 (fully chiral) — should see edge localization for ω ≠ 0.5
- P_edge / P_bulk ratio vs D_r (Fig 3a) — should decrease (less edge localization) as D_r increases
- P_edge / P_bulk vs ω (Fig 3f) — should be independent of ω (!)

**Figures to reproduce:**
- Fig 2(a): P(X,Y) heatmap for ω = 0.5 (chiral walker, OBC)
- Fig 2(f): P(X,Y) for ω = 0 (achiral walker, OBC)
- Fig 3(a): P_edge/P_bulk vs D_r for different L (confirm size-independence)
- Fig 3(f): P_edge/P_bulk vs ω (confirm ω-independence)

**Key insight to internalize:** Edge localization (measured by P_edge/P_bulk) depends on D_r but NOT on ω. Chirality affects the current, not the localization. This is non-obvious and important — it means you can't distinguish the topological phase from the trivial phase just by looking at P(x,y).

**Code notes:**
- For P(x,y): single very long trajectory (ergodic), or equivalently solve the transition matrix steady state exactly via linear algebra (eigenvalue λ=1 of the full L²×4 transition matrix). Paper says both agree.
- Edge sites = boundary layer (sites adjacent to walls). Bulk = everything else.

---

### PHASE 3: Edge currents and their decomposition (Days 5–7)
**Goal:** Measure the chiral edge current. Understand its topological origin.

**What to compute:**
- Total current J(x,y) = net outflow vector from each site. This is the probability flux.
- Decompose: J = J_Dr + J_ω
  - J_Dr: current from steps immediately after a noise step (rotational diffusion event followed by a chiral translation)
  - J_ω: current from all other chiral steps
- Current magnitude and direction along the left wall (Figs 3c-e, 3h-j)
- Current angle θ_{J_Dr} and θ_{J_ω} vs D_r and vs ω

**Figures to reproduce:**
- Fig 2(c)-(e): J, J_ω, J_Dr vector fields for ω = 0.5 (chiral)
- Fig 2(h)-(j): same for ω = 0 (achiral) — should see NO chiral edge current
- Fig 2(k)-(o): with defects — edge current robust, runs along internal boundaries too
- Fig 3(b): |J_Dr|/|J_ω| vs D_r
- Fig 3(c)-(e): current orientation along left edge vs D_r
- Fig 3(g)-(j): same vs ω

**The physics:** The edge current runs OPPOSITE to the chirality of the walker — like a quantum Hall skipping orbit. The chiral move pushes the walker into the wall → it can't translate → waits for a noise rotation → noise rotates it with opposite chirality → this moves it along the wall. The mechanism is explained in Fig 7(c) and the Methods section.

**Key result:** For ω = 0.5 (achiral), J_Dr = 0 exactly. This is because the achiral walker has no preferred rotation direction, so noise steps don't generate net directed motion along the edge.

**Implementation detail:** To decompose J into J_Dr and J_ω, you need to track whether the previous step was a noise step or a chiral step. If the current step is a chiral translation AND the previous step was a noise step, that displacement contributes to J_Dr. Otherwise it goes to J_ω.

---

### PHASE 4: Transition matrix spectrum — the heart of the topology (Days 8–11)
**Goal:** Compute the spectrum of the non-Hermitian transition matrix. See the gap closing. Understand why edge states are topologically protected.

**4A: PBC spectrum in Fourier space**

The bulk transition matrix in Fourier space is a 4×4 matrix P(k) (Eq. 1 of the paper):
```
P(k) = | 0            R1+C1*e^{ikx}  0              R2+C2*e^{-ikx} |
       | R2+C2*e^{iky} 0              R1+C1*e^{-iky}  0              |
       | 0            R2+C2*e^{ikx}  0              R1+C1*e^{-ikx} |
       | R1+C1*e^{iky} 0              R2+C2*e^{-iky}  0              |
```
where C1 = (1−ω)(1−D_r), C2 = ω(1−D_r), R1 = ωD_r, R2 = (1−ω)D_r.

- Diagonalize P(k) for k = (k_x, k_y) on a grid in [-π, π]²
- Plot Re(λ) vs (k_x, k_y) — band structure (Fig 4b)
- Confirm gap closing at ω = 0.5 (the topological phase transition)

**Figures to reproduce:**
- Fig 4(b): Re(λ) band structure for ω = 0.35, 0.5, 0.65 — see the gap close and reopen

**4B: OBC spectrum in real space**

Build the full transition matrix for an L×L grid with OBC. This is a (4L²) × (4L²) sparse matrix (4 internal states per site). Diagonalize it numerically.
- Plot eigenvalues in the complex plane
- Color by edge-localization weight of the eigenvector

**Figures to reproduce:**
- Fig 4(c): spectrum in complex plane for small L (L=2), colored by edge localization, decomposed into Re and Im parts vs (D_r or ω)
- Fig 4(f)-(g): spectrum in complex plane for L=10, varying D_r and ω. In the topological phase (D_r < 0.5), the spectrum spreads in the complex plane and an oval-shaped edge band appears.
- Fig 10: PBC vs OBC vs 1D edge model spectra overlaid — the oval shape in OBC comes from the edge

**The physics:** In PBC, all eigenvalues are real (due to sublattice symmetry). In OBC, eigenvalues become complex — the imaginary parts encode the edge current. The gap closing at ω = 0.5 in PBC signals the bulk topological transition. In OBC, the gapless edge mode appears inside the bulk gap for D_r < 0.5 — this is the bulk-boundary correspondence.

**Code:** numpy.linalg.eig for the 4×4 Fourier-space matrix (swept over k). scipy.sparse.linalg.eigs for the large OBC matrix. For L=10, the matrix is 400×400 — still exact diag is fine.

---

### PHASE 5: Topological invariant — the 2D Zak phase (Days 12–14)
**Goal:** Compute the vectorized 2D Zak phase Φ = (Φ_x, Φ_y) and reproduce the phase diagram.

**The math:**
The system is non-Hermitian, so you need BIORTHOGONAL eigenvectors: left eigenvectors ⟨φ_n| and right eigenvectors |ψ_n⟩ of P(k).

Berry connection for band n:
  A_n(k_x, k_y) = i ⟨φ_n | ∂_k | ψ_n⟩

Vectorized 2D Zak phase:
  Φ_n = (1/2π) ∫_BZ A_n(k_x, k_y) dk_x dk_y

In practice, discretize k-space and use the Wilson-loop / discrete Berry phase formula to avoid gauge issues.

**Figures to reproduce:**
- Fig 8(a): Φ_x in the (ω, D_r) plane — should be 0 for D_r > 0.5 and π for D_r < 0.5 (away from ω = 0.5)
- Fig 8(b): Φ_y — same phase boundary
- Fig 8(c): gapless boundary overlaid — confirms Φ transitions exactly where the gap closes

**What this proves:** The edge current isn't an accident of the dynamics — it's a topological invariant of the bulk. The Zak phase Φ = (π, π) in the non-trivial phase is quantized and protected by the symmetries of P(k). Small perturbations can't destroy it.

**Code:** This requires careful numerics. The key steps are:
1. Compute left and right eigenvectors of P(k) at each k-point (use scipy.linalg.eig with left=True)
2. Normalize biorthogonally: ⟨φ_n|ψ_m⟩ = δ_{nm}
3. Compute Berry phase using discrete Wilson loop: product of ⟨φ_n(k_i)|ψ_n(k_{i+1})⟩ around a closed path
4. Sum over one direction to get Φ_x, Φ_y

**Warning:** Band labeling (which eigenvalue is which band) can be tricky when bands cross. Use continuity in k-space or the Wilson loop method which avoids this issue.

---

### PHASE 6: The 1D effective edge model + applications (Days 15–17)
**Goal:** Understand the analytical edge model. Reproduce maze-solving and self-assembly results.

**6A: 1D edge model**

The effective 1D model (Fig 9) reduces the edge dynamics to a 2-state Markov chain (←, →) with absorption into bulk. The rate matrix A (Eq. 3) gives analytical eigenvalues (Eq. 4). This explains the timescale separation: for D_r < 0.2, the edge band in the complex plane splits into two loops (Fig 10).

- Compute eigenvalues of A(k) analytically and plot them in the complex plane
- Overlay with the full OBC spectrum — the 1D model captures the edge band perfectly

**Figures to reproduce:**
- Fig 7(e): P(τ_edge) distribution, confirm exponential decay with D_r²
- Fig 10: 1D edge model eigenvalues overlaid on PBC and OBC spectra

**6B: Maze solving**

Generate random mazes using Prim's algorithm on an L×L grid. Run TCRW walkers with different (ω, D_r). Measure mean first-passage time (MFPT) τ_M to reach the exit.

**Figures to reproduce:**
- Fig 5(a): sample trajectories through maze for ω = 0, 0.5, 1
- Fig 5(b): MFPT vs starting position and director
- Fig 5(c): MFPT vs ω for different D_r (and rescaled by D_r)
- Fig 5(d): MFPT vs L — confirm τ_M ~ L² for chiral, L³ for achiral

**6C: Self-assembly**

This is the most complex application. Patchy tiles with Hebbian interaction rules, growing a target structure from a seed. The key result is that chiral tiles find the boundary of the growing seed and walk along it, dramatically reducing assembly time.

- Fig 6(d): τ_SA vs D_r for different ω
- Fig 6(e): τ_SA vs ω for different L and D_r

---

## Recommended code structure

```
tcrw_project/
├── tcrw_core.py          # Phase 1: Walker class, simulation engine
├── tcrw_pbc.py            # Phase 1: PBC simulations (MSD, diffusion coeff)
├── tcrw_obc.py            # Phase 2: OBC simulations (P(x,y), edge localization)
├── tcrw_currents.py       # Phase 3: Edge current measurement and decomposition
├── tcrw_spectrum_pbc.py   # Phase 4A: Fourier-space band structure
├── tcrw_spectrum_obc.py   # Phase 4B: Real-space OBC spectrum
├── tcrw_zak_phase.py      # Phase 5: Topological invariant computation
├── tcrw_1d_edge.py        # Phase 6A: Analytical 1D edge model
├── tcrw_maze.py           # Phase 6B: Maze generation and solving
├── tcrw_selfassembly.py   # Phase 6C: Self-assembly simulation
└── tcrw_sim.f90           # (Optional) Fortran version of core MC for speed
```

## Language choice

Python throughout. Here's why:
- The simulation is a discrete Markov chain (~10 lines of logic per step). Even in pure Python with NumPy vectorization over trajectories, 10^8 steps takes seconds.
- The spectral analysis (matrix diag, eigenvectors, Berry phase) requires scipy/numpy linear algebra — no clean Fortran equivalent without LAPACK wrappers.
- Plotting is integral (heatmaps, complex plane, band structures, trajectory coloring) — matplotlib is essential.
- If the MC becomes a bottleneck for very long runs (10^10 steps for P(x,y) convergence), we can write a Fortran version of just the walker loop and call it from Python via f2py, or write a standalone Fortran MC. But start in Python.

## What makes each phase non-trivial

| Phase | Easy part | Hard part |
|-------|-----------|-----------|
| 1 | Writing the walker | Getting boundary conditions exactly right for OBC; vectorizing over trajectories |
| 2 | Running long trajectories | Convergence of P(x,y) to steady state; alternatively, building and solving the exact transition matrix |
| 3 | Measuring total current | Correctly decomposing J into J_Dr and J_ω (need to track previous step type) |
| 4 | Diagonalizing P(k) | Understanding why PBC eigenvalues are real but OBC are complex; interpreting the gap |
| 5 | Computing eigenvectors | Biorthogonal normalization for non-Hermitian matrices; gauge-invariant Berry phase via Wilson loop |
| 6 | Maze generation | Getting MFPT statistics converged; self-assembly interaction rules |

## Timeline reality check

- Phases 1–3 (simulation + edge currents): ~1 week. This is your comfort zone — MC simulation and data analysis. You've done harder things in DP1.
- Phase 4 (spectrum): ~3–4 days. Matrix diag is straightforward, but interpreting the results requires understanding non-Hermitian spectral theory.
- Phase 5 (Zak phase): ~3 days. This is the intellectually hardest part. The biorthogonal Berry phase is a real calculation. But the code is compact once you understand the math.
- Phase 6 (applications): ~3 days. Maze solving is fun and fast. Self-assembly is more involved but optional.

Total: ~2.5–3 weeks for a thorough reproduction.

## What this sets up for DP2

Once you've reproduced this paper, you have a complete toolkit for studying topology in discrete active systems. Natural extensions:
1. **Add memory/persistence** — what if the walker has a persistence length (doesn't rotate every step)? How does the phase diagram change?
2. **Add jerk (third-order dynamics) to the lattice** — this connects directly to Stephy's suggestion 4 (jerk in Tailleur-type lattice models). The TCRW is exactly this kind of lattice model.
3. **Different lattice geometries** — triangular, honeycomb, kagome. Each has different symmetries → different topological classification.
4. **Disorder** — random D_r or ω at each site. Does the topological protection survive?
5. **Many-particle TCRW** — exclusion interactions. Does topology affect MIPS on a lattice?

The most natural DP2 connection: **adding jerk/inertia to the TCRW**. You know jerk from DP1, you know TCRW from this paper. The question "does higher-order dynamics change the topological phase?" is original and connects both halves of your work.
