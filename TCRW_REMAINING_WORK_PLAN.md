# TCRW: Complete Remaining Work Plan

Every missing figure panel, no exceptions. Organized by implementation order (dependencies first).

---

## PHASE A: Generalized Geometry Engine (FOUNDATION — do first)

Everything involving defects, internal walls, hybrid BC depends on this.

### A1. Create `tcrw_geometry.py` — flexible lattice mask

New file. Core abstraction: a mask that says which sites exist and which transitions are allowed.

```
Classes:
  RectangleMask(L)                        — simple L×L box (current behavior)
  RectangleWithDefects(L, blocked_sites)   — L×L with removed sites (internal walls)
  RectangleWithHoles(L, hole_regions)      — disconnected internal boundaries
  HybridBCMask(L, pbc_x, pbc_y)           — mixed PBC/OBC per axis
  
Methods:
  .is_valid(x, y) → bool
  .neighbors(x, y, d) → (nx, ny) or None if blocked
  .valid_sites() → list of (x,y)
  .n_sites → int
  .site_to_index / index_to_site — mapping for variable-size state space
```

Complexity: Medium. No dependencies.

### A2. Generalize `build_transition_matrix` in `tcrw_obc.py`

Add `build_transition_matrix_generic(omega, D_r, mask)` that uses mask.is_valid() instead of hard-coded bounds. Keep old function as wrapper for backward compatibility. Also generalize `state_index` to work with non-rectangular grids (use mask.site_to_index).

### A3. Generalize `exact_currents` in `tcrw_currents.py`

Modify `build_split_matrices` and `exact_currents` to accept a mask. Add boundary-type tagging: for each edge site, tag whether it's external boundary or internal (defect) boundary. This enables decomposing currents by boundary type for Fig 11/12.

### A4. Generalize `simulate_tcrw_obc` in `tcrw_core.py`

Add `simulate_tcrw_geometry(omega, D_r, mask, ...)` — same OBC walker but uses mask for blocking. Keep old function intact.

### A5. Add `obc_spectrum_generic` in `tcrw_spectrum.py`

Same as existing `obc_spectrum(omega, D_r, L)` but takes a mask. Returns eigenvalues + edge weights decomposed by boundary type.

---

## PHASE B: Quick-win spectrum panels (parallel with Phase A)

These just need parameter sweeps of existing code. No new infrastructure.

### B1. Fig 4(d) — Re(λ) vs D_r for ω=1

Add to `tcrw_spectrum.py`. OBC spectrum (L=10), sweep D_r from ~0 to 0.5. Plot Re(λ) of all eigenvalues vs D_r, colored by edge weight. Shows eigenvalue coalescence as D_r→0.

Parameters: ω=1, L=10, D_r ∈ [0.01, 0.5] (say 30 values).

### B2. Fig 4(e) — Re(λ) vs ω for fixed D_r

Same approach, sweep ω at fixed D_r. Shows gap closing at ω=0.5.

Parameters: D_r=0.1 (or read from paper), L=10, ω ∈ [0, 1].

### B3. Fig 4(h) — Hybrid PBC/OBC spectrum

Needs new matrix builder: periodic in y, open in x. Use HybridBCMask from Phase A. Then: sweep k_y ∈ [-π,π), at each k_y build a (4L)×(4L) matrix (only x is real-space), diagonalize. Plot band structure vs k_y.

Key function:
```python
def build_hpbc_matrix(omega, D_r, L, ky):
    """
    Hybrid: x-direction OBC (real space, L sites), y-direction PBC (Fourier, momentum ky).
    Returns (4L)×(4L) matrix.
    """
```

Parameters: ω=1, D_r=0.1, L=10, N_ky=200.

### B4. Fig 4(i) — Band structure on (cos kx, cos ky) circle

PBC spectrum. Sweep full BZ, plot eigenvalues in the (cos kx, cos ky) plane colored by |λ| or Re(λ). Straightforward visualization of existing `pbc_full_bz` data.

### B5. Fig 8(d)-(l) — Nine spectrum panels

3×3 grid of OBC complex-plane spectra at specific (ω, D_r) points. The paper picks points in topological / trivial / transition regions.

Read exact parameter values from paper's Fig 8 caption. For each point: compute OBC spectrum (L=10), plot in complex plane colored by edge weight.

Add as `fig8_spectrum_grid()` in `tcrw_zak_phase.py` or new file.

---

## PHASE C: Missing main-text figure panels

### C1. Fig 2(b) — Single long OBC trajectory

Run existing `simulate_tcrw_obc` with ω=1, D_r=1e-3, L=10, T=10^6 steps, record_traj=True. Plot trajectory on lattice grid, color by time. Use LineCollection or scatter with alpha.

Complexity: Easy. 30 min.

### C2. Fig 2(i)-(o) — Defect boundary panels

Depends on Phase A (generalized geometry).

Paper shows two types of defects:
- **Edge deformation**: notch or bump on the outer boundary (modifies shape but no new internal boundary)
- **Internal defect**: removed block of sites in the bulk (creates a hole with its own internal boundary)

Need to figure out exact geometries from the paper figures. Create 2-3 defect masks, then for each:
- Compute P(X,Y) via exact steady state
- Compute J, J_ω, J_Dr via exact currents
- Show current flows around defects, internal boundary currents

Parameters: ω=1 and ω=0, D_r=1e-3, L=10 or 15.

Panels:
- (i) J_total with edge defects, ω=0
- (j) J_chiral with defects, ω=0
- (k) P(X,Y) with defects
- (l) trajectory with defects showing edge current along boundaries
- (m)-(o) similar with internal defects showing opposite chirality on internal vs external boundary

Complexity: Medium-Hard. Depends on A1-A3.

### C3. Fig 5(b) — MFPT heatmap by starting position

Generate maze (L=20). For each passage cell, for each initial director d ∈ {↑,→,↓,←}, measure MFPT to exit averaged over many trials. Display as 4 heatmaps (one per director) or single heatmap averaged over directors.

Parameters: ω=1, D_r=1e-3, L=20, N_trials_per_cell ~50-100.

Complexity: Medium (compute-heavy, may need Fortran for L=20). Could start with L=10.

### C4. Fig 5(e) — Disconnected mazes

Paper says: "hand-on-the-wall rule fails for non-simply-connected mazes." This means mazes with loops/islands — the walker can get trapped circling an island forever if D_r is too small.

Implementation:
- Generate maze with Prim's algorithm, then remove extra walls to create loops (break simply-connectedness)
- Or: generate maze, add a "bridge" that creates a cycle
- Show MFPT vs D_r: at very small D_r, chiral walker gets stuck on loops → MFPT diverges. At moderate D_r, noise helps escape → MFPT drops.

Complexity: Medium.

### C5. Fig 5(f) — Wide-passage mazes

Modify maze generator to use passage width = 3 (instead of 1). Each passage is 3 cells wide, walls are 1 cell thick. The (2L+1)×(2L+1) grid becomes bigger.

Implementation:
```python
def generate_wide_maze_prim(L, passage_width=3, seed=42):
    # Scale up: each logical cell becomes passage_width × passage_width
    # Walls remain 1 cell thick between passages
```

Parameters: L=8 (logical), passage_width=3, ω=0,0.5,1, D_r=1e-3.

Complexity: Easy-Medium.

---

## PHASE D: Self-Assembly (Fig 6) — Biggest new feature

This is a complete multi-particle simulation. The paper uses a specific Hebbian-like interaction rule.

### D1. Understand the assembly model

From paper: 25 unique patchy tiles on a 5×5 target grid. Each tile has 4 edges. Two tiles interact specifically if they are neighbors in the target structure — "Hebbian-like learning rule." Interaction matrix I[i,j] = 1 if tiles i and j are neighbors in the target.

Assembly process:
1. Start with a seed tile placed at the center of a larger arena
2. Free tiles perform TCRW random walks in the arena
3. When a free tile reaches a position adjacent to the growing cluster AND the interaction is compatible, it binds
4. Measure time until all 25 tiles are assembled = τ_SA

Key physics: Chiral tiles follow the edge of the growing cluster → find binding sites faster → τ_SA is smaller.

### D2. Create `tcrw_assembly.py`

New file. Functions:
```python
def create_target_structure(L_target=5):
    """Define 5×5 target with 25 unique tiles. Return adjacency/interaction matrix."""

def simulate_assembly(omega, D_r, L_arena, L_target, max_steps, seed):
    """
    Full assembly simulation.
    - Place seed tile at center of L_arena × L_arena arena
    - Release one free tile at random position
    - Walk until it binds or timeout
    - Repeat for next tile
    - Return τ_SA = total steps to complete assembly
    """

def fig6a_tiles_and_target():
    """Show 25 tiles, target structure, interaction matrix."""

def fig6b_trajectories():
    """Achiral vs chiral tile trajectory during assembly."""

def fig6c_sample_trajectories():
    """Sample tile trajectories for ω=0.5 and ω=1."""

def fig6d_tau_vs_Dr():
    """τ_SA vs D_r for different ω."""

def fig6e_tau_vs_omega():
    """τ_SA vs ω for different L and D_r."""
```

Parameters from paper:
- L_target = 5 (5×5 = 25 tiles)
- L_arena = ? (probably 20-50, much larger than target)
- D_r range: 10^-3 to 10^-1
- ω values: 0.2, 0.5, 0.8, 1.0
- N_trials: 50+ per parameter point

Complexity: **Hard**. This is the single biggest piece of new work. Estimate 10-15 hours.

---

## PHASE E: Extended data figures (Figs 11, 12)

Both depend on Phase A (generalized geometry).

### E1. Fig 11 — Effect of boundaries on spectrum

- (a) Visualization: L×L OBC lattice with sites colored
- (b) OBC spectrum vs D_r for ω=1, eigenvalues colored by edge weight (similar to existing Fig 4 but as waterfall)
- (c) PBC lattice with internal defect (hole in the middle)
- (d) Spectrum of PBC+defect system showing contributions from internal boundary

For (c)-(d): need PBC matrix with internal defect. The walker wraps at outer boundary (PBC) but is blocked by the internal defect (OBC-like). This is a new boundary condition type.

Implementation: Use `build_transition_matrix_generic` with a mask that has PBC on outer edges and OBC on defect edges.

Parameters: L=10-15, ω=1, D_r sweep.

Complexity: Hard (PBC + internal defect is a new matrix type).

### E2. Fig 12 — Disconnected boundaries and nested spectrum

- (a) Lattice with multiple separated defects (e.g., two holes)
- (b) Spectrum showing external boundary contribution
- (c) Different defect configuration
- (d) Nested spectrum: separate ovals for external vs internal boundaries

This requires:
- Multiple defect regions in one lattice
- Decompose eigenvector weights by WHICH boundary each mode is localized on
- Plot spectra colored/separated by boundary identity

Parameters: L=15-20 with 2-3 defect blocks, ω=1, D_r=0.1.

Complexity: Hard.

---

## IMPLEMENTATION ORDER (dependency-respecting)

```
Week 1:  A1 → A2 → A3,A4,A5 (geometry engine + generalize all builders)
         B1,B2,B4 in parallel (quick spectrum panels, no new infrastructure)
         C1 (Fig 2b trajectory — trivial)

Week 2:  B3 (hybrid BC — needs A1)
         B5 (Fig 8 spectrum grid)
         C2 (Fig 2 defects — needs A1-A3)
         C3 (Fig 5b MFPT heatmap)

Week 3:  C4,C5 (maze variants)
         D1,D2 (self-assembly — biggest single chunk)
         E1 (Fig 11)

Week 4:  D2 continued (assembly tuning + parameter sweeps)
         E2 (Fig 12)
         Final integration + verification
```

---

## EVERY PANEL CHECKLIST

Check off as completed:

### Fig 1 (4 panels)
- [x] (a) schematic — skip
- [x] (b) trajectories
- [x] (c) MSD
- [x] (d) D vs ω

### Fig 2 (15 panels)
- [x] (a) P(X,Y) ω=1
- [ ] (b) trajectory 10^6 steps ← **C1**
- [x] (c) J for ω=1
- [x] (d) J_ω for ω=1
- [x] (e) J_Dr for ω=1
- [x] (f) P(X,Y) ω=0.5
- [x] (g) J for ω=0.5
- [x] (h) J_Dr for ω=0.5
- [ ] (i) J with edge defects ω=0 ← **C2**
- [ ] (j) J_chiral with defects ω=0 ← **C2**
- [ ] (k) P(X,Y) with defects ← **C2**
- [ ] (l) trajectory with defects ← **C2**
- [ ] (m) J with internal defect ← **C2**
- [ ] (n) J_chiral with internal defect ← **C2**
- [ ] (o) internal vs external boundary currents ← **C2**

### Fig 3 (10 panels) — ALL DONE
- [x] (a) through (j)

### Fig 4 (9 panels)
- [x] (a) state labeling — skip (diagram)
- [x] (b) PBC band structure
- [x] (c) OBC complex plane L=2
- [ ] (d) Re(λ) vs D_r ← **B1**
- [ ] (e) Re(λ) vs ω ← **B2**
- [x] (f) OBC complex plane L=10 varying D_r
- [x] (g) OBC complex plane L=10 varying ω
- [ ] (h) hybrid PBC/OBC spectrum ← **B3**
- [ ] (i) band on (cos kx, cos ky) ← **B4**

### Fig 5 (6 panels)
- [x] (a) maze trajectories
- [ ] (b) MFPT heatmap by position ← **C3**
- [x] (c) MFPT vs ω
- [x] (d) MFPT vs L
- [ ] (e) disconnected mazes ← **C4**
- [ ] (f) wide-passage mazes ← **C5**

### Fig 6 (5 panels) — ALL MISSING
- [ ] (a) tile catalog + interaction matrix ← **D2**
- [ ] (b) random vs directed trajectories ← **D2**
- [ ] (c) sample assembly trajectories ← **D2**
- [ ] (d) τ_SA vs D_r ← **D2**
- [ ] (e) τ_SA vs ω ← **D2**

### Fig 7 (5 panels)
- [x] (a)-(d) schematics — skip
- [x] (e) P(τ_edge)

### Fig 8 (12 panels)
- [x] (a) Φ_x phase diagram
- [x] (b) Φ_y phase diagram
- [x] (c) gap boundary
- [ ] (d)-(l) nine spectrum panels ← **B5**

### Fig 9 (2 panels) — schematics, skip

### Fig 10 — DONE
- [x] PBC vs OBC vs 1D edge

### Fig 11 (4 panels) — ALL MISSING
- [ ] (a) OBC lattice visualization ← **E1**
- [ ] (b) spectrum vs D_r colored by edge ← **E1**
- [ ] (c) PBC + internal defect lattice ← **E1**
- [ ] (d) spectrum with internal boundary ← **E1**

### Fig 12 (4 panels) — ALL MISSING
- [ ] (a) multi-defect lattice ← **E2**
- [ ] (b) external boundary spectrum ← **E2**
- [ ] (c) internal boundary lattice ← **E2**
- [ ] (d) nested spectrum ← **E2**

---

## TOTAL COUNT

Done: ~40 panels
Remaining: ~30 panels
Schematics skipped: ~8 panels

After this plan: **every single data panel in the paper reproduced.**
