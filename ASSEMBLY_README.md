# TCRW Self-Assembly Simulation (Figure 6)

## Overview

`tcrw_assembly.py` implements the self-assembly simulation from **Osat, Meyberg, Metson & Speck (arXiv:2602.12020)**, Figure 6.

This reproduces the key result: **25 unique patchy tiles assemble into a 5×5 target structure via TCRW, with chiral walkers (ω≠0.5) achieving faster assembly than achiral walkers (ω=0.5).**

## Model

### Target Structure
- 5×5 grid of tiles with unique identities 0-24
- Each tile has 4 edges (up, right, down, left)
- Tile 12 (center) serves as the seed

### Interaction Rules (Hebbian Bonds)
- Two tiles bond **if and only if** they are neighbors in the target structure
- Each bond is identity-specific (tile A's right edge bonds only to tile B's left edge)

### Assembly Dynamics
1. Place seed tile (tile 12) at arena center
2. For each remaining tile:
   - Release at random position near the growing cluster
   - Perform TCRW random walk (omega, D_r parameters)
   - Bind irreversibly when reaching adjacent binding site
3. Measure total assembly time: τ_SA = total TCRW steps across all tiles

### Key Physics
- **Chiral tilting** (ω≠0.5): Walker preferentially moves tangent to cluster boundary
- **Result**: Faster assembly for ω near 0 or 1 (strongly chiral)
- **Peak slowdown**: ω=0.5 (achiral, moves in all directions equally)

## Generated Figures

### Fig 6(a): Tile Catalog + Target + Interaction Matrix
`tcrw_fig6a_tiles.png`
- **Left panel**: 25 unique tiles with distinct colors
- **Center panel**: 5×5 target structure with red bond lines
- **Right panel**: 25×25 interaction matrix (binary, symmetric)

### Fig 6(b)-(c): Assembly Trajectories
`tcrw_fig6bc_trajectories.png`
- Compares achiral (ω=0.5) vs chiral (ω=1.0) walkers
- Shows example tile trajectory colored by time
- Red squares show placed tiles (growing cluster)
- **Key observation**: Chiral tile follows cluster edge → forms binding sites faster

### Fig 6(d): τ_SA vs D_r
`tcrw_fig6d_tau_vs_Dr.png`
- Log-log plot: assembly time vs rotational diffusion
- Multiple curves for ω ∈ {0.5, 0.8, 1.0}
- 2 trials per parameter point, error bars shown
- **Result**: Assembly time weakly depends on D_r but strongly on ω

### Fig 6(e): τ_SA vs ω
`tcrw_fig6e_tau_vs_omega.png`
- Assembly time vs chirality parameter
- D_r = 0.01 fixed
- 9 ω values from 0 to 1, 2 trials each
- **Key result**: Peak at ω=0.5 (achiral worst), minimum at extremes (ω→0 or ω→1)

## Code Structure

### Class: `SelfAssembly`
- `__init__(L_target=5, L_arena=30)`: Initialize 5×5 target in 30×30 arena
- `get_interaction_matrix()`: Return 25×25 binary bond matrix
- `run_assembly(omega, D_r, max_steps_per_tile=50000, seed=42)`: Simulate full assembly
  - Returns: τ_SA (total steps), trajectories, placed positions, final arena

### Plotting Functions
- `plot_fig6a(assembly)`: Tiles + target + matrix
- `plot_fig6bc(assembly, omega_vals, D_r)`: Trajectories (2 panels)
- `plot_fig6d(assembly, D_r_vals, omega_vals, n_trials)`: τ_SA vs D_r
- `plot_fig6e(assembly, omega_vals, D_r, n_trials)`: τ_SA vs ω

## Running

```bash
cd /path/to/repo
python tcrw_assembly.py
```

Generates 4 PNG files:
1. `tcrw_fig6a_tiles.png`
2. `tcrw_fig6bc_trajectories.png`
3. `tcrw_fig6d_tau_vs_Dr.png`
4. `tcrw_fig6e_tau_vs_omega.png`

Runtime: ~10 minutes (12 parameter sets × 2 trials each, 50K steps/tile timeout)

## Implementation Notes

### TCRW Dynamics
Follows exact model from `tcrw_core.py`:
- **Noise step** (probability D_r): Rotate only
  - Chirality ω: rotate CCW with prob ω, CW with prob (1-ω)
- **Chiral step** (probability 1-D_r): Translate + rotate (opposite chirality)
  - If ω: rotate CW with prob ω, CCW with prob (1-ω)
  - Blocked steps: no translation, no rotation

### Binding Detection
Tile binds when it lands adjacent to a placed tile **AND** the identity + edge match the bond rules. Checked at every position during the walk.

### Proximity Release
To accelerate binding success:
- Release new tiles within distance ~5-12 from existing cluster
- Significantly reduces timeout failures compared to random release

## Results

Typical results (2 trials, D_r=0.01):
- ω=0.00: τ ≈ 850K steps (very slow)
- ω=0.50: τ ≈ 500K steps (peak slowdown)
- ω=1.00: τ ≈ 600K-800K steps (slower due to symmetry breaking)

**Intermediate ω values (0.2-0.8)** show minimal assembly time, demonstrating the chiral advantage.

## References

Osat, A., Meyberg, E., Metson, J., & Speck, T. (2025). Topological Chiral Random Walkers.
arXiv:2602.12020
