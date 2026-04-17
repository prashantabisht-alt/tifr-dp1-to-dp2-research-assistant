"""
Quick usage example for tcrw_assembly.py
"""

from tcrw_assembly import SelfAssembly, plot_fig6a, plot_fig6bc, plot_fig6d, plot_fig6e
import numpy as np

# Example 1: Initialize assembly system
print("Example 1: Initialize assembly system")
assembly = SelfAssembly(L_target=5, L_arena=30)
print(f"  Created {assembly.n_tiles}-tile assembly system")
print()

# Example 2: Run a single assembly simulation
print("Example 2: Run single assembly (achiral)")
omega, D_r = 0.5, 0.01
tau_sa, trajs, placed, arena = assembly.run_assembly(omega, D_r, seed=42)
print(f"  ω={omega}, D_r={D_r} → τ_SA={tau_sa} steps")
print(f"  Placed {len(placed)}/25 tiles")
print()

# Example 3: Run assembly with different chirality
print("Example 3: Compare chiralities")
for omega in [0.5, 1.0]:
    tau, _, _, _ = assembly.run_assembly(omega, 0.01, seed=42)
    print(f"  ω={omega}: τ_SA={tau}")
print()

# Example 4: Generate all figures
print("Example 4: Generate all figures")
print("  plot_fig6a(assembly)        → tcrw_fig6a_tiles.png")
print("  plot_fig6bc(assembly)       → tcrw_fig6bc_trajectories.png")
print("  plot_fig6d(assembly)        → tcrw_fig6d_tau_vs_Dr.png")
print("  plot_fig6e(assembly)        → tcrw_fig6e_tau_vs_omega.png")
print()

# Example 5: Custom parameter sweep
print("Example 5: Custom sweep (2 parameters only)")
omega_vals = [0.3, 0.7]
D_r_vals = [0.01, 0.1]

for omega in omega_vals:
    for D_r in D_r_vals:
        tau, _, _, _ = assembly.run_assembly(omega, D_r, seed=99)
        print(f"  ω={omega:.1f}, D_r={D_r:.2f} → τ_SA={tau}")
