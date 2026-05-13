# Legacy Debug Archive

This folder stores the exploratory scripts and plots from the k-grid/sign-error hunt.
They are not deleted because they document how the result was found, but they
are not the source of truth for future calculations.

Use the root-level files instead:

- `../triangular_jmvr_corrected.py`
- `../fig11_final_hex.py`
- `../export_fig11_final_hex_gnuplot_data.py`
- `../fig11_final_hex.gnu`
- `../kmc_triangular_jmvr.f90`

## Subfolders

| Folder | Contents |
|---|---|
| `scripts/` | Old postmortem, KMC verification, plotting, and forensic scripts |
| `figures/` | Old generated figures from intermediate attempts |
| `gnuplot_attempts/` | Superseded gnuplot scripts |
| `outputs/` | Old text tables and extracted draft-page images |
| `wrong_turns/` | Reasoning paths that were later corrected |

## Important Warning

Some scripts here may have stale relative imports because they were moved after
the debugging phase. Treat them as historical evidence, not active workflows.
