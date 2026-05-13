# Wrong-Turn Diagnostics

These files are kept to document reasoning paths that were later corrected.
Do not use them for the final PI story.

`forensic_pbc_is_fine_WRONG.py` claimed the rectangular and sheared k-grids were equivalent. That was wrong for the triangular torus used by the KMC. The correct diagnostic is now `../../forensic_two_bugs.py`.
