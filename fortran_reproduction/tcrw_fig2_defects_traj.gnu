#=====================================================================
# tcrw_fig2_defects_traj.gnu — Fig 2 row (l): 3D rainbow trajectory
#                               with internal plus-sign defect,
#                               side-by-side with clean box at ω = 0.
#
# Two panels (ω = 0.0 for both):
#   left  : CLEAN box       — helical ribbon hugging outer edge only
#   right : DEFECTS box     — helical ribbon circulates AROUND the
#                             interior plus-sign as well; the defect
#                             carves a visible hole in the path
#                             (walker never visits those 5 cells).
#
# 3D geometry (same as tcrw_fig2_traj.gnu):
#   (x, y) = lattice coords
#   z      = log10(t)
#   color  = log10(t)
#   every_step = 50 — thin to ~2 × 10^4 plotted points per panel so
#                     helical turns at high t are visually resolvable.
#
# On the defects panel, we stamp a column of filled dark squares at the
# base (z = 0) of each defect cell.  Combined with the hole punched in
# the ribbon itself, this makes the plus-sign geometry unambiguous.
#
# Reads :
#   tcrw_fig2_traj_w0.0.txt         (clean, from tcrw_fig2_clean)
#   tcrw_fig2_traj_defects.txt      (defects, from tcrw_fig2_defects)
#   tcrw_fig2_defects_layout.txt    (x, y of each defect cell)
#   t  x  y   (10^6 rows each)
#
# Output: tcrw_fig2_defects_traj.pdf + interactive qt
#
# Usage :
#   bash run.sh tcrw-plot-fig2-defects-traj       # PDF then interactive qt
#   bash run.sh tcrw-plot-fig2-defects-traj-pdf   # PDF only (headless)
#   bash run.sh tcrw-plot-fig2-defects-traj-qt    # interactive qt only
#=====================================================================

if (!exists("mode")) mode = "both"

f_clean   = 'tcrw_fig2_traj_w0.0.txt'
f_defects = 'tcrw_fig2_traj_defects.txt'
f_layout  = 'tcrw_fig2_defects_layout.txt'

# ---- palette (same time-lift rainbow as tcrw_fig2_traj.gnu) ----
set palette defined (0 '#1a0c00', 0.4 '#663300', 0.8 '#ff9933', 1 '#ffe680')
set cbrange [0:6]                       # log10(t), t ∈ [10⁰, 10⁶]
set format cb "10^{%g}"
set cbtics 1
set cblabel "t"

# ---- 3D geometry ----
set xrange [-0.5 : 9.5]
set yrange [-0.5 : 9.5]
set xtics 0, 2, 9
set ytics 0, 2, 9
set xlabel 'X'
set ylabel 'Y'

set zrange [0:6]
set ztics 0, 1, 6
set format z "10^{%g}"
set zlabel "t" rotate parallel

set view 60, 30
set view equal xy
set ticslevel 0
unset grid

every_step = 50

traj_clean   = "'" . f_clean   . "' every " . sprintf("%d", every_step) \
             . " u 2:3:(log10($1)) with lines lw 0.6 lc palette notitle"
traj_defects = "'" . f_defects . "' every " . sprintf("%d", every_step) \
             . " u 2:3:(log10($1)) with lines lw 0.6 lc palette notitle, '" \
             . f_layout . "' u 1:2:(0) with points pt 5 ps 2.8 lc rgb '#222222' notitle"

#--------------------------------------------------------------------
# PDF  (headless first)
#--------------------------------------------------------------------
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 22cm,11cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig2_defects_traj.pdf'

    set multiplot layout 1,2 title \
        "TCRW Fig 2 trajectory (first 10^6 steps, lifted by log_{10} t): clean vs plus-sign defect   |   ω = 0.0,  D_r = 10^{-3}" \
        font 'Helvetica,11'

    set title "clean box" font 'Helvetica,11'
    eval("splot " . traj_clean)

    set title "plus-sign defect at (4, 5)" font 'Helvetica,11'
    eval("splot " . traj_defects)

    unset multiplot
    unset output
    print "Wrote tcrw_fig2_defects_traj.pdf"
}

#--------------------------------------------------------------------
# INTERACTIVE qt
#--------------------------------------------------------------------
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 1200,600 enhanced font 'Helvetica,11'

    set multiplot layout 1,2 title \
        "TCRW Fig 2 trajectory (first 10^6 steps, lifted by log_{10} t): clean vs plus-sign defect" \
        font 'Helvetica,11'

    set title "clean box"
    eval("splot " . traj_clean)

    set title "plus-sign defect at (4, 5)"
    eval("splot " . traj_defects)

    unset multiplot
    pause mouse close
}
