#=====================================================================
# tcrw_fig2_traj.gnu — Fig 2 panels (b), (g) + symmetry check ω=1.0:
#                      3D rainbow trajectories
#
# Paper Fig 2(b)/(g) plots the walker's path as a 3D curve:
#    x-axis, y-axis = lattice coords (0 .. L-1)
#    z-axis         = log10(t)              (time lifted vertically)
#    color          = same log10(t)         (rainbow in time)
#
# 10^6 steps in 2D is an unreadable tangle (every cell visited ~10^4 times).
# Lifting along log10(t) resolves the chronology — CCW/CW circulation
# appears as a helical ribbon with a definite handedness.
#
# Three panels (left→right):
#   ω = 0.0   → left-handed (CCW when viewed from above) helix around edge
#   ω = 0.5   → edge-hugging but no preferred handedness
#   ω = 1.0   → right-handed (CW) helix — MIRROR of ω=0 (symmetry test)
#
# Reads : tcrw_fig2_traj_w0.0.txt
#         tcrw_fig2_traj_w0.5.txt
#         tcrw_fig2_traj_w1.0.txt
#         columns: t  x  y        (10^6 rows each)
#
# Output: tcrw_fig2_traj.pdf + interactive qt
#
# Usage :
#   bash run.sh tcrw-plot-fig2-traj          # PDF then interactive qt
#   bash run.sh tcrw-plot-fig2-traj-pdf      # PDF only (headless)
#   bash run.sh tcrw-plot-fig2-traj-qt       # interactive qt only
#
# Rendering notes:
#   * 3 × 10^6 segments: qt takes ~6 s.  Bump every_step to 5 or 10 if laggy.
#   * `set view 60, 30` is a moderate tilt.  Drag with mouse in qt to rotate.
#=====================================================================

if (!exists("mode")) mode = "both"

f0 = 'tcrw_fig2_traj_w0.0.txt'
f5 = 'tcrw_fig2_traj_w0.5.txt'
f1 = 'tcrw_fig2_traj_w1.0.txt'

# ---- palette ----
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

# subsample factor — was 1, bumped to 50 so individual helical turns at
# high t become visually resolvable.  10^6 segments / 50 = 20,000 plotted
# points per panel; enough density for the high-t band but thin enough
# that the CCW/CW handedness is clearly visible in 3D.
every_step = 50

plot_0 = "'" . f0 . "' every " . sprintf("%d", every_step) . " u 2:3:(log10($1)) with lines lw 0.6 lc palette notitle"
plot_5 = "'" . f5 . "' every " . sprintf("%d", every_step) . " u 2:3:(log10($1)) with lines lw 0.6 lc palette notitle"
plot_1 = "'" . f1 . "' every " . sprintf("%d", every_step) . " u 2:3:(log10($1)) with lines lw 0.6 lc palette notitle"

#--------------------------------------------------------------------
# PDF  (headless first)
#--------------------------------------------------------------------
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 30cm,10cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig2_traj.pdf'

    set multiplot layout 1,3 title \
        "TCRW Fig 2 trajectories (first 10^6 steps, lifted by log_{10} t)   |   D_r = 10^{-3}" \
        font 'Helvetica,11'

    set title "ω = 0.0  (expect CCW helix)"  font 'Helvetica,11'
    eval("splot " . plot_0)

    set title "ω = 0.5  (no handedness)"      font 'Helvetica,11'
    eval("splot " . plot_5)

    set title "ω = 1.0  (expect CW helix)"    font 'Helvetica,11'
    eval("splot " . plot_1)

    unset multiplot
    unset output
    print "Wrote tcrw_fig2_traj.pdf"
}

#--------------------------------------------------------------------
# INTERACTIVE qt
#--------------------------------------------------------------------
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 1600,560 enhanced font 'Helvetica,11'

    set multiplot layout 1,3 title \
        "TCRW Fig 2 trajectories (first 10^6 steps, lifted by log_{10} t)   |   D_r = 10^{-3}" \
        font 'Helvetica,11'

    set title "ω = 0.0  (expect CCW helix)"
    eval("splot " . plot_0)

    set title "ω = 0.5  (no handedness)"
    eval("splot " . plot_5)

    set title "ω = 1.0  (expect CW helix)"
    eval("splot " . plot_1)

    unset multiplot
    pause mouse close
}
