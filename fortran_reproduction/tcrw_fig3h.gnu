#=====================================================================
# tcrw_fig3h.gnu — Fig 3(h):  quiver of J_Dr on the left wall vs ω
#                             (L = 10, D_r = 10^-3)
#
# ω-sweep twin of Fig 3(c).  Matches Osat et al. Fig 3(h):
#   - x-axis:  ω   LINEAR scale, [0, 1]  (ω is naturally linear;
#              no log-x trick is needed here — unlike Fig 3c/d)
#   - y-axis:  wall site index y ∈ {1, ..., 8}  (interior edge sites,
#              corners y = 0 and y = L-1 excluded since they're
#              geometrically anomalous)
#   - arrows:  unit-length vectors in the direction of J_Dr(y) at
#              each (ω, y) cell, COLORED by log|J_Dr|
#
# Paper claim (page 3 + caption):
#   Arrow ROTATES continuously from θ = +π/2 (UP) at ω = 0,
#   through θ = 0 at ω = 0.5 (achirality),
#   to θ = -π/2 (DOWN) at ω = 1.  At ω = 0.5 the noise-driven
#   Jx baseline dominates and arrows point into the bulk (→).
#
# Reads:  tcrw_fig3hij_summary.txt
#         columns:  iW  ω  y  Jx_Dr  Jy_Dr  Jx_om  Jy_om
#                   |J_Dr|  θ_Dr  |J_ω|  θ_ω
#
# Output: tcrw_fig3h.pdf   and/or interactive qt window
#
# Usage:
#   bash run.sh tcrw-plot-fig3h          # PDF then interactive qt
#   bash run.sh tcrw-plot-fig3h-pdf      # PDF only (headless)
#   bash run.sh tcrw-plot-fig3h-qt       # interactive qt only
#=====================================================================

if (!exists("mode")) mode = "both"

f = 'tcrw_fig3hij_summary.txt'

# ---- arrow geometry ----
# ω grid has 21 points spaced Δω = 0.05.  Half-spacing → clean
# separation between adjacent columns, but big enough to read.
ax = 0.025      # x-direction arrow scale (in ω units)
ay = 0.40       # y-direction arrow scale (in site-index units; same as Fig 3c/d)

# ---- axes (linear ω; no log trick needed) ----
set xrange [-0.03 : 1.03]
set yrange [0.5 : 8.5]        # show interior sites y = 1..8 for L = 10
set xtics 0, 0.1, 1.0
set mxtics 2
set ytics 1, 1, 8
set mytics 2
set xlabel '{/Symbol w}'              font ',14'
set ylabel 'left-edge site  y'        font ',14'
set tics scale 0.8
set grid front lc rgb '#dddddd' lw 0.4
set border lw 1.0
set title "TCRW Fig 3(h) — J_{D_r}(y) on left wall  (L = 10, D_r = 10^{-3})"  font 'Helvetica,11'

# visual guide at ω = 0.5 (achirality) so the reader sees the
# mid-point where the arrows flip sense
set arrow from 0.5, 0.5 to 0.5, 8.5 nohead lc rgb '#bbbbbb' dt 2 lw 0.8 front

# ---- colorbar for log|J_Dr| ----
set cblabel '|J_{D_r}|'                font ',12'
set logscale cb
set cbrange [1e1 : 1e6]        # per-site counts of post-noise wall departures
set palette defined ( \
   0 '#440154', \
   1 '#3b528b', \
   2 '#21918c', \
   3 '#5ec962', \
   4 '#fde725' )              # viridis, matches Fig 3(c/d)

# ---- plot command ----
# col 2 = ω, col 3 = y, col 8 = |J_Dr|, col 9 = θ_Dr (radians)
# center-tail arrows: tail at (ω - 0.5·ax·cosθ,  y - 0.5·ay·sinθ),
#                     delta = (ax·cosθ, ay·sinθ).
plot_cmd = \
  "'" . f . "' u ($2 - 0.5*ax*cos($9)):($3 - 0.5*ay*sin($9))" . \
               ":(ax*cos($9)):(ay*sin($9)):8 " . \
               "w vectors head filled size 0.015,20 lw 1.6 lc palette " . \
               "title ''"

#=====================================================================
# PDF  (headless — runs FIRST so no qt window is opened during PDF)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,9cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig3h.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig3h.pdf"
}

#=====================================================================
# INTERACTIVE qt
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 900,580 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
