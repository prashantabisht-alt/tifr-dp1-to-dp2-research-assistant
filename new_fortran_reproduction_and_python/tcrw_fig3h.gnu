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
# separation between adjacent columns.
ax = 0.025      # max x-direction arrow scale (in ω units)
ay = 0.40       # max y-direction arrow scale (in site-index units)

# ---- arrow length scaling (Python exact reference) ----
# Fig 3(h) is the one vector panel where arrow length shrinks with |J_Dr|;
# Python uses s = 0.05 + 0.95 * |J_Dr| / max_panel(|J_Dr|).
s_floor = 0.05
stats f u (($3 >= 1 && $3 <= 8) ? $13 : 1/0) nooutput
J_max = (STATS_max > 0 ? STATS_max : 1.0)
slen(J) = (J > 0 ? s_floor + (1.0 - s_floor) * J / J_max : s_floor)

# ---- axes (linear ω; no log trick needed) ----
set xrange [-0.03 : 1.03]
set yrange [0.5 : 8.5]        # paper shows y = 1..8 (skip y=9 corner-adjacent)
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

# ---- colorbar for log10|J_Dr| (Python exact reference) ----
set cblabel 'log_{10}|J_{D_r}|'        font ',12'
unset logscale cb
set cbrange [-7 : -5]
set format cb "%g"
set palette defined ( \
   0.000 '#053061', \
   0.125 '#2166ac', \
   0.250 '#4393c3', \
   0.375 '#92c5de', \
   0.500 '#f7f7f7', \
   0.625 '#f4a582', \
   0.750 '#d6604d', \
   0.875 '#b2182b', \
   1.000 '#67001f' )          # RdBu_r-like, matches Python exact

# ---- plot command ----
# After fig3hij rerun: col 13 = |J_Dr|/T_use and col 9 = θ_Dr.
# Length scales linearly with col 13; color is log10(col 13).
plot_cmd = \
  "'" . f . "' u ($2 - 0.5*ax*slen($13)*cos($9)):($3 - 0.5*ay*slen($13)*sin($9))" . \
               ":(ax*slen($13)*cos($9)):(ay*slen($13)*sin($9)):(log10($13)) " . \
               "w vectors head filled size 0.012,20 lw 1.6 lc palette " . \
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
