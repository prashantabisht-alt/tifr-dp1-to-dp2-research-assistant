#=====================================================================
# tcrw_fig3d.gnu — Fig 3(d):  quiver of J_ω on the left wall vs D_r
#                             (L = 10, ω = 1)
#
# Matches Osat et al. Fig 3(d):
#   - x-axis:  D_r    log scale, [10^-4, 10^0]
#   - y-axis:  wall site index y  ∈ {1, ..., 8}
#   - arrows:  unit-length vectors in the direction of J_ω(y) at each
#              (D_r, y) cell, COLORED by log|J_ω|
#
# Paper claim (page 3 + caption):
#   "Note, that the chiral current J_ω is consistently directed
#    towards the bulk with π/4."
#   → arrows point at θ = +π/4  ( +x , +y )  uniformly across the
#     full (D_r, y) grid.  Magnitude drops with increasing D_r (just
#     fewer chiral-continuation events per unit time).  In the plot
#     we should see a rigid sea of identical-angle arrows with the
#     color gradient encoding |J_ω|.
#
# Same trick for log-x as Fig 3(c) (plot log10(D_r) on a linear axis).
#
# Reads:  tcrw_fig3cde_summary.txt
#         columns:  iD  D_r  y  Jx_Dr  Jy_Dr  Jx_om  Jy_om  |J_Dr|  θ_Dr  |J_ω|  θ_ω
#
# Output: tcrw_fig3d.pdf   and/or interactive qt window
#=====================================================================

if (!exists("mode")) mode = "both"

f = 'tcrw_fig3cde_summary.txt'

# ---- arrow geometry (same as Fig 3c for consistency) ----
ax = 0.20
ay = 0.40

# ---- axes ----
set xrange [-4.3 : 0.3]
set yrange [0.5 : 8.5]
set xtics ('10^{-4}' -4, '10^{-3}' -3, '10^{-2}' -2, '10^{-1}' -1, '10^{0}' 0)
set mxtics 2
set ytics 1, 1, 8
set mytics 2
set xlabel 'D_r'                       font ',14'
set ylabel 'left-edge site  y'         font ',14'
set tics scale 0.8
set grid front lc rgb '#dddddd' lw 0.4
set border lw 1.0
set title "TCRW Fig 3(d) — J_{/Symbol w}(y) on left wall  (L = 10, ω = 1)"  font 'Helvetica,11'

# ---- colorbar for log|J_ω| ----
set cblabel '|J_{/Symbol w}|'          font ',12'
set logscale cb
set cbrange [1e1 : 1e6]
set palette defined ( \
   0 '#440154', \
   1 '#3b528b', \
   2 '#21918c', \
   3 '#5ec962', \
   4 '#fde725' )

# ---- plot command (col 10 = |J_ω|, col 11 = θ_ω) ----
plot_cmd = \
  "'" . f . "' u (log10($2) - 0.5*ax*cos($11)):($3 - 0.5*ay*sin($11))" . \
               ":(ax*cos($11)):(ay*sin($11)):10 " . \
               "w vectors head filled size 0.08,20 lw 1.6 lc palette " . \
               "title ''"

#=====================================================================
# PDF  (runs FIRST)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,9cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig3d.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig3d.pdf"
}

#=====================================================================
# INTERACTIVE qt
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 900,580 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
