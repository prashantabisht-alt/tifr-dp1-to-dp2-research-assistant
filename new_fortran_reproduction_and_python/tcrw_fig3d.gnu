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

# ---- arrow length scaling (Python exact reference) ----
# Fig 3(d) is an orientation plot: all nonzero arrows have fixed length;
# magnitude is encoded by colour only.
s_floor   = 0.05
slen(J) = (J > 0 ? 1.0 : s_floor)

# ---- axes ----
set xrange [-4.3 : 0.3]
set yrange [0.5 : 8.5]      # paper convention: y = 1..8 (skip y=9)
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

# ---- colorbar for log10|J_ω| (Python exact reference) ----
set cblabel 'log_{10}|J_{/Symbol w}|'  font ',12'
unset logscale cb
set cbrange [-5 : -3]
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
   1.000 '#67001f' )

# ---- plot command ----
# After fig3cde rerun: col 14 = |J_om|/T_use, col 11 = θ_ω.
# Length is fixed for nonzero arrows; color is log10(col 14).
plot_cmd = \
  "'" . f . "' u (log10($2) - 0.5*ax*slen($14)*cos($11)):($3 - 0.5*ay*slen($14)*sin($11))" . \
               ":(ax*slen($14)*cos($11)):(ay*slen($14)*sin($11)):(log10($14)) " . \
               "w vectors head filled size 0.06,20 lw 1.6 lc palette " . \
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
