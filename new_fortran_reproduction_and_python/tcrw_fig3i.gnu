#=====================================================================
# tcrw_fig3i.gnu — Fig 3(i):  quiver of J_ω on the left wall vs ω
#                             (L = 10, D_r = 10^-3)
#
# ω-sweep twin of Fig 3(d).  Matches Osat et al. Fig 3(i):
#   - x-axis:  ω   LINEAR scale, [0, 1]
#   - y-axis:  wall site index y ∈ {1, ..., 8}
#   - arrows:  unit-length vectors in the direction of J_ω(y) at each
#              (ω, y) cell, COLORED by log|J_ω|
#
# Paper claim (page 3 + caption):
#   Arrow ROTATES continuously from θ = -π/4 (DOWN-RIGHT) at ω = 0,
#   through θ = 0 at ω = 0.5 (achirality),
#   to θ = +π/4 (UP-RIGHT) at ω = 1.
#   Magnitude is roughly constant across ω (unlike |J_Dr|): the
#   chiral current J_ω comes from the chiral-continuation orbit,
#   whose RATE barely depends on ω at fixed D_r.
#
# Reads:  tcrw_fig3hij_summary.txt
#         columns:  iW  ω  y  Jx_Dr  Jy_Dr  Jx_om  Jy_om
#                   |J_Dr|  θ_Dr  |J_ω|  θ_ω
#
# Output: tcrw_fig3i.pdf   and/or interactive qt window
#=====================================================================

if (!exists("mode")) mode = "both"

f = 'tcrw_fig3hij_summary.txt'

# ---- arrow geometry (same as Fig 3h for consistency) ----
ax = 0.025
ay = 0.40

# ---- arrow length scaling (Python exact reference) ----
# Fig 3(i) is an orientation plot: all nonzero arrows have fixed length;
# magnitude is encoded by colour only.
s_floor   = 0.05
slen(J) = (J > 0 ? 1.0 : s_floor)

# ---- axes ----
set xrange [-0.03 : 1.03]
set yrange [0.5 : 8.5]        # paper convention: y = 1..8 (skip y=9)
set xtics 0, 0.1, 1.0
set mxtics 2
set ytics 1, 1, 8
set mytics 2
set xlabel '{/Symbol w}'              font ',14'
set ylabel 'left-edge site  y'        font ',14'
set tics scale 0.8
set grid front lc rgb '#dddddd' lw 0.4
set border lw 1.0
set title "TCRW Fig 3(i) — J_{/Symbol w}(y) on left wall  (L = 10, D_r = 10^{-3})"  font 'Helvetica,11'

# visual guide at ω = 0.5 (achirality)
set arrow from 0.5, 0.5 to 0.5, 9.5 nohead lc rgb '#bbbbbb' dt 2 lw 0.8 front

# ---- colorbar for log10|J_ω| (Python exact reference) ----
set cblabel 'log_{10}|J_{/Symbol w}|'  font ',12'
unset logscale cb
set cbrange [-4.65 : -4.50]
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
# After fig3hij rerun: col 14 = |J_om|/T_use, col 11 = θ_ω.
# Length is fixed for nonzero arrows; color is log10(col 14).
plot_cmd = \
  "'" . f . "' u ($2 - 0.5*ax*slen($14)*cos($11)):($3 - 0.5*ay*slen($14)*sin($11))" . \
               ":(ax*slen($14)*cos($11)):(ay*slen($14)*sin($11)):(log10($14)) " . \
               "w vectors head filled size 0.012,20 lw 1.6 lc palette " . \
               "title ''"

#=====================================================================
# PDF  (runs FIRST)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,9cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig3i.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig3i.pdf"
}

#=====================================================================
# INTERACTIVE qt
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 900,580 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
