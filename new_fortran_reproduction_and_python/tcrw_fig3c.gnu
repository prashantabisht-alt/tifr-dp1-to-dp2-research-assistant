#=====================================================================
# tcrw_fig3c.gnu — Fig 3(c):  quiver of J_Dr on the left wall vs D_r
#                             (L = 10, ω = 1)
#
# Matches Osat et al. Fig 3(c):
#   - x-axis:  D_r    log scale, [10^-4, 10^0]
#   - y-axis:  wall site index y  ∈ {1, ..., 8}  (interior edge sites,
#              excluding the two corners y = 0 and y = L-1 which
#              are geometrically anomalous)
#   - arrows:  unit-length vectors in the direction of J_Dr(y) at
#              each (D_r, y) cell, COLORED by log|J_Dr|
#
# Paper claim (page 3):
#   "Initially J_Dr points downwards as expected from the counter-
#    clockwise edge current. Increasing D_r leads to scattering of
#    the current towards the bulk."  → arrows rotate from θ = −π/2
#    (straight DOWN) at small D_r to θ ≈ 0 (straight RIGHT, into the
#    bulk) at large D_r.  Fig 3(e) plots this angle as a scalar.
#
# Trick for quiver on a log-x axis
# --------------------------------
#   `with vectors` draws (x, y) → (x+dx, y+dy) in plot units.  On a
#   log-x axis, a fixed dx makes arrow length depend on x.  We side-
#   step that by plotting log10(D_r) on a LINEAR x-axis and putting
#   the 10^k tick labels back manually.  Then the arrow dx is in
#   "decade-fraction" units and has a uniform visual length.
#
# Reads:  tcrw_fig3cde_summary.txt
#         columns:  iD  D_r  y  Jx_Dr  Jy_Dr  Jx_om  Jy_om  |J_Dr|  θ_Dr  |J_ω|  θ_ω
#
# Output: tcrw_fig3c.pdf   and/or interactive qt window
#
# Usage:
#   bash run.sh tcrw-plot-fig3c          # PDF then interactive qt
#   bash run.sh tcrw-plot-fig3c-pdf      # PDF only (headless)
#   bash run.sh tcrw-plot-fig3c-qt       # interactive qt only
#=====================================================================

if (!exists("mode")) mode = "both"

f = 'tcrw_fig3cde_summary.txt'

# ---- arrow geometry  (tuned for 14 cm × 9 cm panel, L = 10) ----
ax = 0.20       # max x-direction arrow scale (decade fractions)
ay = 0.40       # max y-direction arrow scale (site-index units)

# ---- arrow length scaling (Python exact reference) ----
# Fig 3(c) is an orientation plot: all nonzero arrows have fixed length;
# magnitude is encoded by colour only.
s_floor   = 0.05
slen(J) = (J > 0 ? 1.0 : s_floor)

# ---- axes (log-x rendered via log10(D_r) on linear axis) ----
set xrange [-4.3 : 0.3]
set yrange [0.5 : 8.5]      # paper shows y = 1..8 (skip y=9 corner-adjacent)
set xtics ('10^{-4}' -4, '10^{-3}' -3, '10^{-2}' -2, '10^{-1}' -1, '10^{0}' 0)
set mxtics 2
set ytics 1, 1, 8
set mytics 2
set xlabel 'D_r'                       font ',14'
set ylabel 'left-edge site  y'         font ',14'
set tics scale 0.8
set grid front lc rgb '#dddddd' lw 0.4
set border lw 1.0
set title "TCRW Fig 3(c) — J_{D_r}(y) on left wall  (L = 10, ω = 1)"  font 'Helvetica,11'

# ---- colorbar for log10|J_Dr| (Python exact reference) ----
set cblabel 'log_{10}|J_{D_r}|'        font ',12'
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
   1.000 '#67001f' )          # RdBu_r-like, matches Python exact

# ---- plot command ----
# After fig3cde rerun, col 13 = |J_Dr|/T_use and col 9 = θ_Dr.
# Length is fixed for nonzero arrows; color is log10(col 13).
plot_cmd = \
  "'" . f . "' u (log10($2) - 0.5*ax*slen($13)*cos($9)):($3 - 0.5*ay*slen($13)*sin($9))" . \
               ":(ax*slen($13)*cos($9)):(ay*slen($13)*sin($9)):(log10($13)) " . \
               "w vectors head filled size 0.06,20 lw 1.6 lc palette " . \
               "title ''"

#=====================================================================
# PDF  (headless — runs FIRST so no qt window is opened during PDF)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,9cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig3c.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig3c.pdf"
}

#=====================================================================
# INTERACTIVE qt
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 900,580 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
