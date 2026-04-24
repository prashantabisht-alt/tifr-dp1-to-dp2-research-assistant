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
ax = 0.20       # x-direction arrow scale (decade fractions)
ay = 0.40       # y-direction arrow scale (site-index units)

# ---- axes (log-x rendered via log10(D_r) on linear axis) ----
set xrange [-4.3 : 0.3]
set yrange [0.5 : 9.5]      # show interior sites y = 1..9 for L_paper=10 (11×11 grid, indices 0..10)
set xtics ('10^{-4}' -4, '10^{-3}' -3, '10^{-2}' -2, '10^{-1}' -1, '10^{0}' 0)
set mxtics 2
set ytics 1, 1, 9
set mytics 2
set xlabel 'D_r'                       font ',14'
set ylabel 'left-edge site  y'         font ',14'
set tics scale 0.8
set grid front lc rgb '#dddddd' lw 0.4
set border lw 1.0
set title "TCRW Fig 3(c) — J_{D_r}(y) on left wall  (L = 10, ω = 1)"  font 'Helvetica,11'

# ---- colorbar for log|J_Dr| ----
set cblabel '|J_{D_r}|'                font ',12'
set logscale cb
set cbrange [1e1 : 1e6]        # per-site counts of post-noise wall departures
set palette defined ( \
   0 '#440154', \
   1 '#3b528b', \
   2 '#21918c', \
   3 '#5ec962', \
   4 '#fde725' )              # viridis, matches other Fig 3 panels

# ---- plot command ----
# col 2 = D_r, col 3 = y, col 8 = |J_Dr|, col 9 = θ_Dr (radians)
# center-tail arrows: start at (log10(D_r) - 0.5·dx,  y - 0.5·dy),
#                     delta = (dx, dy)     with dx = ax·cos θ,  dy = ay·sin θ
plot_cmd = \
  "'" . f . "' u (log10($2) - 0.5*ax*cos($9)):($3 - 0.5*ay*sin($9))" . \
               ":(ax*cos($9)):(ay*sin($9)):8 " . \
               "w vectors head filled size 0.08,20 lw 1.6 lc palette " . \
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
