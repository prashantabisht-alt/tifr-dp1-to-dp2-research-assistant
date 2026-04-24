#=====================================================================
# tcrw_fig2_defects_occ.gnu — Fig 2 row (k): P(X,Y) with internal
#                              plus-sign defect, compared side-by-side
#                              with the clean box at the same ω.
#
# Two panels (ω = 0.0 for both):
#   left  : CLEAN box        — single localization ring along outer edge
#   right : DEFECTS box      — outer edge ring  +  SECOND localization
#                              ring around the interior plus-sign defect
#
# Physics:
#   The noise (D_r) rotates without translating, and does so at every
#   step with probability D_r.  Any immediate subsequent chiral step
#   that would drive the walker into a wall — outer OR interior —
#   is blocked, so the walker accumulates against every wall segment.
#   Adding an internal wall therefore ADDS a new localization ring at
#   the site of the defect; the outer ring is essentially unchanged.
#
# Shared log colorbar [1e-5 : 1e-1] → direct visual comparison of the
# two rings against the bulk.
#
# The defect cells themselves show up as blank squares (P = 0, filtered
# to NaN) and are additionally overlaid with dark grey filled squares
# from tcrw_fig2_defects_layout.txt so the plus-sign shape is obvious.
#
# Reads :
#   tcrw_fig2_occ_w0.0.txt          (clean  ω=0, from tcrw_fig2_clean)
#   tcrw_fig2_occ_defects.txt       (defects ω=0, from tcrw_fig2_defects)
#   tcrw_fig2_defects_layout.txt    (x, y of each defect cell)
#
# Output: tcrw_fig2_defects_occ.pdf + interactive qt
#
# Usage :
#   bash run.sh tcrw-plot-fig2-defects-occ       # PDF then interactive qt
#   bash run.sh tcrw-plot-fig2-defects-occ-pdf   # PDF only (headless)
#   bash run.sh tcrw-plot-fig2-defects-occ-qt    # interactive qt only
#=====================================================================

if (!exists("mode")) mode = "both"

f_clean   = 'tcrw_fig2_occ_w0.0.txt'
f_defects = 'tcrw_fig2_occ_defects.txt'
f_layout  = 'tcrw_fig2_defects_layout.txt'

# ---- palette (same as clean fig for comparison) ----
set palette defined (0 '#4de2ff', 0.5 '#ffffff', 1 '#ff4da6')
set logscale cb
set cbrange [1e-5:1e-1]
set format cb "10^{%L}"

# ---- common per-panel setup ----
set size ratio -1
set xrange [-0.5 : 9.5]
set yrange [-0.5 : 9.5]
set xtics 0, 2, 9
set ytics 0, 2, 9
set xlabel 'X'
set ylabel 'Y'
unset grid

# NaN filter: walls AND defect cells (P=0) rendered blank.
# Second 'plot' layer uses boxxyerror (x:y:xdelta:ydelta) to draw a FILLED
# CELL-SIZED dark grey square on each defect cell — so the plus-sign shape
# is unambiguous and distinguishable from the outer walls (which are just
# the blank ring around the edge).
set style fill solid 1.0 noborder
plot_clean   = "'" . f_clean   . "' u 1:2:($3 > 0 ? $3 : NaN) with image notitle"
plot_defects = "'" . f_defects . "' u 1:2:($3 > 0 ? $3 : NaN) with image notitle, '" \
             . f_layout  . "' u 1:2:(0.5):(0.5) with boxxyerror fs solid 1.0 noborder fc rgb '#222222' notitle"

#--------------------------------------------------------------------
# PDF  (headless first)
#--------------------------------------------------------------------
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 20cm,9cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig2_defects_occ.pdf'

    set multiplot layout 1,2 title \
        "TCRW Fig 2 occupancy P(X,Y): clean vs plus-sign defect   |   ω = 0.0,  D_r = 10^{-3},  T = 10^{10}" \
        font 'Helvetica,11'

    set title "clean box" font 'Helvetica,11'
    eval("plot " . plot_clean)

    set title "plus-sign defect at (4, 5)" font 'Helvetica,11'
    eval("plot " . plot_defects)

    unset multiplot
    unset output
    print "Wrote tcrw_fig2_defects_occ.pdf"
}

#--------------------------------------------------------------------
# INTERACTIVE qt
#--------------------------------------------------------------------
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 1100,500 enhanced font 'Helvetica,11'

    set multiplot layout 1,2 title \
        "TCRW Fig 2 occupancy P(X,Y): clean vs plus-sign defect" \
        font 'Helvetica,11'

    set title "clean box"
    eval("plot " . plot_clean)

    set title "plus-sign defect at (4, 5)"
    eval("plot " . plot_defects)

    unset multiplot
    pause mouse close
}
