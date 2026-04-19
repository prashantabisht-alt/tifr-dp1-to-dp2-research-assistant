#=====================================================================
# tcrw_fig2_occ.gnu — Fig 2 panels (a), (f) + symmetry check ω=1.0:
#                      P(X,Y) heatmaps
#
# Three panels side-by-side:
#   left   : ω = 0.0  (fully chiral, CCW)   — boundary-localized
#   middle : ω = 0.5  (achiral)             — boundary-localized, symmetric
#   right  : ω = 1.0  (fully chiral, CW)    — mirror of ω=0; same density
#                                             (symmetry check, not a paper panel)
#
# Paper expectation: P(x,y) depends on |ω−1/2| only (by parity symmetry
# of the step rule), so ω=0 and ω=1 should give the IDENTICAL density
# up to statistical noise.  We use this as a sanity test.
#
# Shared log colorbar [1e-5 : 1e-1] for direct comparison.
#
# Reads : tcrw_fig2_occ_w0.0.txt
#         tcrw_fig2_occ_w0.5.txt
#         tcrw_fig2_occ_w1.0.txt
#         columns: x  y  P(x,y)
#         ALL L×L sites written; wall sites have P = 0.
#
# Output: tcrw_fig2_occ.pdf + interactive qt
#
# Usage :
#   bash run.sh tcrw-plot-fig2-occ          # PDF then interactive qt
#   bash run.sh tcrw-plot-fig2-occ-pdf      # PDF only (headless)
#   bash run.sh tcrw-plot-fig2-occ-qt       # interactive qt only
#=====================================================================

if (!exists("mode")) mode = "both"

f0 = 'tcrw_fig2_occ_w0.0.txt'
f5 = 'tcrw_fig2_occ_w0.5.txt'
f1 = 'tcrw_fig2_occ_w1.0.txt'

# ---- palette ----
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

# Filter: walls (P=0) rendered as NaN → blank.
plot_cmd_0 = "'" . f0 . "' u 1:2:($3 > 0 ? $3 : NaN) with image notitle"
plot_cmd_5 = "'" . f5 . "' u 1:2:($3 > 0 ? $3 : NaN) with image notitle"
plot_cmd_1 = "'" . f1 . "' u 1:2:($3 > 0 ? $3 : NaN) with image notitle"

#--------------------------------------------------------------------
# PDF  (headless first)
#--------------------------------------------------------------------
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 28cm,9cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig2_occ.pdf'

    set multiplot layout 1,3 title \
        "TCRW Fig 2 occupancy P(X,Y)   |   L = 10,  D_r = 10^{-3},  T = 10^{10}" \
        font 'Helvetica,11'

    set title "ω = 0.0  (fully chiral, CCW)"   font 'Helvetica,11'
    eval("plot " . plot_cmd_0)

    set title "ω = 0.5  (achiral)"             font 'Helvetica,11'
    eval("plot " . plot_cmd_5)

    set title "ω = 1.0  (fully chiral, CW)"    font 'Helvetica,11'
    eval("plot " . plot_cmd_1)

    unset multiplot
    unset output
    print "Wrote tcrw_fig2_occ.pdf"
}

#--------------------------------------------------------------------
# INTERACTIVE qt  (LAST — blocks on `pause mouse close`)
#--------------------------------------------------------------------
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 1500,520 enhanced font 'Helvetica,11'

    set multiplot layout 1,3 title \
        "TCRW Fig 2 occupancy P(X,Y)   |   L = 10,  D_r = 10^{-3},  T = 10^{10}" \
        font 'Helvetica,11'

    set title "ω = 0.0  (fully chiral, CCW)"
    eval("plot " . plot_cmd_0)

    set title "ω = 0.5  (achiral)"
    eval("plot " . plot_cmd_5)

    set title "ω = 1.0  (fully chiral, CW)"
    eval("plot " . plot_cmd_1)

    unset multiplot
    pause mouse close
}
