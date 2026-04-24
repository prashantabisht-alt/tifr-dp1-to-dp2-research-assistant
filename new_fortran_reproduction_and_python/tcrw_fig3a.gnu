#=====================================================================
# tcrw_fig3a.gnu — Fig 3(a): P_edge/P_bulk vs D_r  (ω = 1, various L)
#
# Matches Osat et al. Fig 3(a):
#   - x-axis: D_r   log scale, ~[10^-4, 1]
#   - y-axis: per-site edge/bulk ratio,  log scale
#   - one curve per L ∈ {4, 9, 19, 49},  markers + joining line
#   - curves collapse onto a common envelope at small D_r, diverge
#     only near the small-L saturation cap (L = 4)
#
# Reads:   tcrw_fig3a_summary.txt
#          columns:  L  D_r  ratio  P_edge_norm  P_bulk_norm  n_edge  n_bulk
#          (we plot col 2 vs col 3, split by col 1 == L)
#
# Output:  tcrw_fig3a.pdf   and/or interactive qt window
#
# Usage:
#   bash run.sh tcrw-plot-fig3a          # PDF then interactive qt
#   bash run.sh tcrw-plot-fig3a-pdf      # PDF only (headless)
#   bash run.sh tcrw-plot-fig3a-qt       # interactive qt only
#
# Notes
# -----
#   - We filter rows by L via the `($1==L ? $x : 1/0)` trick; gnuplot
#     drops 1/0 points silently so each L curve only draws its own 25
#     points.
#   - Palette: viridis-like 4-shade ramp, deep blue (small L) → bright
#     yellow-green (large L), so the eye can follow the L → ∞ trend.
#   - No error bars in the paper; we follow the paper.
#=====================================================================

if (!exists("mode")) mode = "both"

f = 'tcrw_fig3a_summary.txt'

# ---- styles:  4 L curves, viridis-like ramp ----
set style line 1 lc rgb '#440154' pt 7  ps 1.0 lw 1.8    # L = 4  (deep purple)
set style line 2 lc rgb '#3b528b' pt 5  ps 1.0 lw 1.8    # L = 9  (blue)
set style line 3 lc rgb '#21918c' pt 9  ps 1.0 lw 1.8    # L = 19 (teal)
set style line 4 lc rgb '#5ec962' pt 13 ps 1.0 lw 1.8    # L = 49 (green)
set style line 99 lc rgb '#bbbbbb' dt 2 lw 0.8            # ratio=1 guide

# ---- axes ----
set logscale x 10
set logscale y 10
set xrange [5e-5 : 2.0]
set yrange [0.5  : 1e6]
set xlabel 'D_r'                      font ',14'
set ylabel 'P_{edge} / P_{bulk}'      font ',14'
set format x '10^{%T}'
set format y '10^{%T}'
set tics scale 0.8
set grid   lc rgb '#dddddd' lw 0.4
set border lw 1.0
set key    right top box opaque samplen 1.8 spacing 1.2 font ',11'
set title  "TCRW Fig 3(a) — edge localization vs D_r  (ω = 1, T = max(10^8, 100·max(L^2, 1/D_r)/D_r))" \
           font 'Helvetica,11'

# horizontal reference line at ratio = 1 (well-mixed)
set arrow from 5e-5, 1 to 2.0, 1 nohead ls 99

# ---- common plot command (shared across qt and pdf) ----
# col 1 = L, col 2 = D_r, col 3 = ratio
plot_cmd = \
  "'" . f . "' u ($1==4  ? $2 : 1/0):($1==4  ? $3 : 1/0) w lp ls 1 title 'L = 4', " .  \
  "'" . f . "' u ($1==9  ? $2 : 1/0):($1==9  ? $3 : 1/0) w lp ls 2 title 'L = 9', " .  \
  "'" . f . "' u ($1==19 ? $2 : 1/0):($1==19 ? $3 : 1/0) w lp ls 3 title 'L = 19', " . \
  "'" . f . "' u ($1==49 ? $2 : 1/0):($1==49 ? $3 : 1/0) w lp ls 4 title 'L = 49'"

#=====================================================================
# PDF  (headless — runs FIRST so no qt window is created during PDF)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,10cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig3a.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig3a.pdf"
}

#=====================================================================
# INTERACTIVE qt  (LAST — blocks on `pause mouse close`)
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 820,640 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
