#=====================================================================
# tcrw_fig1d.gnu — Fig 1(d): D(ω) for two D_r values (paper match)
#
# Matches Osat et al. Fig 1(d):
#   - x-axis: ω ∈ [0, 1]  linear
#   - y-axis: D         linear, range [0, ~0.3]
#   - multiple D_r curves collapse onto the same linear "tent"
#     → caption: "D decreases linearly with chirality ω independent
#                 of the value of D_r"
#   - markers + joining line, NO error bars (paper style)
#
# We plot two of the paper's three D_r values (10⁻³, 10⁻²); that's
# already enough to demonstrate the collapse (factor of 10 apart).
#
# Reads:   tcrw_fig1d_D_dr3.txt   (D_r = 10⁻³)
#          tcrw_fig1d_D_dr2.txt   (D_r = 10⁻²)
#          columns:  omega   D_fit   D_err   slope   slope_err   D_end
#   (we plot column 2, D_fit — the log-log regression estimate)
#
# Output:  tcrw_fig1d.pdf + interactive qt
#
# Usage:
#   bash run.sh tcrw-plot-fig1d          # PDF then interactive qt
#   bash run.sh tcrw-plot-fig1d-pdf      # PDF only (headless)
#   bash run.sh tcrw-plot-fig1d-qt       # interactive qt only
#
# Notes:
#   - The error bars stored in the datafile (columns 3, 5) are NOT
#     plotted, matching the paper. If you want to see them, use the
#     companion `tcrw_fig1d_logy.gnu` (not yet written) or just
#     replace `w lp` with `w yerrorlines u 1:2:3` below.
#=====================================================================

if (!exists("mode")) mode = "both"

f3 = 'tcrw_fig1d_D_dr3.txt'   # D_r = 10⁻³
f2 = 'tcrw_fig1d_D_dr2.txt'   # D_r = 10⁻²

# ---- styles (paper-like palette: two shades of red-orange) ----
set style line 1 lc rgb '#d6604d' pt 5 ps 1.2 lw 1.8   # filled square, D_r=1e-3
set style line 2 lc rgb '#8b0000' pt 7 ps 1.2 lw 1.8   # filled circle, D_r=1e-2
set style line 99 lc rgb '#bbbbbb' dt 2 lw 0.8         # ω=0.5 symmetry guide

# ---- axes ----
unset logscale
set xrange [-0.02 : 1.02]
set yrange [0    : 0.30]
set xlabel 'ω' font ',14'
set ylabel 'D'  font ',14'
set xtics 0, 0.2, 1.0
set ytics 0, 0.05, 0.3
set mxtics 2
set mytics 5
set tics scale 0.8
set grid   lc rgb '#dddddd' lw 0.4
set border lw 1.0
set key    right top box opaque samplen 1.8 spacing 1.2 font ',11'
set title  "TCRW Fig 1(d) — D(ω),  T = 10^6,  N_{traj} = 500" \
           font 'Helvetica,12'

# ω=0.5 symmetry axis (visual guide — paper also shows curves peaked here)
set arrow from 0.5, 0 to 0.5, 0.30 nohead ls 99

# ---- common plot command (shared across qt and pdf) ----
plot_cmd = \
  "'" . f3 . "' u 1:2 w lp ls 1 title 'D_r = 10^{-3}', " . \
  "'" . f2 . "' u 1:2 w lp ls 2 title 'D_r = 10^{-2}'"

#=====================================================================
# PDF  (headless, runs FIRST so no window is created during PDF pass)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,10cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig1d.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig1d.pdf"
}

#=====================================================================
# INTERACTIVE qt  (LAST — script blocks on `pause mouse close` so
# zoom / pan / toolbar stay responsive until you close the window)
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 820,640 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
