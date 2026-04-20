#=====================================================================
# tcrw_fig3g.gnu — Fig 3(g): |J_Dr|_wall / |J_ω|_wall  vs  ω
#                            (D_r = 10^-3, various L)
#
# Matches Osat et al. Fig 3(g):
#   - x-axis: ω ∈ [0, 1]            (linear)
#   - y-axis: wall-ratio            (log, ~[10^-3, 10^1])
#   - one curve per L ∈ {10, 19, 49},  markers + joining line
#   - deep symmetric DIP at ω = 0.5;  peaks at ω = 0 and ω = 1.
#
# Physics punchline (full derivation in tcrw_fig3g.f90 header):
#   - At ω ≈ 0 or ω ≈ 1 (fully chiral) the noise-triggered edge
#     current |J_Dr|_wall is robust — each noise event drifts the
#     trapped walker by one lattice site along the wall (skipping
#     orbit).  The ratio ~ O(1).
#   - At ω = 0.5 (achiral) the walker is an ordinary persistent RW:
#     signed (Δx, Δy) sums at the wall have ZERO MEAN and both
#     |J_Dr| and |J_ω| drop by several decades, but J_Dr drops faster
#     (paper page 3: "|J_Dr| decreases several orders of magnitude
#     ... however, |J_ω| is almost constant").  Ratio → 10^-2 or so.
#   - Curve is SYMMETRIC about ω = 0.5 (CW ↔ CCW relabelling
#     invariance of the chain).
#
# Reads:   tcrw_fig3g_summary.txt
#          columns:  L  ω  ratio  |J_Dr|_wall  |J_ω|_wall  Jx_w_Dr  Jy_w_Dr  Jx_w_om  Jy_w_om
#          (we plot col 2 vs col 3, split by col 1 == L)
#
# Output:  tcrw_fig3g.pdf   and/or interactive qt window
#
# Usage:
#   bash run.sh tcrw-plot-fig3g          # PDF then interactive qt
#   bash run.sh tcrw-plot-fig3g-pdf      # PDF only (headless)
#   bash run.sh tcrw-plot-fig3g-qt       # interactive qt only
#
# Notes
# -----
#   - Palette: 3-shade viridis, SAME RGBs as Fig 3(f) (blue, teal,
#     green) so the two ω-sweep panels read as one pair.
#   - Faint vertical dashed guide at ω = 0.5 marks the achiral point
#     (where the V-shape has its minimum).
#   - MC scatter is strongest near ω = 0.5 (both numerator and
#     denominator are noise-level) and at ω = 0, 1 (wall-trap mixing
#     time is 1/D_r² = 10^6, so T_use = 10^8 means ~100 effective
#     samples per cell).  Factor-of-2 scatter is normal there.
#   - No error bars in the paper; we follow the paper.
#=====================================================================

if (!exists("mode")) mode = "both"

f = 'tcrw_fig3g_summary.txt'

# ---- styles:  3 L curves, SAME viridis as Fig 3(f) ----
set style line 1 lc rgb '#3b528b' pt 7  ps 1.0 lw 1.8    # L = 10  (blue)
set style line 2 lc rgb '#21918c' pt 5  ps 1.0 lw 1.8    # L = 19  (teal)
set style line 3 lc rgb '#5ec962' pt 9  ps 1.0 lw 1.8    # L = 49  (green)
set style line 99 lc rgb '#bbbbbb' dt 2 lw 0.8            # ω = 0.5 guide

# ---- axes ----
unset logscale x
set logscale y 10
set xrange [-0.02 : 1.02]
set yrange [1e-3 : 1e1]
set xlabel '{/Symbol w}'                                     font ',14'
set ylabel '|J_{D_r}|_{wall} / |J_{/Symbol w}|_{wall}'       font ',14'
# explicit y-tic labels (avoid the gnuplot '10^{%T}' quirk we hit
# in Fig 3(f) where every minor tic in a single decade shows the
# same floor(log10) label).
set ytics ( "10^{-3}" 0.001, "10^{-2}" 0.01, "10^{-1}" 0.1, \
            "10^{0}"  1.0,   "10^{1}"  10.0 )
set mytics 10
set xtics 0, 0.1, 1.0
set mxtics 2
set tics scale 0.8
set grid   lc rgb '#dddddd' lw 0.4
set border lw 1.0
set key    right top box opaque samplen 1.8 spacing 1.2 font ',11'
set title  "TCRW Fig 3(g) — |J_{D_r}|/|J_{/Symbol w}| on left wall vs ω   (D_r = 10^{-3}, OBC)" \
           font 'Helvetica,11'

# vertical guide at ω = 0.5 (achiral, V-minimum)
set arrow from 0.5, 1e-3 to 0.5, 1e1 nohead ls 99
set label "ω = 0.5  (achiral)" at 0.52, 3e-3 left textcolor rgb '#777777' font ',10'

# ---- common plot command (shared across qt and pdf) ----
# col 1 = L, col 2 = ω, col 3 = ratio
plot_cmd = \
  "'" . f . "' u ($1==9 ? $2 : 1/0):($1==9 ? $3 : 1/0) w lp ls 1 title 'L = 9', "   . \
  "'" . f . "' u ($1==19 ? $2 : 1/0):($1==19 ? $3 : 1/0) w lp ls 2 title 'L = 19', " . \
  "'" . f . "' u ($1==49 ? $2 : 1/0):($1==49 ? $3 : 1/0) w lp ls 3 title 'L = 49'"

#=====================================================================
# PDF  (headless — runs FIRST so no qt window is created during PDF)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,10cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig3g.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig3g.pdf"
}

#=====================================================================
# INTERACTIVE qt  (LAST — blocks on `pause mouse close`)
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 820,640 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
