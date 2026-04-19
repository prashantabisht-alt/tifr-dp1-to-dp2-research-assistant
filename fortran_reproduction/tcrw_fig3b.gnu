#=====================================================================
# tcrw_fig3b.gnu — Fig 3(b): |J_Dr|_wall / |J_ω|_wall  vs  D_r
#                            (ω = 1, various L)
#
# Matches Osat et al. Fig 3(b):
#   - x-axis: D_r          log scale, ~[10^-4, 1]
#   - y-axis: wall-ratio   log scale, ~[10^-3, 10^2]
#   - one curve per L ∈ {4, 9, 19, 49},  markers + joining line
#   - curves overlap in the plateau region and all diverge together
#     as D_r → 1.
#
# Physics punchline (full derivation in tcrw_fig3b.f90 header):
#   Two sources of the ratio:
#     (a) the small-D_r PLATEAU is the real physics — the noise-
#         triggered current J_Dr on the left wall has no strong
#         directional preference, so even though the left wall carries
#         a large CCW edge current |J_ω|, there is still some leakage
#         to J_Dr whose magnitude grows linearly with the rate of
#         noise events, giving a D_r-independent ratio.
#     (b) the D_r → 1 DIVERGENCE is a BOOKKEEPING artifact — at high
#         D_r almost every translation is preceded by a noise step,
#         so by construction it gets attributed to J_Dr, and the ratio
#         blows up without any physical edge-current change.
#   Trivial counting gives ratio ≈ D_r / (1 - D_r) for the divergence
#   at D_r → 1 (dashed line); measured ratio tracks this slope.
#
# Reads:   tcrw_fig3b_summary.txt
#          columns:  L  D_r  ratio  |J_Dr|_wall  |J_ω|_wall  Jx_w_Dr  Jy_w_Dr  Jx_w_om  Jy_w_om
#          (we plot col 2 vs col 3, split by col 1 == L)
#
# Output:  tcrw_fig3b.pdf   and/or interactive qt window
#
# Usage:
#   bash run.sh tcrw-plot-fig3b          # PDF then interactive qt
#   bash run.sh tcrw-plot-fig3b-pdf      # PDF only (headless)
#   bash run.sh tcrw-plot-fig3b-qt       # interactive qt only
#
# Notes
# -----
#   - Palette: viridis-like 4-shade ramp, SAME RGBs as Fig 3(a) so
#     the two panels read as one panel pair visually.
#   - The dashed reference line r(D_r) = D_r/(1-D_r) is the "naive
#     counting" prediction — observed data should asymptote to it at
#     D_r → 1 but sit ABOVE it (on a log plot) in the small-D_r
#     plateau region (edge-current physics adds a constant offset).
#   - No error bars in the paper; we follow the paper.
#=====================================================================

if (!exists("mode")) mode = "both"

f = 'tcrw_fig3b_summary.txt'

# ---- styles:  4 L curves, SAME viridis ramp as Fig 3(a) ----
set style line 1 lc rgb '#440154' pt 7  ps 1.0 lw 1.8    # L = 4  (deep purple)
set style line 2 lc rgb '#3b528b' pt 5  ps 1.0 lw 1.8    # L = 9  (blue)
set style line 3 lc rgb '#21918c' pt 9  ps 1.0 lw 1.8    # L = 19 (teal)
set style line 4 lc rgb '#5ec962' pt 13 ps 1.0 lw 1.8    # L = 49 (green)
set style line 99 lc rgb '#bbbbbb' dt 2 lw 0.8            # D_r/(1-D_r) guide

# ---- axes ----
set logscale x 10
set logscale y 10
set xrange [5e-5 : 2.0]
set yrange [1e-3 : 1e2]
set xlabel 'D_r'                              font ',14'
set ylabel '|J_{D_r}|_{wall} / |J_{/Symbol w}|_{wall}'   font ',14'
set format x '10^{%T}'
set format y '10^{%T}'
set tics scale 0.8
set grid   lc rgb '#dddddd' lw 0.4
set border lw 1.0
set key    left top box opaque samplen 1.8 spacing 1.2 font ',11'
set title  "TCRW Fig 3(b) — |J_{D_r}|/|J_{/Symbol w}| on left wall vs D_r   (ω = 1, OBC)" \
           font 'Helvetica,11'

# Reference:  trivial-counting asymptote y = D_r / (1 - D_r)
# (valid as a large-D_r divergence; at small D_r it goes like D_r
#  itself, so the plateau of the data sitting above it is the real
#  edge-current physics.)
guide(x) = x / (1.0 - x)

# ---- common plot command (shared across qt and pdf) ----
# col 1 = L, col 2 = D_r, col 3 = ratio
plot_cmd = \
  "'" . f . "' u ($1==4  ? $2 : 1/0):($1==4  ? $3 : 1/0) w lp ls 1 title 'L = 4', "  . \
  "'" . f . "' u ($1==9  ? $2 : 1/0):($1==9  ? $3 : 1/0) w lp ls 2 title 'L = 9', "  . \
  "'" . f . "' u ($1==19 ? $2 : 1/0):($1==19 ? $3 : 1/0) w lp ls 3 title 'L = 19', " . \
  "'" . f . "' u ($1==49 ? $2 : 1/0):($1==49 ? $3 : 1/0) w lp ls 4 title 'L = 49', " . \
  "guide(x)                                             w l  ls 99 title 'D_r/(1-D_r)'"

#=====================================================================
# PDF  (headless — runs FIRST so no qt window is created during PDF)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,10cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig3b.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig3b.pdf"
}

#=====================================================================
# INTERACTIVE qt  (LAST — blocks on `pause mouse close`)
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 820,640 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
