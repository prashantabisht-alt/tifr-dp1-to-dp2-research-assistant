#=====================================================================
# tcrw_fig3f.gnu — Fig 3(f): P_edge/P_bulk vs ω  at fixed D_r = 10^-3
#
# Matches Osat et al. Fig 3(f):
#   - x-axis: ω ∈ [0, 1]                (linear)
#   - y-axis: per-site edge/bulk ratio  (log; values of order 10^2–10^3)
#   - one curve per L ∈ {10, 19, 49}
#   - markers + joining line, viridis-like palette consistent with 3(a)
#
# Physics punchline (see tcrw_fig3f.f90 header):
#   The ratio is APPROXIMATELY INDEPENDENT of ω.  Edge localization is
#   a topological band-structure effect (non-trivial Chern number), not
#   a chirality effect.  A dashed horizontal reference at r_ref = 700
#   marks the expected Fig 3(a) master-curve value at D_r = 10^-3;
#   all three L curves should hover near that line with mild MC noise.
#
# Reads:   tcrw_fig3f_summary.txt
#          columns:  L   ω   ratio   P_edge_norm   P_bulk_norm   n_edge   n_bulk
#          we plot  col 2  vs  col 3, filtered by col 1 == L.
#
# Output:  tcrw_fig3f.pdf   and/or interactive qt window
#
# Usage:
#   bash run.sh tcrw-plot-fig3f          # PDF then interactive qt
#   bash run.sh tcrw-plot-fig3f-pdf      # PDF only (headless)
#   bash run.sh tcrw-plot-fig3f-qt       # interactive qt only
#
# Notes
# -----
#   - Viridis ramp uses only 3 shades (L = 10, 19, 49), drawn from the
#     same palette as Fig 3(a) for visual continuity.
#   - Reference line at r_ref = 700 is an estimate; if the Fig 3(a)
#     summary file is present, one can read its D_r = 10^-3 row to
#     cross-check, but we do not require it here.
#   - No error bars in the paper; we follow the paper.
#=====================================================================

if (!exists("mode")) mode = "both"

f = 'tcrw_fig3f_summary.txt'

# expected Fig 3(a) master-curve value at D_r = 10^-3  (horizontal guide)
# ratio at ω = 1 for L ∈ {19, 49} from tcrw_fig3a_summary.txt is 200-300;
# we put the guideline at the mid-range ≈ 260 so a flat Fig 3(f) curve
# sits right on it.
r_ref = 260.0

# ---- styles:  3 L curves, viridis-like ramp ----
set style line 1 lc rgb '#3b528b' pt 7  ps 1.0 lw 1.8    # L = 10  (blue)
set style line 2 lc rgb '#21918c' pt 5  ps 1.0 lw 1.8    # L = 19  (teal)
set style line 3 lc rgb '#5ec962' pt 9  ps 1.0 lw 1.8    # L = 49  (green)
set style line 99 lc rgb '#bbbbbb' dt 2 lw 0.8            # r_ref guide

# ---- axes ----
unset logscale x
set logscale y 10
set xrange [-0.02 : 1.02]
set yrange [100 : 1000]
set xlabel '{/Symbol w}'                   font ',14'
set ylabel 'P_{edge} / P_{bulk}'           font ',14'
# explicit ytic positions with numeric labels — avoids the gnuplot quirk
# where `set format y "10^{%T}"` labels every minor tick in [10^2, 10^3]
# as "10^2" because %T = floor(log10(y)).
set ytics ( "100" 100, "150" 150, "200" 200, "300" 300, "500" 500, \
            "700" 700, "1000" 1000 )
set mytics 5
set xtics 0, 0.1, 1.0
set mxtics 2
set tics scale 0.8
set grid   lc rgb '#dddddd' lw 0.4
set border lw 1.0
set key    right bottom box opaque samplen 1.8 spacing 1.2 font ',11'
set title  "TCRW Fig 3(f) — edge/bulk ratio vs ω   (D_r = 10^{-3}, flat ⇒ topological origin)" \
           font 'Helvetica,11'

# horizontal reference at Fig 3(a) master-curve value
set arrow from -0.02, r_ref to 1.02, r_ref nohead ls 99
set label sprintf("r_{ref} ≈ %.0f  (Fig 3a master curve at D_r = 10^{-3})", r_ref) \
    at 0.05, r_ref*1.10 left textcolor rgb '#777777' font ',10'

# ---- common plot command (shared across qt and pdf) ----
# col 1 = L, col 2 = ω, col 3 = ratio
plot_cmd = \
  "'" . f . "' u ($1==10 ? $2 : 1/0):($1==10 ? $3 : 1/0) w lp ls 1 title 'L = 10', "  . \
  "'" . f . "' u ($1==19 ? $2 : 1/0):($1==19 ? $3 : 1/0) w lp ls 2 title 'L = 19', "  . \
  "'" . f . "' u ($1==49 ? $2 : 1/0):($1==49 ? $3 : 1/0) w lp ls 3 title 'L = 49'"

#=====================================================================
# PDF  (headless — runs FIRST so no qt window is created during PDF)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,10cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig3f.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig3f.pdf"
}

#=====================================================================
# INTERACTIVE qt  (LAST — blocks on `pause mouse close`)
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 820,640 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
