#=====================================================================
# tcrw_fig1c.gnu — Fig 1(c): MSD(t) on log–log axes
#
# Matches Osat et al. Fig 1(c):
#   - x-axis: t,   log scale  [10^0, 10^6]
#   - y-axis: MSD, log scale  [10^0, 10^6]
#   - dashed reference line:  MSD = t  (slope 1 → normal diffusion)
#   - 4 curves for ω ∈ {0.5, 0.7, 0.9, 1.0}
#
# Reads:   tcrw_fig1c_msd_w0.5.txt
#          tcrw_fig1c_msd_w0.7.txt
#          tcrw_fig1c_msd_w0.9.txt
#          tcrw_fig1c_msd_w1.0.txt    (columns: t   <|r|^2>)
# Output:  interactive qt window and/or tcrw_fig1c.pdf
#
# Usage:
#   bash run.sh tcrw-plot-fig1c          # PDF then interactive qt
#   bash run.sh tcrw-plot-fig1c-pdf      # PDF only (headless)
#   bash run.sh tcrw-plot-fig1c-qt       # interactive qt only (fastest feedback)
#=====================================================================

if (!exists("mode")) mode = "both"

# ---- one data file per chirality -----------------------------------
f05 = 'tcrw_fig1c_msd_w0.5.txt'
f07 = 'tcrw_fig1c_msd_w0.7.txt'
f09 = 'tcrw_fig1c_msd_w0.9.txt'
f10 = 'tcrw_fig1c_msd_w1.0.txt'

# ---- color choices roughly mimic paper Fig 1(c) ----
# (the paper uses a qualitative 4-color scheme — teal / red / blue / orange)
set style line 1 lc rgb '#2ca02c' pt  7 ps 0.8 lw 1.8   # green, ω = 0.5
set style line 2 lc rgb '#d62728' pt  5 ps 0.8 lw 1.8   # red,   ω = 0.7
set style line 3 lc rgb '#1f77b4' pt  9 ps 0.9 lw 1.8   # blue,  ω = 0.9
set style line 4 lc rgb '#ff7f0e' pt 11 ps 0.9 lw 1.8   # orange,ω = 1.0
set style line 99 lc rgb '#555555' dt 2       lw 1.3   # dashed reference

# ---- axes --------------------------------------------------------------
set logscale xy
set xrange [1    : 1e6]
set yrange [1e-1 : 1e7]
set xlabel 't'  font ',12'
set ylabel 'MSD ⟨|r(t)|²⟩' font ',12'
set format x "10^{%L}"
set format y "10^{%L}"
set tics   scale 0.8
set grid   lc rgb '#dddddd' lw 0.4
set border lw 1.0
set key    left top box opaque samplen 1.8 spacing 1.2 font ',10'
set title  "TCRW Fig 1(c) — MSD vs t,  D_r = 10^{-3},  N_{traj} = 1000" \
           font 'Helvetica,12'

# ---- common plot command as a string so qt and pdf don't diverge -----
plot_cmd = \
  "x  w l ls 99 title 'MSD ∝ t (slope 1)', " . \
  "f05 u 1:2 w lp ls 1 title 'ω = 0.5', " . \
  "f07 u 1:2 w lp ls 2 title 'ω = 0.7', " . \
  "f09 u 1:2 w lp ls 3 title 'ω = 0.9', " . \
  "f10 u 1:2 w lp ls 4 title 'ω = 1.0'"

#=====================================================================
# PDF  (headless — run FIRST so no window is opened during the PDF pass)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,12cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig1c.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig1c.pdf"
}

#=====================================================================
# INTERACTIVE qt  (LAST — so the qt terminal is the live one at script
# end, and `pause mouse close` holds the window until you close it,
# keeping zoom / pan / toolbar responsive on macOS Homebrew builds
# where -persist alone is flaky after a mid-script terminal switch.)
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 820,720 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
