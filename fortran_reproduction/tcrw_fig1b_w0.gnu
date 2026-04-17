#=====================================================================
# tcrw_fig1b_w0.gnu — single-panel preview: ω = 0.00, paper style
#
# Matches Osat et al. Fig 1(b) visual choices:
#   - INFERNO colormap  (black → purple → red → orange → yellow)
#   - LOG-SCALE colorbar  t ∈ [10⁰, 10⁶]
#     (this is the key — linear scale would squash all structure into
#      the late-time (yellow) end of the palette)
#   - trajectory drawn with lines, not points
#   - green disc  = Start  (t = 0)
#   - red disc    = End    (t = T)
#
# Reads:    tcrw_fig1b_traj_w0.00.txt      (columns: step  x  y)
#           "step" is the REAL step number (0, 100, 200, ..., 10⁶).
# Output:   qt window and/or tcrw_fig1b_w0.pdf
#=====================================================================

if (!exists("mode")) mode = "both"

# ---- count data rows ONCE up front; used by both qt and pdf blocks -
stats 'tcrw_fig1b_traj_w0.00.txt' using 0 nooutput
Nlast = STATS_records - 1          # index of the last data row

# ---- inferno palette (9-stop matplotlib reference) -----------------
set palette defined ( \
    0.000 '#000004', \
    0.125 '#1b0c41', \
    0.250 '#4a0c6b', \
    0.375 '#781c6d', \
    0.500 '#a52c60', \
    0.625 '#cf4446', \
    0.750 '#ed6925', \
    0.875 '#fb9a06', \
    1.000 '#fcffa4' )

# ---- log-scale colorbar, 10^0 .. 10^6 ------------------------------
set logscale cb
set cbrange  [1 : 1000000]
set cblabel  "t" font ',11'
set cbtics   ( "10^{0}" 1, "10^{1}" 10, "10^{2}" 100, \
               "10^{3}" 1000, "10^{4}" 10000, \
               "10^{5}" 100000, "10^{6}" 1000000 )

# ---- start / end markers -------------------------------------------
set style line 101 lc rgb '#2ca02c' pt 7 ps 1.6 lw 1.5   # green disc — Start
set style line 102 lc rgb '#d62728' pt 7 ps 1.6 lw 1.5   # red   disc — End

# ---- aesthetics ----------------------------------------------------
set grid   lc rgb '#e0e0e0' lw 0.4
set border lw 1.0
set key    right bottom reverse Left samplen 1.0 spacing 1.1 font ',10'

set size square
set xlabel "x"
set ylabel "y"
set title  "TCRW Fig 1(b),  ω = 0.00   (D_r = 10^{-3},  T = 10^{6})" \
           font 'Helvetica,12'

# The "(column(1) + 1)" expression shifts step 0 up to 1 so that log-cb
# has a finite color to assign to the very first point (log(0) is −∞).
# Downstream samples (step 100, 200, …) are unaffected at 4-sig-fig
# precision because 100+1 ≈ 100.

#=====================================================================
# INTERACTIVE qt
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 760,720 enhanced font 'Helvetica,11' persist
    plot \
      'tcrw_fig1b_traj_w0.00.txt' using 2:3:(column(1)+1) \
                                    w l palette lw 0.9 notitle, \
      ''                          u 2:3 every ::0::0 \
                                    w p ls 101 title 'Start', \
      ''                          u 2:3 every ::Nlast::Nlast \
                                    w p ls 102 title 'End'
}

#=====================================================================
# PDF
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 13cm,12cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig1b_w0.pdf'
    plot \
      'tcrw_fig1b_traj_w0.00.txt' using 2:3:(column(1)+1) \
                                    w l palette lw 0.8 notitle, \
      ''                          u 2:3 every ::0::0 \
                                    w p ls 101 title 'Start', \
      ''                          u 2:3 every ::Nlast::Nlast \
                                    w p ls 102 title 'End'
    unset output
    print "Wrote tcrw_fig1b_w0.pdf"
}
