#=====================================================================
# tcrw_fig3e.gnu — Fig 3(e):  θ_JDr of the TOTAL left-wall current
#                             vs  D_r   (L = 10, ω = 1)
#
# Matches Osat et al. Fig 3(e):
#   - x-axis: D_r   log scale, [10^-4, 10^0]
#   - y-axis: θ_JDr ∈ [-π/2, 0]
#   - one curve: θ of the wall-summed J_Dr vector
#
# Paper claim (page 3 + caption):
#   "The angle of the current J_Dr along the left edge goes from
#    −π/2 to zero as we vary D_r from zero to one."
#   → monotone sigmoid-like curve from −π/2 at D_r → 0+ to ~0 at
#    D_r → 1.  The transition happens where the chiral-orbit time
#    (1/D_r) crosses the bulk-diffusion time, roughly D_r ≈ 1/L.
#
# We also plot θ_Jω (same data, for Fig 3(j) completeness / sanity —
# paper says θ_Jω is flat near +π/4 across the full D_r range).
#
# Reads:  tcrw_fig3e_summary.txt
#         columns:  iD  D_r  Jx_Dr_tot  Jy_Dr_tot  Jx_om_tot  Jy_om_tot  θ_Dr_tot  θ_ω_tot
#
# Output: tcrw_fig3e.pdf   and/or interactive qt window
#=====================================================================

if (!exists("mode")) mode = "both"

f = 'tcrw_fig3e_summary.txt'

# ---- styles ----
set style line 1 lc rgb '#440154' pt 7 ps 0.9 lw 1.8    # θ_JDr  (deep purple)
set style line 2 lc rgb '#5ec962' pt 5 ps 0.9 lw 1.8    # θ_Jω   (green, reference)
set style line 99 lc rgb '#bbbbbb' dt 2 lw 0.8           # guides at 0, ±π/2, ±π/4

# ---- axes ----
set logscale x 10
set xrange [5e-5 : 2.0]
set yrange [-1.7 : 1.2]
set xlabel 'D_r'                              font ',14'
set ylabel '{/Symbol q}  (radians)'           font ',14'
set format x '10^{%T}'
# custom y-tics at multiples of π/4
set ytics ('-π/2' -1.5708, '-π/4' -0.7854, '0' 0.0, 'π/4' 0.7854, 'π/2' 1.5708)
set mytics 2
set tics scale 0.8
set grid   lc rgb '#dddddd' lw 0.4
set border lw 1.0
set key    right top box opaque samplen 1.8 spacing 1.2 font ',11'
set title  "TCRW Fig 3(e) — θ of total left-wall currents vs D_r  (L = 10, ω = 1)" \
           font 'Helvetica,11'

# reference horizontal guides at 0, ±π/4, ±π/2
set arrow from 5e-5,  0.0000 to 2.0,  0.0000 nohead ls 99
set arrow from 5e-5, -1.5708 to 2.0, -1.5708 nohead ls 99
set arrow from 5e-5,  0.7854 to 2.0,  0.7854 nohead ls 99

# ---- plot command ----
plot_cmd = \
  "'" . f . "' u 2:7 w lp ls 1 title 'θ_{J_{D_r}}', " . \
  "'" . f . "' u 2:8 w lp ls 2 title 'θ_{J_{/Symbol w}}  (ref: +π/4)'"

#=====================================================================
# PDF  (runs FIRST)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,10cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig3e.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig3e.pdf"
}

#=====================================================================
# INTERACTIVE qt
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 820,640 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
