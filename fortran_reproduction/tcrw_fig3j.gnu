#=====================================================================
# tcrw_fig3j.gnu — Fig 3(j):  θ_JDr(ω) and θ_Jω(ω) of the TOTAL
#                             left-wall current at fixed D_r = 10^-3
#                             (L = 10)
#
# ω-sweep twin of Fig 3(e).  Matches Osat et al. Fig 3(j):
#   - x-axis: ω   LINEAR scale, [0, 1]
#   - y-axis: θ   ∈ [-π/2, π/2],  custom ticks at 0, ±π/4, ±π/2
#   - two curves:
#       θ_JDr(ω)  sigmoid  +π/2  → 0 → -π/2
#       θ_Jω(ω)   sigmoid  -π/4  → 0 → +π/4
#     crossing zero together at ω = 0.5 (achirality).
#
# Physical explanation (derived, matches paper)
# ---------------------------------------------
#   At ω = 1: spatial orbit CCW → wall escape DOWN (Jy_Dr<0, θ≈-π/2)
#             and chiral continuation departs UP-RIGHT (θ_Jω≈+π/4).
#   At ω = 0: spatial orbit  CW → wall escape   UP (Jy_Dr>0, θ≈+π/2)
#             and chiral continuation departs DOWN-RIGHT (θ_Jω≈-π/4).
#   At ω = ½: both signed components vanish on average; Jx still
#             positive (noise-scatter bias into the bulk), so θ ≈ 0.
#
#   Note: even at ω = 0 and ω = 1, |θ_JDr| saturates slightly below
#   π/2 (paper shows ~1.46 rad ≈ 83.6°): there's always a small
#   positive Jx baseline from noise scatter, so atan2(Jy, Jx) can't
#   reach exactly ±π/2.  θ_Jω, by contrast, is bounded exactly at
#   ±π/4 by construction of the chiral-orbit geometry.
#
# Reads:  tcrw_fig3j_summary.txt
#         columns:  iW  ω  Jx_Dr_tot  Jy_Dr_tot  Jx_om_tot  Jy_om_tot
#                   θ_JDr_tot  θ_Jω_tot
#
# Output: tcrw_fig3j.pdf   and/or interactive qt window
#=====================================================================

if (!exists("mode")) mode = "both"

f = 'tcrw_fig3j_summary.txt'

# ---- styles ----
set style line 1 lc rgb '#440154' pt 7 ps 0.9 lw 1.8    # θ_JDr  (deep purple)
set style line 2 lc rgb '#5ec962' pt 5 ps 0.9 lw 1.8    # θ_Jω   (green)
set style line 99 lc rgb '#bbbbbb' dt 2 lw 0.8          # guides at 0, ±π/4, ±π/2

# ---- axes (linear ω) ----
set xrange [-0.02 : 1.02]
set yrange [-1.7 : 1.7]
set xtics 0, 0.1, 1.0
set mxtics 2
set xlabel '{/Symbol w}'                    font ',14'
set ylabel '{/Symbol q}  (radians)'         font ',14'
# custom y-tics at multiples of π/4
set ytics ('-π/2' -1.5708, '-π/4' -0.7854, '0' 0.0, 'π/4' 0.7854, 'π/2' 1.5708)
set mytics 2
set tics scale 0.8
set grid   lc rgb '#dddddd' lw 0.4
set border lw 1.0
set key    right top box opaque samplen 1.8 spacing 1.2 font ',11'
set title  "TCRW Fig 3(j) — θ of total left-wall currents vs ω  (L = 10, D_r = 10^{-3})" \
           font 'Helvetica,11'

# horizontal guides at 0, ±π/4, ±π/2
set arrow from -0.02,  0.0000 to 1.02,  0.0000 nohead ls 99
set arrow from -0.02, -1.5708 to 1.02, -1.5708 nohead ls 99
set arrow from -0.02,  1.5708 to 1.02,  1.5708 nohead ls 99
set arrow from -0.02, -0.7854 to 1.02, -0.7854 nohead ls 99
set arrow from -0.02,  0.7854 to 1.02,  0.7854 nohead ls 99
# vertical guide at ω = 0.5 (achirality)
set arrow from 0.5, -1.7 to 0.5, 1.7 nohead ls 99

# ---- plot command ----
#   col 2 = ω, col 7 = θ_JDr_tot, col 8 = θ_Jω_tot
plot_cmd = \
  "'" . f . "' u 2:7 w lp ls 1 title 'θ_{J_{D_r}}', " . \
  "'" . f . "' u 2:8 w lp ls 2 title 'θ_{J_{/Symbol w}}'"

#=====================================================================
# PDF  (runs FIRST)
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 14cm,10cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig3j.pdf'
    eval("plot " . plot_cmd)
    unset output
    print "Wrote tcrw_fig3j.pdf"
}

#=====================================================================
# INTERACTIVE qt
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 820,640 enhanced font 'Helvetica,11'
    eval("plot " . plot_cmd)
    pause mouse close
}
