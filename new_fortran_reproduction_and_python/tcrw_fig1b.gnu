#=====================================================================
# tcrw_fig1b.gnu — Fig 1(b): sample trajectories, paper style
#
# Matches Osat et al. Fig 1(b) (same treatment as single-panel preview
# tcrw_fig1b_w0.gnu, tiled 2×2 across the 4 chiralities):
#   - INFERNO colormap  (black → purple → red → orange → yellow)
#   - LOG-SCALE colorbar  t ∈ [10⁰, 10⁶]
#     (linear scale squashes all structure into the yellow end — the
#      log scale is what lets the early CCW/CW loops stay visible.)
#   - trajectories drawn with lines (lw 0.9, same as the approved
#     tcrw_fig1b_w0.gnu preview). For full-resolution skip=1 runs with
#     10⁶ segments per panel, drop LW to 0.2–0.3 to avoid blobbing.
#   - green disc  = Start  (t = 0)
#   - red   disc  = End    (t = T = 10⁶)
#
# Layout:   2×2     ω = 0.00   top-left      ω = 0.50   top-right
#                   ω = 0.75   bottom-left   ω = 1.00   bottom-right
#
# Reads:    tcrw_fig1b_traj_w{0.00,0.50,0.75,1.00}.txt
#           columns: step  x  y    (subsampled trajectory records)
# Output:   qt window, tcrw_fig1b.pdf, and/or tcrw_fig1b.png
#
# Usage:
#   gnuplot -persist tcrw_fig1b.gnu                 # qt + PDF (default)
#   gnuplot -e "mode='qt'"  -persist tcrw_fig1b.gnu # interactive only
#   gnuplot -e "mode='pdf'"         tcrw_fig1b.gnu  # PDF only (headless)
#   gnuplot -e "mode='png'"         tcrw_fig1b.gnu  # PNG only
#=====================================================================

if (!exists("mode")) mode = "both"

# ---- one data file per chirality -----------------------------------
f00 = 'tcrw_fig1b_traj_w0.00.txt'
f05 = 'tcrw_fig1b_traj_w0.50.txt'
f07 = 'tcrw_fig1b_traj_w0.75.txt'
f10 = 'tcrw_fig1b_traj_w1.00.txt'

# ---- count rows in each file (needed for 'End' marker) -------------
# Files are the same length for matched T and skip, but we stat each
# one so the script stays robust if someone tweaks a single ω.
stats f00 using 0 nooutput;  N00 = STATS_records - 1
stats f05 using 0 nooutput;  N05 = STATS_records - 1
stats f07 using 0 nooutput;  N07 = STATS_records - 1
stats f10 using 0 nooutput;  N10 = STATS_records - 1

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

# ---- start / end markers (matches tcrw_fig1b_w0.gnu) ---------------
set style line 101 lc rgb '#2ca02c' pt 7 ps 1.6 lw 1.5   # green disc — Start
set style line 102 lc rgb '#d62728' pt 7 ps 1.6 lw 1.5   # red   disc — End

# ---- aesthetics ----------------------------------------------------
set grid   lc rgb '#e8e8e8' lw 0.3
set border lw 1.0
set tics   scale 0.7
set key    off

# Line width sized for full-resolution skip=1 runs (10⁶ segments per
# panel). At lw 0.9 the small ω=0 / ω=1 cages overdraw into a solid
# blob; lw 0.3 keeps the woven unit-square mesh visible.
# For the earlier skip=100 data (~10⁴ segments), bump this to LW = 0.9.
LW = 0.3

# The `(column(1) + 1)` shifts step 0 up to 1 so that log-cb has a
# finite color for the very first point (log(0) is −∞). Downstream
# samples are unaffected at 4 sig-fig precision.

#=====================================================================
# INTERACTIVE qt
#=====================================================================
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 1040,1040 enhanced font 'Helvetica,11' persist
    set multiplot layout 2,2 \
        title "TCRW Fig 1(b) — sample trajectories,  D_r = 10^{-3},  T = 10^{6}" \
        font 'Helvetica,13'

    set size square
    set xlabel 'x'; set ylabel 'y'
    unset colorbox                      # suppress per-panel colorbars

    # ---- ω = 0.00 (top-left) ----
    set title 'ω = 0.00'
    plot f00 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 notitle,   \
         '< tail -n 1 '.f00               u 2:3 w p ls 102 notitle

    # ---- ω = 0.50 (top-right) ----
    set title 'ω = 0.50'
    plot f05 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 notitle,   \
         '< tail -n 1 '.f05               u 2:3 w p ls 102 notitle

    # ---- ω = 0.75 (bottom-left) ----
    set title 'ω = 0.75'
    plot f07 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 notitle,   \
         '< tail -n 1 '.f07               u 2:3 w p ls 102 notitle

    # ---- ω = 1.00 (bottom-right) — also carry the shared colorbar --
    set colorbox vertical user origin 0.93,0.10 size 0.018,0.80 \
        front noinvert
    set title 'ω = 1.00'
    plot f10 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 title 'Start', \
         '< tail -n 1 '.f10               u 2:3 w p ls 102 title 'End'

    unset multiplot
}

#=====================================================================
# PDF
#=====================================================================
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 22cm,22cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig1b.pdf'
    set multiplot layout 2,2 \
        title "TCRW Fig 1(b) — sample trajectories,  D_r = 10^{-3},  T = 10^{6}" \
        font 'Helvetica,12'

    set size square
    set xlabel 'x'; set ylabel 'y'
    unset colorbox

    set title 'ω = 0.00'
    plot f00 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 notitle,   \
         '< tail -n 1 '.f00               u 2:3 w p ls 102 notitle

    set title 'ω = 0.50'
    plot f05 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 notitle,   \
         '< tail -n 1 '.f05               u 2:3 w p ls 102 notitle

    set title 'ω = 0.75'
    plot f07 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 notitle,   \
         '< tail -n 1 '.f07               u 2:3 w p ls 102 notitle

    set colorbox vertical user origin 0.93,0.10 size 0.018,0.80 \
        front noinvert
    set title 'ω = 1.00'
    plot f10 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 title 'Start', \
         '< tail -n 1 '.f10               u 2:3 w p ls 102 title 'End'

    unset multiplot
    unset output
    print "Wrote tcrw_fig1b.pdf"
}

#=====================================================================
# PNG
#=====================================================================
if (mode eq "png") {
    set terminal pngcairo size 1800,1800 enhanced font 'Helvetica,11'
    set output 'tcrw_fig1b.png'
    set multiplot layout 2,2 \
        title "TCRW Fig 1(b) — sample trajectories,  D_r = 10^{-3},  T = 10^{6}" \
        font 'Helvetica,14'

    set size square
    set xlabel 'x'; set ylabel 'y'
    unset colorbox

    set title 'ω = 0.00'
    plot f00 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 notitle,   \
         '< tail -n 1 '.f00               u 2:3 w p ls 102 notitle

    set title 'ω = 0.50'
    plot f05 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 notitle,   \
         '< tail -n 1 '.f05               u 2:3 w p ls 102 notitle

    set title 'ω = 0.75'
    plot f07 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 notitle,   \
         '< tail -n 1 '.f07               u 2:3 w p ls 102 notitle

    set colorbox vertical user origin 0.93,0.10 size 0.018,0.80 \
        front noinvert
    set title 'ω = 1.00'
    plot f10 using 2:3:(column(1)+1) w l palette lw LW notitle, \
         ''  u 2:3 every ::0::0           w p ls 101 title 'Start', \
         '< tail -n 1 '.f10               u 2:3 w p ls 102 title 'End'

    unset multiplot
    unset output
    print "Wrote tcrw_fig1b.png"
}
