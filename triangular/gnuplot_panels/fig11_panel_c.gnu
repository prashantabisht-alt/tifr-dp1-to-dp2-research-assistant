# Panel (c) — Horizontal cross-section through the walker at y=0
# Usage:
#   gnuplot fig11_panel_c.gnu                  # PNG + PDF
#   gnuplot -e "mode='qt'"  fig11_panel_c.gnu  # interactive qt (mouse zoom)
#   gnuplot -e "mode='png'" fig11_panel_c.gnu  # PNG only
#   gnuplot -e "mode='pdf'" fig11_panel_c.gnu  # PDF only

if (!exists("mode")) mode = "both"

png_out = "fig11_panel_c.png"
pdf_out = "fig11_panel_c.pdf"
data_file = "outputs/fig11_final_hex_cross_section.txt"

set encoding utf8
set grid
set key top center horizontal
set xlabel "x at y=0 (minimum-image coordinates)"
set ylabel "P(x, y=0, t)"
set xrange [-33:33]
set title "(c) Horizontal cross-section   gamma=0.01, epsilon=0.15, t=50, N=1e8"

do_plot = "plot data_file using 1:2 with points pointtype 7 pointsize 0.9 linecolor rgb 'black' title 'KMC',     '' using 1:3 with lines linewidth 2.5 linecolor rgb '#d62728' title 'buggy theory (RMS=1.43e-4)',     '' using 1:4 with lines linewidth 2.5 linecolor rgb '#0057ff' title 'corrected theory (RMS=3.43e-6)'"

if (mode eq "png" || mode eq "both") {
    set terminal pngcairo enhanced color font "Helvetica,12" size 1600,900
    set output png_out
    eval(do_plot)
    set output
    print "saved ".png_out
}

if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo enhanced color font "Helvetica,12" size 11in,6in
    set output pdf_out
    eval(do_plot)
    set output
    print "saved ".pdf_out
}

if (mode eq "qt") {
    set terminal qt enhanced font "Helvetica,12" size 1500,850 \
        title "Fig 11 panel (c) — cross-section  [right-drag: zoom | u: unzoom | a: full | L: log toggle]"
    eval(do_plot)
    pause mouse close
}
