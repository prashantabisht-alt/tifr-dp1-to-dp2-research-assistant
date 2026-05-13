# Panel (a) — corrected exact theory, P(x,y,t) on triangular lattice
# Usage:
#   gnuplot fig11_panel_a.gnu                  # PNG + PDF
#   gnuplot -e "mode='qt'"  fig11_panel_a.gnu  # interactive qt (mouse zoom)
#   gnuplot -e "mode='png'" fig11_panel_a.gnu  # PNG only
#   gnuplot -e "mode='pdf'" fig11_panel_a.gnu  # PDF only

if (!exists("mode")) mode = "both"

png_out = "fig11_panel_a.png"
pdf_out = "fig11_panel_a.pdf"
data_file = "outputs/fig11_final_hex_theory_points.txt"

set encoding utf8
set palette defined (0 '#050505', 0.22 '#32105f', 0.42 '#b42bd6', 0.60 '#ed1c24', 0.78 '#ff8c00', 1.00 '#fff200')
set cbrange [0:0.0018]
unset key
set colorbox
set size ratio -1
set xrange [-34:34]
set yrange [-39:39]
set xlabel "x"
set ylabel "y"
set title "(a) Corrected exact theory  P(x,y,t)   gamma=0.01, epsilon=0.15, t=50"

do_plot = "plot data_file using 1:2:3 with points pointtype 7 pointsize 1.25 palette notitle"

if (mode eq "png" || mode eq "both") {
    set terminal pngcairo enhanced color font "Helvetica,12" size 1200,1300
    set output png_out
    eval(do_plot)
    set output
    print "saved ".png_out
}

if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo enhanced color font "Helvetica,12" size 8in,8.5in
    set output pdf_out
    eval(do_plot)
    set output
    print "saved ".pdf_out
}

if (mode eq "qt") {
    set terminal qt enhanced font "Helvetica,12" size 950,1000 \
        title "Fig 11 panel (a) — corrected theory  [right-drag: zoom | u: unzoom | a: full]"
    eval(do_plot)
    pause mouse close
}
