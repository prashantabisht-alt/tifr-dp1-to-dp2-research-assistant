# Final gnuplot version of Fig. 11.
#
# Run from the triangular folder:
#   python3 export_fig11_final_hex_gnuplot_data.py
#   gnuplot fig11_final_hex.gnu
#
# Interactive only:
#   gnuplot -persist -e "mode='qt'" fig11_final_hex.gnu
#
# PNG only:
#   gnuplot -e "mode='png'" fig11_final_hex.gnu
#
# PDF only:
#   gnuplot -e "mode='pdf'" fig11_final_hex.gnu

if (!exists("mode")) mode = "both"

png_out = "fig11_final_hex_gnuplot.png"
pdf_out = "fig11_final_hex_gnuplot.pdf"

data_dir = "outputs"
theory_file = data_dir . "/fig11_final_hex_theory_points.txt"
kmc_file = data_dir . "/fig11_final_hex_kmc_points.txt"
cross_file = data_dir . "/fig11_final_hex_cross_section.txt"

common_setup = "set encoding utf8; set palette defined (0 '#050505', 0.22 '#32105f', 0.42 '#b42bd6', 0.60 '#ed1c24', 0.78 '#ff8c00', 1.00 '#fff200'); set cbrange [0:0.0018]; set border linewidth 1.0; set tics scale 0.65; set tics font ',14'; set xlabel font ',15'; set ylabel font ',15'; set title font ',16'"

do_plot = " \
@common_setup; \
set multiplot title 'Triangular active walker: corrected theory matches KMC' font ',18'; \
set label 99 'Minimum-image triangular display.  MC noise floor approx. 4.22e-6; corrected theory is at the noise floor.  The old c_3 sign gives the red curve.' at screen 0.5,0.035 center font ',13'; \
unset key; set colorbox; set size ratio -1; set xrange [-28:28]; set yrange [-25:25]; set xtics 10; set ytics 10; set cbtics 0.0004; set format cb '%.4f'; \
set xlabel 'x'; set ylabel 'y'; \
set origin 0.065,0.535; set size 0.405,0.385; \
set title '(a) Corrected exact theory  (gamma=0.01, epsilon=0.15, t=50)'; \
plot theory_file using 1:2:(0.78):3 with circles lc palette fill solid 1.0 noborder notitle; \
set origin 0.555,0.535; set size 0.405,0.385; \
set title '(b) Kinetic Monte Carlo  (N=100,000,000)'; \
plot kmc_file using 1:2:(0.78):3 with circles lc palette fill solid 1.0 noborder notitle; \
unset colorbox; set size noratio; set origin 0.08,0.115; set size 0.84,0.32; \
set autoscale; set xrange [-33:33]; set yrange [0.0001:0.00185]; set ytics 0.00025; set xtics autofreq; set format y '%.5f'; \
set xlabel 'x at y=0 (minimum-image coordinates)'; set ylabel 'P(x,y=0,t)'; \
set title '(c) Horizontal cross-section through the walker'; set grid linecolor rgb '#dddddd'; set key top center horizontal font ',13' samplen 2.0; \
plot cross_file using 1:2 with points pointtype 7 pointsize 1.15 linecolor rgb 'black' title 'KMC', \
     cross_file using 1:3 with lines linewidth 2.5 linecolor rgb '#d62728' title 'buggy theory (RMS=1.43e-4)', \
     cross_file using 1:4 with lines linewidth 2.5 linecolor rgb '#0057ff' title 'corrected theory (RMS=3.43e-6)'; \
unset label 99; unset multiplot"

if (mode eq "png" || mode eq "both") {
    set terminal pngcairo enhanced color font "Helvetica,16" size 2686,1901
    set output png_out
    eval(do_plot)
    set output
}

if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo enhanced color font "Helvetica,16" size 13.5in,9.55in
    set output pdf_out
    eval(do_plot)
    set output
}

if (mode eq "qt") {
    set terminal qt enhanced font "Helvetica,14" size 1400,990
    eval(do_plot)
    pause mouse close
}
