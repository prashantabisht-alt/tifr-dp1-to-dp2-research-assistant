# Original-Confinement-style gnuplot Fig. 11.
#
# Build data:
#   python3 export_fig11_original_style_gnuplot_data.py
#
# If KMC data is present:
#   gfortran -O3 -fno-range-check -ffree-line-length-none kmc_triangular_jmvr_L60.f90 -o kmc_triangular_L60
#   ./kmc_triangular_L60 > kmc_L60_run.log
#   python3 export_fig11_original_style_gnuplot_data.py
#
# Plot:
#   gnuplot fig11_original_style.gnu

if (!exists("mode")) mode = "both"
if (!exists("cmax")) cmax = 0.001

png_out = "fig11_original_style_gnuplot.png"
pdf_out = "fig11_original_style_gnuplot.pdf"

data_dir = "outputs"
theory_file = data_dir . "/fig11_original_style_theory_points.txt"
kmc_file = data_dir . "/fig11_original_style_kmc_points.txt"
cross_file = data_dir . "/fig11_original_style_cross_section.txt"

if (!exists("has_kmc")) has_kmc = system(sprintf("test -f '%s' && echo 1 || echo 0", kmc_file)) + 0
panel_b_file = has_kmc ? kmc_file : theory_file
panel_b_title = has_kmc ? "(b) KMC: gamma=0.01, epsilon=0.15, t=50" : "(b) corrected theory again: no L60 KMC file"

common_setup = "set encoding utf8; set palette defined (0 '#050505', 0.22 '#32105f', 0.42 '#b42bd6', 0.60 '#ed1c24', 0.78 '#ff8c00', 1.00 '#fff200'); set cbrange [0:cmax]; set border linewidth 1.0; set tics scale 0.65; set tics font ',14'; set xlabel font ',16'; set ylabel font ',16'; set title font ',16'"

cross_plot_with_kmc = "plot cross_file using 1:2 with points pointtype 7 pointsize 0.75 linecolor rgb 'black' title 'KMC', cross_file using 1:3 with lines dashtype 2 linewidth 2.0 linecolor rgb '#7e2f8e' title 'full original buggy', cross_file using 1:4 with lines linewidth 2.0 linecolor rgb '#d62728' title 'old c_3 sign', cross_file using 1:5 with lines linewidth 2.0 linecolor rgb '#0057ff' title 'corrected'"
cross_plot_no_kmc = "plot cross_file using 1:2 with lines dashtype 2 linewidth 2.0 linecolor rgb '#7e2f8e' title 'full original buggy', cross_file using 1:3 with lines linewidth 2.0 linecolor rgb '#d62728' title 'old c_3 sign', cross_file using 1:4 with lines linewidth 2.0 linecolor rgb '#0057ff' title 'corrected'"
cross_plot = has_kmc ? cross_plot_with_kmc : cross_plot_no_kmc

do_plot = " \
@common_setup; \
set multiplot title 'Original-style Fig. 11 reproduction, L=60' font ',18'; \
unset key; set colorbox; set size ratio -1; set xrange [-1:60]; set yrange [-1:52]; set xtics 10; set ytics 10; set cbtics 0.0001; set format cb '%.4f'; \
set xlabel 'x'; set ylabel 'Y'; \
set origin 0.070,0.535; set size 0.405,0.385; \
set title '(a) corrected theory: gamma=0.01, epsilon=0.15, t=50'; \
plot theory_file using 1:2:(0.44):3 with circles lc palette fill solid 1.0 noborder notitle; \
set origin 0.555,0.535; set size 0.405,0.385; \
set title panel_b_title; \
plot panel_b_file using 1:2:(0.44):3 with circles lc palette fill solid 1.0 noborder notitle; \
unset colorbox; set size noratio; set origin 0.31,0.115; set size 0.38,0.31; \
set autoscale; set xrange [-1:60]; set yrange [0:0.00105]; set ytics 0.0001; set format y '%.4f'; \
set xlabel 'x'; set ylabel 'Probability'; set title 'off-center slice: gamma=0.010, epsilon=0.15, t=50'; \
set key bottom center horizontal font ',10' samplen 1.6; \
eval(cross_plot); \
unset multiplot"

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
