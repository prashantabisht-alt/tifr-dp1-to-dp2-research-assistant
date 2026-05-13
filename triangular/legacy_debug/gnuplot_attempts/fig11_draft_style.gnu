# Draft-style Fig. 11 plot.
#
# Run from the triangular folder:
#   python3 export_fig11_draft_style_data.py
#   gnuplot fig11_draft_style.gnu

if (!exists("mode")) mode = "both"

pdf_out = "fig11_corrected_draft_style.pdf"
png_out = "fig11_corrected_draft_style.png"

data_dir = "outputs"
theory_file = data_dir . "/fig11_draft_theory_points.txt"
kmc_file = data_dir . "/fig11_draft_kmc_points.txt"
cross_file = data_dir . "/fig11_draft_cross_section.txt"

plot_body = "set palette defined (0 'black', 0.18 'purple', 0.45 'red', 0.75 'orange', 1 'yellow'); set cbrange [0:0.0018]; set multiplot layout 2,2 title 'Triangular active walker: corrected Fig. 11, draft-style coordinates'; set size ratio -1; set xrange [0:60]; set yrange [0:30]; set xlabel 'x'; set ylabel 'y'; set colorbox; unset key; set title '(a) Theory, corrected c_3'; plot theory_file using 1:2:3 with points pointtype 7 pointsize 0.42 palette notitle; set title '(b) Kinetic Monte Carlo'; plot kmc_file using 1:2:3 with points pointtype 7 pointsize 0.42 palette notitle; unset colorbox; set size noratio; set xrange [0:60]; set yrange [0:*]; set xlabel 'x'; set ylabel 'Probability'; set grid; set key top right; set title '(c) Cross-section through y=L/2'; plot cross_file using 1:2 with points pointtype 7 pointsize 0.7 linecolor rgb 'black' title 'KMC', cross_file using 1:3 with lines linewidth 2 linecolor rgb 'red' title 'buggy theory', cross_file using 1:4 with lines linewidth 2 linecolor rgb 'blue' title 'corrected theory'; unset grid; unset key; unset border; unset xtics; unset ytics; unset xlabel; unset ylabel; set xrange [0:1]; set yrange [0:1]; set title ''; plot 0.5 with lines linecolor rgb 'white' notitle; unset multiplot"

if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo enhanced color font "Helvetica,10" size 11in,8.5in
    set output pdf_out
    @plot_body
    set output
}

if (mode eq "png" || mode eq "both") {
    set terminal pngcairo enhanced color font "Helvetica,10" size 1600,1200
    set output png_out
    set border
    set xtics
    set ytics
    @plot_body
    set output
}

if (mode eq "qt") {
    set terminal qt enhanced font "Helvetica,10" size 1200,900
    set border
    set xtics
    set ytics
    @plot_body
    pause mouse close
}
