# Run from the triangular folder:
#   python3 export_fig11_gnuplot_data.py
#   gnuplot fig11_corrected.gnu
#
# Interactive only:
#   gnuplot -persist -e "mode='qt'" fig11_corrected.gnu
#
# PDF only:
#   gnuplot -e "mode='pdf'" fig11_corrected.gnu
#
# PNG only:
#   gnuplot -e "mode='png'" fig11_corrected.gnu

if (!exists("mode")) mode = "both"

pdf_out = "fig11_corrected_gnuplot.pdf"
png_out = "fig11_corrected_gnuplot.png"

data_dir = "outputs"
theory_file = data_dir . "/fig11_theory_corrected_xyz.txt"
kmc_file = data_dir . "/fig11_kmc_xyz.txt"
cross_file = data_dir . "/fig11_cross_section.txt"

load_settings = "set palette defined (0 'black', 0.25 'red', 0.55 'orange', 0.8 'yellow', 1 'white'); set cbrange [0:*]"

if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo enhanced color font "Helvetica,10" size 11in,8.5in
    set output pdf_out

    @load_settings
    set autoscale
    set multiplot layout 2,2 title "Triangular active walker: corrected theory vs KMC"

    set view map
    set size ratio -1
    set pm3d interpolate 1,1 corners2color c1
    unset key
    set colorbox
    set title "(a) Corrected exact theory"
    set xlabel "x = 2 n_1 + n_2"
    set ylabel "y = sqrt(3) n_2"
    splot theory_file using 1:2:3 with pm3d notitle

    set title "(b) Kinetic Monte Carlo"
    set xlabel "x = 2 n_1 + n_2"
    set ylabel "y = sqrt(3) n_2"
    splot kmc_file using 1:2:3 with pm3d notitle

    unset pm3d
    unset colorbox
    set size noratio
    set title "(c) Cross-section through n_2 = L/2"
    set xlabel "n_1"
    set ylabel "P(n_1, n_2=L/2, t)"
    set grid
    set key top right
    plot cross_file using 1:2 with points pointtype 7 pointsize 0.8 linecolor rgb "black" title "KMC", \
         cross_file using 1:3 with lines linewidth 2 linecolor rgb "red" title "buggy theory", \
         cross_file using 1:4 with lines linewidth 2 linecolor rgb "blue" title "corrected theory"

    unset grid
    unset key
    unset border
    unset xtics
    unset ytics
    unset xlabel
    unset ylabel
    set xrange [0:1]
    set yrange [0:1]
    set title ""
    plot 0.5 with lines linecolor rgb "white" notitle

    unset multiplot
    set output
}

if (mode eq "png" || mode eq "both") {
    set terminal pngcairo enhanced color font "Helvetica,10" size 1600,1200
    set output png_out

    @load_settings
    set autoscale
    set border
    set xtics
    set ytics
    set multiplot layout 2,2 title "Triangular active walker: corrected theory vs KMC"

    set view map
    set size ratio -1
    set pm3d interpolate 1,1 corners2color c1
    unset key
    set colorbox
    set title "(a) Corrected exact theory"
    set xlabel "x = 2 n_1 + n_2"
    set ylabel "y = sqrt(3) n_2"
    splot theory_file using 1:2:3 with pm3d notitle

    set title "(b) Kinetic Monte Carlo"
    set xlabel "x = 2 n_1 + n_2"
    set ylabel "y = sqrt(3) n_2"
    splot kmc_file using 1:2:3 with pm3d notitle

    unset pm3d
    unset colorbox
    set size noratio
    set title "(c) Cross-section through n_2 = L/2"
    set xlabel "n_1"
    set ylabel "P(n_1, n_2=L/2, t)"
    set grid
    set key top right
    plot cross_file using 1:2 with points pointtype 7 pointsize 0.8 linecolor rgb "black" title "KMC", \
         cross_file using 1:3 with lines linewidth 2 linecolor rgb "red" title "buggy theory", \
         cross_file using 1:4 with lines linewidth 2 linecolor rgb "blue" title "corrected theory"

    unset grid
    unset key
    unset border
    unset xtics
    unset ytics
    unset xlabel
    unset ylabel
    set xrange [0:1]
    set yrange [0:1]
    set title ""
    plot 0.5 with lines linecolor rgb "white" notitle

    unset multiplot
    set output
}

if (mode eq "qt") {
    set terminal qt enhanced font "Helvetica,10" size 1200,900

    @load_settings
    set autoscale
    set border
    set xtics
    set ytics
    set multiplot layout 2,2 title "Triangular active walker: corrected theory vs KMC"

    set view map
    set size ratio -1
    set pm3d interpolate 1,1 corners2color c1
    unset key
    set colorbox
    set title "(a) Corrected exact theory"
    set xlabel "x = 2 n_1 + n_2"
    set ylabel "y = sqrt(3) n_2"
    splot theory_file using 1:2:3 with pm3d notitle

    set title "(b) Kinetic Monte Carlo"
    set xlabel "x = 2 n_1 + n_2"
    set ylabel "y = sqrt(3) n_2"
    splot kmc_file using 1:2:3 with pm3d notitle

    unset pm3d
    unset colorbox
    set size noratio
    set title "(c) Cross-section through n_2 = L/2"
    set xlabel "n_1"
    set ylabel "P(n_1, n_2=L/2, t)"
    set grid
    set key top right
    plot cross_file using 1:2 with points pointtype 7 pointsize 0.8 linecolor rgb "black" title "KMC", \
         cross_file using 1:3 with lines linewidth 2 linecolor rgb "red" title "buggy theory", \
         cross_file using 1:4 with lines linewidth 2 linecolor rgb "blue" title "corrected theory"

    unset grid
    unset key
    unset border
    unset xtics
    unset ytics
    unset xlabel
    unset ylabel
    set xrange [0:1]
    set yrange [0:1]
    set title ""
    plot 0.5 with lines linecolor rgb "white" notitle

    unset multiplot
    pause mouse close
}
