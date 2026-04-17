set terminal qt

# ============================================================
# Plot 1: MSD vs t
# ============================================================
set xlabel "t (steps)"
set ylabel "MSD(t)"
set title "TCRW: Mean Square Displacement (PBC)"
set key top left
set logscale xy
set grid

plot "tcrw_pbc_achiral.dat" u 1:2 every ::1 w l lw 2 lc rgb "#2166ac" \
     title "{/Symbol w}=0.5, D_r=0.001 (achiral)", \
     "tcrw_pbc_chiral.dat" u 1:2 every ::1 w l lw 2 lc rgb "#b2182b" \
     title "{/Symbol w}=0.9, D_r=0.1 (chiral)", \
     x lw 1 dt 2 lc rgb "gray50" title "~ t (diffusive)"

pause -1 "Press ENTER for next plot..."
unset logscale

# ============================================================
# Plot 2: D(omega)
# ============================================================
set xlabel "{/Symbol w}"
set ylabel "D_{eff}"
set title "TCRW: Effective Diffusion Coefficient vs {/Symbol w}  (D_r=0.1)"
set key top right
set grid
set xrange [0:1]

plot "tcrw_D_vs_omega.dat" u 1:2 w lp lw 2 pt 7 ps 1.2 lc rgb "#d6604d" \
     title "D_{eff}({/Symbol w})"

pause -1 "Press ENTER for next plot..."
unset xrange

# ============================================================
# Plot 3: Single trajectory
# ============================================================
set xlabel "x"
set ylabel "y"
set title "TCRW: Single Trajectory ({/Symbol w}=0.9, D_r=0.05)"
set key off
set grid
set size ratio -1

plot "tcrw_trajectory.dat" u 2:3 w l lw 0.5 lc rgb "#4393c3"

pause -1 "Press ENTER for next plot..."
unset size

# ============================================================
# Plot 4: OBC density heatmap
# ============================================================
set xlabel "x"
set ylabel "y"
set title "TCRW OBC: Steady-State Density P(x,y)  ({/Symbol w}=0.9, D_r=0.1)"
set key off
set size ratio 1
set view map
set pm3d interpolate 2,2
set palette defined (0 "white", 0.5 "#fee090", 1 "#d73027")
set cblabel "P(x,y)"

splot "tcrw_obc_visits.dat" u 1:2:3 w pm3d

pause -1 "Press ENTER for next plot..."
unset view
unset pm3d

# ============================================================
# Plot 5: OBC current field
# ============================================================
set xlabel "x"
set ylabel "y"
set title "TCRW OBC: Current Field J(x,y)  ({/Symbol w}=0.9, D_r=0.1)"
set key off
set size ratio 1
set xrange [-0.5:19.5]
set yrange [-0.5:19.5]

scale = 5000.0
plot "tcrw_obc_currents.dat" u 1:2:(scale*$3):(scale*$4) w vectors \
     head size 0.3,20 filled lw 1.5 lc rgb "#1a9850"

pause -1 "Press ENTER to exit..."
