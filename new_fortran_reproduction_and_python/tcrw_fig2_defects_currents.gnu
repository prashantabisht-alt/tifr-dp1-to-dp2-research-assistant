#=====================================================================
# tcrw_fig2_defects_currents.gnu — Fig 2 rows (m), (n), (o):
#                                   J_tot, J_ω, J_Dr with plus-sign
#                                   defect, compared side-by-side
#                                   with clean box at ω = 0.
#
# 2 × 3 multiplot:
#                      J_tot              J_ω                J_Dr
#   clean box        weak CCW residual  strong CCW edge   strong CW tangential
#   defects box      circulation            CCW flow            JDr picks up
#                    around outer +         around outer +      NORMAL component
#                    around inner +         around inner +      at defect walls
#                    (hidden chirality      (chiral edge J      (this is the
#                     cancellation mostly    also wraps the      paper's "J_Dr
#                     survives + adapts      defect boundary)    perpendicular to
#                     to defect geometry)                        the edge" signal)
#
# Paper punchline — the ENTIRE reason for adding defects:
#   On a straight outer edge, reflection symmetry forbids J_Dr from
#   having a normal component, so it must run tangentially.  Interior
#   walls have a well-defined outward normal ⇒ J_Dr is no longer
#   symmetry-locked to the tangential direction, and acquires a real
#   radial (perpendicular-to-wall) component near the defect ring.
#
#   That normal J_Dr is the topological signature being advertised by
#   the paper — it appears whenever the geometry admits one, and is
#   invisible in a clean box.
#
# The bottom (defects) row has the defect cells overlaid as dark grey
# filled squares so you can see which arrows are next to a wall and
# which are in the bulk.
#
# Reads :
#   clean   : tcrw_fig2_Jtot_w0.0.txt,  tcrw_fig2_Jomega_w0.0.txt,  tcrw_fig2_JDr_w0.0.txt
#   defects : tcrw_fig2_Jtot_defects.txt, tcrw_fig2_Jomega_defects.txt, tcrw_fig2_JDr_defects.txt
#   layout  : tcrw_fig2_defects_layout.txt     (x, y of each defect cell)
#   each current file: x   y   Jx   Jy   |J|
#
# Output: tcrw_fig2_defects_currents.pdf + interactive qt
#
# Usage :
#   bash run.sh tcrw-plot-fig2-defects-currents       # PDF + qt
#   bash run.sh tcrw-plot-fig2-defects-currents-pdf   # PDF only (headless)
#   bash run.sh tcrw-plot-fig2-defects-currents-qt    # interactive qt only
#=====================================================================

if (!exists("mode")) mode = "both"

# ---- filenames ----
Jtot_c = 'tcrw_fig2_Jtot_w0.0.txt'
Jom_c  = 'tcrw_fig2_Jomega_w0.0.txt'
Jdr_c  = 'tcrw_fig2_JDr_w0.0.txt'

Jtot_d = 'tcrw_fig2_Jtot_defects.txt'
Jom_d  = 'tcrw_fig2_Jomega_defects.txt'
Jdr_d  = 'tcrw_fig2_JDr_defects.txt'

f_layout = 'tcrw_fig2_defects_layout.txt'

# ---- vector plot settings ----
arrow_len = 0.45

# Clean-row vector expression: only the vector field.
vec_clean(f) = sprintf("'%s' u ($3*$3+$4*$4 > 1e-30 ? $1 : NaN) : 2 : ($3/sqrt($3*$3+$4*$4)*%f) : ($4/sqrt($3*$3+$4*$4)*%f) : (sqrt($3*$3+$4*$4)) with vectors head filled size 0.12,20,60 lw 1.2 lc palette notitle", f, arrow_len, arrow_len)

# Defects-row vector expression: vector field + defect-cell overlay.
# The overlay uses `boxxyerror` (x:y:xdelta:ydelta) so each dark grey box
# covers EXACTLY one lattice cell (±0.5 around the cell centre), making the
# plus-sign geometry immediately readable on the vector-field panels.
vec_def(f) = sprintf("'%s' u ($3*$3+$4*$4 > 1e-30 ? $1 : NaN) : 2 : ($3/sqrt($3*$3+$4*$4)*%f) : ($4/sqrt($3*$3+$4*$4)*%f) : (sqrt($3*$3+$4*$4)) with vectors head filled size 0.12,20,60 lw 1.2 lc palette notitle, '%s' u 1:2:(0.5):(0.5) with boxxyerror fs solid 1.0 noborder fc rgb '#222222' notitle", f, arrow_len, arrow_len, f_layout)

# fill style for the boxxyerror overlay (gnuplot requires `set style fill`)
set style fill solid 1.0 noborder

# ---- palette ----
set palette defined (0 '#fff4d6', 0.3 '#fab66a', 0.7 '#c0381a', 1 '#3d0a00')
set logscale cb
set format cb "10^{%L}"

# ---- common per-panel geometry ----
set size ratio -1
set xrange [-0.5 : 9.5]
set yrange [-0.5 : 9.5]
set xtics 0, 2, 9
set ytics 0, 2, 9
set xlabel 'X'
set ylabel 'Y'
unset grid

#--------------------------------------------------------------------
# PDF  (headless first)
#--------------------------------------------------------------------
if (mode eq "pdf" || mode eq "both") {
    set terminal pdfcairo size 26cm,17cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig2_defects_currents.pdf'

    set multiplot layout 2,3 title \
        "TCRW Fig 2 currents: clean vs plus-sign defect   |   J = J_ω + J_{D_r}   |   ω = 0.0,  D_r = 10^{-3},  T = 10^{10}" \
        font 'Helvetica,11'

    # --- row 1: clean box ---
    set title "J_{tot}   (clean)"     font 'Helvetica,11'
    eval("plot " . vec_clean(Jtot_c))

    set title "J_ω     (clean)"       font 'Helvetica,11'
    eval("plot " . vec_clean(Jom_c))

    set title "J_{D_r}    (clean)"    font 'Helvetica,11'
    eval("plot " . vec_clean(Jdr_c))

    # --- row 2: defects box ---
    set title "J_{tot}   (+ defect)"  font 'Helvetica,11'
    eval("plot " . vec_def(Jtot_d))

    set title "J_ω     (+ defect)"    font 'Helvetica,11'
    eval("plot " . vec_def(Jom_d))

    set title "J_{D_r}    (+ defect)" font 'Helvetica,11'
    eval("plot " . vec_def(Jdr_d))

    unset multiplot
    unset output
    print "Wrote tcrw_fig2_defects_currents.pdf"
}

#--------------------------------------------------------------------
# INTERACTIVE qt
#--------------------------------------------------------------------
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 1400,900 enhanced font 'Helvetica,11'

    set multiplot layout 2,3 title \
        "TCRW Fig 2 currents: clean vs plus-sign defect   |   J = J_ω + J_{D_r}" \
        font 'Helvetica,11'

    set title "J_{tot}   (clean)"
    eval("plot " . vec_clean(Jtot_c))

    set title "J_ω     (clean)"
    eval("plot " . vec_clean(Jom_c))

    set title "J_{D_r}    (clean)"
    eval("plot " . vec_clean(Jdr_c))

    set title "J_{tot}   (+ defect)"
    eval("plot " . vec_def(Jtot_d))

    set title "J_ω     (+ defect)"
    eval("plot " . vec_def(Jom_d))

    set title "J_{D_r}    (+ defect)"
    eval("plot " . vec_def(Jdr_d))

    unset multiplot
    pause mouse close
}
