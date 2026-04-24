#=====================================================================
# tcrw_fig2_currents.gnu — Fig 2 panels (c-e), (h-j) + symmetry-check
#                          row ω=1.0:   J_tot, J_ω, J_Dr vector fields
#
# 3 × 3 multiplot:
#                 J_tot               J_ω                  J_Dr
#   ω = 0.0     weak CCW residual    strong CCW edge     strong CW tangential
#                (= J_ω + J_Dr ≪ each) (chiral edge J)     (anti-chiral post-noise J)
#   ω = 0.5     ~zero                ~zero               ~zero
#                (symmetry kills any net circulation)
#   ω = 1.0     weak CW residual     strong CW edge      strong CCW tangential
#                (mirror of ω=0.0 row — sanity check)
#
# Paper punchline:
#   * J_tot = J_ω + J_Dr   (decomposition by prev-step type).
#   * For ω = 0, J_ω circulates CCW; J_Dr runs tangentially CW with
#     essentially the same magnitude (|radial| ~ 3e-8 vs |tangential| ~ 3e-5).
#     The two cancel to leave J_tot one decade smaller — "hidden chirality".
#   * For ω = 1, everything flips: J_ω CW, J_Dr CCW, J_tot a weak CW residual.
#     Pixel-by-pixel this row should be the vector negative of the ω=0 row.
#   * NOTE: the paper describes J_Dr as "perpendicular to the edge" for the
#     DEFECTIVE box (separate driver), where interior walls give rebound a
#     normal direction.  In a CLEAN box the only symmetry-compatible J_Dr
#     on a straight edge is tangential.
#   * For ω = 0.5, both J_ω and J_Dr are small and symmetric — no net current.
#
# Reads :
#   tcrw_fig2_Jtot_w0.0.txt    tcrw_fig2_Jomega_w0.0.txt    tcrw_fig2_JDr_w0.0.txt
#   tcrw_fig2_Jtot_w0.5.txt    tcrw_fig2_Jomega_w0.5.txt    tcrw_fig2_JDr_w0.5.txt
#   tcrw_fig2_Jtot_w1.0.txt    tcrw_fig2_Jomega_w1.0.txt    tcrw_fig2_JDr_w1.0.txt
#   each with columns:   x   y   Jx   Jy   |J|
#   (allowed inner sites only, 64 rows per file for L = 10)
#
# Output: tcrw_fig2_currents.pdf + interactive qt
#
# Usage:
#   bash run.sh tcrw-plot-fig2-currents         # PDF then interactive qt
#   bash run.sh tcrw-plot-fig2-currents-pdf     # PDF only (headless)
#   bash run.sh tcrw-plot-fig2-currents-qt      # interactive qt only
#=====================================================================

if (!exists("mode")) mode = "both"

# ---- filenames (3 ω × 3 components) ----
Jtot0 = 'tcrw_fig2_Jtot_w0.0.txt'
Jom0  = 'tcrw_fig2_Jomega_w0.0.txt'
Jdr0  = 'tcrw_fig2_JDr_w0.0.txt'
Jtot5 = 'tcrw_fig2_Jtot_w0.5.txt'
Jom5  = 'tcrw_fig2_Jomega_w0.5.txt'
Jdr5  = 'tcrw_fig2_JDr_w0.5.txt'
Jtot1 = 'tcrw_fig2_Jtot_w1.0.txt'
Jom1  = 'tcrw_fig2_Jomega_w1.0.txt'
Jdr1  = 'tcrw_fig2_JDr_w1.0.txt'

# ---- vector plot settings shared by all nine panels ----
arrow_len = 0.45                        # lattice units; keeps arrows < half a cell

# One-line sprintf (no backslash-inside-string — avoids parser version quirks).
vec_expr(f) = sprintf("'%s' u ($3*$3+$4*$4 > 1e-30 ? $1 : NaN) : 2 : ($3/sqrt($3*$3+$4*$4)*%f) : ($4/sqrt($3*$3+$4*$4)*%f) : (sqrt($3*$3+$4*$4)) with vectors head filled size 0.12,20,60 lw 1.2 lc palette notitle", f, arrow_len, arrow_len)

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
    set terminal pdfcairo size 26cm,24cm enhanced font 'Helvetica,10'
    set output 'tcrw_fig2_currents.pdf'

    set multiplot layout 3,3 title \
        "TCRW Fig 2 currents   |   J = J_ω + J_{D_r}   |   L = 10,  D_r = 10^{-3},  T = 10^{10}" \
        font 'Helvetica,11'

    # --- row 1: ω = 0.0 ---
    set title "J_{tot}   (ω = 0.0)"     font 'Helvetica,11'
    eval("plot " . vec_expr(Jtot0))

    set title "J_ω     (ω = 0.0)"       font 'Helvetica,11'
    eval("plot " . vec_expr(Jom0))

    set title "J_{D_r}    (ω = 0.0)"    font 'Helvetica,11'
    eval("plot " . vec_expr(Jdr0))

    # --- row 2: ω = 0.5 ---
    set title "J_{tot}   (ω = 0.5)"     font 'Helvetica,11'
    eval("plot " . vec_expr(Jtot5))

    set title "J_ω     (ω = 0.5)"       font 'Helvetica,11'
    eval("plot " . vec_expr(Jom5))

    set title "J_{D_r}    (ω = 0.5)"    font 'Helvetica,11'
    eval("plot " . vec_expr(Jdr5))

    # --- row 3: ω = 1.0 ---
    set title "J_{tot}   (ω = 1.0)"     font 'Helvetica,11'
    eval("plot " . vec_expr(Jtot1))

    set title "J_ω     (ω = 1.0)"       font 'Helvetica,11'
    eval("plot " . vec_expr(Jom1))

    set title "J_{D_r}    (ω = 1.0)"    font 'Helvetica,11'
    eval("plot " . vec_expr(Jdr1))

    unset multiplot
    unset output
    print "Wrote tcrw_fig2_currents.pdf"
}

#--------------------------------------------------------------------
# INTERACTIVE qt
#--------------------------------------------------------------------
if (mode eq "qt" || mode eq "both") {
    set terminal qt size 1400,1300 enhanced font 'Helvetica,11'

    set multiplot layout 3,3 title \
        "TCRW Fig 2 currents   |   J = J_ω + J_{D_r}" \
        font 'Helvetica,11'

    set title "J_{tot}   (ω = 0.0)"
    eval("plot " . vec_expr(Jtot0))

    set title "J_ω     (ω = 0.0)"
    eval("plot " . vec_expr(Jom0))

    set title "J_{D_r}    (ω = 0.0)"
    eval("plot " . vec_expr(Jdr0))

    set title "J_{tot}   (ω = 0.5)"
    eval("plot " . vec_expr(Jtot5))

    set title "J_ω     (ω = 0.5)"
    eval("plot " . vec_expr(Jom5))

    set title "J_{D_r}    (ω = 0.5)"
    eval("plot " . vec_expr(Jdr5))

    set title "J_{tot}   (ω = 1.0)"
    eval("plot " . vec_expr(Jtot1))

    set title "J_ω     (ω = 1.0)"
    eval("plot " . vec_expr(Jom1))

    set title "J_{D_r}    (ω = 1.0)"
    eval("plot " . vec_expr(Jdr1))

    unset multiplot
    pause mouse close
}
