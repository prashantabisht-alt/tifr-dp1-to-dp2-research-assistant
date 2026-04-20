#!/usr/bin/env bash
#=====================================================================
# run.sh — build + run + plot driver for the TCRW Fortran reproduction.
#
# Style follows the user's own run.sh conventions:
#   - strict bash (`set -euo pipefail`)
#   - .build/ holds object files / executables
#   - named targets, per figure AND per subpanel (for fine-grained review)
#
# Usage:
#   bash run.sh <target>
#
# Available targets (we'll add rows here figure-by-figure):
#   tcrw-check              probe gfortran + gnuplot + LAPACK link
#   tcrw-fig1b              compile + run tcrw_fig1b.f90
#   tcrw-plot-fig1b-w0      plot ONLY the ω=0 trajectory (single panel preview)
#   tcrw-plot-fig1b         plot Fig 1(b) (qt + PDF, 2×2 paper layout)
#   tcrw-plot-fig1b-pdf     plot Fig 1(b) (PDF only; headless)
#   tcrw-plot-fig1b-png     plot Fig 1(b) (PNG only; easy visual diff)
#   tcrw-fig1b-all          tcrw-fig1b → tcrw-plot-fig1b
#   tcrw-fig1c              compile + run tcrw_fig1c.f90  (MSD ensemble, ~10 min)
#   tcrw-plot-fig1c         plot Fig 1(c) MSD (qt + PDF)
#   tcrw-plot-fig1c-pdf     plot Fig 1(c) MSD (PDF only; headless)
#   tcrw-plot-fig1c-qt      plot Fig 1(c) MSD (qt only; fastest feedback)
#   tcrw-fig1c-all          tcrw-fig1c → tcrw-plot-fig1c
#   tcrw-fig1d              compile + run tcrw_fig1d.f90  (D(ω) sweep, ~7–15 min)
#   tcrw-plot-fig1d         plot Fig 1(d) D(ω) (PDF then interactive qt)
#   tcrw-plot-fig1d-pdf     plot Fig 1(d) D(ω) (PDF only; headless)
#   tcrw-plot-fig1d-qt      plot Fig 1(d) D(ω) (interactive qt only)
#   tcrw-fig1d-all          tcrw-fig1d → tcrw-plot-fig1d
#   tcrw-fig2-clean             compile + run tcrw_fig2_clean.f90 (2 ω × 10^9 steps, ~1–2 min)
#   tcrw-plot-fig2-occ          plot Fig 2 P(X,Y) heatmap pair   (PDF + qt)
#   tcrw-plot-fig2-occ-pdf      plot Fig 2 P(X,Y)                (PDF only)
#   tcrw-plot-fig2-occ-qt       plot Fig 2 P(X,Y)                (qt only)
#   tcrw-plot-fig2-traj         plot Fig 2 trajectory pair       (PDF + qt)
#   tcrw-plot-fig2-traj-pdf     plot Fig 2 trajectory pair       (PDF only)
#   tcrw-plot-fig2-traj-qt      plot Fig 2 trajectory pair       (qt only)
#   tcrw-plot-fig2-currents     plot Fig 2 2×3 vector fields     (PDF + qt)
#   tcrw-plot-fig2-currents-pdf plot Fig 2 2×3 vector fields     (PDF only)
#   tcrw-plot-fig2-currents-qt  plot Fig 2 2×3 vector fields     (qt only)
#   tcrw-fig2-clean-all         tcrw-fig2-clean → all three Fig 2 plots
#   tcrw-fig2-defects                    compile + run tcrw_fig2_defects.f90 (1 ω × 10^10 steps, ~5 min)
#   tcrw-plot-fig2-defects-occ           plot Fig 2 defects P(X,Y) pair   (PDF + qt)
#   tcrw-plot-fig2-defects-occ-pdf       plot Fig 2 defects P(X,Y)        (PDF only)
#   tcrw-plot-fig2-defects-occ-qt        plot Fig 2 defects P(X,Y)        (qt only)
#   tcrw-plot-fig2-defects-traj          plot Fig 2 defects trajectory    (PDF + qt)
#   tcrw-plot-fig2-defects-traj-pdf      plot Fig 2 defects trajectory    (PDF only)
#   tcrw-plot-fig2-defects-traj-qt       plot Fig 2 defects trajectory    (qt only)
#   tcrw-plot-fig2-defects-currents      plot Fig 2 defects 2×3 currents  (PDF + qt)
#   tcrw-plot-fig2-defects-currents-pdf  plot Fig 2 defects 2×3 currents  (PDF only)
#   tcrw-plot-fig2-defects-currents-qt   plot Fig 2 defects 2×3 currents  (qt only)
#   tcrw-fig2-defects-all                tcrw-fig2-defects → all three defects plots
#   tcrw-fig3a              compile + run tcrw_fig3a.f90  (P_edge/P_bulk vs D_r, ~3–5 min)
#   tcrw-plot-fig3a         plot Fig 3(a) (PDF then interactive qt)
#   tcrw-plot-fig3a-pdf     plot Fig 3(a) (PDF only; headless)
#   tcrw-plot-fig3a-qt      plot Fig 3(a) (interactive qt only)
#   tcrw-fig3a-all          tcrw-fig3a → tcrw-plot-fig3a
#   tcrw-fig3b              compile + run tcrw_fig3b.f90  (|J_Dr|/|J_ω| wall vs D_r, ~25 min)
#   tcrw-plot-fig3b         plot Fig 3(b) (PDF then interactive qt)
#   tcrw-plot-fig3b-pdf     plot Fig 3(b) (PDF only; headless)
#   tcrw-plot-fig3b-qt      plot Fig 3(b) (interactive qt only)
#   tcrw-fig3b-all          tcrw-fig3b → tcrw-plot-fig3b
#   tcrw-fig3cde            compile + run tcrw_fig3cde.f90 (per-site left-wall currents, ~15–20 min)
#   tcrw-plot-fig3c         plot Fig 3(c) quiver J_Dr   (PDF + qt)
#   tcrw-plot-fig3c-pdf     plot Fig 3(c) quiver J_Dr   (PDF only)
#   tcrw-plot-fig3c-qt      plot Fig 3(c) quiver J_Dr   (qt only)
#   tcrw-plot-fig3d         plot Fig 3(d) quiver J_ω    (PDF + qt)
#   tcrw-plot-fig3d-pdf     plot Fig 3(d) quiver J_ω    (PDF only)
#   tcrw-plot-fig3d-qt      plot Fig 3(d) quiver J_ω    (qt only)
#   tcrw-plot-fig3e         plot Fig 3(e) θ_JDr vs D_r  (PDF + qt)
#   tcrw-plot-fig3e-pdf     plot Fig 3(e)                (PDF only)
#   tcrw-plot-fig3e-qt      plot Fig 3(e)                (qt only)
#   tcrw-fig3cde-all        tcrw-fig3cde → all three (c/d/e) plots
#   tcrw-fig3f              compile + run tcrw_fig3f.f90  (P_edge/P_bulk vs ω,  ~5–10 min)
#   tcrw-plot-fig3f         plot Fig 3(f) (PDF then interactive qt)
#   tcrw-plot-fig3f-pdf     plot Fig 3(f) (PDF only; headless)
#   tcrw-plot-fig3f-qt      plot Fig 3(f) (interactive qt only)
#   tcrw-fig3f-all          tcrw-fig3f → tcrw-plot-fig3f
#   tcrw-fig3g              compile + run tcrw_fig3g.f90  (|J_Dr|/|J_ω| wall vs ω, ~10–12 min)
#   tcrw-plot-fig3g         plot Fig 3(g) (PDF then interactive qt)
#   tcrw-plot-fig3g-pdf     plot Fig 3(g) (PDF only; headless)
#   tcrw-plot-fig3g-qt      plot Fig 3(g) (interactive qt only)
#   tcrw-fig3g-all          tcrw-fig3g → tcrw-plot-fig3g
#   tcrw-clean              remove .build/ and generated .txt / .pdf
#
# (Fig 3, ... will get more target rows as we add panels.)
#
# Requires:  gfortran (Homebrew gcc on macOS), gnuplot.
#            Apple Accelerate is auto-detected here for later figures
#            that need LAPACK (diagonalization of the Bloch matrix).
#=====================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD="$ROOT/.build"
mkdir -p "$BUILD"

# ------------ link flag for LAPACK (used later; harmless now) ---------
case "$(uname)" in
  Darwin) LAPACK_LIBS="-framework Accelerate" ;;
  Linux)  LAPACK_LIBS="-llapack -lblas"       ;;
  *)      LAPACK_LIBS=""                      ;;
esac

FFLAGS_BASE="-O2 -fno-range-check -ffree-line-length-none"

# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------
need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "✗ required command '$1' not found in PATH" >&2
    exit 1
  }
}

# compile_and_run <tag> <src.f90> [extra link flags...]
compile_and_run() {
  local tag="$1"; shift
  local src="$1"; shift
  local exe="$BUILD/$tag"

  need_cmd gfortran
  echo "==> gfortran $FFLAGS_BASE $src $* -o $exe"
  ( cd "$ROOT" && gfortran $FFLAGS_BASE "$src" "$@" -o "$exe" )

  # macOS may kill unsigned binaries (Killed: 9); codesign if available
  if [[ "$(uname)" == "Darwin" ]] && command -v codesign >/dev/null 2>&1; then
    codesign --sign - "$exe" 2>/dev/null || true
  fi

  echo "==> running $tag"
  ( cd "$ROOT" && "$exe" )
  echo "==> $tag done"
}

# plot_gnu <gnufile> [mode]     mode = qt | pdf | both   (default: both)
plot_gnu() {
  local gnu="$1"
  local mode="${2:-both}"
  need_cmd gnuplot
  echo "==> gnuplot (mode=$mode) $gnu"
  ( cd "$ROOT" && gnuplot -persist -e "mode='$mode'" "$gnu" )
}

# ---------------------------------------------------------------------
# targets
# ---------------------------------------------------------------------

tcrw_check() {
  echo "--- environment check ---"
  need_cmd gfortran
  gfortran --version | head -1
  need_cmd gnuplot
  gnuplot --version
  echo
  echo "Planned LAPACK link flag: $LAPACK_LIBS"
  echo "(run ../lapack_check.sh at repo root to verify end-to-end.)"
}

# ---- Fig 1(b): sample trajectories ----
tcrw_fig1b()          { compile_and_run tcrw_fig1b tcrw_fig1b.f90; }
tcrw_plot_fig1b_w0()  { plot_gnu tcrw_fig1b_w0.gnu both; }   # single ω=0 preview
tcrw_plot_fig1b()     { plot_gnu tcrw_fig1b.gnu both; }
tcrw_plot_fig1b_pdf() { plot_gnu tcrw_fig1b.gnu pdf;  }
tcrw_plot_fig1b_png() { plot_gnu tcrw_fig1b.gnu png;  }
tcrw_fig1b_all()      { tcrw_fig1b; tcrw_plot_fig1b; }

# ---- Fig 1(c): MSD ensemble ----
tcrw_fig1c()          { compile_and_run tcrw_fig1c tcrw_fig1c.f90; }
tcrw_plot_fig1c()     { plot_gnu tcrw_fig1c.gnu both; }
tcrw_plot_fig1c_pdf() { plot_gnu tcrw_fig1c.gnu pdf;  }
tcrw_plot_fig1c_qt()  { plot_gnu tcrw_fig1c.gnu qt;   }
tcrw_fig1c_all()      { tcrw_fig1c; tcrw_plot_fig1c; }

# ---- Fig 1(d): D(ω) sweep ----
tcrw_fig1d()          { compile_and_run tcrw_fig1d tcrw_fig1d.f90; }
tcrw_plot_fig1d()     { plot_gnu tcrw_fig1d.gnu both; }
tcrw_plot_fig1d_pdf() { plot_gnu tcrw_fig1d.gnu pdf;  }
tcrw_plot_fig1d_qt()  { plot_gnu tcrw_fig1d.gnu qt;   }
tcrw_fig1d_all()      { tcrw_fig1d; tcrw_plot_fig1d; }

# ---- Fig 2 (clean OBC box): 2 ω × 10^9 steps ----
tcrw_fig2_clean()               { compile_and_run tcrw_fig2_clean tcrw_fig2_clean.f90; }
tcrw_plot_fig2_occ()            { plot_gnu tcrw_fig2_occ.gnu      both; }
tcrw_plot_fig2_occ_pdf()        { plot_gnu tcrw_fig2_occ.gnu      pdf;  }
tcrw_plot_fig2_occ_qt()         { plot_gnu tcrw_fig2_occ.gnu      qt;   }
tcrw_plot_fig2_traj()           { plot_gnu tcrw_fig2_traj.gnu     both; }
tcrw_plot_fig2_traj_pdf()       { plot_gnu tcrw_fig2_traj.gnu     pdf;  }
tcrw_plot_fig2_traj_qt()        { plot_gnu tcrw_fig2_traj.gnu     qt;   }
tcrw_plot_fig2_currents()       { plot_gnu tcrw_fig2_currents.gnu both; }
tcrw_plot_fig2_currents_pdf()   { plot_gnu tcrw_fig2_currents.gnu pdf;  }
tcrw_plot_fig2_currents_qt()    { plot_gnu tcrw_fig2_currents.gnu qt;   }
tcrw_fig2_clean_all() {
  tcrw_fig2_clean
  tcrw_plot_fig2_occ_pdf
  tcrw_plot_fig2_traj_pdf
  tcrw_plot_fig2_currents_pdf
}

# ---- Fig 2 (plus-sign defect): 1 ω × 10^10 steps ----
tcrw_fig2_defects()                      { compile_and_run tcrw_fig2_defects tcrw_fig2_defects.f90; }
tcrw_plot_fig2_defects_occ()             { plot_gnu tcrw_fig2_defects_occ.gnu      both; }
tcrw_plot_fig2_defects_occ_pdf()         { plot_gnu tcrw_fig2_defects_occ.gnu      pdf;  }
tcrw_plot_fig2_defects_occ_qt()          { plot_gnu tcrw_fig2_defects_occ.gnu      qt;   }
tcrw_plot_fig2_defects_traj()            { plot_gnu tcrw_fig2_defects_traj.gnu     both; }
tcrw_plot_fig2_defects_traj_pdf()        { plot_gnu tcrw_fig2_defects_traj.gnu     pdf;  }
tcrw_plot_fig2_defects_traj_qt()         { plot_gnu tcrw_fig2_defects_traj.gnu     qt;   }
tcrw_plot_fig2_defects_currents()        { plot_gnu tcrw_fig2_defects_currents.gnu both; }
tcrw_plot_fig2_defects_currents_pdf()    { plot_gnu tcrw_fig2_defects_currents.gnu pdf;  }
tcrw_plot_fig2_defects_currents_qt()     { plot_gnu tcrw_fig2_defects_currents.gnu qt;   }
tcrw_fig2_defects_all() {
  tcrw_fig2_defects
  tcrw_plot_fig2_defects_occ_pdf
  tcrw_plot_fig2_defects_traj_pdf
  tcrw_plot_fig2_defects_currents_pdf
}

# ---- Fig 3(a): P_edge/P_bulk vs D_r (ω = 1, L ∈ {4,9,19,49}) ----
tcrw_fig3a()          { compile_and_run tcrw_fig3a tcrw_fig3a.f90; }
tcrw_plot_fig3a()     { plot_gnu tcrw_fig3a.gnu both; }
tcrw_plot_fig3a_pdf() { plot_gnu tcrw_fig3a.gnu pdf;  }
tcrw_plot_fig3a_qt()  { plot_gnu tcrw_fig3a.gnu qt;   }
tcrw_fig3a_all()      { tcrw_fig3a; tcrw_plot_fig3a; }

# ---- Fig 3(b): |J_Dr|_wall / |J_ω|_wall vs D_r  (ω = 1, L ∈ {4,9,19,49}) ----
tcrw_fig3b()          { compile_and_run tcrw_fig3b tcrw_fig3b.f90; }
tcrw_plot_fig3b()     { plot_gnu tcrw_fig3b.gnu both; }
tcrw_plot_fig3b_pdf() { plot_gnu tcrw_fig3b.gnu pdf;  }
tcrw_plot_fig3b_qt()  { plot_gnu tcrw_fig3b.gnu qt;   }
tcrw_fig3b_all()      { tcrw_fig3b; tcrw_plot_fig3b; }

# ---- Fig 3(c), 3(d), 3(e):  per-site left-wall currents (L = 10, ω = 1) ----
tcrw_fig3cde()        { compile_and_run tcrw_fig3cde tcrw_fig3cde.f90; }
tcrw_plot_fig3c()     { plot_gnu tcrw_fig3c.gnu both; }
tcrw_plot_fig3c_pdf() { plot_gnu tcrw_fig3c.gnu pdf;  }
tcrw_plot_fig3c_qt()  { plot_gnu tcrw_fig3c.gnu qt;   }
tcrw_plot_fig3d()     { plot_gnu tcrw_fig3d.gnu both; }
tcrw_plot_fig3d_pdf() { plot_gnu tcrw_fig3d.gnu pdf;  }
tcrw_plot_fig3d_qt()  { plot_gnu tcrw_fig3d.gnu qt;   }
tcrw_plot_fig3e()     { plot_gnu tcrw_fig3e.gnu both; }
tcrw_plot_fig3e_pdf() { plot_gnu tcrw_fig3e.gnu pdf;  }
tcrw_plot_fig3e_qt()  { plot_gnu tcrw_fig3e.gnu qt;   }
tcrw_fig3cde_all() {
  tcrw_fig3cde
  tcrw_plot_fig3c_pdf
  tcrw_plot_fig3d_pdf
  tcrw_plot_fig3e_pdf
}

# ---- Fig 3(f): P_edge/P_bulk vs ω  (D_r = 10^-3, L ∈ {10,19,49}) ----
tcrw_fig3f()          { compile_and_run tcrw_fig3f tcrw_fig3f.f90; }
tcrw_plot_fig3f()     { plot_gnu tcrw_fig3f.gnu both; }
tcrw_plot_fig3f_pdf() { plot_gnu tcrw_fig3f.gnu pdf;  }
tcrw_plot_fig3f_qt()  { plot_gnu tcrw_fig3f.gnu qt;   }
tcrw_fig3f_all()      { tcrw_fig3f; tcrw_plot_fig3f; }

# ---- Fig 3(g): |J_Dr|_wall / |J_ω|_wall vs ω  (D_r = 10^-3, L ∈ {10,19,49}) ----
tcrw_fig3g()          { compile_and_run tcrw_fig3g tcrw_fig3g.f90; }
tcrw_plot_fig3g()     { plot_gnu tcrw_fig3g.gnu both; }
tcrw_plot_fig3g_pdf() { plot_gnu tcrw_fig3g.gnu pdf;  }
tcrw_plot_fig3g_qt()  { plot_gnu tcrw_fig3g.gnu qt;   }
tcrw_fig3g_all()      { tcrw_fig3g; tcrw_plot_fig3g; }

tcrw_clean() {
  echo "==> cleaning .build and generated outputs"
  rm -rf "$BUILD"
  rm -f "$ROOT"/tcrw_fig1b_traj_w*.txt
  rm -f "$ROOT"/tcrw_fig1b.pdf
  rm -f "$ROOT"/tcrw_fig1b.png
  rm -f "$ROOT"/tcrw_fig1b_w0.pdf
  rm -f "$ROOT"/tcrw_fig1c_msd_w*.txt
  rm -f "$ROOT"/tcrw_fig1c.pdf
  rm -f "$ROOT"/tcrw_fig1d_msd_dr*_w*.txt
  rm -f "$ROOT"/tcrw_fig1d_D_dr*.txt
  rm -f "$ROOT"/tcrw_fig1d.pdf
  rm -f "$ROOT"/tcrw_fig2_occ_w*.txt
  rm -f "$ROOT"/tcrw_fig2_traj_w*.txt
  rm -f "$ROOT"/tcrw_fig2_Jtot_w*.txt
  rm -f "$ROOT"/tcrw_fig2_Jomega_w*.txt
  rm -f "$ROOT"/tcrw_fig2_JDr_w*.txt
  rm -f "$ROOT"/tcrw_fig2_occ.pdf
  rm -f "$ROOT"/tcrw_fig2_traj.pdf
  rm -f "$ROOT"/tcrw_fig2_currents.pdf
  rm -f "$ROOT"/tcrw_fig2_occ_defects.txt
  rm -f "$ROOT"/tcrw_fig2_traj_defects.txt
  rm -f "$ROOT"/tcrw_fig2_Jtot_defects.txt
  rm -f "$ROOT"/tcrw_fig2_Jomega_defects.txt
  rm -f "$ROOT"/tcrw_fig2_JDr_defects.txt
  rm -f "$ROOT"/tcrw_fig2_defects_layout.txt
  rm -f "$ROOT"/tcrw_fig2_defects_occ.pdf
  rm -f "$ROOT"/tcrw_fig2_defects_traj.pdf
  rm -f "$ROOT"/tcrw_fig2_defects_currents.pdf
  rm -f "$ROOT"/tcrw_fig3a_summary.txt
  rm -f "$ROOT"/tcrw_fig3a.pdf
  rm -f "$ROOT"/tcrw_fig3b_summary.txt
  rm -f "$ROOT"/tcrw_fig3b.pdf
  rm -f "$ROOT"/tcrw_fig3cde_summary.txt
  rm -f "$ROOT"/tcrw_fig3e_summary.txt
  rm -f "$ROOT"/tcrw_fig3c.pdf
  rm -f "$ROOT"/tcrw_fig3d.pdf
  rm -f "$ROOT"/tcrw_fig3e.pdf
  rm -f "$ROOT"/tcrw_fig3f_summary.txt
  rm -f "$ROOT"/tcrw_fig3f.pdf
  rm -f "$ROOT"/tcrw_fig3g_summary.txt
  rm -f "$ROOT"/tcrw_fig3g.pdf
  # (more cleanup rows will be added as we create more figures)
  echo "   done."
}

# ---------------------------------------------------------------------
# dispatch
# ---------------------------------------------------------------------
if [[ $# -lt 1 ]]; then
  echo "usage: bash run.sh <target>"
  echo "targets:"
  grep -E '^#   tcrw-' "$0" | sed 's/^# //'
  exit 1
fi

case "$1" in
  tcrw-check)           tcrw_check            ;;
  tcrw-fig1b)           tcrw_fig1b            ;;
  tcrw-plot-fig1b-w0)   tcrw_plot_fig1b_w0    ;;
  tcrw-plot-fig1b)      tcrw_plot_fig1b       ;;
  tcrw-plot-fig1b-pdf)  tcrw_plot_fig1b_pdf   ;;
  tcrw-plot-fig1b-png)  tcrw_plot_fig1b_png   ;;
  tcrw-fig1b-all)       tcrw_fig1b_all        ;;
  tcrw-fig1c)           tcrw_fig1c            ;;
  tcrw-plot-fig1c)      tcrw_plot_fig1c       ;;
  tcrw-plot-fig1c-pdf)  tcrw_plot_fig1c_pdf   ;;
  tcrw-plot-fig1c-qt)   tcrw_plot_fig1c_qt    ;;
  tcrw-fig1c-all)       tcrw_fig1c_all        ;;
  tcrw-fig1d)           tcrw_fig1d            ;;
  tcrw-plot-fig1d)      tcrw_plot_fig1d       ;;
  tcrw-plot-fig1d-pdf)  tcrw_plot_fig1d_pdf   ;;
  tcrw-plot-fig1d-qt)   tcrw_plot_fig1d_qt    ;;
  tcrw-fig1d-all)       tcrw_fig1d_all        ;;
  tcrw-fig2-clean)              tcrw_fig2_clean              ;;
  tcrw-plot-fig2-occ)           tcrw_plot_fig2_occ           ;;
  tcrw-plot-fig2-occ-pdf)       tcrw_plot_fig2_occ_pdf       ;;
  tcrw-plot-fig2-occ-qt)        tcrw_plot_fig2_occ_qt        ;;
  tcrw-plot-fig2-traj)          tcrw_plot_fig2_traj          ;;
  tcrw-plot-fig2-traj-pdf)      tcrw_plot_fig2_traj_pdf      ;;
  tcrw-plot-fig2-traj-qt)       tcrw_plot_fig2_traj_qt       ;;
  tcrw-plot-fig2-currents)      tcrw_plot_fig2_currents      ;;
  tcrw-plot-fig2-currents-pdf)  tcrw_plot_fig2_currents_pdf  ;;
  tcrw-plot-fig2-currents-qt)   tcrw_plot_fig2_currents_qt   ;;
  tcrw-fig2-clean-all)          tcrw_fig2_clean_all          ;;
  tcrw-fig2-defects)                    tcrw_fig2_defects                    ;;
  tcrw-plot-fig2-defects-occ)           tcrw_plot_fig2_defects_occ           ;;
  tcrw-plot-fig2-defects-occ-pdf)       tcrw_plot_fig2_defects_occ_pdf       ;;
  tcrw-plot-fig2-defects-occ-qt)        tcrw_plot_fig2_defects_occ_qt        ;;
  tcrw-plot-fig2-defects-traj)          tcrw_plot_fig2_defects_traj          ;;
  tcrw-plot-fig2-defects-traj-pdf)      tcrw_plot_fig2_defects_traj_pdf      ;;
  tcrw-plot-fig2-defects-traj-qt)       tcrw_plot_fig2_defects_traj_qt       ;;
  tcrw-plot-fig2-defects-currents)      tcrw_plot_fig2_defects_currents      ;;
  tcrw-plot-fig2-defects-currents-pdf)  tcrw_plot_fig2_defects_currents_pdf  ;;
  tcrw-plot-fig2-defects-currents-qt)   tcrw_plot_fig2_defects_currents_qt   ;;
  tcrw-fig2-defects-all)                tcrw_fig2_defects_all                ;;
  tcrw-fig3a)           tcrw_fig3a            ;;
  tcrw-plot-fig3a)      tcrw_plot_fig3a       ;;
  tcrw-plot-fig3a-pdf)  tcrw_plot_fig3a_pdf   ;;
  tcrw-plot-fig3a-qt)   tcrw_plot_fig3a_qt    ;;
  tcrw-fig3a-all)       tcrw_fig3a_all        ;;
  tcrw-fig3b)           tcrw_fig3b            ;;
  tcrw-plot-fig3b)      tcrw_plot_fig3b       ;;
  tcrw-plot-fig3b-pdf)  tcrw_plot_fig3b_pdf   ;;
  tcrw-plot-fig3b-qt)   tcrw_plot_fig3b_qt    ;;
  tcrw-fig3b-all)       tcrw_fig3b_all        ;;
  tcrw-fig3cde)         tcrw_fig3cde          ;;
  tcrw-plot-fig3c)      tcrw_plot_fig3c       ;;
  tcrw-plot-fig3c-pdf)  tcrw_plot_fig3c_pdf   ;;
  tcrw-plot-fig3c-qt)   tcrw_plot_fig3c_qt    ;;
  tcrw-plot-fig3d)      tcrw_plot_fig3d       ;;
  tcrw-plot-fig3d-pdf)  tcrw_plot_fig3d_pdf   ;;
  tcrw-plot-fig3d-qt)   tcrw_plot_fig3d_qt    ;;
  tcrw-plot-fig3e)      tcrw_plot_fig3e       ;;
  tcrw-plot-fig3e-pdf)  tcrw_plot_fig3e_pdf   ;;
  tcrw-plot-fig3e-qt)   tcrw_plot_fig3e_qt    ;;
  tcrw-fig3cde-all)     tcrw_fig3cde_all      ;;
  tcrw-fig3f)           tcrw_fig3f            ;;
  tcrw-plot-fig3f)      tcrw_plot_fig3f       ;;
  tcrw-plot-fig3f-pdf)  tcrw_plot_fig3f_pdf   ;;
  tcrw-plot-fig3f-qt)   tcrw_plot_fig3f_qt    ;;
  tcrw-fig3f-all)       tcrw_fig3f_all        ;;
  tcrw-fig3g)           tcrw_fig3g            ;;
  tcrw-plot-fig3g)      tcrw_plot_fig3g       ;;
  tcrw-plot-fig3g-pdf)  tcrw_plot_fig3g_pdf   ;;
  tcrw-plot-fig3g-qt)   tcrw_plot_fig3g_qt    ;;
  tcrw-fig3g-all)       tcrw_fig3g_all        ;;
  tcrw-clean)           tcrw_clean            ;;
  *)
    echo "unknown target: $1"
    echo "available targets:"
    grep -E '^#   tcrw-' "$0" | sed 's/^# //'
    exit 1
    ;;
esac
