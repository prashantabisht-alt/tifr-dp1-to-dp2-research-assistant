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
#   tcrw-clean              remove .build/ and generated .txt / .pdf
#
# (Fig 2, 3, ... will get their own target rows as we go.)
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
  tcrw-clean)           tcrw_clean            ;;
  *)
    echo "unknown target: $1"
    echo "available targets:"
    grep -E '^#   tcrw-' "$0" | sed 's/^# //'
    exit 1
    ;;
esac
