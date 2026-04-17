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
#   tcrw-plot-fig1b         plot Fig 1(b) (qt + PDF, all 4 ω in 1×4 row)
#   tcrw-plot-fig1b-pdf     plot Fig 1(b) (PDF only; headless)
#   tcrw-plot-fig1b-png     plot Fig 1(b) (PNG only; easy visual diff)
#   tcrw-fig1b-all          tcrw-fig1b → tcrw-plot-fig1b
#   tcrw-clean              remove .build/ and generated .txt / .pdf
#
# (Figs 1c, 1d, 2, 3, ... will get their own target rows as we go.)
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

tcrw_clean() {
  echo "==> cleaning .build and generated outputs"
  rm -rf "$BUILD"
  rm -f "$ROOT"/tcrw_fig1b_traj_w*.txt
  rm -f "$ROOT"/tcrw_fig1b.pdf
  rm -f "$ROOT"/tcrw_fig1b.png
  rm -f "$ROOT"/tcrw_fig1b_w0.pdf
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
  tcrw-clean)           tcrw_clean            ;;
  *)
    echo "unknown target: $1"
    echo "available targets:"
    grep -E '^#   tcrw-' "$0" | sed 's/^# //'
    exit 1
    ;;
esac
